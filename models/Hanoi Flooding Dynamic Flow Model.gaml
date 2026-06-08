/**
* Name: Hanoi Flooding Dynamic Flow Model
* Author: Thanh-Do Nguyen / ARCHIVES Project
*
* Technique: MASS-CONSERVING CELLULAR FLOW (depth routing), 1 step = 1 HOUR.
*   Every cell carries its OWN water DEPTH; the water-surface elevation is
*   WSE = ground + depth. Each pass, water physically MOVES from higher-WSE cells to
*   lower-WSE neighbours (proportional to head, diagonals weighted), conserving volume -
*   nothing is created by "filling to a level". The river feeds the system at its stage.
*
*   So the surface is genuinely NON-UNIFORM and volume-limited: deep near the river and
*   the breaches, shallower as you move away, and each basin settles at the level the
*   water that actually reached it can support - not a blanket fill to the peak stage.
*   Spreading is paced by how fast water can flow there, so far basins fill late and
*   stay shallow until enough water arrives.
*
* Dykes (recorded breaks hold strict; other dykes overtop):
*   A BREAK="YES" dyke (Dykes.shp) is a STRICT wall until its recorded DATE (DD-MM 1926):
*   the river cannot overtop it even where the 50 m DEM made it too short, so the recorded
*   breaches drive when water crosses the MAIN dykes. On its date it is carved to the
*   channel floor and opens. Every OTHER dyke holds to its native DEM crest and is
*   OVERTOPPED when water rises above it from ANY side - so the river spills over low
*   dykes as it rises, and water filling a basin behind a breach overtops the NEXT dyke
*   and cascades on, basin to basin. (allow_overtopping OFF makes every dyke a strict wall.)
*
* Source / boundary: the river (RedRiver1925) is the ONLY source; its channel is
*   carved flat and its cells are pinned each hour to the river stage, which rises
*   from the start level to the 1926 Hanoi (Long Bien) peak of 11.93 m, then holds.
*   Lakes (Lakes1925) are DISPLAY ONLY.
*
* Terrain: mnt-gz50.asc (50 m, Vietnam TM + matching .prj; shares the shapefiles'
*          CRS so river/lakes/dykes align with the grid).
*/

model HanoiFloodingDynamicFlowModel

global {

    // === TERRAIN (50 m DEM, georeferenced) ===
    file dem_file <- grid_file("../includes/mnt-gz50.asc");
    geometry shape <- envelope(dem_file);
    field elevation_map <- field(dem_file);

    // === SHAPEFILES ===
    file river_shapefile <- file("../includes/RedRiver1925.shp"); // SOURCE
    file lakes_shapefile <- file("../includes/Lakes1925.shp");    // display only
    file dykes_shape_file <- shape_file("../includes/Dykes.shp"); // breaches + display

    // === TIME: 1 cycle = 1 hour, anchored to the 1926 flood ===
    date starting_date <- date([1926, 7, 22, 0, 0, 0]);
    float step <- 1 #h;

    // === SOURCE (river stage hydrograph) PARAMETERS ===
    float start_level  <- 0.0   min: 0.0 max: 15.0;   // initial river-stage level (m)
    float target_level <- 11.93 min: 0.0 max: 15.0;   // peak stage = 1926 Hanoi (m)
    float rise_rate    <- 0.0775 min: 0.0 max: 5.0 step: 0.005; // level rise per hour (m/h)

    // === FLOW PARAMETERS (mass-conserving depth routing) ===
    int   spread_iters <- 3 min: 1 max: 30;     // flow passes per hour (front speed)
    int   init_spread_passes <- 6;              // passes run once at t=0 to wet the channel
    // flow_relax must stay <= ~0.146: per pass a low cell receives up to ~6.83*flow_relax*head
    // (diagonal-weighted 8-neighbour), so beyond that a pit overshoots into a "water tower".
    float flow_relax        <- 0.12;            // share of head moved per pass (stable)
    float min_visible_depth <- 0.01;            // m, render / "is flooded" threshold
    float min_flow_depth    <- 0.001;           // m, a cell shallower than this is not a flow source

    // === DYKE / BREACH PARAMETERS ===
    // BREAK="YES" dykes are STRICT walls until their recorded break date (the river cannot
    // overtop them, so the recorded breaches drive when water crosses the main dykes).
    // Every OTHER dyke holds to its native DEM crest and is OVERTOPPED when water rises
    // above it (from any side: the river rising, or a basin filling behind a breach then
    // cascading over the next dyke). Set allow_overtopping OFF to make every dyke a strict
    // wall instead. Native DEM crest heights are used as-is (nothing is raised).
    bool  allow_overtopping <- true;
    bool  enable_breaches   <- true;
    int   breach_foot_radius <- 1 min: 0 max: 5;        // extra cell rings opened around a breach

    // === STOPPING ===
    float equilibrium_tol <- 0.02;   // % flooded growth in an hour below which it counts as "no progress"
    int   settle_patience <- 48;     // stop only after this many CONSECUTIVE no-progress hours, at peak
                                     // and with no pending breach. A breach restarts the clock.
    int   max_hours       <- 600;    // safety cap (25 days) so a run always terminates

    // === MONITORING POINTS (interactive arrival-time logging) ===
    bool   marking_mode <- false;
    int    marker_seq   <- 0;
    string markers_csv  <- "../exported_results/flood_arrival_times_flow.csv";
    bool   load_fixed_points <- true;
    file   fixed_points_file <- file("../includes/5_arrival_time.shp");

    // === STATE ===
    float source_level <- start_level;   // current river stage (set in init, rises each hour)
    bool  flooding     <- true;
    bool  peak_logged  <- false;
    bool  done_logged  <- false;
    int   last_progress_hour <- 0;       // last hour the flood grew by >= equilibrium_tol

    // === DERIVED ===
    field water_field;        // the NON-FLAT water surface (ground + depth on wet cells)
    field inflow_field;       // per-pass inflow accumulator (synchronous mass transfer)
    list<cell> active_cells;
    list<cell> river_cells;
    list<cell> wet_cells;     // maintained working set: cells with water (wet = true)
    int   n_active;           // cached length(active_cells)
    float flooded_pct      <- 0.0;   // cached once per cycle
    float prev_flooded_pct <- 0.0;
    float max_altitude;
    float min_altitude;

    init {
        write "=== Hanoi Flooding Dynamic Flow Model (mass-conserving flow) ===";
        write "DEM: " + elevation_map.columns + " x " + elevation_map.rows
              + " cells @ 50 m (Vietnam TM, matches shapefiles)";
        write "Time: 1 step = 1 hour from " + starting_date;
        write "Stage: start " + start_level + " m -> peak " + target_level + " m"
              + " | flow " + spread_iters + " passes/hour, relax " + flow_relax + ", rise " + rise_rate + " m/h";

        water_field  <- field(elevation_map.columns, elevation_map.rows);
        inflow_field <- field(elevation_map.columns, elevation_map.rows);

        // Read terrain; flag NODATA (<= -9000) cells as inactive.
        ask cell {
            ground <- elevation_map[grid_x, grid_y];
            is_inactive <- ground <= -9000.0;
            if (is_inactive) { ground <- 0.0; }
            ground0 <- ground;     // native DEM elevation (before any dyke raising / carving)
        }
        active_cells <- cell where !(each.is_inactive);
        ask active_cells {
            flow_neighbours <- (self neighbors_at 1) where !(each.is_inactive);
        }

        // SOURCE = river only.
        create river from: river_shapefile { cells_concerned <- cell overlapping self; }
        river_cells <- remove_duplicates(river accumulate (each.cells_concerned))
                         where !(each.is_inactive);
        ask river_cells { is_river <- true; }

        // Lakes = display only.
        create lake from: lakes_shapefile { cells_concerned <- cell overlapping self; }

        // Dykes = HARD barriers (baked crest + blocked) + recorded breach metadata.
        create digue from: dykes_shape_file with: [
            will_break::(string(read("BREAK")) = "YES"),
            break_date_str::string(read("DATE")),
            commune::string(read("Commune"))
        ];
        list<cell> dyke_cells <- remove_duplicates(digue accumulate (each.cells_concerned))
                                   where (!each.is_inactive and !each.is_river);
        ask dyke_cells { is_dyke <- true; }
        // BREAK="YES" dykes are STRICT walls until their recorded date: the river cannot
        // overtop them even where the 50 m DEM made them too short, so the recorded
        // breaches drive when water crosses the MAIN dykes. Other dykes overtop normally.
        list<cell> breaking_cells <- remove_duplicates((digue where each.will_break) accumulate (each.cells_concerned))
                                       where (!each.is_inactive and !each.is_river);
        ask breaking_cells { blocked <- true; }
        if (!allow_overtopping) { ask dyke_cells { blocked <- true; } }   // full-strict: every dyke a wall
        int nd <- length(dyke_cells);
        int dhold <- length(dyke_cells where (each.ground0 >= target_level));
        write "Dyke cells: " + nd + " (" + length(breaking_cells)
              + " on breaking dykes, held strict until their date) | native crest (raw DEM): min "
              + ((dyke_cells min_of each.ground0) with_precision 1) + " / mean "
              + ((dyke_cells mean_of each.ground0) with_precision 1) + " / max "
              + ((dyke_cells max_of each.ground0) with_precision 1) + " m | "
              + dhold + "/" + nd + " (" + ((100.0 * dhold / nd) with_precision 0) + "%) sit above the stage";
        write "Dyke mode: " + (allow_overtopping
              ? "breaking dykes STRICT until date; others OVERTOP (spill above crest) + cascade"
              : "STRICT (every dyke a wall; cross only where broken)");

        // Hydro-condition the river channel: burn it to one flat low bed so the
        // source reads as a single smooth water body (the polygon catches high banks).
        float river_floor <- empty(river_cells) ? start_level
                              : ((river_cells min_of (each.ground)) + 2.0);
        ask river_cells {
            ground <- min(ground, river_floor);
            elevation_map[grid_x, grid_y] <- ground;
        }
        source_level <- max(start_level, river_floor);

        max_altitude <- active_cells max_of (each.ground);
        min_altitude <- active_cells min_of (each.ground);
        n_active <- length(active_cells);
        write "Elevation range (after baking dykes / carving channel): "
              + (min_altitude with_precision 2) + " .. " + (max_altitude with_precision 2) + " m";
        write "Active cells: " + length(active_cells)
              + " | River (source) cells: " + length(river_cells)
              + " | Dyke cells: " + length(dyke_cells)
              + " | Dyke segments: " + length(digue)
              + " (recorded breaches BREAK=YES: " + length(digue where each.will_break) + ")";
        write "River carved to flat bed " + (river_floor with_precision 2)
              + " m; initial stage " + (source_level with_precision 2) + " m";

        // Seed the source (river depth = stage above its bed) and let it spread.
        ask river_cells { depth <- max(0.0, source_level - ground); }
        wet_cells <- active_cells where (each.depth > min_flow_depth);
        loop times: init_spread_passes { do spread_step; }
        do update_flooded_pct;
        do update_water_field;
        write "Initial flooded: " + flooded_percent() + " %";

        // Fixed comparison points (same locations for every model).
        if (load_fixed_points and fixed_points_file != nil) {
            create marker from: fixed_points_file with: [id::int(read("id")), from_file::true];
            marker_seq <- empty(marker) ? 0 : (marker max_of (each.id));
            write "Loaded " + length(marker) + " fixed monitoring points from 5_arrival_time.shp";
        }
    }

    // Fire each historical breach once its recorded date arrives (runs first).
    reflex trigger_dyke_breaches when: enable_breaches {
        int broken_before <- length(digue where each.has_broken);
        ask digue where (each.will_break and !each.has_broken and each.breach_date != nil) {
            if (current_date >= breach_date) { do break_dyke; }
        }
        if (length(digue where each.has_broken) > broken_before) {
            last_progress_hour <- cycle;   // restart the settle clock: a breach just opened
        }
    }

    // One simulated hour: raise the river stage, then advance the mass-conserving flow.
    reflex flood when: flooding {
        source_level <- min(target_level, source_level + rise_rate);
        if (!peak_logged and source_level >= target_level) {
            peak_logged <- true;
            write "River stage reached peak " + (target_level with_precision 2)
                  + " m at " + current_date + " (now holding; waiting on breaches)";
        }
        // feed the source at the new stage and keep the river in the working set
        ask river_cells { depth <- max(0.0, source_level - ground); }
        wet_cells <- remove_duplicates(wet_cells + river_cells);

        prev_flooded_pct <- flooded_pct;
        loop times: spread_iters { do spread_step; }
        do update_flooded_pct;
        do update_water_field;

        if ((flooded_pct - prev_flooded_pct) >= equilibrium_tol) { last_progress_hour <- cycle; }
        bool pending <- enable_breaches and
            !empty(digue where (each.will_break and !each.has_broken and each.breach_date != nil));
        if (source_level >= target_level and !pending
                and (cycle - last_progress_hour) >= settle_patience) {
            flooding <- false;
            write "Flood settled at peak stage (" + current_date + ") -> "
                  + flooded_percent() + " % flooded (no progress for " + settle_patience + " h)";
            do report_connectivity;
        }
        if (cycle >= max_hours) {
            flooding <- false;
            write "Reached max_hours safety (" + max_hours + " h) at " + current_date + ". Pausing.";
            do report_connectivity;
        }
    }

    reflex pause_at_end when: !flooding {
        if (!done_logged) {
            done_logged <- true;
            write "Done at hour " + cycle + " | flooded " + flooded_percent() + " %.";
        }
        do pause;
    }

    // ----- ONE mass-conserving flow pass (synchronous) --------------------------
    // PHASE 1: each wet, non-blocked cell pushes water DOWN-gradient (toward lower
    //   water-surface neighbours) in proportion to head, diagonals weighted by 1/sqrt(2);
    //   the total it gives is capped by its own depth and by flow_relax*head, so volume
    //   is conserved and no pit is over-filled. Contributions land in inflow_field.
    // PHASE 2: depth += inflow - outflow; the river is re-fed to the stage (the source);
    //   blocked (unbroken breaking) dykes never change.
    action spread_step {
        ask wet_cells {
            if (!blocked) {
                float my_wse <- ground + depth;
                float total_head <- 0.0;
                loop nb over: flow_neighbours {
                    if (!nb.blocked) {
                        float hd <- my_wse - (nb.ground + nb.depth);
                        if (hd > 0.0) {
                            float w <- ((nb.grid_x != grid_x) and (nb.grid_y != grid_y)) ? 0.70710678 : 1.0;
                            total_head <- total_head + w * hd;
                        }
                    }
                }
                if (total_head > 0.0) {
                    float give <- min(depth, flow_relax * total_head);
                    float kk <- give / total_head;
                    loop nb over: flow_neighbours {
                        if (!nb.blocked) {
                            float hd <- my_wse - (nb.ground + nb.depth);
                            if (hd > 0.0) {
                                float w <- ((nb.grid_x != grid_x) and (nb.grid_y != grid_y)) ? 0.70710678 : 1.0;
                                inflow_field[nb.grid_x, nb.grid_y] <-
                                    inflow_field[nb.grid_x, nb.grid_y] + kk * w * hd;
                            }
                        }
                    }
                    outflow <- give;
                }
            }
        }
        // apply over the wet front + the neighbours that may have received water
        list<cell> touched <- remove_duplicates(wet_cells + (wet_cells accumulate (each.flow_neighbours)));
        ask touched {
            if (is_river) {
                depth <- max(0.0, source_level - ground);   // infinite feed at the river stage
            } else if (!blocked) {
                depth <- max(0.0, depth - outflow + inflow_field[grid_x, grid_y]);
            }
            outflow <- 0.0;
            inflow_field[grid_x, grid_y] <- 0.0;
        }
        wet_cells <- touched where (each.depth > min_flow_depth);
    }

    // Render the NON-FLAT surface: wet cells show their level; dry cells park just
    // below THEIR OWN terrain so they hide under the terrain mesh.
    action update_water_field {
        ask active_cells {
            water_field[grid_x, grid_y] <- (depth > min_visible_depth)
                                           ? (ground + depth) : (ground - 1.0);
        }
    }

    // Cached once per cycle, counted within the (flood-sized) wet set.
    action update_flooded_pct {
        flooded_pct <- ((length(wet_cells where (each.depth > min_visible_depth))) / n_active)
                       * 100.0;
    }
    float flooded_percent {
        return flooded_pct with_precision 2;
    }

    // Diagnostic: how much below-stage area is reachable from the river, ignoring
    // the dyke walls vs. respecting them (current breach state), vs. what actually
    // flooded. Pinpoints whether the cap is terrain, dyke walls, or the flow.
    action report_connectivity {
        int ceiling <- length(active_cells where (each.ground < target_level));
        // BFS 1: ignore dyke walls (pure terrain connectivity below the stage)
        ask active_cells { bfs_mark <- false; }
        ask river_cells  { bfs_mark <- true; }
        list<cell> fr <- river_cells;
        loop while: !empty(fr) {
            list<cell> nx <- remove_duplicates(fr accumulate (each.flow_neighbours))
                              where (!each.bfs_mark and (each.ground < target_level));
            ask nx { bfs_mark <- true; }
            fr <- nx;
        }
        int ri <- length(active_cells where each.bfs_mark);
        // BFS 2: respect current dyke walls (blocked cells impassable)
        ask active_cells { bfs_mark <- false; }
        ask river_cells  { bfs_mark <- true; }
        list<cell> fr2 <- river_cells;
        loop while: !empty(fr2) {
            list<cell> nx <- remove_duplicates(fr2 accumulate (each.flow_neighbours))
                              where (!each.bfs_mark and !each.blocked and (each.ground < target_level));
            ask nx { bfs_mark <- true; }
            fr2 <- nx;
        }
        int rw <- length(active_cells where each.bfs_mark);
        // BFS 3: treat every dyke RIDGE as an open gap (passable). Shows how much
        // floodplain sits behind the ridges (what a breach could reach) vs. what
        // is cut off by natural high terrain regardless of dykes.
        ask active_cells { bfs_mark <- false; }
        ask river_cells  { bfs_mark <- true; }
        list<cell> fr3 <- river_cells;
        loop while: !empty(fr3) {
            list<cell> nx <- remove_duplicates(fr3 accumulate (each.flow_neighbours))
                              where (!each.bfs_mark and (each.is_dyke or (each.ground < target_level)));
            ask nx { bfs_mark <- true; }
            fr3 <- nx;
        }
        int ro <- length(active_cells where each.bfs_mark);
        int fl <- length(wet_cells where (each.depth > min_visible_depth));
        write "=== CONNECTIVITY @ stage " + (target_level with_precision 2)
              + " m (active " + n_active + ") ===";
        write "  below stage (ceiling)          : " + ceiling
              + " = " + ((100.0 * ceiling / n_active) with_precision 1) + " %";
        write "  reachable IGNORING dyke walls   : " + ri
              + " = " + ((100.0 * ri / n_active) with_precision 1) + " %";
        write "  reachable WITH dyke walls (now) : " + rw
              + " = " + ((100.0 * rw / n_active) with_precision 1) + " %";
        write "  reachable if dyke RIDGES open    : " + ro
              + " = " + ((100.0 * ro / n_active) with_precision 1) + " %";
        write "  actually flooded                : " + fl
              + " = " + ((100.0 * fl / n_active) with_precision 1) + " %";
        write "  => behind dyke ridges: " + (ro - rw)
              + " cells; reachable-but-unfilled " + (rw - fl) + " cells";
    }

    // --- interactive monitoring points ---------------------------------------
    action toggle_marking {
        marking_mode <- !marking_mode;
        write marking_mode
            ? "MARKING ON  -> click on the map to drop monitoring points"
            : "MARKING OFF";
    }

    action clear_markers {
        ask marker { do die; }
        marker_seq <- 0;
        write "All monitoring points cleared.";
    }

    action report_markers {
        list<marker> ms <- marker sort_by each.id;
        write "=== FLOOD ARRIVAL TIMES (" + length(ms) + " points) ===";
        string csv <- "id,source,easting,northing,reached,arrival_hour,arrival_date\n";
        loop m over: ms {
            write "  Point #" + m.id + (m.from_file ? " (fixed)" : " (click)")
                + " (" + int(m.location.x) + ", " + int(m.location.y) + "): "
                + (m.reached ? ("hour " + m.reached_hour + "  (" + m.reached_date + ")")
                             : "NOT reached");
            csv <- csv + ("" + m.id + "," + (m.from_file ? "fixed" : "click") + ","
                  + int(m.location.x) + "," + int(m.location.y) + ","
                  + m.reached + "," + m.reached_hour + ","
                  + (m.reached ? string(m.reached_date) : "") + "\n");
        }
        save csv to: markers_csv rewrite: true;
        write "Saved -> " + markers_csv;
    }
}

// Terrain + level-pool grid. Each cell carries its own water level, so the water
// surface is non-flat.
grid cell file: dem_file neighbors: 8 frequency: 0
     use_regular_agents: false use_individual_shapes: false {
    float ground;                  // terrain elevation (carved channel baked in)
    float ground0     <- 0.0;      // native DEM elevation (before carving)
    bool  is_inactive <- false;    // NODATA
    bool  is_river    <- false;    // perennial source (fed at the stage)
    bool  is_dyke     <- false;    // dyke footprint
    bool  blocked     <- false;    // impermeable (unbroken breaking dyke); cleared on breach
    float depth       <- 0.0;      // water depth (mass-conserving); WSE = ground + depth
    float outflow     <- 0.0;      // water given away this pass (synchronous-update scratch)
    bool  bfs_mark    <- false;    // scratch for the connectivity diagnostic
    list<cell> flow_neighbours;    // active (non-NODATA) 8-neighbours, cached
}

// River = water source
species river {
    list<cell> cells_concerned;
    aspect geometry { draw shape color: rgb("blue") depth: 1.0; }
}

// Lake = display only
species lake {
    list<cell> cells_concerned;
    aspect geometry { draw shape color: rgb("green") depth: 0.5; }
}

// Dyke = hard barrier (baked crest + blocked) that opens on its recorded date
species digue {
    list<cell> cells_concerned;
    bool   will_break <- false;     // BREAK = "YES"
    string break_date_str;          // DATE = "DD-MM"
    string commune;
    date   breach_date;
    bool   has_broken <- false;

    init {
        cells_concerned <- cell overlapping self;
        if (will_break and break_date_str != nil and break_date_str != "") {
            list<string> p <- break_date_str split_with "-";
            if (length(p) >= 2) {
                breach_date <- date([1926, int(p[1]), int(p[0]), 6, 0, 0]);
            }
        }
    }

    // Open the gap: footprint (+ foot rings) carved down to the adjacent channel
    // floor and unblocked, so the protected basin fills through it from this hour.
    action break_dyke {
        has_broken <- true;
        list<cell> breach_cells <- cells_concerned where !(each.is_inactive);
        loop times: world.breach_foot_radius {
            breach_cells <- remove_duplicates(
                breach_cells + (breach_cells accumulate (each.flow_neighbours)));
        }
        breach_cells <- breach_cells where !(each.is_inactive);
        list<cell> ring <- remove_duplicates(breach_cells accumulate (each.flow_neighbours))
                             where !(breach_cells contains each);
        float floor_alt <- empty(ring) ? (breach_cells min_of each.ground)
                                        : (ring min_of each.ground);
        ask breach_cells {
            ground  <- min(ground, floor_alt);
            blocked <- false;                              // now permeable
            world.elevation_map[grid_x, grid_y] <- ground;
        }
        write "DYKE BREACH @ " + current_date + " | " + commune
              + " | opened " + length(breach_cells) + " cells down to "
              + (floor_alt with_precision 2) + " m";
    }

    aspect geometry {
        draw shape color: (has_broken ? rgb("orange") : rgb("red")) depth: 2.0;
    }
}

// Interactive monitoring point: logs the FIRST hour the cell under it gets wet.
species marker {
    int  id;
    bool from_file <- false;        // true = loaded from 5_arrival_time.shp (fixed)
    bool reached <- false;
    int  reached_hour <- -1;
    date reached_date;

    reflex detect when: !reached {
        cell c <- cell closest_to self;
        if (c != nil and (c.depth > world.min_visible_depth)) {
            reached <- true;
            reached_hour <- cycle;
            reached_date <- current_date;
            write "FLOOD reached point #" + id + (from_file ? " (fixed)" : "")
                  + " at hour " + cycle + "  (" + current_date + ")";
        }
    }

    aspect default {
        cell c <- cell closest_to self;
        float z <- (c != nil ? c.ground : 0.0);
        rgb col <- reached ? #red : (from_file ? #cyan : #yellow);
        draw line([{location.x, location.y, z}, {location.x, location.y, z + 150}]) color: #black;
        draw sphere(110) at: {location.x, location.y, z + 150} color: col;
        draw (string(id) + (reached ? (" h" + reached_hour) : ""))
             at: {location.x, location.y, z + 320} color: #white font: font("Arial", 12, #bold);
    }
}

experiment HanoiFloodSimulation type: gui {

    parameter "Start stage (m)"          var: start_level      category: "Source";
    parameter "Peak stage (m)"           var: target_level     category: "Source";
    parameter "Level rise per hour (m/h)" var: rise_rate       category: "Source";
    parameter "Flow passes per hour"     var: spread_iters     category: "Flow";
    parameter "Flow relax (share/pass)"  var: flow_relax       category: "Flow";
    parameter "Settle patience (h)"      var: settle_patience  category: "Flow";
    parameter "Allow overtopping (off = strict walls)" var: allow_overtopping category: "Dyke";
    parameter "Enable dyke breaches"     var: enable_breaches  category: "Dyke";
    parameter "Breach foot radius"       var: breach_foot_radius category: "Dyke";
    parameter "Load fixed points (5_arrival_time)" var: load_fixed_points category: "Points";

    output {

        // Main 3D view: terrain field + non-flat flood-water field
        display "Flood 3D" type: opengl {
            mesh elevation_map scale: 1 grayscale: true smooth: false triangulation: true
                 refresh: true;     // refresh so breach notches show
            mesh water_field   scale: 1 color: rgb(100, 150, 255, 180) smooth: false
                 triangulation: true refresh: true;
            species river aspect: geometry;
            species lake aspect: geometry;
            species digue aspect: geometry;
            species marker aspect: default;
            light #ambient intensity: 130;

            event #mouse_down {
                point clicked <- #user_location;
                ask world {
                    if (marking_mode) {
                        marker_seq <- marker_seq + 1;
                        create marker with: [location::clicked, id::marker_seq];
                        write "Monitoring point #" + marker_seq + " placed at ("
                              + int(clicked.x) + ", " + int(clicked.y) + ").";
                    }
                }
            }

            overlay position: {5, 5} size: {360, 150} background: #black
                    transparency: 0.25 border: #white {
                draw "HANOI 1926 — Level-pool (non-flat)" at: {10, 18}
                     color: #yellow font: font("Arial", 13, #bold);
                draw ("" + current_date) at: {10, 40} color: #white;
                draw ("River stage: " + (world.source_level with_precision 2) + " m")
                     at: {10, 56} color: #cyan;
                draw ("Flooded: " + world.flooded_percent() + " %")
                     at: {10, 72} color: #deepskyblue;
                draw (flooding ? "SPREADING" : "DONE") at: {10, 88}
                     color: (flooding ? #orange : #lime);
                draw ("Breaches opened: " + length(digue where each.has_broken)
                     + " / " + length(digue where each.will_break)) at: {10, 104}
                     color: #tomato;
                draw ("Marking: " + (marking_mode ? "ON (click to add)" : "off")
                     + "   points: " + length(marker)
                     + " (" + length(marker where each.reached) + " reached)")
                     at: {10, 124} color: (marking_mode ? #lime : #lightgray);
            }
        }

        monitor "Date"             value: current_date;
        monitor "Hour"             value: cycle;
        monitor "River stage (m)"  value: source_level with_precision 2;
        monitor "Flooded %"        value: world.flooded_percent();
        monitor "Breaches opened"  value: length(digue where each.has_broken);
        monitor "Marking mode"     value: marking_mode;
        monitor "Points placed"    value: length(marker);
        monitor "Points reached"   value: length(marker where each.reached);
    }

    user_command "Toggle marking points" action: {ask world {do toggle_marking;}};
    user_command "Clear points"          action: {ask world {do clear_markers;}};
    user_command "Report arrival times"  action: {ask world {do report_markers;}};
    user_command "Report connectivity"   action: {ask world {do report_connectivity;}};
}
