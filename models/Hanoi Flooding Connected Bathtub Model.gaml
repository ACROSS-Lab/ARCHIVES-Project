/**
* Name: Hanoi Flooding Connected Bathtub Model
* Author: Thanh-Do Nguyen / ARCHIVES Project
*
* Technique: CONNECTED BATHTUB with SPREAD-THEN-RISE dynamics (1 step = 1 HOUR).
*   Each hour the flood first SPREADS to equilibrium: from the RIVER source it
*   fills every connected cell whose terrain is below the current water level.
*   Only once it can spread no further does the level RISE - and it rises to the
*   next spill point (lowest dry lip touching the water), so a basin fills,
*   overflows its lowest rim, fills the next basin, and so on. Low pockets not
*   connected to the river stay dry (unlike a naive bathtub).
*
* Recorded dyke breaches (date-triggered, from Dykes.shp):
*   A dyke with BREAK = "YES" is a solid barrier (its crest is baked into the
*   DEM) UNTIL its recorded DATE (DD-MM, 1926). On that date its footprint is
*   carved down to the adjacent channel floor and water can pass through it.
*   Other dykes never open. While any dated breach is still pending, the water
*   holds at the peak level instead of finishing, so the east bank floods only
*   after the historical breaks (28-29 Jul 1926).
*
* Conventions follow the repo's 3D models:
*   - River (RedRiver1925) is the ONLY water source.
*   - Lakes (Lakes1925) are loaded for DISPLAY ONLY (no water seeded).
*   - Dykes baked into the DEM crest (dem_contains_dykes); shapefile drives the
*     breaches and the red overlay.
*   - elevation_map / water_field fields drive the OpenGL mesh display.
*
* Terrain: DEM_Hanoi_50m.asc (50 m, Vietnam TM + matching .prj; 5x downsample of
*          DEM_Hanoi.tif). The .asc/.prj shares the shapefiles' CRS so they align.
* Context: 1926 Red-River flood; 11.93 m is the recorded Hanoi (Long Bien) peak
*          stage at which the left-bank dykes broke (28-29 Jul 1926).
*/

model HanoiFloodingConnectedBathtubModel

global {

    // === TERRAIN (new 50 m DEM, georeferenced) ===
    // Use the .asc (+ matching .prj) so the DEM shares the shapefiles' CRS exactly
    // (Vietnam TM, k=0.9999). Loading the .tif tagged EPSG:32648 instead makes GAMA
    // reproject the river/lakes off the grid, leaving almost no source cells.
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

    // === FLOOD (connected bathtub) PARAMETERS ===
    float start_level  <- 0.0   min: 0.0 max: 15.0;   // initial water-surface level (m)
    float target_level <- 11.93 min: 0.0 max: 15.0;   // peak level = 1926 Hanoi stage (m)
    float rise_rate    <- 0.0775  min: 0.0  max: 5.0 step: 0.01;    // MIN rise when the level does rise (m/h)
    float min_visible_depth <- 0.01;                   // m, render threshold

    // === BREACH PARAMETERS (date-triggered) ===
    bool enable_breaches <- true;
    int  breach_foot_radius <- 1 min: 0 max: 5;        // cells around the breach also opened

    // === STATE ===
    float water_level <- start_level;
    bool  flooding    <- true;

    // === DERIVED ===
    field water_field;
    list<cell> active_cells;
    list<cell> river_cells;
    float max_altitude;
    float min_altitude;
    int neighbours_type <- 8;

    init {
        write "=== Hanoi Flooding Connected Bathtub Model ===";
        write "DEM: " + elevation_map.columns + " x " + elevation_map.rows
              + " cells @ 50 m (Vietnam TM, matches shapefiles)";
        write "Time: 1 step = 1 hour from " + starting_date
              + " | spread-then-rise until peak reached and all breaches fired";
        write "Level: start " + start_level + " m -> target " + target_level + " m";

        water_field <- field(elevation_map.columns, elevation_map.rows);

        ask cell {
            altitude2 <- elevation_map[grid_x, grid_y];
            is_inactive <- altitude2 <= -9000.0;
            if (is_inactive) { altitude2 <- 0.0; }
        }
        active_cells <- cell where !(each.is_inactive);
        ask active_cells {
            active_neighbours <- (self neighbors_at 1) where !(each.is_inactive);
        }
        max_altitude <- active_cells max_of (each.altitude2);
        min_altitude <- active_cells min_of (each.altitude2);
        write "Elevation range: " + (min_altitude with_precision 2) + " .. "
              + (max_altitude with_precision 2) + " m";

        // Park the water surface below terrain everywhere (hidden on dry cells)
        ask cell { water_field[grid_x, grid_y] <- min_altitude - 1.0; }

        // SOURCE = river only
        create river from: river_shapefile { cells_concerned <- cell overlapping self; }
        river_cells <- remove_duplicates(river accumulate (each.cells_concerned))
                         where !(each.is_inactive);

        // Lakes = display only
        create lake from: lakes_shapefile { cells_concerned <- cell overlapping self; }

        // Dykes = barriers (in DEM) + recorded breach metadata
        create digue from: dykes_shape_file with: [
            will_break::(string(read("BREAK")) = "YES"),
            break_date_str::string(read("DATE")),
            commune::string(read("Commune"))
        ];

        write "Active cells: " + length(active_cells)
              + " | River (source) cells: " + length(river_cells)
              + " | Dyke segments: " + length(digue)
              + " (recorded breaches BREAK=YES: " + length(digue where each.will_break) + ")";

        // Hydro-condition the river channel. The DEM's river polygon spans beds
        // from ~4 m up to ~19 m (sandbars / banks caught in the polygon), so a
        // single water surface would show holes wherever a bed pokes above the
        // stage. Burn the whole channel down to one flat low stage so the river
        // reads as a single smooth, connected water body - the smooth flat source
        // the BFS model seeds. The river is the perennial source: wet from t = 0.
        float river_stage <- empty(river_cells) ? start_level
                              : ((river_cells min_of (each.altitude2)) + 2.0);
        water_level <- max(start_level, river_stage);
        ask river_cells {
            altitude2 <- min(altitude2, river_stage);     // burn banks to the channel
            elevation_map[grid_x, grid_y] <- altitude2;   // keep the terrain mesh consistent
            is_river <- true;
            flooded <- true;
        }
        write "River carved to flat stage " + (river_stage with_precision 2)
              + " m; initial water level " + (water_level with_precision 2) + " m";

        // Spread the initial pool onto any land already below the start stage
        do grow_flood;
        do update_water_field;
        write "Initial flooded: " + flooded_percent() + " %";
    }

    // Fire each historical breach once its recorded date arrives (runs first).
    reflex trigger_dyke_breaches when: enable_breaches {
        ask digue where (each.will_break and !each.has_broken and each.breach_date != nil) {
            if (current_date >= breach_date) { do break_dyke; }
        }
    }

    // SPREAD-THEN-RISE: advance the flood by ONE ring per step; the level only
    // rises on a step where the water could spread no further at the current level.
    reflex flood when: flooding {
        // 1) SPREAD one ring: flood the immediate dry neighbours of the current
        //    flood (the river is already wet) that lie below the current level.
        list<cell> ring <- active_cells where (
            !each.flooded and each.altitude2 < water_level
            and !empty(each.active_neighbours where (each.flooded)));
        ask ring { flooded <- true; }

        // 2) Only when nothing spread this step is the flood at equilibrium -> RISE
        if (empty(ring)) {
            if (water_level < target_level) {
                water_level <- min(target_level, water_level + rise_rate);
            } else {
                bool pending <- enable_breaches and
                    !empty(digue where (each.will_break and !each.has_broken and each.breach_date != nil));
                if (!pending) {
                    flooding <- false;
                    write "Flood fully spread at " + (water_level with_precision 2)
                          + " m (" + current_date + ") -> " + flooded_percent() + " %";
                }
                // else: hold at peak and wait for the recorded breach dates
            }
        }

        do update_water_field;
    }

    reflex pause_at_end when: !flooding {
        write "Done at hour " + cycle + " | flooded " + flooded_percent() + " %. Pausing.";
        do pause;
    }

    // Connectivity-constrained level fill: region-grow from the river source
    // through neighbours lying below the current water level.
    action grow_flood {
        ask river_cells where (!each.flooded and each.altitude2 < water_level) {
            flooded <- true;
        }
        list<cell> frontier <- active_cells where (each.flooded);
        loop while: !empty(frontier) {
            list<cell> next <- remove_duplicates(frontier accumulate (each.active_neighbours))
                                 where (!each.flooded and each.altitude2 < water_level);
            ask next { flooded <- true; }
            frontier <- next;
        }
    }

    action update_water_field {
        // Flooded cells render as one flat sheet at water_level (carved river +
        // submerged land all sit below it). Dry cells are parked just below THEIR
        // OWN terrain (altitude2 - 1.0, like the BFS model) so they hide under the
        // terrain mesh and blend in, instead of plunging to a global low and
        // showing as dark pits at the flood edge.
        ask active_cells {
            water_depth <- flooded ? (water_level - altitude2) : 0.0;
            water_field[grid_x, grid_y] <- flooded ? water_level : (altitude2 - 1.0);
        }
    }

    float flooded_percent {
        return (((length(active_cells where (each.flooded))) / length(active_cells)) * 100.0)
               with_precision 2;
    }
}

// Terrain + flood grid (elevation comes from the DEM via grid_value)
grid cell file: dem_file neighbors: 8 frequency: 0
     use_regular_agents: false use_individual_shapes: false {
    float altitude2;
    bool  is_inactive <- false;
    bool  is_river <- false;
    bool  flooded <- false;
    float water_depth <- 0.0;
    list<cell> active_neighbours;
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

// Dyke = barrier baked into DEM + recorded breach (BREAK=YES on its DATE)
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

    // Open the gap: footprint (+ foot rings) dropped to the channel floor.
    action break_dyke {
        has_broken <- true;
        list<cell> breach_cells <- cells_concerned where !(each.is_inactive);
        loop times: world.breach_foot_radius {
            breach_cells <- remove_duplicates(
                breach_cells + (breach_cells accumulate (each.active_neighbours)));
        }
        breach_cells <- breach_cells where !(each.is_inactive);
        list<cell> ring <- remove_duplicates(breach_cells accumulate (each.active_neighbours))
                             where !(breach_cells contains each);
        float floor_alt <- empty(ring) ? (breach_cells min_of each.altitude2)
                                        : (ring min_of each.altitude2);
        ask breach_cells {
            altitude2 <- min(altitude2, floor_alt);
            world.elevation_map[grid_x, grid_y] <- altitude2;
        }
        write "DYKE BREACH @ " + current_date + " | " + commune
              + " | opened " + length(breach_cells) + " cells down to "
              + (floor_alt with_precision 2) + " m";
    }

    aspect geometry {
        draw shape color: (has_broken ? rgb("orange") : rgb("red")) depth: 2.0;
    }
}

experiment HanoiFloodSimulation type: gui {

    parameter "Start level (m)"         var: start_level     category: "Flood";
    parameter "Target level (m)"        var: target_level    category: "Flood";
    parameter "Min rise per step (m/h)" var: rise_rate        category: "Flood";
    parameter "Enable dyke breaches"    var: enable_breaches  category: "Breach";
    parameter "Breach foot radius"      var: breach_foot_radius category: "Breach";

    output {

        // Main 3D view: terrain field + flood-water field
        display "Flood 3D" type: opengl {
            mesh elevation_map scale: 1 grayscale: true smooth: false triangulation: true
                 refresh: true;     // refresh so breach notches show
            mesh water_field   scale: 1 color: rgb(100, 150, 255, 180) smooth: false
                 triangulation: true refresh: true;
            species river aspect: geometry;
            species lake aspect: geometry;
            species digue aspect: geometry;
            light #ambient intensity: 130;

            overlay position: {5, 5} size: {330, 128} background: #black
                    transparency: 0.25 border: #white {
                draw "HANOI 1926 — Connected Bathtub" at: {10, 18}
                     color: #yellow font: font("Arial", 13, #bold);
                draw ("" + current_date) at: {10, 40} color: #white;
                draw ("Water level: " + (water_level with_precision 2) + " m")
                     at: {10, 56} color: #cyan;
                draw ("Flooded: " + world.flooded_percent() + " %")
                     at: {10, 72} color: #deepskyblue;
                draw (flooding ? "RISING / SPREADING" : "DONE") at: {10, 88}
                     color: (flooding ? #orange : #lime);
                draw ("Breaches opened: " + length(digue where each.has_broken)
                     + " / " + length(digue where each.will_break)) at: {10, 104}
                     color: #tomato;
            }
        }

        monitor "Date"            value: current_date;
        monitor "Hour"            value: cycle;
        monitor "Water level (m)" value: water_level with_precision 2;
        monitor "Flooded %"       value: world.flooded_percent();
        monitor "Breaches opened" value: length(digue where each.has_broken);
    }
}
