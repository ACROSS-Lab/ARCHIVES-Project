/**
* Name: Hanoi FastBathtub Flood Model 1926 (Discharge-driven)
* Author: Thành Đô Nguyễn (2026-06-25)
*
* A fork of "Hanoi FastBathtub Flood Model 1926.gaml" that KEEPS the breach
* engine but DROPS the river-stage-file mechanism (the synthetic gaussian stage
* hydrograph), so the stage is driven SOLELY by the observed discharge through
* the Manning rating curve - exactly as the Đồng Hới and Phú Thọ models do.
*
* WHAT CHANGED vs the original Hanoi FastBathtub model, and WHY.
*   The dyke + breach scenario is IDENTICAL: dykes are hard barriers, water
*   reaches the plain ONLY through a breach corridor cut into a dyke on the
*   Dykes.shp BREAK/DATE dates, and the spill + front fields are rebuilt when a
*   breach opens. The ENGINE (spread-limited, connectivity-constrained,
*   hysteretic level-pool) is unchanged.
*
*   ONLY the forcing changed. The original model had TWO stage paths:
*     (1) observed discharge (WaterDischarge.csv) -> rating curve -> stage, and
*     (2) a SYNTHETIC GAUSSIAN stage hydrograph fallback (peak_date,
*         sigma_rise_days, sigma_fall_days) used when no discharge was loaded.
*   Path (2) - the "river stage file" / fabricated-stage-series mechanism - is
*   REMOVED here. The Đồng Hới and Phú Thọ models have no such synthetic stage
*   curve; their stage is the discharge pushed through the rating curve (with a
*   measured gauge as the only alternative). This model follows that pattern:
*   river_stage = rating(WaterDischarge.csv), the SOLE driver. If the discharge
*   record is ever missing the stage simply holds at base_stage (no synthetic
*   hydrograph is fabricated).
*
* Inputs (../includes/): DEM selectable via the "DEM file" parameter; default
*   mnt-gz50.asc. Plus RedRiver1925.shp, Buildings1925.shp, Lakes1925.shp,
*   5_arrival_time.shp, Dykes.shp (BREAK/DATE/Commune), WaterDischarge.csv.
*   NOTE: RedRiverStage1926_hourly.csv is deliberately NOT used.
*/
model HanoiFastBathtubFlood1926Discharge

global {

	// ------------------------------------------------------------------ input
	// DEM: default mnt-gz50.asc. mnt-gz40.asc = the ARCHIVAL 1926-datum MNT the ABM
	// runs on (so the two are directly comparable); mnt-gz50.asc = modern DEM used by
	// the V5 LISEM model (different vertical datum).
	string dem_name     <- "mnt-gz50.asc";
	file dem_file       <- grid_file("../includes/" + dem_name);
	file river_file     <- shape_file("../includes/RedRiver1925.shp");
	file buildings_file <- shape_file("../includes/Buildings1925.shp");
	file lakes_file     <- shape_file("../includes/Lakes1925.shp");
	file points_file    <- shape_file("../includes/5_arrival_time.shp");
	file dykes_file     <- shape_file("../includes/Dykes.shp");
	file discharge_file <- csv_file("../includes/WaterDischarge.csv", ",", true);
	geometry shape <- envelope(dem_file);

	// ------------------------------------------------------------------ time
	date starting_date <- date("1926-07-20 00:00:00");
	date end_date      <- date("1926-08-08 00:00:00");
	float step <- 1 #h;                       // ONE CYCLE = ONE HOUR

	// ------------------------------------------------------------------ river stage forcing
	// DISCHARGE-DRIVEN ONLY (no synthetic stage hydrograph). The observed daily
	// discharge (WaterDischarge.csv) is pushed through a Manning power-law rating
	// curve to a stage; the two anchors q_min -> base_stage and q_max -> peak_stage
	// set the curve. There is NO gaussian fallback: if the discharge record is
	// missing the stage holds at base_stage.
	bool  use_discharge_csv <- true;
	float rating_exponent <- 0.6;   // Manning h ~ Q^(3/5)
	float base_stage <- 7.0;        // m, stage at the lowest recorded discharge (~22-07)
	float peak_stage <- 16.00;      // m, observed 1926 peak at Hanoi (Gourou fig. 9)
	float river_stage <- base_stage;
	float river_discharge <- 0.0;
	list<date>  q_dates  <- [];
	list<float> q_values <- [];
	float q_min <- 1.0; float q_max <- 2.0;

	// ------------------------------------------------------------------ breach scenario
	// (identical geometry to the original model: invert from the lowest
	//  protected-side ground, corridor cut through the embankment width)
	int   breach_hour <- 6;
	float breach_floor_min <- 2.0;       // m, breach invert never below this
	float breach_freeboard <- 0.2;       // m, invert sits this above the land side
	float breach_search_radius <- 300.0; // m, search radius for the invert ground level
	float breach_cut_halfwidth <- 60.0;  // m, half-width of the corridor cut through the embankment

	// ------------------------------------------------------------------ engine / thresholds
	float flood_threshold <- 0.05;  // m, depth (above initial) counted as flooded
	float wet_thr   <- 0.1;         // m, building WET threshold (matches the ABM)
	float flood_thr <- 1.0;         // m, building FLOODED threshold (matches the ABM)
	int   max_spill_sweeps <- 120;  // fast-sweeping rounds for the spill + front-distance fields
	float lake_initial_depth <- 0.5;// m, initial standing water in 1925 lakes (as the ABM seeds)
	float datum_offset <- 0.0;      // m, gauge datum minus DEM datum
	bool  auto_pause <- true;

	// ------------------------------------------------------------------ spread front (the datum-robust knob)
	// The flood front advances from the breaches at a finite celerity, so at any
	// moment only cells within geodesic distance R(t) = front_celerity * (time
	// since the first breach) are wetted - even if the level field (S <= stage)
	// would allow more. This makes the extent SPREAD-limited (like the ABM) and
	// datum-robust. front_limit = false disables the front and recovers the pure
	// level-pool. front_celerity is the one knob to calibrate to the ABM peak area.
	bool  front_limit <- true;
	float front_celerity <- 0.02;   // m/s, flood-front spreading speed (CALIBRATE to the ABM)

	// ------------------------------------------------------------------ bookkeeping
	int grid_cols; int grid_rows;
	float cell_dx <- 50.0;
	float z_min <- 0.0; float z_max <- 1.0;
	float SPILL_BIG <- 1e9;                 // unreachable spill level
	list<cell> river_cells <- [];
	list<cell> lake_cells  <- [];
	list<cell> floodable   <- [];           // cells a breach can reach at SOME level (finite S)
	list<cell> metric_cells<- [];           // floodable + lakes (the only cells that can be wet)
	bool spill_dirty <- false;              // a breach changed the terrain -> recompute S
	int  breaches_open <- 0;
	float flooded_area_km2 <- 0.0;
	float flood_volume_mm3 <- 0.0;
	float peak_flooded_km2 <- 0.0;
	int   n_bldg_wet <- 0;
	int   n_bldg_flooded <- 0;
	date  first_break_time <- nil;
	float front_reach_m <- 0.0;             // current front radius R(t) (m), for the monitor
	bool  sim_finished <- false;

	// ------------------------------------------------------------------ geotiff export
	// One georeferenced GeoTIFF per hour holding the water DEPTH (m) of every cell.
	bool   export_geotiff <- true;
	string export_dir <- "../exported_results/hanoi_fastbathtub_1926_discharge/";

	init {
		write "=== Hanoi FastBathtub Flood 1926 (breach engine, DISCHARGE-driven stage) ===";
		grid_cols <- 1 + max(cell collect each.grid_x);
		grid_rows <- 1 + max(cell collect each.grid_y);
		cell_dx <- first(cell).shape.width;

		// --- terrain + neighbour references + inline pit fill -----------------
		ask cell {
			z <- grid_value;
			nE <- grid_x < grid_cols - 1 ? cell[grid_x + 1, grid_y] : nil;
			nW <- grid_x > 0             ? cell[grid_x - 1, grid_y] : nil;
			nS <- grid_y < grid_rows - 1 ? cell[grid_x, grid_y + 1] : nil;
			nN <- grid_y > 0             ? cell[grid_x, grid_y - 1] : nil;
		}
		ask cell {
			float ze <- nE = nil ? z : nE.z;
			float zw <- nW = nil ? z : nW.z;
			float zn <- nN = nil ? z : nN.z;
			float zs <- nS = nil ? z : nS.z;
			float pit <- min(max(0.0, ze - z), min(max(0.0, zw - z), min(max(0.0, zn - z), max(0.0, zs - z))));
			z_dyn <- z + pit;          // pit-filled terrain used by the level-pool
		}
		z_min <- cell min_of each.z;
		z_max <- cell max_of each.z;
		float z_mean <- cell mean_of each.z;
		write "Grid: " + grid_cols + " x " + grid_rows + " cells of " + (cell_dx with_precision 2)
			+ " m, z " + (z_min with_precision 1) + ".." + (z_max with_precision 1) + " m (mean " + (z_mean with_precision 2) + ")";

		// --- vector layers ----------------------------------------------------
		create river_poly from: river_file;
		create lake from: lakes_file;
		create building from: buildings_file;
		create dyke from: dykes_file with: [
			break_s::string(read("BREAK")),
			date_s::string(read("DATE")),
			commune::string(read("Commune"))
		];
		create observation_point from: points_file with: [pid::int(read("id"))];

		// dykes are HARD barriers (water enters only through a breach corridor),
		// exactly like the ABM strict-dyke rule
		ask dyke {
			will_break <- break_s = "YES";
			if will_break and length(date_s) >= 5 {
				int dd <- int(copy_between(date_s, 0, 2));
				int mm <- int(copy_between(date_s, 3, 5));
				breach_date <- date([1926, mm, dd, breach_hour, 0, 0]);
			} else {
				will_break <- false;
			}
			my_cells <- cell overlapping self;
			ask my_cells { is_dyke <- true; }
		}
		first_break_time <- (dyke count each.will_break) > 0 ? (dyke where each.will_break) min_of each.breach_date : nil;

		ask river_poly { ask cell overlapping self { is_river <- true; } }
		ask cell where each.is_dyke { is_river <- false; }   // a crest cell is a barrier, not river
		river_cells <- cell where each.is_river;

		// 1925 lakes: initial standing water (perched ponds, like the ABM seeds)
		ask lake { ask cell overlapping self where (!each.is_river and !each.is_dyke) { is_lake <- true; } }
		lake_cells <- cell where each.is_lake;
		ask lake_cells { h <- lake_initial_depth; wsl <- z_dyn + h; }
		metric_cells <- list(lake_cells);   // lakes are wet before any breach; grows when the spill field is built

		ask building { my_cell <- first(cell overlapping location); }
		ask observation_point { my_cell <- first(cell overlapping self); }

		// baseline (flood metrics are measured ABOVE this, as the ABM does)
		ask cell { h0 <- h; }

		// --- observed discharge record (the SOLE driver) ---------------------
		matrix qm <- matrix(discharge_file);
		loop r over: rows_list(qm) {
			string ds <- string(r[0]);
			float qv <- float(r[1]);
			if length(ds) > 0 and qv > 0.0 {
				list<string> parts <- ds split_with " ";
				list<string> dmy <- first(parts) split_with "/";
				q_dates  <+ date([int(dmy[2]), int(dmy[0]), int(dmy[1]), 0, 0, 0]);
				q_values <+ qv;
			}
		}
		if empty(q_values) {
			use_discharge_csv <- false;
			write "WaterDischarge.csv empty -> stage holds at base_stage (no synthetic fallback)";
		} else {
			q_min <- min(q_values); q_max <- max(q_values);
			write "Discharge (driver): " + length(q_values) + " values, " + first(q_dates) + ".." + last(q_dates)
				+ ", " + q_min + "-" + q_max + " m3/s";
		}

		// --- spill + front fields: the river bank seeds the river side from t=0 (so
		//     the river fills up to the dyke as the stage rises); a breach re-runs
		//     compute_fields to add the protected side through the corridor ----------
		do compute_fields;

		// --- colours ----------------------------------------------------------
		ask cell {
			float shade <- (z - z_min) / max(0.001, z_max - z_min);
			terrain_color <- is_dyke ? rgb(110, 80, 60)
				: rgb(70 + int(150 * shade), 80 + int(130 * shade), 60 + int(110 * shade));
			color <- terrain_color;
		}

		river_stage <- stage_at(starting_date);
		ask river_cells { h <- max(0.0, river_stage + datum_offset - z_dyn); }
		do refresh_colors;

		write "river cells: " + length(river_cells) + " | dyke cells: " + (cell count each.is_dyke)
			+ " | lake cells: " + length(lake_cells);
		write "river-side fills to the dyke from t=0; protected plain dry until the first breach (" + first_break_time + ").";
		write "Init done. Simulation: " + starting_date + " -> " + end_date;
	}

	// ====================================================================== filename helpers
	string pad2 (int v) { return (v < 10 ? "0" : "") + v; }
	string pad4 (int v) { string s <- "" + v; loop while: (length(s) < 4) { s <- "0" + s; } return s; }

	// ====================================================================== forcing
	float discharge_at (date d) {
		if empty(q_values) { return 0.0; }
		if d <= first(q_dates) { return first(q_values); }
		if d >= last(q_dates)  { return last(q_values); }
		loop i from: 1 to: length(q_dates) - 1 {
			if d <= q_dates[i] {
				float f <- (d - q_dates[i - 1]) / max(1.0, q_dates[i] - q_dates[i - 1]);
				return q_values[i - 1] + f * (q_values[i] - q_values[i - 1]);
			}
		}
		return last(q_values);
	}

	// discharge -> power-law rating curve -> stage: the SOLE forcing here. Anchors
	// q_min -> base_stage and q_max -> peak_stage set the curve; there is NO
	// synthetic gaussian stage hydrograph. If the discharge record is missing, the
	// stage simply holds at base_stage (nothing is fabricated).
	float stage_at (date d) {
		if !use_discharge_csv or empty(q_values) { return base_stage; }
		float q <- discharge_at(d);
		float fq <- (q ^ rating_exponent - q_min ^ rating_exponent)
		          / max(1e-6, q_max ^ rating_exponent - q_min ^ rating_exponent);
		return base_stage + (peak_stage - base_stage) * min(1.0, max(0.0, fq));
	}

	reflex update_stage {
		river_discharge <- discharge_at(current_date);
		river_stage <- stage_at(current_date);
	}

	// ====================================================================== breaching
	reflex open_breaches {
		ask dyke where (each.will_break and !each.opened and current_date >= each.breach_date) {
			do open_breach;
		}
		breaches_open <- dyke count (each.opened);
		if spill_dirty {
			do compute_fields;                // terrain changed -> recompute connectivity + front distance
			spill_dirty <- false;
		}
	}

	// ====================================================================== spill + front fields
	// Two static fields of the terrain + open breaches, both solved by the same
	// fast sweeping (4 directional orders), recomputed only when a breach opens.
	// Non-breached dyke cells and river cells are barriers, so water can ONLY
	// enter through a breach corridor - the ABM strict-dyke rule.
	action compute_fields {
		float inv_cel <- 1.0 / max(1e-6, front_celerity);   // s per metre of front travel
		ask cell {
			passable <- (!is_river and !is_dyke) or is_corridor;
			spill_lvl <- SPILL_BIG;
			front_dist <- SPILL_BIG;
			front_arrival <- SPILL_BIG;
		}
		// breach entry: each corridor cell seeds the arrival clock at ITS OWN breach
		// open time (breach_open_s), so a breach that opens a day later does not
		// inherit any spreading credit from an earlier breach (no burst).
		ask cell where each.is_corridor {
			spill_lvl <- z_dyn;
			front_dist <- 0.0;
			front_arrival <- breach_open_s;
		}

		// river-bank entry: every passable LAND cell touching the river floods from
		// the RIVER SIDE as the stage rises (front clock = 0 at the bank), BEFORE any
		// breach. Non-breached dyke cells are barriers, so this fills only up to the
		// dyke on the river side; the protected plain stays dry until a corridor opens.
		// (skip corridor cells: those are the breach gate, seeded at breach_open_s
		//  above - the river side must NOT stamp them with t=0, or that credit leaks
		//  through the breach into the protected plain as a burst.)
		ask cell where each.is_river {
			loop nb over: [nE, nW, nN, nS] {
				if nb != nil and nb.passable and !nb.is_corridor {
					nb.spill_lvl <- min(nb.spill_lvl, nb.z_dyn);
					nb.front_dist <- 0.0;
					nb.front_arrival <- 0.0;
				}
			}
		}

		list<list<cell>> sweep_orders <- [
			cell sort_by (float(each.grid_y * grid_cols + each.grid_x)),
			cell sort_by (float(each.grid_y * grid_cols - each.grid_x)),
			cell sort_by (float(-(each.grid_y * grid_cols) + each.grid_x)),
			cell sort_by (float(-(each.grid_y * grid_cols) - each.grid_x))
		];
		bool changed <- true;
		int rounds <- 0;
		loop while: (changed and rounds < max_spill_sweeps) {
			changed <- false;
			rounds <- rounds + 1;
			loop ord over: sweep_orders {
				loop ce over: ord {
					if ce.passable {
						float bestS <- ce.spill_lvl;
						float bestD <- ce.front_dist;
						float bestA <- ce.front_arrival;
						loop nb over: [ce.nE, ce.nW, ce.nN, ce.nS] {
							if nb != nil and nb.passable {
								// spill: the path to ce through nb must clear max(sill@nb, ground@ce)
								float candS <- max(nb.spill_lvl, ce.z_dyn);
								if candS < bestS { bestS <- candS; }
								// front distance: geodesic, one cell step (diagnostic)
								float candD <- nb.front_dist + cell_dx;
								if candD < bestD { bestD <- candD; }
								// front arrival TIME: neighbour's arrival + travel time of one cell.
								// A corridor cell is a TIME GATE: water cannot cross the breach
								// before it opens, so floor any candidate at this cell's own
								// breach_open_s. (Without this the river-side seed's t=0 credit
								// would leak through the corridor the instant the breach opens -> burst.)
								float candA <- nb.front_arrival + cell_dx * inv_cel;
								if ce.is_corridor and candA < ce.breach_open_s { candA <- ce.breach_open_s; }
								if candA < bestA { bestA <- candA; }
							}
						}
						if bestS < ce.spill_lvl - 1e-6   { ce.spill_lvl <- bestS;    changed <- true; }
						if bestD < ce.front_dist - 1e-3  { ce.front_dist <- bestD;   changed <- true; }
						if bestA < ce.front_arrival - 1.0 { ce.front_arrival <- bestA; changed <- true; }
					}
				}
			}
		}
		floodable    <- cell where (each.passable and each.spill_lvl < 0.5 * SPILL_BIG);
		metric_cells <- remove_duplicates(floodable + lake_cells);
		write "" + current_date + "  spill + front fields rebuilt in " + rounds + " sweep rounds; "
			+ length(floodable) + " cells reachable by a breach (S "
			+ ((empty(floodable) ? 0.0 : floodable min_of each.spill_lvl) with_precision 2) + ".."
			+ ((empty(floodable) ? 0.0 : floodable max_of each.spill_lvl) with_precision 2) + " m, D up to "
			+ ((empty(floodable) ? 0.0 : floodable max_of each.front_dist) with_precision 0) + " m).";
	}

	// ====================================================================== dynamic flood (the engine)
	// Spread-limited, connectivity-constrained, hysteretic level-pool. One
	// filtered parallel ask per hour, no iteration. A cell is wetted when BOTH the
	// flood front has arrived (front_arrival[c] <= elapsed) AND it is
	// river-connected at the stage (S[c] <= stage); then it tracks the stage up and
	// down. Cells that fall out of (S <= L) keep their water -> trapped casiers.
	reflex dynamic_flood {
		float L <- river_stage + datum_offset;
		ask river_cells { h <- max(0.0, L - z_dyn); }       // river held at the stage (source/boundary)
		if !empty(floodable) {
			float elapsed_s <- front_limit ? (current_date - starting_date) : SPILL_BIG;
			front_reach_m <- front_limit ? front_celerity * (current_date - starting_date) : SPILL_BIG;
			ask (floodable where (each.spill_lvl <= L and each.front_arrival <= elapsed_s)) parallel: true {
				h <- L - z_dyn;                              // >= 0 since z_dyn <= S <= L
				wsl <- z_dyn + h;
			}
			// cells the front has not reached (D > R): stay dry until it arrives;
			// perched cells (S > L): retain h from the last connected hour = trapped casier
		}
	}

	// ====================================================================== bookkeeping
	reflex bookkeeping {
		ask metric_cells parallel: true {
			float exc <- h - h0;                             // depth above the initial baseline
			if exc > flood_threshold {
				if arrival_h < 0.0 { arrival_h <- (current_date - starting_date) / 3600.0; }
				if h > h_peak { h_peak <- h; }
			}
		}
		list<cell> wet <- metric_cells where ((each.h - each.h0) > flood_threshold and !each.is_river);
		flooded_area_km2 <- length(wet) * cell_dx * cell_dx / 1e6;
		flood_volume_mm3 <- (wet sum_of (each.h - each.h0)) * cell_dx * cell_dx / 1e6;
		peak_flooded_km2 <- max(peak_flooded_km2, flooded_area_km2);

		ask building {
			depth_w <- my_cell = nil ? 0.0 : max(0.0, my_cell.h - my_cell.h0);
			status <- depth_w > flood_thr ? 2 : (depth_w > wet_thr ? 1 : 0);
		}
		n_bldg_wet     <- building count (each.status = 1);
		n_bldg_flooded <- building count (each.status = 2);

		ask observation_point where (each.arrival_h < 0.0) {
			if my_cell != nil and (my_cell.h - my_cell.h0) > flood_threshold {
				arrival_h <- (current_date - starting_date) / 3600.0;
				write "Observation point " + pid + " reached on " + current_date
					+ " (depth " + ((my_cell.h - my_cell.h0) with_precision 2) + " m)";
			}
		}

		do refresh_colors;
		if current_date.hour mod 6 = 0 {
			write "" + current_date + " | stage " + (river_stage with_precision 2) + " m | breaches " + breaches_open
				+ " | flooded " + (flooded_area_km2 with_precision 2) + " km2 | vol "
				+ (flood_volume_mm3 with_precision 1) + " Mm3 | bldg flooded " + n_bldg_flooded;
		}
	}

	// ====================================================================== geotiff export (every step)
	reflex export_water_height when: export_geotiff {
		// flooded footprint = same rule as flooded_area (excess above baseline, not river)
		list<cell> wet_land <- metric_cells where ((each.h - each.h0) > flood_threshold and !each.is_river);
		float peak_h <- empty(wet_land) ? 0.0 : wet_land max_of each.h;
		// band value = water depth: river channel depth on river cells, flood depth on
		// flooded land, 0 on dry land
		ask cell { grid_value <- (is_river or ((h - h0) > flood_threshold)) ? h : 0.0; }
		string ts <- "" + current_date.year + pad2(current_date.month) + pad2(current_date.day)
			+ "-" + pad2(current_date.hour) + pad2(current_date.minute);
		string fname <- export_dir + "wh_step" + pad4(cycle) + "_" + ts
			+ "_area" + (flooded_area_km2 with_precision 2) + "km2"
			+ "_peak" + (peak_h with_precision 2) + "m.tif";
		save cell to: fname format: "geotiff";
	}

	action refresh_colors {
		ask metric_cells + river_cells {
			if h > 0.02 {
				float f <- min(1.0, h / 4.0);
				color <- rgb(int(150 * (1 - f)), int(190 * (1 - f) + 30), int(180 + 75 * f));
				was_wet_color <- true;
			} else if was_wet_color {
				color <- terrain_color;
				was_wet_color <- false;
			}
		}
	}

	reflex stop_simulation when: current_date >= end_date {
		sim_finished <- true;
		write "End of event. Flooded area: " + (flooded_area_km2 with_precision 2)
			+ " km2 (peak " + (peak_flooded_km2 with_precision 2) + " km2).";
		ask observation_point {
			write "Point " + pid + " arrival: " + (arrival_h < 0.0 ? "never" : string(arrival_h with_precision 1) + " h");
		}
		if auto_pause { do pause; }
	}
}

// ==========================================================================
//  Raster domain (minimal grid agents, no scheduler stepping)
// ==========================================================================
grid cell file: dem_file neighbors: 4
	use_regular_agents: false use_individual_shapes: false use_neighbors_cache: false schedules: [] {
	float z;                 // raw DEM elevation (dyke crests included)
	float z_dyn;             // pit-filled terrain; lowered in a breach corridor
	float h <- 0.0;          // water depth (m)
	float h0 <- 0.0;         // initial depth (lakes) - flood metrics are measured above this
	float wsl <- 0.0;        // water-surface level z_dyn + h
	float spill_lvl <- 1e9;  // S[c]: lowest stage at which the cell connects to a breach
	float front_dist <- 1e9; // D[c]: geodesic distance (m) from the nearest breach corridor
	float front_arrival <- 1e9; // T[c]: sim-seconds when the flood front first reaches the cell
	float breach_open_s <- 1e9; // seconds since sim start when a breach first cut this cell (corridor seed)
	bool  passable <- false; // traversable by the spill / front sweep (not a barrier)
	bool  is_river <- false;
	bool  is_lake  <- false;
	bool  is_dyke  <- false;
	bool  is_corridor <- false; // crest cell opened by a breach (passable entry)
	float h_peak <- 0.0;
	float arrival_h <- -1.0;
	cell nE; cell nW; cell nN; cell nS;
	rgb terrain_color <- #gray;
	bool was_wet_color <- false;
}

// ==========================================================================
//  Vector species (no reflexes; behaviour is driven by the global asks)
// ==========================================================================
species dyke schedules: [] {
	string break_s; string date_s; string commune;
	bool will_break <- false;
	bool opened <- false;
	date breach_date;
	list<cell> my_cells;

	// open the breach: invert from the lowest protected-side ground within
	// breach_search_radius, cut as a corridor through the full embankment width;
	// the cut cells become a passable entry, then the spill-level field is rebuilt
	action open_breach {
		opened <- true;
		list<cell> search_zone <- cell overlapping (shape + breach_search_radius);
		list<cell> ground <- search_zone where (!each.is_dyke and !each.is_river);
		float target <- empty(ground)
			? (my_cells min_of each.z_dyn) - 5.0
			: (ground min_of each.z_dyn) + breach_freeboard;
		target <- max(breach_floor_min, target);
		list<cell> corridor <- (cell overlapping (shape + breach_cut_halfwidth)) where (!each.is_river);
		float open_s <- current_date - starting_date;        // seconds since sim start, this breach
		ask corridor {
			z_dyn <- min(z_dyn, target);
			is_corridor <- true;
			breach_open_s <- min(breach_open_s, open_s);     // earliest breach to cut this cell wins
		}
		spill_dirty <- true;
		write "BREACH at " + commune + " on " + current_date + " (invert lowered to "
			+ (target with_precision 2) + " m, " + length(corridor) + " cells cut)";
	}

	aspect default {
		draw shape color: opened ? #red : (will_break ? #orange : #darkgreen) width: 3;
	}
}

species river_poly schedules: [] { aspect default { draw shape color: rgb(70, 130, 180, 120) border: #steelblue; } }
species lake       schedules: [] { aspect default { draw shape color: rgb(150, 200, 230, 130); } }
species building   schedules: [] {
	cell my_cell;
	float depth_w <- 0.0;
	int status <- 0;          // 0 dry, 1 wet, 2 flooded
	aspect default { draw shape color: status = 2 ? #red : (status = 1 ? #orange : rgb(90, 90, 90)); }
}
species observation_point schedules: [] {
	int pid;
	cell my_cell;
	float arrival_h <- -1.0;
	aspect default {
		draw circle(120) color: arrival_h < 0.0 ? #white : #red border: #black;
		draw string(pid) + (arrival_h < 0.0 ? "" : (" : " + (arrival_h with_precision 1) + " h"))
			at: location + {150, -100} color: #black font: font("SansSerif", 14, #bold);
	}
}

// ==========================================================================
//  Experiments
// ==========================================================================
experiment fastbathtub_1926_discharge type: gui {
	parameter "DEM file" var: dem_name among: ["DEM_Hanoi_Asc.asc", "mnt-gz10.asc", "mnt-gz40.asc", "mnt-gz25.asc", "mnt-gz50.asc", "mnt-gz50-1926.asc"] category: "Terrain";
	parameter "Use observed discharge (WaterDischarge.csv)" var: use_discharge_csv category: "Forcing";
	parameter "Rating curve exponent" var: rating_exponent category: "Forcing";
	parameter "Base river stage (m)" var: base_stage category: "Forcing";
	parameter "Peak river stage (m)" var: peak_stage category: "Forcing";
	parameter "Gauge datum offset (m)" var: datum_offset min: -3.0 max: 3.0 category: "Forcing";
	parameter "Breach hour of day" var: breach_hour category: "Breaching";
	parameter "Breach invert search radius (m)" var: breach_search_radius category: "Breaching";
	parameter "Breach cut half-width (m)" var: breach_cut_halfwidth category: "Breaching";
	parameter "Initial lake depth (m)" var: lake_initial_depth min: 0.0 max: 2.0 category: "Initial state";
	parameter "Flood threshold (m)" var: flood_threshold min: 0.01 max: 0.5 category: "Engine";
	parameter "Limit spread by front (ABM-like / datum-robust)" var: front_limit category: "Spread front";
	parameter "Front celerity (m/s) - CALIBRATE to ABM" var: front_celerity min: 0.005 max: 2.0 category: "Spread front";
	parameter "Export water-height GeoTIFF each step" var: export_geotiff category: "Export";
	parameter "Auto pause at end" var: auto_pause category: "Engine";

	output {
		layout #split;
		display "Flood simulation" type: 2d background: #black {
			grid cell;
			graphics "static landscape" refresh: false {
				loop rp over: river_poly { draw rp.shape color: rgb(70, 130, 180, 120) border: #steelblue; }
				loop lk over: lake { draw lk.shape color: rgb(150, 200, 230, 90); }
			}
			species building;
			species dyke;
			species observation_point;
			graphics "info" {
				draw string(current_date) + "   stage: " + (river_stage with_precision 2)
					+ " m   flooded: " + (flooded_area_km2 with_precision 1) + " km2"
					at: {world.shape.width * 0.02, world.shape.height * 0.03}
					color: #white font: font("SansSerif", 16, #bold);
			}
		}
		display "Time series" type: 2d {
			chart "1926 flood event" type: series x_label: "hours since 20-07 00:00" {
				data "River stage (m)" value: river_stage color: #blue marker: false;
				data "Discharge (1000 m3/s)" value: river_discharge / 1000.0 color: #darkblue marker: false;
				data "Flooded area (km2)" value: flooded_area_km2 color: #red marker: false;
				data "Flood volume (10^6 m3)" value: flood_volume_mm3 color: #darkorange marker: false;
			}
		}
		monitor "Date" value: current_date;
		monitor "River discharge (m3/s)" value: river_discharge with_precision 1;
		monitor "River stage (m)" value: river_stage with_precision 2;
		monitor "Breaches open" value: breaches_open;
		monitor "Reachable cells" value: length(floodable);
		monitor "Front reach (m)" value: front_limit ? int(front_reach_m) : -1;
		monitor "Flooded area (km2)" value: flooded_area_km2 with_precision 2;
		monitor "Peak flooded (km2)" value: peak_flooded_km2 with_precision 2;
		monitor "Flood volume (10^6 m3)" value: flood_volume_mm3 with_precision 2;
		monitor "Buildings wet / flooded" value: "" + n_bldg_wet + " / " + n_bldg_flooded;
	}
}

// Display-free production / calibration experiment (numbers only, fastest).
experiment fastbathtub_1926_discharge_fast type: gui {
	parameter "DEM file" var: dem_name among: ["mnt-gz40.asc", "mnt-gz25.asc", "mnt-gz50.asc", "mnt-gz50-1926.asc"] category: "Terrain";
	parameter "Peak river stage (m)" var: peak_stage category: "Forcing";
	parameter "Breach cut half-width (m)" var: breach_cut_halfwidth category: "Breaching";
	parameter "Limit spread by front (ABM-like / datum-robust)" var: front_limit category: "Spread front";
	parameter "Front celerity (m/s) - CALIBRATE to ABM" var: front_celerity min: 0.005 max: 2.0 category: "Spread front";
	parameter "Export water-height GeoTIFF each step" var: export_geotiff category: "Export";

	output {
		display "Time series" type: 2d {
			chart "1926 flood event" type: series x_label: "hours since 20-07 00:00" {
				data "River stage (m)" value: river_stage color: #blue marker: false;
				data "Flooded area (km2)" value: flooded_area_km2 color: #red marker: false;
				data "Flood volume (10^6 m3)" value: flood_volume_mm3 color: #darkorange marker: false;
			}
		}
		monitor "Date" value: current_date;
		monitor "River discharge (m3/s)" value: river_discharge with_precision 1;
		monitor "River stage (m)" value: river_stage with_precision 2;
		monitor "Breaches open" value: breaches_open;
		monitor "Front reach (m)" value: front_limit ? int(front_reach_m) : -1;
		monitor "Flooded area (km2)" value: flooded_area_km2 with_precision 2;
		monitor "Peak flooded (km2)" value: peak_flooded_km2 with_precision 2;
		monitor "Flood volume (10^6 m3)" value: flood_volume_mm3 with_precision 2;
	}
}
