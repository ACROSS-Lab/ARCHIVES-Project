/**
* Name: Phu Tho FastBathtub Flood Model
* Author: Thành Đô Nguyễn (2026-06-17)
*
* A fork of "Dong Hoi FastBathtub Flood Model.gaml" for the Thao/Red-River
* floodplain at Phú Thọ (Typhoon Yagi, Sep-2024 event).
*
* WHAT CHANGED vs the Đồng Hới model, and WHY.
*   Same fast engine - a spread-limited, connectivity-constrained, hysteretic
*   level-pool solved by one fast-sweeping pass, SEEDED FROM THE RIVER BANK
*   (no dyke ring; the plain floods by the river spilling over its banks). Only
*   the site, the datum and the forcing change:
*     - DEM is the ~30 m SRTM tile dem-phutho.tif. Unlike the Đồng Hới model
*       (which used a hand-reprojected EPSG:3857 .asc), EVERY Phú Thọ layer (this
*       DEM and both shapefiles) is in EPSG:4326 (degrees), so they are mutually
*       consistent and GAMA's default metric auto-projection reprojects them
*       together -> cell width comes out ~30 m (check the init log: if it prints a
*       tiny cell_dx ~0.0003, GAMA did NOT project and the DEM must be reprojected
*       to a metric CRS first).
*     - Forcing is STAGE-DRIVEN by default (use_discharge_csv = false). The
*       measured forcing we trust is the Phú Thọ stage record PhuThoStage2024_hourly.csv
*       (daily levels Sep 7-22, two verified points: crest 18.34 m on Sep 11,
*       17.42 m on Sep 12, hourly-interpolated). The discharge file
*       WaterDischarge_PhuTho.csv is a TIME-SYNCED PLACEHOLDER (a synthetic
*       0->3000 ramp re-dated onto the event window); it is loaded for wiring but
*       must NOT drive the model until a real Sep-2024 hydrograph replaces it.
*
* INPUTS (../includes/ and ../includes/phu-tho/).
*   DEM: dem-phutho.tif - ~30 m SRTM (EPSG:4326, z -10..363 m, mean 31.5 m),
*        nodata -32767. GAMA auto-projects it to metres on load.
*   River: water_phutho.shp - the Thao/Red-River water polygon (stage boundary).
*   Buildings: building_phutho.shp (clipped to the domain in init).
*   Forcing: PhuThoStage2024_hourly.csv (measured stage, DRIVING) with
*        WaterDischarge_PhuTho.csv as a (placeholder) discharge alternative.
*        DATUM: the Phú Thọ gauge and the SRTM/EGM96 DEM are nearly co-datum here
*        (inland reach: stage 13.3..18.34 m ~ absolute elevation), so datum_offset
*        starts at 0. Calibrate it so the flood extent matches the observed event.
*/
model PhuThoFastBathtubFlood

global {

	// ------------------------------------------------------------------ input
	// DEM in EPSG:4326 (degrees); GAMA reprojects to a metric CRS on load.
	string dem_name     <- "dem-phutho-epsg-3857.tif";
	file dem_file       <- grid_file("../includes/phu-tho/" + dem_name);
	file river_file     <- shape_file("../includes/phu-tho/water_phutho.shp");
	file buildings_file <- shape_file("../includes/phu-tho/building_phutho.shp");
	file stage_file     <- csv_file("../includes/PhuThoStage2024_hourly.csv", ",", true);
	file discharge_file <- csv_file("../includes/WaterDischarge_PhuTho_hourly.csv", ",", true);
	geometry shape <- envelope(dem_file);

	// ------------------------------------------------------------------ time
	// the stage record is hourly, 7-22 Sep 2024 (360 h)
	date starting_date <- date("2024-09-07 00:00:00");
	date end_date      <- date("2024-09-22 00:00:00");
	float step <- 1 #h;                       // ONE CYCLE = ONE HOUR

	// ------------------------------------------------------------------ river stage forcing
	// STAGE-DRIVEN by default: the measured Phú Thọ gauge (PhuThoStage2024_hourly.csv)
	// drives the level directly. use_discharge_csv = true would instead push the
	// (PLACEHOLDER) discharge WaterDischarge_PhuTho.csv through the Manning rating
	// curve - do NOT enable until a real Sep-2024 hydrograph is supplied.
	bool  use_discharge_csv <- true;
	float rating_exponent <- 0.6;             // Manning h ~ Q^(3/5)
	float base_stage <- 13.3;                 // m, stage at the lowest recorded discharge (gauge min)
	float peak_stage <- 25.00;                // m, stage at the peak discharge (verified crest, Sep 11)
	float stage_cap  <- 30.0;                 // m, HARD CEILING on river stage (uncapped rating curve can overshoot)
	float river_baseflow <- 1.5;              // m, minimum water depth kept in the channel (bed + baseflow),
	                                          // so the river stays CONTINUOUSLY filled along its sloping bed
	                                          // even upstream where the flat stage L is below the local bed.
	float river_stage <- 13.3;
	float river_discharge <- 0.0;
	list<date>  q_dates  <- [];
	list<float> q_values <- [];
	float q_min <- 1.0; float q_max <- 2.0;
	// measured-gauge series (the DEFAULT driver here, since use_discharge_csv = false)
	list<float> stage_series <- [];
	int   n_stage <- 0;

	// datum_offset = (gauge zero) - (DEM/SRTM zero), in metres. Water surface used
	// by the engine is  L = river_stage + datum_offset.  This is THE calibration knob.
	// STARTING ANCHOR (verified by sampling the DEM under water_phutho.shp):
	//   This is a MOUNTAINOUS reach - the river channel itself sits HIGH in the SRTM
	//   (z 13..41 m, median ~20 m), so the gauge "zero" is at a high absolute
	//   elevation. With offset 0 the gauge band (13.3..18.34 m) falls BELOW the river
	//   bed -> almost nothing floods. A POSITIVE offset lifts the whole band onto the
	//   ~20 m floodplain WITHOUT changing the ~5 m flood amplitude:
	//     base L ~ 15-16 m (channel wet, banks dry) and crest L ~ 21-22 m (spills onto
	//     the floodplain) are both met near offset +3..+4 m.
	//   This is the DATUM (position) knob - to change HOW MUCH the river rises, edit
	//   peak_stage instead; to lift WHERE it sits, edit this.
	// CALIBRATED TO EXTENT via the DEM hypsometry (area below level L):
	//   L=22 -> ~113 km2, L=23 -> ~128 km2, L=24 -> ~141 km2 (UNCONSTRAINED bathtub;
	//   connectivity trims the model below this). Target HEC-RAS ~128 km2 -> crest
	//   L ~24-25 m. crest L = peak_stage(18.34) + datum_offset, so offset 6 -> L 24.3.
	//   Sweep +5..+7 watching the "Flooded area" monitor to land on 128 km2.
	float datum_offset <- 0.0;

	// ------------------------------------------------------------------ engine / thresholds
	float flood_threshold <- 0.05;  // m, depth above baseline counted as flooded
	float wet_thr   <- 0.1;         // m, building WET threshold
	float flood_thr <- 1.0;         // m, building FLOODED threshold
	int   max_spill_sweeps <- 200;  // fast-sweeping rounds for the spill + front fields
	bool  auto_pause <- true;

	// ------------------------------------------------------------------ spread front (datum-robust knob)
	// The flood front advances from the river bank at a finite celerity, so at any
	// moment only cells within geodesic travel-time of the bank are wet, even if
	// the level field (S <= stage) would allow more. This makes the extent
	// spread-limited and datum-robust. front_limit = false recovers the pure
	// level-pool. NOTE: over this 360 h event 0.05 m/s reaches ~65 km - i.e. the
	// front barely limits a ~20 km domain; lower it to make spreading visible.
	bool  front_limit <- true;
	float front_celerity <- 0.05;   // m/s, flood-front spreading speed (CALIBRATE)

	// ------------------------------------------------------------------ bookkeeping
	int grid_cols; int grid_rows;
	float cell_dx <- 30.0;
	float z_min <- 0.0; float z_max <- 1.0;
	float SPILL_BIG <- 1e9;
	list<cell> river_cells <- [];
	list<cell> floodable   <- [];           // cells the river can reach at SOME stage (finite S)
	list<cell> metric_cells<- [];
	float flooded_area_km2 <- 0.0;
	float flood_volume_mm3 <- 0.0;
	float peak_flooded_km2 <- 0.0;
	int   n_bldg_wet <- 0;
	int   n_bldg_flooded <- 0;
	float front_reach_m <- 0.0;
	bool  sim_finished <- false;

	init {
		write "=== Phu Tho FastBathtub Flood (river-seeded, dyke-free level-pool) ===";
		grid_cols <- 1 + max(cell collect each.grid_x);
		grid_rows <- 1 + max(cell collect each.grid_y);
		cell_dx <- first(cell).shape.width;
		if cell_dx < 1.0 {
			write "WARNING: cell_dx = " + cell_dx + " (looks like DEGREES, not metres). GAMA did NOT "
				+ "project the DEM to a metric CRS - reproject dem-phutho.tif before trusting any area/celerity.";
		}

		// --- terrain + neighbour references + inline pit fill -----------------
		ask cell {
			z <- grid_value;
			is_nodata <- z < -1000.0;           // -32767 SRTM sentinel
			nE <- grid_x < grid_cols - 1 ? cell[grid_x + 1, grid_y] : nil;
			nW <- grid_x > 0             ? cell[grid_x - 1, grid_y] : nil;
			nS <- grid_y < grid_rows - 1 ? cell[grid_x, grid_y + 1] : nil;
			nN <- grid_y > 0             ? cell[grid_x, grid_y - 1] : nil;
		}
		ask cell where (!each.is_nodata) {
			float ze <- (nE = nil or nE.is_nodata) ? z : nE.z;
			float zw <- (nW = nil or nW.is_nodata) ? z : nW.z;
			float zn <- (nN = nil or nN.is_nodata) ? z : nN.z;
			float zs <- (nS = nil or nS.is_nodata) ? z : nS.z;
			float pit <- min(max(0.0, ze - z), min(max(0.0, zw - z), min(max(0.0, zn - z), max(0.0, zs - z))));
			z_dyn <- z + pit;
		}
		list<cell> valid_cells <- cell where (!each.is_nodata);
		z_min <- valid_cells min_of each.z;
		z_max <- valid_cells max_of each.z;
		float z_mean <- valid_cells mean_of each.z;
		write "Grid: " + grid_cols + " x " + grid_rows + " cells of " + (cell_dx with_precision 2)
			+ " m, z " + (z_min with_precision 1) + ".." + (z_max with_precision 1) + " m (mean " + (z_mean with_precision 2) + ")";

		// --- vector layers ----------------------------------------------------
		create river_poly from: river_file;
		create building from: buildings_file;

		// river/water cells = stage boundary (held at the stage every hour)
		ask river_poly { ask cell overlapping self where (!each.is_nodata) { is_river <- true; } }
		river_cells <- cell where each.is_river;

		// buildings: bind to a cell, drop those outside the DEM domain
		ask building { my_cell <- first(cell overlapping location); }
		ask building where (each.my_cell = nil) { do die; }

		// baseline (flood metrics measured ABOVE this; no lakes here so h0 = 0)
		ask cell { h0 <- h; }

		// --- measured stage record (the DEFAULT driver) ----------------------
		matrix sm <- matrix(stage_file);
		loop r over: rows_list(sm) {
			string ds <- string(r[0]);
			if length(ds) > 0 and ds != "datetime" {
				stage_series <+ float(r[1]);
			}
		}
		n_stage <- length(stage_series);
		if n_stage = 0 {
			write "PhuThoStage2024_hourly.csv empty -> measured driver unavailable";
		} else {
			write "Stage record (driver): " + n_stage + " h, " + (min(stage_series) with_precision 2)
				+ ".." + (max(stage_series) with_precision 2) + " m (gauge datum)";
		}

		// --- discharge record (placeholder; dates are M/D/YYYY H:MM) ----------
		matrix qm <- matrix(discharge_file);
		loop r over: rows_list(qm) {
			string ds <- string(r[0]);
			float qv <- float(r[1]);
			if length(ds) > 0 and qv > 0.0 {
				list<string> parts <- ds split_with " ";
				list<string> dmy <- first(parts) split_with "/";
				int hh <- 0;
				if length(parts) > 1 { hh <- int(first(parts[1] split_with ":")); }
				q_dates  <+ date([int(dmy[2]), int(dmy[0]), int(dmy[1]), hh, 0, 0]);
				q_values <+ qv;
			}
		}
		if empty(q_values) {
			use_discharge_csv <- false;
			write "WaterDischarge_PhuTho.csv empty -> using measured stage gauge";
		} else {
			q_min <- min(q_values); q_max <- max(q_values);
			write "Discharge (PLACEHOLDER): " + length(q_values) + " values, " + first(q_dates) + ".." + last(q_dates)
				+ ", " + (q_min with_precision 1) + "-" + (q_max with_precision 1) + " m3/s"
				+ (use_discharge_csv ? "  (DRIVING via rating curve - WARNING: placeholder data)" : "  (loaded; measured gauge is driving)");
		}

		// --- spill + front fields (static: river bank is the seed) -----------
		do compute_fields;

		// --- colours ----------------------------------------------------------
		ask cell where (!each.is_nodata) {
			float shade <- (z - z_min) / max(0.001, z_max - z_min);
			terrain_color <- is_river ? rgb(70, 130, 180)
				: rgb(70 + int(150 * shade), 80 + int(130 * shade), 60 + int(110 * shade));
			color <- terrain_color;
		}

		river_stage <- stage_at(starting_date);
		// channel = local bed + baseflow, OR the flat stage where that is higher
		ask river_cells { h <- max(river_baseflow, river_stage + datum_offset - z_dyn); }
		do refresh_colors;

		write "river cells: " + length(river_cells) + " | floodable cells: " + length(floodable)
			+ " | buildings in domain: " + length(building);
		write "Init done. Simulation: " + starting_date + " -> " + end_date;
	}

	// ====================================================================== forcing
	// placeholder discharge interpolated in time (Hanoi-style)
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

	// measured-gauge stage (the DEFAULT driver, used when use_discharge_csv = false)
	float measured_stage_at (date d) {
		if n_stage = 0 { return base_stage; }
		int hr <- int((d - starting_date) / 3600.0);
		if hr < 0        { return first(stage_series); }
		if hr >= n_stage { return last(stage_series); }
		return stage_series[hr];
	}

	// measured gauge; or discharge -> power-law rating curve -> stage.
	// The two anchors (q_min -> base_stage, q_max -> peak_stage) only set the
	// curve's SLOPE - there is NO upper cap, so above the reference discharge the
	// stage keeps rising (discharge fully drives the height). Only the lower bound
	// is kept (stage never drops below base_stage for Q < q_min).
	float stage_at (date d) {
		float s;
		if use_discharge_csv and !empty(q_values) {
			float q <- discharge_at(d);
			float fq <- (q ^ rating_exponent - q_min ^ rating_exponent)
			          / max(1e-6, q_max ^ rating_exponent - q_min ^ rating_exponent);
			s <- base_stage + (peak_stage - base_stage) * max(0.0, fq);
		} else {
			s <- measured_stage_at(d);
		}
		return min(stage_cap, s);   // hard 30 m ceiling
	}

	reflex update_stage {
		river_discharge <- discharge_at(current_date);
		river_stage <- stage_at(current_date);
	}

	// ====================================================================== spill + front fields
	// Two static fields, both solved by the same fast sweeping (4 directional
	// orders). The river is the boundary: every passable LAND cell adjacent to a
	// river cell is a seed (it connects to the river once the stage clears its own
	// ground), with the front clock starting at t=0. River and no-data cells are
	// barriers (the front travels through land, not through the channel).
	//   S[c] (spill_lvl) = lowest stage at which c connects to the river.
	//   T[c] (front_arrival) = sim-seconds when the front first reaches c.
	action compute_fields {
		float inv_cel <- 1.0 / max(1e-6, front_celerity);   // s per metre of front travel
		ask cell {
			passable <- !is_river and !is_nodata;
			spill_lvl <- SPILL_BIG;
			front_dist <- SPILL_BIG;
			front_arrival <- SPILL_BIG;
		}
		// seed: land cells touching the river enter at their own ground level,
		// connected from the first hour (front clock = 0 at the bank)
		ask cell where each.is_river {
			loop nb over: [nE, nW, nN, nS] {
				if nb != nil and nb.passable {
					if nb.z_dyn < nb.spill_lvl { nb.spill_lvl <- nb.z_dyn; }
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
								float candS <- max(nb.spill_lvl, ce.z_dyn);
								if candS < bestS { bestS <- candS; }
								float candD <- nb.front_dist + cell_dx;
								if candD < bestD { bestD <- candD; }
								float candA <- nb.front_arrival + cell_dx * inv_cel;
								if candA < bestA { bestA <- candA; }
							}
						}
						if bestS < ce.spill_lvl - 1e-6   { ce.spill_lvl <- bestS;     changed <- true; }
						if bestD < ce.front_dist - 1e-3  { ce.front_dist <- bestD;    changed <- true; }
						if bestA < ce.front_arrival - 1.0 { ce.front_arrival <- bestA; changed <- true; }
					}
				}
			}
		}
		floodable    <- cell where (each.passable and each.spill_lvl < 0.5 * SPILL_BIG);
		metric_cells <- floodable;
		write "spill + front fields built in " + rounds + " sweep rounds; "
			+ length(floodable) + " cells reachable from the river (S "
			+ ((empty(floodable) ? 0.0 : floodable min_of each.spill_lvl) with_precision 2) + ".."
			+ ((empty(floodable) ? 0.0 : floodable max_of each.spill_lvl) with_precision 2) + " m, D up to "
			+ ((empty(floodable) ? 0.0 : floodable max_of each.front_dist) with_precision 0) + " m).";
	}

	// ====================================================================== dynamic flood (the engine)
	// Spread-limited, connectivity-constrained, hysteretic level-pool. One
	// filtered parallel ask per hour, no iteration. A cell is wet when BOTH the
	// front has arrived (T[c] <= elapsed) AND it is river-connected at the stage
	// (S[c] <= L); then it tracks the stage up and down. Cells that fall out of
	// S <= L on the receding limb keep their water -> trapped ponds.
	reflex dynamic_flood {
		float L <- river_stage + datum_offset;
		ask river_cells { h <- max(river_baseflow, L - z_dyn); }   // channel = bed+baseflow, or stage if higher
		if !empty(floodable) {
			float elapsed_s <- front_limit ? (current_date - starting_date) : SPILL_BIG;
			front_reach_m <- front_limit ? front_celerity * (current_date - starting_date) : SPILL_BIG;
			ask (floodable where (each.spill_lvl <= L and each.front_arrival <= elapsed_s)) parallel: true {
				h <- L - z_dyn;                              // >= 0 since z_dyn <= S <= L
				wsl <- z_dyn + h;
			}
		}
	}

	// ====================================================================== bookkeeping
	reflex bookkeeping {
		ask metric_cells parallel: true {
			float exc <- h - h0;
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

		do refresh_colors;
		if current_date.hour mod 6 = 0 {
			write "" + current_date + " | stage " + (river_stage with_precision 2) + " m (L "
				+ ((river_stage + datum_offset) with_precision 2) + ") | flooded "
				+ (flooded_area_km2 with_precision 2) + " km2 | vol "
				+ (flood_volume_mm3 with_precision 1) + " Mm3 | bldg flooded " + n_bldg_flooded;
		}
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
		if auto_pause { do pause; }
	}
}

// ==========================================================================
//  Raster domain (minimal grid agents, no scheduler stepping)
// ==========================================================================
grid cell file: dem_file neighbors: 4
	use_regular_agents: false use_individual_shapes: false use_neighbors_cache: false schedules: [] {
	float z;                 // raw DEM elevation (m, SRTM/EGM96)
	float z_dyn;             // pit-filled terrain used by the level-pool
	float h <- 0.0;          // water depth (m)
	float h0 <- 0.0;         // initial depth (0 here; metrics measured above this)
	float wsl <- 0.0;        // water-surface level z_dyn + h
	float spill_lvl <- 1e9;  // S[c]: lowest stage at which the cell connects to the river
	float front_dist <- 1e9; // D[c]: geodesic distance (m) from the nearest river bank
	float front_arrival <- 1e9; // T[c]: sim-seconds when the front first reaches the cell
	bool  passable <- false;
	bool  is_river  <- false;
	bool  is_nodata <- false;
	float h_peak <- 0.0;
	float arrival_h <- -1.0;
	cell nE; cell nW; cell nN; cell nS;
	rgb terrain_color <- #gray;
	bool was_wet_color <- false;
}

// ==========================================================================
//  Vector species (no reflexes; behaviour is driven by the global asks)
// ==========================================================================
species river_poly schedules: [] { aspect default { draw shape color: rgb(70, 130, 180, 120) border: #steelblue; } }
species building   schedules: [] {
	cell my_cell;
	float depth_w <- 0.0;
	int status <- 0;          // 0 dry, 1 wet, 2 flooded
	aspect default { draw shape color: status = 2 ? #red : (status = 1 ? #orange : rgb(90, 90, 90)); }
}

// ==========================================================================
//  Experiments
// ==========================================================================
experiment phutho_fastbathtub type: gui {
	parameter "Use placeholder discharge (WaterDischarge_PhuTho.csv)" var: use_discharge_csv category: "Forcing";
	parameter "Rating curve exponent" var: rating_exponent category: "Forcing";
	parameter "Base river stage (m)" var: base_stage category: "Forcing";
	parameter "Peak river stage (m)" var: peak_stage category: "Forcing";
	parameter "Stage cap (m) - hard ceiling" var: stage_cap min: 5.0 max: 60.0 category: "Forcing";
	parameter "Channel baseflow depth (m)" var: river_baseflow min: 0.0 max: 5.0 category: "Forcing";
	parameter "Gauge datum offset (m) - CALIBRATE" var: datum_offset min: -15.0 max: 15.0 category: "Forcing";
	parameter "Flood threshold (m)" var: flood_threshold min: 0.01 max: 0.5 category: "Engine";
	parameter "Limit spread by front (datum-robust)" var: front_limit category: "Spread front";
	parameter "Front celerity (m/s) - CALIBRATE" var: front_celerity min: 0.005 max: 2.0 category: "Spread front";
	parameter "Auto pause at end" var: auto_pause category: "Engine";

	output {
		layout #split;
		display "Flood simulation" type: 2d background: #black {
			grid cell;
			graphics "static landscape" refresh: false {
				loop rp over: river_poly { draw rp.shape color: rgb(70, 130, 180, 120) border: #steelblue; }
			}
			species building;
			graphics "info" {
				draw string(current_date) + "   stage: " + (river_stage with_precision 2)
					+ " m   flooded: " + (flooded_area_km2 with_precision 1) + " km2"
					at: {world.shape.width * 0.02, world.shape.height * 0.03}
					color: #white font: font("SansSerif", 16, #bold);
			}
		}
		display "Time series" type: 2d {
			chart "Phu Tho flood event" type: series x_label: "hours since 07-09 00:00" {
				data "River stage (m)" value: river_stage color: #blue marker: false;
				data "Discharge (1000 m3/s)" value: river_discharge / 1000.0 color: #darkblue marker: false;
				data "Flooded area (km2)" value: flooded_area_km2 color: #red marker: false;
				data "Flood volume (10^6 m3)" value: flood_volume_mm3 color: #darkorange marker: false;
			}
		}
		monitor "Date" value: current_date;
		monitor "River discharge (m3/s)" value: river_discharge with_precision 1;
		monitor "River stage (m)" value: river_stage with_precision 2;
		monitor "Water surface L (m)" value: (river_stage + datum_offset) with_precision 2;
		monitor "Reachable cells" value: length(floodable);
		monitor "Front reach (m)" value: front_limit ? int(front_reach_m) : -1;
		monitor "Flooded area (km2)" value: flooded_area_km2 with_precision 2;
		monitor "Peak flooded (km2)" value: peak_flooded_km2 with_precision 2;
		monitor "Flood volume (10^6 m3)" value: flood_volume_mm3 with_precision 2;
		monitor "Buildings wet / flooded" value: "" + n_bldg_wet + " / " + n_bldg_flooded;
	}
}

// Display-free production / calibration experiment (numbers only, fastest).
experiment phutho_fastbathtub_fast type: gui {
	parameter "Use placeholder discharge (WaterDischarge_PhuTho.csv)" var: use_discharge_csv category: "Forcing";
	parameter "Rating curve exponent" var: rating_exponent category: "Forcing";
	parameter "Peak river stage (m)" var: peak_stage category: "Forcing";
	parameter "Stage cap (m) - hard ceiling" var: stage_cap min: 5.0 max: 60.0 category: "Forcing";
	parameter "Gauge datum offset (m) - CALIBRATE" var: datum_offset min: -15.0 max: 15.0 category: "Forcing";
	parameter "Limit spread by front (datum-robust)" var: front_limit category: "Spread front";
	parameter "Front celerity (m/s) - CALIBRATE" var: front_celerity min: 0.005 max: 2.0 category: "Spread front";

	output {
		display "Time series" type: 2d {
			chart "Phu Tho flood event" type: series x_label: "hours since 07-09 00:00" {
				data "River stage (m)" value: river_stage color: #blue marker: false;
				data "Discharge (1000 m3/s)" value: river_discharge / 1000.0 color: #darkblue marker: false;
				data "Flooded area (km2)" value: flooded_area_km2 color: #red marker: false;
				data "Flood volume (10^6 m3)" value: flood_volume_mm3 color: #darkorange marker: false;
			}
		}
		monitor "Date" value: current_date;
		monitor "River discharge (m3/s)" value: river_discharge with_precision 1;
		monitor "River stage (m)" value: river_stage with_precision 2;
		monitor "Water surface L (m)" value: (river_stage + datum_offset) with_precision 2;
		monitor "Front reach (m)" value: front_limit ? int(front_reach_m) : -1;
		monitor "Flooded area (km2)" value: flooded_area_km2 with_precision 2;
		monitor "Peak flooded (km2)" value: peak_flooded_km2 with_precision 2;
		monitor "Flood volume (10^6 m3)" value: flood_volume_mm3 with_precision 2;
	}
}
