/**
* Name: Dong Hoi FastBathtub Flood Model (Discharge-driven)
* Author: Thành Đô Nguyễn (2026-06-25)
*
* A fork of "Dong Hoi FastBathtub Flood Model.gaml" that DROPS the river-stage-file
* mechanism entirely: there is no measured-gauge CSV path (no stage_file,
* stage_series, n_stage or measured_stage_at). The river stage is driven SOLELY by
* the observed discharge (WaterDischarge_HamNinh.csv) through the Manning rating
* curve. If the discharge record is missing the stage simply holds at base_stage.
*
* Everything else is IDENTICAL to the parent: the dyke-free, river-seeded,
* spread-limited, connectivity-constrained, hysteretic level-pool engine, the
* DEM (dong-hoi_utm48n.tif), the river polygon stage boundary, the buildings,
* the observation points and the per-step GeoTIFF export.
*
* INPUTS (../includes/ and ../includes/dong-hoi/).
*   DEM: dong-hoi_utm48n.tif (EPSG:32648 metres).
*   River: water_donghoi.shp (Nhật Lệ estuary polygon) - the stage boundary.
*   Buildings: building_multipolygon.shp (clipped to the domain in init).
*   Forcing: WaterDischarge_HamNinh.csv -> Manning rating curve -> river stage.
*   NOTE: DongHoiStage2020_hourly.csv is deliberately NOT used.
*/
model DongHoiFastBathtubFloodDischarge

global {

	// ------------------------------------------------------------------ input
	// DEM in UTM zone 48N (EPSG:32648) metres (see the parent model for why 3857 is broken here).
	string dem_name     <- "dong-hoi_utm48n.tif";
	file dem_file       <- grid_file("../includes/dong-hoi/" + dem_name);
	file river_file     <- shape_file("../includes/dong-hoi/water_donghoi.shp");
	file buildings_file <- shape_file("../includes/dong-hoi/building_multipolygon.shp");
	file points_file    <- shape_file("../includes/dong-hoi/depth_arrival_time.shp");  // observation points (id), reprojected to UTM 48N
	file discharge_file <- csv_file("../includes/WaterDischarge_HamNinh.csv", ",", true);
	geometry shape <- envelope(dem_file);

	// ------------------------------------------------------------------ time
	date starting_date <- date("2020-10-07 00:00:00");
	date end_date      <- date("2020-10-08 23:00:00");
	float step <- 1 #h;                       // ONE CYCLE = ONE HOUR

	// ------------------------------------------------------------------ river stage forcing
	// DISCHARGE-DRIVEN ONLY (no measured-gauge stage file). Observed discharge
	// (WaterDischarge_HamNinh.csv) -> Manning rating curve -> stage. Anchors:
	// q_min -> base_stage, q_max -> peak_stage. If the discharge record is missing
	// the stage holds at base_stage (no gauge/stage-file fallback).
	bool  use_discharge_csv <- true;
	float rating_exponent <- 0.6;             // Manning h ~ Q^(3/5)
	float base_stage <- 1.0;                  // m, stage at the lowest recorded discharge
	float peak_stage <- 7.0;                  // m, stage at the peak discharge (≈ observed WL_DongHoi peak)
	float river_stage <- 7.0;
	float river_discharge <- 0.0;
	list<date>  q_dates  <- [];
	list<float> q_values <- [];
	float q_min <- 1.0; float q_max <- 2.0;

	// datum_offset = (gauge zero) - (DEM/SRTM zero), in metres. THE calibration knob:
	// water surface used by the engine is  L = river_stage + datum_offset.
	float datum_offset <- 0;

	// ------------------------------------------------------------------ engine / thresholds
	float flood_threshold <- 0.05;  // m, depth above baseline counted as flooded
	float wet_thr   <- 0.1;         // m, building WET threshold
	float flood_thr <- 1.0;         // m, building FLOODED threshold
	int   max_spill_sweeps <- 200;  // fast-sweeping rounds for the spill + front fields
	bool  auto_pause <- true;

	// ------------------------------------------------------------------ spread front (datum-robust knob)
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

	// ------------------------------------------------------------------ geotiff export
	bool   export_geotiff <- true;
	string export_dir <- "../exported_results/donghoi_fastbathtub_discharge/";

	init {
		write "=== Dong Hoi FastBathtub Flood (river-seeded, DISCHARGE-driven stage) ===";
		grid_cols <- 1 + max(cell collect each.grid_x);
		grid_rows <- 1 + max(cell collect each.grid_y);
		cell_dx <- first(cell).shape.width;

		// --- terrain + neighbour references + inline pit fill -----------------
		ask cell {
			z <- grid_value;
			is_nodata <- z < -1000.0;           // -32767 sentinel (none in this DEM, but guard)
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

		// river/estuary cells = stage boundary (held at the stage every hour)
		ask river_poly { ask cell overlapping self where (!each.is_nodata) { is_river <- true; } }
		river_cells <- cell where each.is_river;

		// buildings: bind to a cell, drop those outside the DEM domain
		ask building { my_cell <- first(cell overlapping location); }
		ask building where (each.my_cell = nil) { do die; }

		// observation points: bind each to its cell (arrival time recorded as the flood reaches it)
		create observation_point from: points_file with: [pid::int(read("id"))];
		ask observation_point { my_cell <- first(cell overlapping self); }

		// baseline (flood metrics measured ABOVE this; no lakes here so h0 = 0)
		ask cell { h0 <- h; }

		// --- observed discharge record (the SOLE driver; dates are M/D/YYYY H:MM) -
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
			write "WaterDischarge_HamNinh.csv empty -> stage holds at base_stage (no stage-file fallback)";
		} else {
			q_min <- min(q_values); q_max <- max(q_values);
			write "Discharge (driver): " + length(q_values) + " values, " + first(q_dates) + ".." + last(q_dates)
				+ ", " + (q_min with_precision 1) + "-" + (q_max with_precision 1) + " m3/s";
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
		ask river_cells { h <- max(0.0, river_stage + datum_offset - z_dyn); }
		do refresh_colors;

		write "river cells: " + length(river_cells) + " | floodable cells: " + length(floodable)
			+ " | buildings in domain: " + length(building);
		write "Init done. Simulation: " + starting_date + " -> " + end_date;
	}

	// ====================================================================== filename helpers
	string pad2 (int v) { return (v < 10 ? "0" : "") + v; }
	string pad4 (int v) { string s <- "" + v; loop while: (length(s) < 4) { s <- "0" + s; } return s; }

	// ====================================================================== forcing
	// observed discharge interpolated in time
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

	// discharge -> Manning rating curve -> stage: the SOLE forcing. No measured-gauge
	// stage file. If the discharge record is missing the stage holds at base_stage.
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

	// ====================================================================== spill + front fields
	// Two static fields, both solved by the same fast sweeping (4 directional
	// orders). The river is the boundary: every passable LAND cell adjacent to a
	// river cell is a seed. River and no-data cells are barriers.
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
	reflex dynamic_flood {
		float L <- river_stage + datum_offset;
		ask river_cells { h <- max(0.0, L - z_dyn); }       // river held at the stage (boundary)
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

		// record the flood ARRIVAL TIME at each observation point (first hour its cell goes wet)
		ask observation_point where (each.arrival_h < 0.0) {
			if my_cell != nil and (my_cell.h - my_cell.h0) > flood_threshold {
				arrival_h <- (current_date - starting_date) / 3600.0;
				write "Observation point " + pid + " reached on " + current_date
					+ " (depth " + ((my_cell.h - my_cell.h0) with_precision 2) + " m)";
			}
		}

		do refresh_colors;
		if current_date.hour mod 3 = 0 {
			write "" + current_date + " | stage " + (river_stage with_precision 2) + " m (L "
				+ ((river_stage + datum_offset) with_precision 2) + ") | flooded "
				+ (flooded_area_km2 with_precision 2) + " km2 | vol "
				+ (flood_volume_mm3 with_precision 1) + " Mm3 | bldg flooded " + n_bldg_flooded;
		}
	}

	// ====================================================================== geotiff export (every step)
	reflex export_water_height when: export_geotiff {
		list<cell> wet_land <- metric_cells where ((each.h - each.h0) > flood_threshold and !each.is_river);
		float peak_h <- empty(wet_land) ? 0.0 : wet_land max_of each.h;
		// band value = water depth: river channel depth on river cells, flood depth on
		// flooded land, 0 on dry land, -9999 on no-data
		ask cell { grid_value <- is_nodata ? -9999.0 : ((is_river or ((h - h0) > flood_threshold)) ? h : 0.0); }
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
	float z;                 // raw DEM elevation (m, UTM 48N)
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
species observation_point schedules: [] {
	int pid;
	cell my_cell;
	float arrival_h <- -1.0;          // hours since start when the flood first reached this point (-1 = never)
	aspect default {
		draw circle(120) color: arrival_h < 0.0 ? #white : #red border: #black;
		draw string(pid) + (arrival_h < 0.0 ? "" : (" : " + (arrival_h with_precision 1) + " h"))
			at: location + {150, -100} color: #black font: font("SansSerif", 14, #bold);
	}
}

// ==========================================================================
//  Experiments
// ==========================================================================
experiment donghoi_fastbathtub_discharge type: gui {
	parameter "Use observed discharge (WaterDischarge_HamNinh.csv)" var: use_discharge_csv category: "Forcing";
	parameter "Rating curve exponent" var: rating_exponent category: "Forcing";
	parameter "Base river stage (m)" var: base_stage category: "Forcing";
	parameter "Peak river stage (m)" var: peak_stage category: "Forcing";
	parameter "Gauge datum offset (m) - CALIBRATE" var: datum_offset min: -15.0 max: 5.0 category: "Forcing";
	parameter "Flood threshold (m)" var: flood_threshold min: 0.01 max: 0.5 category: "Engine";
	parameter "Limit spread by front (datum-robust)" var: front_limit category: "Spread front";
	parameter "Front celerity (m/s) - CALIBRATE" var: front_celerity min: 0.005 max: 2.0 category: "Spread front";
	parameter "Export water-height GeoTIFF each step" var: export_geotiff category: "Export";
	parameter "Auto pause at end" var: auto_pause category: "Engine";

	output {
		layout #split;
		display "Flood simulation" type: 2d background: #black {
			grid cell;
			graphics "static landscape" refresh: false {
				loop rp over: river_poly { draw rp.shape color: rgb(70, 130, 180, 120) border: #steelblue; }
			}
			species building;
			species observation_point;
			graphics "info" {
				draw string(current_date) + "   stage: " + (river_stage with_precision 2)
					+ " m   flooded: " + (flooded_area_km2 with_precision 1) + " km2"
					at: {world.shape.width * 0.02, world.shape.height * 0.03}
					color: #white font: font("SansSerif", 16, #bold);
			}
		}
		display "Time series" type: 2d {
			chart "Dong Hoi flood event" type: series x_label: "hours since 07-10 00:00" {
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
experiment donghoi_fastbathtub_discharge_fast type: gui {
	parameter "Use observed discharge (WaterDischarge_HamNinh.csv)" var: use_discharge_csv category: "Forcing";
	parameter "Rating curve exponent" var: rating_exponent category: "Forcing";
	parameter "Peak river stage (m)" var: peak_stage category: "Forcing";
	parameter "Gauge datum offset (m) - CALIBRATE" var: datum_offset min: -15.0 max: 5.0 category: "Forcing";
	parameter "Limit spread by front (datum-robust)" var: front_limit category: "Spread front";
	parameter "Front celerity (m/s) - CALIBRATE" var: front_celerity min: 0.005 max: 2.0 category: "Spread front";
	parameter "Export water-height GeoTIFF each step" var: export_geotiff category: "Export";

	output {
		display "Time series" type: 2d {
			chart "Dong Hoi flood event" type: series x_label: "hours since 07-10 00:00" {
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
