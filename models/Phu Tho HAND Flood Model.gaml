/**
* Name: Phu Tho HAND Flood Model
* Author: Thành Đô Nguyễn (2026-06-17)
*
* A SLOPING-BASELINE (HAND) sibling of "Phu Tho FastBathtub Flood Model.gaml",
* built to be run side-by-side and COMPARED with the flat-pool version on the
* same Thao/Red-River reach at Phú Thọ (Typhoon Yagi, Sep-2024).
*
* WHY A SECOND MODEL.
*   The FastBathtub model uses ONE flat water surface L for the whole domain
*   (L = river_stage + datum_offset). On a mountainous reach where the channel
*   itself drops ~28 m (z 13..41 m), a flat surface is wrong: it leaves the
*   upstream channel dry and over-deepens the downstream end (11 m of water at
*   the crest). HEC-RAS instead computes a water-surface PROFILE that slopes with
*   the valley. This model approximates that with HAND (Height Above Nearest
*   Drainage):
*     ref_bed[c] = bed elevation of the river cell nearest to c (propagated along
*                  the geodesic shortest path to the river bank).
*     HAND[c]    = z_dyn[c] - ref_bed[c]   (height of the cell above its channel).
*     flow_depth(t) = river_baseflow + max(0, river_stage(t) - base_stage)
*                  = water depth above the local channel bed (baseflow + the
*                    gauge RISE above its base reading).
*   Each cell's LOCAL water level is WL[c] = ref_bed[c] + flow_depth(t); the cell
*   is wet when it is river-connected at that local level (spill_lvl <= WL[c]) and
*   the front has arrived. depth = WL[c] - z_dyn[c], bounded by flow_depth, so a
*   few metres everywhere - a sloping surface, like HEC-RAS, NOT a flat pool.
*
* KEY DIFFERENCE FROM THE FLAT-POOL MODEL.
*   - NO datum_offset. HAND references the DEM's own channel bed, so a uniform
*     vertical DEM shift cancels out - the model is datum-robust by construction.
*   - The global flat L is replaced by the per-cell local level WL[c].
*   Everything else (inputs, discharge forcing, fast-sweep fields, front limit,
*   bookkeeping, displays) is the same so the two are directly comparable.
*
* INPUTS: identical to the FastBathtub model (../includes/ and ../includes/phu-tho/).
*/
model PhuThoHANDFlood

global {

	// ------------------------------------------------------------------ input (same as flat-pool model)
	string dem_name     <- "dem-phutho-epsg-3857.tif";
	file dem_file       <- grid_file("../includes/phu-tho/" + dem_name);
	file river_file     <- shape_file("../includes/phu-tho/water_phutho.shp");
	file buildings_file <- shape_file("../includes/phu-tho/building_phutho.shp");
	file stage_file     <- csv_file("../includes/PhuThoStage2024_hourly.csv", ",", true);
	file discharge_file <- csv_file("../includes/WaterDischarge_PhuTho_hourly.csv", ",", true);
	geometry shape <- envelope(dem_file);

	// ------------------------------------------------------------------ time
	date starting_date <- date("2024-09-07 00:00:00");
	date end_date      <- date("2024-09-22 00:00:00");
	float step <- 1 #h;

	// ------------------------------------------------------------------ forcing (discharge-driven, same as flat-pool)
	bool  use_discharge_csv <- true;
	float rating_exponent <- 0.6;
	float base_stage <- 13.3;                 // m, gauge reading at lowest flow (RISE measured above this)
	float peak_stage <- 18.34;                // m, verified crest (Sep 11)
	float stage_cap  <- 30.0;                 // m, hard ceiling on river stage
	float river_stage <- 13.3;
	float river_discharge <- 0.0;
	list<date>  q_dates  <- [];
	list<float> q_values <- [];
	float q_min <- 1.0; float q_max <- 2.0;
	list<float> stage_series <- [];
	int   n_stage <- 0;

	// flow_depth = water depth ABOVE the local channel bed (baseflow + gauge rise).
	// This REPLACES the flat-pool's (river_stage + datum_offset). NO datum_offset
	// is needed - HAND references the DEM's own bed, so it is datum-robust.
	float river_baseflow <- 1.5;              // m, channel depth at base flow
	float flow_depth <- 1.5;                  // m, current depth above the local bed (set each hour)

	// ------------------------------------------------------------------ engine / thresholds
	float flood_threshold <- 0.05;
	float wet_thr   <- 0.1;
	float flood_thr <- 1.0;
	int   max_spill_sweeps <- 200;
	bool  auto_pause <- true;

	// ------------------------------------------------------------------ spread front
	bool  front_limit <- true;
	float front_celerity <- 0.01;

	// ------------------------------------------------------------------ bookkeeping
	int grid_cols; int grid_rows;
	float cell_dx <- 30.0;
	float z_min <- 0.0; float z_max <- 1.0;
	float SPILL_BIG <- 1e9;
	list<cell> river_cells <- [];
	list<cell> floodable   <- [];
	list<cell> metric_cells<- [];
	float flooded_area_km2 <- 0.0;
	float flood_volume_mm3 <- 0.0;
	float peak_flooded_km2 <- 0.0;
	int   n_bldg_wet <- 0;
	int   n_bldg_flooded <- 0;
	float front_reach_m <- 0.0;
	bool  sim_finished <- false;

	init {
		write "=== Phu Tho HAND Flood (sloping-baseline, datum-robust) ===";
		grid_cols <- 1 + max(cell collect each.grid_x);
		grid_rows <- 1 + max(cell collect each.grid_y);
		cell_dx <- first(cell).shape.width;
		if cell_dx < 1.0 {
			write "WARNING: cell_dx = " + cell_dx + " (looks like DEGREES). Reproject the DEM to a metric CRS.";
		}

		// terrain + neighbours + inline pit fill
		ask cell {
			z <- grid_value;
			is_nodata <- z < -1000.0;
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

		create river_poly from: river_file;
		create building from: buildings_file;
		ask river_poly { ask cell overlapping self where (!each.is_nodata) { is_river <- true; } }
		river_cells <- cell where each.is_river;
		ask building { my_cell <- first(cell overlapping location); }
		ask building where (each.my_cell = nil) { do die; }
		ask cell { h0 <- h; }

		// stage record (fallback) + discharge record (driver) -- same parsing as flat-pool
		matrix sm <- matrix(stage_file);
		loop r over: rows_list(sm) {
			string ds <- string(r[0]);
			if length(ds) > 0 and ds != "datetime" { stage_series <+ float(r[1]); }
		}
		n_stage <- length(stage_series);
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
			write "discharge CSV empty -> using measured stage gauge";
		} else {
			q_min <- min(q_values); q_max <- max(q_values);
			write "Discharge: " + length(q_values) + " values, " + (q_min with_precision 1) + "-" + (q_max with_precision 1)
				+ " m3/s" + (use_discharge_csv ? "  (DRIVING via rating curve)" : "  (loaded; gauge driving)");
		}

		// HAND + spill + front fields (static)
		do compute_fields;

		ask cell where (!each.is_nodata) {
			float shade <- (z - z_min) / max(0.001, z_max - z_min);
			terrain_color <- is_river ? rgb(70, 130, 180)
				: rgb(70 + int(150 * shade), 80 + int(130 * shade), 60 + int(110 * shade));
			color <- terrain_color;
		}

		river_stage <- stage_at(starting_date);
		flow_depth  <- river_baseflow + max(0.0, river_stage - base_stage);
		ask river_cells { h <- flow_depth; }
		do refresh_colors;

		write "river cells: " + length(river_cells) + " | floodable cells: " + length(floodable)
			+ " | buildings in domain: " + length(building);
		write "Init done. Simulation: " + starting_date + " -> " + end_date;
	}

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

	float measured_stage_at (date d) {
		if n_stage = 0 { return base_stage; }
		int hr <- int((d - starting_date) / 3600.0);
		if hr < 0        { return first(stage_series); }
		if hr >= n_stage { return last(stage_series); }
		return stage_series[hr];
	}

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
		return min(stage_cap, s);
	}

	reflex update_stage {
		river_discharge <- discharge_at(current_date);
		river_stage <- stage_at(current_date);
		flow_depth  <- river_baseflow + max(0.0, river_stage - base_stage);   // depth above the LOCAL bed
	}

	// ====================================================================== HAND + spill + front fields
	// Same fast sweeping as the flat-pool model, with ONE addition: alongside the
	// geodesic distance we carry ref_bed = the bed elevation of the river cell at
	// the start of the shortest path (the NEAREST DRAINAGE). HAND = z_dyn - ref_bed.
	action compute_fields {
		float inv_cel <- 1.0 / max(1e-6, front_celerity);
		ask cell {
			passable <- !is_river and !is_nodata;
			spill_lvl <- SPILL_BIG;
			front_dist <- SPILL_BIG;
			front_arrival <- SPILL_BIG;
			ref_bed <- SPILL_BIG;
		}
		// seed: land cells touching the river. ref_bed = lowest adjacent channel bed.
		ask cell where each.is_river {
			loop nb over: [nE, nW, nN, nS] {
				if nb != nil and nb.passable {
					if nb.z_dyn < nb.spill_lvl { nb.spill_lvl <- nb.z_dyn; }
					if nb.front_dist > 0.0 { nb.front_dist <- 0.0; nb.front_arrival <- 0.0; }
					if z_dyn < nb.ref_bed { nb.ref_bed <- z_dyn; }
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
						float bestRef <- ce.ref_bed;
						float bestA <- ce.front_arrival;
						loop nb over: [ce.nE, ce.nW, ce.nN, ce.nS] {
							if nb != nil and nb.passable {
								float candS <- max(nb.spill_lvl, ce.z_dyn);
								if candS < bestS { bestS <- candS; }
								float candD <- nb.front_dist + cell_dx;
								if candD < bestD { bestD <- candD; bestRef <- nb.ref_bed; }   // carry nearest bed
								float candA <- nb.front_arrival + cell_dx * inv_cel;
								if candA < bestA { bestA <- candA; }
							}
						}
						if bestS < ce.spill_lvl - 1e-6   { ce.spill_lvl <- bestS;     changed <- true; }
						if bestD < ce.front_dist - 1e-3  { ce.front_dist <- bestD; ce.ref_bed <- bestRef; changed <- true; }
						if bestA < ce.front_arrival - 1.0 { ce.front_arrival <- bestA; changed <- true; }
					}
				}
			}
		}
		floodable <- cell where (each.passable and each.spill_lvl < 0.5 * SPILL_BIG and each.ref_bed < 0.5 * SPILL_BIG);
		ask floodable { hand <- z_dyn - ref_bed; }
		metric_cells <- floodable;
		write "HAND + spill + front fields built in " + rounds + " rounds; " + length(floodable)
			+ " cells reachable (HAND "
			+ ((empty(floodable) ? 0.0 : floodable min_of each.hand) with_precision 2) + ".."
			+ ((empty(floodable) ? 0.0 : floodable max_of each.hand) with_precision 2) + " m, ref_bed "
			+ ((empty(floodable) ? 0.0 : floodable min_of each.ref_bed) with_precision 1) + ".."
			+ ((empty(floodable) ? 0.0 : floodable max_of each.ref_bed) with_precision 1) + " m).";
	}

	// ====================================================================== dynamic flood (HAND engine)
	// Per-cell LOCAL water level WL[c] = ref_bed[c] + flow_depth(t). A cell is wet
	// when the front has arrived AND it is river-connected at its own local level
	// (spill_lvl <= WL[c]); depth = WL[c] - z_dyn, bounded by flow_depth. This is a
	// SLOPING surface (high upstream where ref_bed is high, low downstream), unlike
	// the flat-pool model's single global L.
	reflex dynamic_flood {
		ask river_cells { h <- flow_depth; wsl <- z_dyn + flow_depth; }   // channel = bed + flow_depth
		if !empty(floodable) {
			float elapsed_s <- front_limit ? (current_date - starting_date) : SPILL_BIG;
			front_reach_m <- front_limit ? front_celerity * (current_date - starting_date) : SPILL_BIG;
			ask (floodable where (each.front_arrival <= elapsed_s and (each.ref_bed + flow_depth) >= each.spill_lvl)) parallel: true {
				float WL <- ref_bed + flow_depth;
				h <- WL - z_dyn;          // >= 0 since spill_lvl >= z_dyn and WL >= spill_lvl
				wsl <- z_dyn + h;
			}
		}
	}

	// ====================================================================== bookkeeping (same as flat-pool)
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
			write "" + current_date + " | stage " + (river_stage with_precision 2) + " m | flow_depth "
				+ (flow_depth with_precision 2) + " m | flooded " + (flooded_area_km2 with_precision 2)
				+ " km2 | vol " + (flood_volume_mm3 with_precision 1) + " Mm3 | bldg flooded " + n_bldg_flooded;
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
//  Raster domain
// ==========================================================================
grid cell file: dem_file neighbors: 4
	use_regular_agents: false use_individual_shapes: false use_neighbors_cache: false schedules: [] {
	float z;
	float z_dyn;
	float h <- 0.0;
	float h0 <- 0.0;
	float wsl <- 0.0;
	float spill_lvl <- 1e9;     // S[c]: lowest stage at which the cell connects to the river
	float ref_bed   <- 1e9;     // bed elevation of the NEAREST river cell (drainage reference)
	float hand      <- 0.0;     // HAND: z_dyn - ref_bed (height above nearest drainage)
	float front_dist <- 1e9;
	float front_arrival <- 1e9;
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
//  Vector species
// ==========================================================================
species river_poly schedules: [] { aspect default { draw shape color: rgb(70, 130, 180, 120) border: #steelblue; } }
species building   schedules: [] {
	cell my_cell;
	float depth_w <- 0.0;
	int status <- 0;
	aspect default { draw shape color: status = 2 ? #red : (status = 1 ? #orange : rgb(90, 90, 90)); }
}

// ==========================================================================
//  Experiments
// ==========================================================================
experiment phutho_hand type: gui {
	parameter "Use discharge (WaterDischarge_PhuTho_hourly.csv)" var: use_discharge_csv category: "Forcing";
	parameter "Rating curve exponent" var: rating_exponent category: "Forcing";
	parameter "Base river stage (m)" var: base_stage category: "Forcing";
	parameter "Peak river stage (m)" var: peak_stage category: "Forcing";
	parameter "Stage cap (m) - hard ceiling" var: stage_cap min: 5.0 max: 60.0 category: "Forcing";
	parameter "Channel baseflow depth (m)" var: river_baseflow min: 0.0 max: 5.0 category: "Forcing";
	parameter "Flood threshold (m)" var: flood_threshold min: 0.01 max: 0.5 category: "Engine";
	parameter "Limit spread by front" var: front_limit category: "Spread front";
	parameter "Front celerity (m/s)" var: front_celerity min: 0.005 max: 2.0 category: "Spread front";
	parameter "Auto pause at end" var: auto_pause category: "Engine";

	output {
		layout #split;
		display "Flood simulation (HAND)" type: 2d background: #black {
			grid cell;
			graphics "static landscape" refresh: false {
				loop rp over: river_poly { draw rp.shape color: rgb(70, 130, 180, 120) border: #steelblue; }
			}
			species building;
			graphics "info" {
				draw string(current_date) + "   stage: " + (river_stage with_precision 2)
					+ " m   flow depth: " + (flow_depth with_precision 2)
					+ " m   flooded: " + (flooded_area_km2 with_precision 1) + " km2"
					at: {world.shape.width * 0.02, world.shape.height * 0.03}
					color: #white font: font("SansSerif", 16, #bold);
			}
		}
		display "Time series" type: 2d {
			chart "Phu Tho flood event (HAND)" type: series x_label: "hours since 07-09 00:00" {
				data "River stage (m)" value: river_stage color: #blue marker: false;
				data "Flow depth above bed (m)" value: flow_depth color: #teal marker: false;
				data "Discharge (1000 m3/s)" value: river_discharge / 1000.0 color: #darkblue marker: false;
				data "Flooded area (km2)" value: flooded_area_km2 color: #red marker: false;
				data "Flood volume (10^6 m3)" value: flood_volume_mm3 color: #darkorange marker: false;
			}
		}
		monitor "Date" value: current_date;
		monitor "River discharge (m3/s)" value: river_discharge with_precision 1;
		monitor "River stage (m)" value: river_stage with_precision 2;
		monitor "Flow depth above bed (m)" value: flow_depth with_precision 2;
		monitor "Reachable cells" value: length(floodable);
		monitor "Front reach (m)" value: front_limit ? int(front_reach_m) : -1;
		monitor "Flooded area (km2)" value: flooded_area_km2 with_precision 2;
		monitor "Peak flooded (km2)" value: peak_flooded_km2 with_precision 2;
		monitor "Flood volume (10^6 m3)" value: flood_volume_mm3 with_precision 2;
		monitor "Buildings wet / flooded" value: "" + n_bldg_wet + " / " + n_bldg_flooded;
	}
}

experiment phutho_hand_fast type: gui {
	parameter "Use discharge (WaterDischarge_PhuTho_hourly.csv)" var: use_discharge_csv category: "Forcing";
	parameter "Peak river stage (m)" var: peak_stage category: "Forcing";
	parameter "Channel baseflow depth (m)" var: river_baseflow min: 0.0 max: 5.0 category: "Forcing";
	parameter "Limit spread by front" var: front_limit category: "Spread front";
	parameter "Front celerity (m/s)" var: front_celerity min: 0.005 max: 2.0 category: "Spread front";

	output {
		display "Time series" type: 2d {
			chart "Phu Tho flood event (HAND)" type: series x_label: "hours since 07-09 00:00" {
				data "River stage (m)" value: river_stage color: #blue marker: false;
				data "Flow depth above bed (m)" value: flow_depth color: #teal marker: false;
				data "Flooded area (km2)" value: flooded_area_km2 color: #red marker: false;
				data "Flood volume (10^6 m3)" value: flood_volume_mm3 color: #darkorange marker: false;
			}
		}
		monitor "Date" value: current_date;
		monitor "River stage (m)" value: river_stage with_precision 2;
		monitor "Flow depth above bed (m)" value: flow_depth with_precision 2;
		monitor "Flooded area (km2)" value: flooded_area_km2 with_precision 2;
		monitor "Peak flooded (km2)" value: peak_flooded_km2 with_precision 2;
		monitor "Flood volume (10^6 m3)" value: flood_volume_mm3 with_precision 2;
	}
}
