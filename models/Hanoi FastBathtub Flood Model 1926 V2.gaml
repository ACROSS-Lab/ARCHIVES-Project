/**
* Name: Hanoi FastBathtub Flood Model 1926 V2 (gate-volume budget)
* Author: Thành Đô Nguyễn (2026-06-13)
*
* SAME GOAL as V1 ("Hanoi FastBathtub Flood Model 1926.gaml"): the speed of the
* FastFlood/LISEM relaxation model with the flood extent of the ABM. V2 is a
* DIFFERENT limiter, built so you can compare it against V1 side by side.
*
*   V1 limits the spread GEOMETRICALLY: a flood front advances from the breaches
*       at a fixed celerity; cells within the front AND below the stage flood.
*       Downside (seen on gz40): given enough time the front covers the whole
*       basin and the depth is still stage - z everywhere, so it relaxes back to
*       the full datum-sensitive bathtub.
*
*   V2 limits the spread by a MASS BUDGET, exactly like the ABM. Each hour the
*       breaches pass a weir discharge
*           Q = weir_coef * width * (stage - invert)^1.5 * f_Villemonte(pool)
*       (the same broad-crested weir + Villemonte submergence as the ABM gate),
*       the volume is accumulated, and that volume is poured into the connected
*       basin as a flat pool whose LEVEL is found by volume conservation
*       (volume(level) = accumulated volume), capped at the river stage. The
*       Villemonte factor shuts the gate as the pool approaches the stage, so the
*       inflow self-throttles - the extent is set by how much water actually came
*       in, not by "everything below the stage". Connected cells track the pool
*       up and down (drain back through the gate); perched cells keep their water
*       -> trapped casiers on the receding limb.
*
*   Still one volume->level solve (a 40-step bisection over the connected cells)
*   plus one filtered parallel ask per hour - no iteration -> LISEM-class speed.
*
* WHY THIS MIGHT (or might NOT) be more datum-robust than V1: the gate volume is
*   anchored to the LOCAL breach invert and self-throttles, so the delivered
*   volume is similar on the gz40 and gz50 DEMs. BUT this is a LUMPED single-pool
*   model: it spreads the volume as one flat sheet over the whole connected
*   basin, so on the flat low gz40 the pool level rises slowly, the gate stays
*   wide open, and it may still fill toward the full bathtub. The ABM avoids that
*   because its water piles up LOCALLY near each gate (diffusion is slow), which
*   raises the level AT the gate and throttles it early. V2 cannot reproduce that
*   local pile-up. Treat V2 as a test of "is the gate mass budget alone enough?"
*   - if it is still datum-sensitive on gz40, the conclusion is that the ABM's
*   robustness needs the spatial/diffusion structure (i.e. accelerate the ABM
*   rather than replace it), or combine V1's front with V2's budget (a future V3).
*
* Inputs, forcing and breach scenario are identical to V1 (DEM selectable; default
* mnt-gz40.asc to match the ABM).
*/
model HanoiFastBathtubFlood1926V2

global {

	// ------------------------------------------------------------------ input
	string dem_name     <- "mnt-gz10.asc";   // default = the ABM's archival DEM; mnt-gz50.asc = modern (LISEM)
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
	float step <- 1 #h;

	// ------------------------------------------------------------------ river stage forcing (identical to V1)
	bool  use_discharge_csv <- true;
	float rating_exponent <- 0.6;
	float base_stage <- 7.0;
	float peak_stage <- 11.93;
	date  peak_date  <- date("1926-07-30 00:00:00");
	float sigma_rise_days <- 3.5;
	float sigma_fall_days <- 3.0;
	float river_stage <- base_stage;
	float river_discharge <- 0.0;
	list<date>  q_dates  <- [];
	list<float> q_values <- [];
	float q_min <- 1.0; float q_max <- 2.0;

	// ------------------------------------------------------------------ breach scenario (identical to V1)
	int   breach_hour <- 6;
	float breach_floor_min <- 2.0;
	float breach_freeboard <- 0.2;
	float breach_search_radius <- 300.0;
	float breach_cut_halfwidth <- 60.0;

	// ------------------------------------------------------------------ gate-volume engine (the V2 knob)
	// weir_coef is the ABM's broad-crested weir coefficient; raise it -> more
	// inflow per hour -> faster fill / larger extent. This is the calibration
	// knob (analogue of V1's front_celerity). The Villemonte submergence factor
	// (gate shuts as the pool nears the stage) is hard-wired, as in the ABM.
	float weir_coef <- 1.5;
	float gate_width_fraction <- 1.0;   // share of each breach crest acting as a gate (ABM: breach_width_fraction)

	// ------------------------------------------------------------------ thresholds / misc
	float flood_threshold <- 0.05;
	float wet_thr   <- 0.1;
	float flood_thr <- 1.0;
	int   max_spill_sweeps <- 120;
	float lake_initial_depth <- 0.5;
	float datum_offset <- 0.0;
	bool  auto_pause <- true;

	// ------------------------------------------------------------------ bookkeeping
	int grid_cols; int grid_rows;
	float cell_dx <- 50.0;
	float cell_a <- 2500.0;                 // cell area (m2)
	float z_min <- 0.0; float z_max <- 1.0;
	float SPILL_BIG <- 1e9;
	list<cell> river_cells <- [];
	list<cell> lake_cells  <- [];
	list<cell> floodable   <- [];
	list<cell> metric_cells<- [];
	bool spill_dirty <- false;
	int  breaches_open <- 0;
	float cumulative_V <- 0.0;              // m3 of water delivered by the gates, net (the budget)
	float fill_level <- -1e9;               // m, current flat-pool level over the connected basin
	float gate_Q_now <- 0.0;                // m3/s, net gate flow this hour (+ = into land)
	float flooded_area_km2 <- 0.0;
	float flood_volume_mm3 <- 0.0;
	float peak_flooded_km2 <- 0.0;
	int   n_bldg_wet <- 0;
	int   n_bldg_flooded <- 0;
	date  first_break_time <- nil;
	bool  sim_finished <- false;

	init {
		write "=== Hanoi FastBathtub Flood 1926 V2 (gate-volume budget) ===";
		grid_cols <- 1 + max(cell collect each.grid_x);
		grid_rows <- 1 + max(cell collect each.grid_y);
		cell_dx <- first(cell).shape.width;
		cell_a  <- cell_dx * cell_dx;

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
			z_dyn <- z + pit;
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
		ask cell where each.is_dyke { is_river <- false; }
		river_cells <- cell where each.is_river;

		ask lake { ask cell overlapping self where (!each.is_river and !each.is_dyke) { is_lake <- true; } }
		lake_cells <- cell where each.is_lake;
		ask lake_cells { h <- lake_initial_depth; wsl <- z_dyn + h; }
		metric_cells <- list(lake_cells);

		ask building { my_cell <- first(cell overlapping location); }
		ask observation_point { my_cell <- first(cell overlapping self); }

		ask cell { h0 <- h; }

		// --- observed discharge record ---------------------------------------
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
			write "WaterDischarge.csv empty -> gaussian hydrograph fallback";
		} else {
			q_min <- min(q_values); q_max <- max(q_values);
			write "Discharge: " + length(q_values) + " values, " + first(q_dates) + ".." + last(q_dates)
				+ ", " + q_min + "-" + q_max + " m3/s";
		}

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
		write "strict-dyke gate-volume pool: dry until the first breach (" + first_break_time + ").";
		write "Init done. Simulation: " + starting_date + " -> " + end_date;
	}

	// ====================================================================== forcing (identical to V1)
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

	float stage_at (date d) {
		if use_discharge_csv and !empty(q_values) {
			float q <- discharge_at(d);
			float fq <- (q ^ rating_exponent - q_min ^ rating_exponent)
			          / max(1e-6, q_max ^ rating_exponent - q_min ^ rating_exponent);
			return base_stage + (peak_stage - base_stage) * min(1.0, max(0.0, fq));
		}
		float t_days <- (d - peak_date) / 86400.0;
		float sigma <- t_days < 0.0 ? sigma_rise_days : sigma_fall_days;
		return base_stage + (peak_stage - base_stage) * exp(-0.5 * (t_days / sigma) ^ 2);
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
			do compute_spill;
			spill_dirty <- false;
		}
	}

	// ====================================================================== spill field (connectivity, V1's S only)
	// S[c] = lowest stage at which c connects to a breach (minimax path level from
	// the corridor seeds; dykes + river are barriers). Recomputed only on breach.
	action compute_spill {
		ask cell {
			passable <- (!is_river and !is_dyke) or is_corridor;
			spill_lvl <- SPILL_BIG;
		}
		ask cell where each.is_corridor { spill_lvl <- z_dyn; }
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
						loop nb over: [ce.nE, ce.nW, ce.nN, ce.nS] {
							if nb != nil and nb.passable {
								float candS <- max(nb.spill_lvl, ce.z_dyn);
								if candS < bestS { bestS <- candS; }
							}
						}
						if bestS < ce.spill_lvl - 1e-6 { ce.spill_lvl <- bestS; changed <- true; }
					}
				}
			}
		}
		floodable    <- cell where (each.passable and each.spill_lvl < 0.5 * SPILL_BIG);
		metric_cells <- remove_duplicates(floodable + lake_cells);
		if fill_level < -1e8 and !empty(floodable) { fill_level <- floodable min_of each.z_dyn; }
		write "" + current_date + "  spill field rebuilt in " + rounds + " sweep rounds; "
			+ length(floodable) + " cells reachable by a breach (S "
			+ ((empty(floodable) ? 0.0 : floodable min_of each.spill_lvl) with_precision 2) + ".."
			+ ((empty(floodable) ? 0.0 : floodable max_of each.spill_lvl) with_precision 2) + " m).";
	}

	// ====================================================================== gate-volume helpers
	// broad-crested weir + Villemonte submergence, identical form to the ABM gate
	float weir_q (float invert, float w, float up_lvl, float dn_lvl) {
		float hu <- up_lvl - invert;
		if hu <= 0.0 { return 0.0; }
		float hd <- max(0.0, dn_lvl - invert);
		float f <- hd >= hu ? 0.0 : (1.0 - (hd / hu) ^ 1.5) ^ 0.385;
		return weir_coef * w * (hu ^ 1.5) * f;        // m3/s
	}

	// volume of water held in the breach-connected basin if filled to a flat level L
	float pool_volume (float L) {
		if empty(floodable) { return 0.0; }
		return (floodable sum_of (each.spill_lvl <= L ? (L - each.z_dyn) : 0.0)) * cell_a;
	}

	// ====================================================================== dynamic flood (gate-volume engine)
	reflex dynamic_flood {
		float L <- river_stage + datum_offset;
		ask river_cells { h <- max(0.0, L - z_dyn); }
		if !empty(floodable) and breaches_open > 0 {
			// 1. net weir inflow this hour, throttled by the current pool level
			float dV <- 0.0;
			loop d over: (dyke where each.opened) {
				if L > fill_level {
					dV <- dV + weir_q(d.breach_invert, d.breach_width, L, fill_level) * step;        // river -> land
				} else {
					dV <- dV - weir_q(d.breach_invert, d.breach_width, fill_level, L) * step;        // land -> river (drain)
				}
			}
			// 2. accumulate the budget, but move MONOTONICALLY toward the stage-equilibrium
			//    volume without crossing it. The explicit weir flow over a 1 h step is far
			//    larger than the head difference warrants, so an unclamped update over-fills
			//    then over-drains -> a fill/drain limit cycle (violent on a small basin like
			//    gz50). Clamping at the equilibrium removes the oscillation.
			float Veq <- pool_volume(L);
			float newV <- cumulative_V + dV;
			newV <- (cumulative_V <= Veq) ? min(newV, Veq) : max(newV, Veq);
			newV <- max(0.0, newV);
			gate_Q_now   <- (newV - cumulative_V) / step;   // actual net flow after the clamp
			cumulative_V <- newV;
			// 3. invert volume -> flat pool level (bisection over the connected basin)
			float lo <- floodable min_of each.z_dyn;
			if cumulative_V <= 0.0 {
				fill_level <- lo;
			} else {
				float hi <- L;
				loop times: 40 {
					float mid <- 0.5 * (lo + hi);
					if pool_volume(mid) < cumulative_V { lo <- mid; } else { hi <- mid; }
				}
				fill_level <- 0.5 * (lo + hi);
			}
			// 4. fill connected cells to the pool level; perched cells keep h (trapped casier)
			ask (floodable where (each.spill_lvl <= fill_level)) parallel: true {
				h <- fill_level - z_dyn;
				wsl <- z_dyn + h;
			}
		}
	}

	// ====================================================================== bookkeeping (identical to V1)
	reflex bookkeeping {
		ask metric_cells parallel: true {
			float exc <- h - h0;
			if exc > flood_threshold {
				if arrival_h < 0.0 { arrival_h <- (current_date - starting_date) / 3600.0; }
				if h > h_peak { h_peak <- h; }
			}
		}
		list<cell> wet <- metric_cells where ((each.h - each.h0) > flood_threshold and !each.is_river);
		flooded_area_km2 <- length(wet) * cell_a / 1e6;
		flood_volume_mm3 <- (wet sum_of (each.h - each.h0)) * cell_a / 1e6;
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
				+ " | gate " + int(gate_Q_now) + " m3/s | pool " + (fill_level with_precision 2)
				+ " m | flooded " + (flooded_area_km2 with_precision 2) + " km2 | vol "
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
			+ " km2 (peak " + (peak_flooded_km2 with_precision 2) + " km2), delivered "
			+ (cumulative_V / 1e6 with_precision 1) + " Mm3 through the gates.";
		ask observation_point {
			write "Point " + pid + " arrival: " + (arrival_h < 0.0 ? "never" : string(arrival_h with_precision 1) + " h");
		}
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
	float spill_lvl <- 1e9;  // S[c]: lowest stage at which the cell connects to a breach
	bool  passable <- false;
	bool  is_river <- false;
	bool  is_lake  <- false;
	bool  is_dyke  <- false;
	bool  is_corridor <- false;
	float h_peak <- 0.0;
	float arrival_h <- -1.0;
	cell nE; cell nW; cell nN; cell nS;
	rgb terrain_color <- #gray;
	bool was_wet_color <- false;
}

// ==========================================================================
//  Vector species
// ==========================================================================
species dyke schedules: [] {
	string break_s; string date_s; string commune;
	bool will_break <- false;
	bool opened <- false;
	date breach_date;
	list<cell> my_cells;
	float breach_invert <- 0.0;   // m, weir invert (lowest protected-side ground + freeboard)
	float breach_width  <- 0.0;   // m, gate width (crest cells x cell size)

	// open the breach: cut a corridor, set the weir invert and width (the gate
	// that feeds the volume budget). Identical geometry to V1.
	action open_breach {
		opened <- true;
		list<cell> search_zone <- cell overlapping (shape + breach_search_radius);
		list<cell> ground <- search_zone where (!each.is_dyke and !each.is_river);
		float target <- empty(ground)
			? (my_cells min_of each.z_dyn) - 5.0
			: (ground min_of each.z_dyn) + breach_freeboard;
		target <- max(breach_floor_min, target);
		breach_invert <- target;
		breach_width  <- length(my_cells) * cell_dx * gate_width_fraction;
		list<cell> corridor <- (cell overlapping (shape + breach_cut_halfwidth)) where (!each.is_river);
		ask corridor {
			z_dyn <- min(z_dyn, target);
			is_corridor <- true;
		}
		spill_dirty <- true;
		write "BREACH at " + commune + " on " + current_date + " (invert " + (target with_precision 2)
			+ " m, gate width " + int(breach_width) + " m, " + length(corridor) + " cells cut)";
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
	int status <- 0;
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
experiment fastbathtub_v2 type: gui {
	parameter "DEM file" var: dem_name among: ["mnt-gz40.asc", "mnt-gz25.asc", "mnt-gz50.asc", "mnt-gz50-1926.asc"] category: "Terrain";
	parameter "Use observed discharge (WaterDischarge.csv)" var: use_discharge_csv category: "Forcing";
	parameter "Peak river stage (m)" var: peak_stage category: "Forcing";
	parameter "Gauge datum offset (m)" var: datum_offset min: -3.0 max: 3.0 category: "Forcing";
	parameter "Breach cut half-width (m)" var: breach_cut_halfwidth category: "Breaching";
	parameter "Weir coefficient - CALIBRATE to ABM" var: weir_coef min: 0.2 max: 3.0 category: "Gate budget";
	parameter "Gate width fraction" var: gate_width_fraction min: 0.1 max: 1.0 category: "Gate budget";
	parameter "Initial lake depth (m)" var: lake_initial_depth min: 0.0 max: 2.0 category: "Initial state";
	parameter "Flood threshold (m)" var: flood_threshold min: 0.01 max: 0.5 category: "Engine";
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
			chart "1926 flood event (V2 gate-volume)" type: series x_label: "hours since 20-07 00:00" {
				data "River stage (m)" value: river_stage color: #blue marker: false;
				data "Gate flow (100 m3/s)" value: gate_Q_now / 100.0 color: #green marker: false;
				data "Flooded area (km2)" value: flooded_area_km2 color: #red marker: false;
				data "Flood volume (10^6 m3)" value: flood_volume_mm3 color: #darkorange marker: false;
			}
		}
		monitor "Date" value: current_date;
		monitor "River stage (m)" value: river_stage with_precision 2;
		monitor "Breaches open" value: breaches_open;
		monitor "Gate flow (m3/s)" value: int(gate_Q_now);
		monitor "Pool level (m)" value: fill_level with_precision 2;
		monitor "Delivered volume (Mm3)" value: (cumulative_V / 1e6) with_precision 2;
		monitor "Flooded area (km2)" value: flooded_area_km2 with_precision 2;
		monitor "Peak flooded (km2)" value: peak_flooded_km2 with_precision 2;
		monitor "Buildings wet / flooded" value: "" + n_bldg_wet + " / " + n_bldg_flooded;
	}
}

experiment fastbathtub_v2_fast type: gui {
	parameter "DEM file" var: dem_name among: ["mnt-gz40.asc", "mnt-gz25.asc", "mnt-gz50.asc", "mnt-gz50-1926.asc"] category: "Terrain";
	parameter "Peak river stage (m)" var: peak_stage category: "Forcing";
	parameter "Weir coefficient - CALIBRATE to ABM" var: weir_coef min: 0.2 max: 3.0 category: "Gate budget";
	parameter "Gate width fraction" var: gate_width_fraction min: 0.1 max: 1.0 category: "Gate budget";

	output {
		display "Time series" type: 2d {
			chart "1926 flood event (V2 gate-volume)" type: series x_label: "hours since 20-07 00:00" {
				data "River stage (m)" value: river_stage color: #blue marker: false;
				data "Gate flow (100 m3/s)" value: gate_Q_now / 100.0 color: #green marker: false;
				data "Flooded area (km2)" value: flooded_area_km2 color: #red marker: false;
				data "Flood volume (10^6 m3)" value: flood_volume_mm3 color: #darkorange marker: false;
			}
		}
		monitor "Date" value: current_date;
		monitor "River stage (m)" value: river_stage with_precision 2;
		monitor "Breaches open" value: breaches_open;
		monitor "Gate flow (m3/s)" value: int(gate_Q_now);
		monitor "Pool level (m)" value: fill_level with_precision 2;
		monitor "Delivered volume (Mm3)" value: (cumulative_V / 1e6) with_precision 2;
		monitor "Flooded area (km2)" value: flooded_area_km2 with_precision 2;
		monitor "Peak flooded (km2)" value: peak_flooded_km2 with_precision 2;
	}
}
