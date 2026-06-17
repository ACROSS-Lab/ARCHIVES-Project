/**
* Name: Dong Hoi ABM Flood Model
* Author: Nguyen Thanh Do (2026-06-16)
*
* A dyke-free, river-seeded AGENT-BASED model of the October 2020 Nhật Lệ / Đồng Hới
* coastal flood. It is the ABM analogue of "Dong Hoi FastBathtub Flood Model.gaml":
* it ports the REAL diffusion-wave engine of "Hanoi ABM Flood Model 1926 O.gaml" to
* Đồng Hới, exactly as the FastBathtub fork ported the fast level-pool surrogate.
*
* WHY THIS ENGINE (and how it differs from the two parents).
*   - The Hanoi ABM moves water between ~69k cell agents with DIFFUSION-WAVE Manning
*     fluxes (the default engine of HEC-RAS 2D): a two-phase, per-face-capped,
*     strictly mass-conservative scheme. The river is a fixed head that can ONLY
*     reach the floodplain through BROAD-CRESTED WEIR GATES opened by a dyke breach.
*   - Đồng Hới has NO dyke ring; it floods by the estuary spilling over its natural
*     banks. So the breach machinery (dyke FSM, weir gates, invert search, auto-seal
*     ring, reseal) is removed and replaced by a RIVER-BANK OVERBANK EXCHANGE: every
*     land cell touching the estuary is a bidirectional spill face whose sill is its
*     OWN bank ground z. Water spills in when the stage L exceeds the bank, drains
*     back when L falls, and is routed across the plain by the SAME diffusion-wave
*     flux as everywhere else - one flux law, no weir coefficient (there is no
*     structure to parameterise, only terrain).
*   - Unlike the FastBathtub fork, the flood front here advances at the PHYSICAL
*     diffusion-wave rate, so there is NO front_celerity knob: this model is the
*     reference the FastBathtub's front_celerity is calibrated against.
*
* MASS LEDGER. The domain edge is closed and the river is the only source/sink, so
*   stored_V (water on the land) == river_in_V - river_out_V at all times; mass_err_%
*   should stay ~0. There are no lakes here, so the baseline (water0) is dry land.
*
* TIME. ONE CYCLE = ONE HOUR, aligned with the hourly records (cycle k <-> CSV row k)
*   so simulated and observed series compare row by row. Inside each cycle the engine
*   runs in hydraulic sub-steps (sub_dt_s, default 60 s) and exits early once flows
*   die out, so quiet hours are nearly free. Every face transfer is capped at cap_frac
*   of the head difference and by donor storage -> unconditionally stable.
*
* FORCING. Default = the MEASURED hourly Đồng Hới gauge (DongHoiStage2020_hourly.csv),
*   used directly as the estuary head - faithful to the Hanoi ABM, which is driven by
*   its measured hourly stage. Set use_discharge_csv = true to drive the stage from
*   observed discharge (WaterDischarge_HamNinh.csv) through a Manning rating curve
*   instead - this reproduces the FastBathtub fork's default forcing exactly, for a
*   like-for-like comparison.
*
* DATUM. The gauge datum is NOT the SRTM datum, so the raw stage is NOT a DEM
*   elevation. datum_offset = (gauge zero) - (DEM zero) is THE calibration knob; the
*   engine uses the water surface  L = river_stage + datum_offset.  The starting
*   value -5.185 m is carried over from the calibrated FastBathtub fork (same DEM,
*   same gauge): L ~ 2.0 m at the start and ~6.8 m at the peak, flooding the low plain
*   (z < ~7 m) while sparing the dunes and hills.
*
* INPUTS (../includes/ and ../includes/dong-hoi/):
*   dong-hoi/dong-hoi_3857.asc      SRTM Nhật Lệ DEM reprojected to EPSG:3857 (metres),
*                                   30 m, z -4..33 m. The geographic dong-hoi.asc is in
*                                   degrees and must NOT be used (cell width ~0.0003).
*   dong-hoi/water_donghoi.shp      Nhật Lệ estuary polygon - the stage boundary.
*   dong-hoi/building_multipolygon.shp  buildings (clipped to the domain in init).
*   DongHoiStage2020_hourly.csv     measured hourly gauge (default forcing).
*   WaterDischarge_HamNinh.csv      observed hourly discharge (rating-curve forcing).
*/
model DongHoiABMFlood

global {

	// ------------------------------------------------------------------ parameters
	// DEM in EPSG:3857 metres (reprojected from the geographic dong-hoi.asc). Fixed:
	// the geographic .asc would make cell sizes/areas/fluxes meaningless.
	string dem_name     <- "dong-hoi_3857.asc";
	file dem_file       <- grid_file("../includes/dong-hoi/" + dem_name);
	file river_file     <- shape_file("../includes/dong-hoi/water_donghoi.shp");
	file buildings_file <- shape_file("../includes/dong-hoi/building_multipolygon.shp");
	file stage_file     <- csv_file("../includes/DongHoiStage2020_hourly.csv", ",", true);
	file discharge_file <- csv_file("../includes/WaterDischarge_HamNinh.csv", ",", true);
	geometry shape <- envelope(dem_file);

	// ------------------------------------------------------------------ engine
	float  sub_dt_s <- 60.0;                  // hydraulic sub-step (s) inside each 1-hour cycle; 30 for final runs, 120-300 for previews
	float  cap_frac <- 0.5;                   // fraction of a head difference equalised per sub-step (0.5 = pair-stable max)
	bool   idle_sleep <- true;                // cells with 2 zero-flux sub-steps sleep until pushed or next hour

	float  n_water <- 0.025;                  // Manning n: open water (estuary face)
	float  n_field <- 0.05;                   // Manning n: floodplain / fields
	float  n_urban <- 0.12;                   // Manning n: cells containing buildings

	// ------------------------------------------------------------------ forcing
	// Default: the MEASURED hourly gauge drives the stage (faithful to the Hanoi ABM).
	// use_discharge_csv = true -> observed discharge through a Manning rating curve
	// drives it instead (the FastBathtub fork's default, for a like-for-like compare).
	bool  use_discharge_csv <- false;
	float rating_exponent <- 0.6;             // Manning h ~ Q^(3/5)
	float base_stage <- 7.0;                  // m, stage at the lowest recorded discharge
	float peak_stage <- 12.2;                 // m, stage at the peak discharge (~ observed peak)
	float river_stage <- 7.0;                 // gauge stage (m), this hour
	float river_discharge <- 0.0;             // discharge (m3/s), this hour (chart / context)
	list<date>  q_dates  <- [];
	list<float> q_values <- [];
	float q_min <- 1.0; float q_max <- 2.0;
	list<float> stage_series <- [];           // measured hourly gauge (index = hours since record_start)
	int   n_stage <- 0;

	// datum_offset = (gauge zero) - (DEM/SRTM zero). THE calibration knob; the engine
	// uses the water surface L = river_stage + datum_offset. -5.185 carried over from
	// the calibrated FastBathtub fork (same DEM + gauge).
	float datum_offset <- 0;
	float water_level <- 0.0;                 // L = river_stage + datum_offset (DEM datum), refreshed hourly

	// ------------------------------------------------------------------ thresholds / control
	float wet_thr   <- 0.1;                   // m, building WET threshold AND flooded-area threshold (matches the ABM)
	float flood_thr <- 1.0;                   // m, building FLOODED threshold
	bool  auto_pause <- true;                 // pause the GUI when the record ends

	// ------------------------------------------------------------------ time
	date record_start  <- date("2020-10-07 00:00:00");   // first row of the hourly records
	date starting_date <- date("2020-10-07 00:00:00");
	date end_date      <- date("2020-10-09 00:00:00");   // 48 hourly rows processed (07-00:00 .. 08-23:00)
	float step <- 3600.0;                                // ONE CYCLE = ONE HOUR

	// ------------------------------------------------------------------ engine state
	float cell_s <- 30.0;                     // cell size (m), measured at init
	float cell_a <- 900.0;                    // cell area (m2), measured at init
	float inv_cell_a <- 1.0 / 900.0;          // 1/cell_a (multiply instead of divide in hot loops)
	float cap_vol <- 450.0;                   // cell_a * cap_frac, refreshed hourly
	float EPS_WET  <- 1e-4;                   // m, minimum depth considered wet
	float EPS_HEAD <- 1e-3;                   // m, minimum head difference that drives flow
	list<cell> active_cells  <- [];           // wet land cells (and cells that ever got water) - metrics
	list<cell> pending_cells <- [];
	list<cell> flow_list     <- [];           // land cells with MOVING water (+ ring) - the engine only visits these
	list<cell> pending_flow  <- [];
	list<cell> river_cells   <- [];           // estuary cells (fixed head, display-only depth, blocked from the land flux)
	list<cell> bank_cells    <- [];           // land cells touching the estuary - the overbank spill boundary

	// ------------------------------------------------------------------ ledger & indicators
	float river_in_V  <- 0.0;                 // m3 spilled over the banks (river -> land)
	float river_out_V <- 0.0;                 // m3 drained back over the banks (land -> river)
	float river_Q_now <- 0.0;                 // m3/s, net over all banks (+ = into land), hourly average
	float prev_net_river_V <- 0.0;            // m3, net bank volume at the previous hourly tick
	float moved_V_sub <- 0.0;                 // m3 moved in the current sub-step (early-exit test)
	float stored_V <- 0.0;                    // m3 currently on the floodplain
	float mass_err_pct <- 0.0;                // ledger closure error, should stay ~0
	float flooded_km2 <- 0.0;                 // land area with > wet_thr of water
	float flood_volume_mm3 <- 0.0;            // 10^6 m3 stored on the land
	int   n_bldg_wet <- 0;
	int   n_bldg_flooded <- 0;
	float peak_flooded_km2 <- 0.0;
	int   peak_bldg_flooded <- 0;
	bool  sim_finished <- false;

	init {
		write "=== Dong Hoi ABM Flood (dyke-free, river-bank overbank spill, diffusion-wave) ===";

		// ---- terrain
		ask cell {
			z <- grid_value;
			is_nodata <- z < -1000.0;           // -32767 sentinel (none in this DEM, but guard)
			neigh <- neighbors;
		}
		cell one <- first(cell);
		cell_a <- one.shape.area;
		cell_s <- sqrt(cell_a);
		inv_cell_a <- 1.0 / cell_a;
		cap_vol <- cell_a * cap_frac;
		list<cell> valid_cells <- cell where (!each.is_nodata);
		float z_mean <- valid_cells mean_of each.z;

		// ---- vector agents
		create river_poly from: river_file;
		create building from: buildings_file;

		// ---- cell flags. Estuary cells and no-data are BLOCKED from the land flux:
		// the river is a fixed head reached only by the overbank exchange (the
		// dyke-free analogue of the Hanoi ABM's "river is reached only through a gate").
		ask river_poly { ask cell overlapping self where (!each.is_nodata) { is_river <- true; } }
		ask building   { ask cell overlapping self { is_urban <- true; } }
		ask cell { blocked <- is_river or is_nodata; }
		river_cells <- cell where each.is_river;

		ask cell {
			n_man <- is_river ? n_water : (is_urban ? n_urban : n_field);
		}

		// hot-loop precomputation: open-neighbour lists (blocked never changes after this)
		// and per-face conductance K = w / (n_face * sqrt(dx)); the sub-step flux is then
		// just heff^(5/3) * sqrt(dW) * K
		float sq_dx <- sqrt(cell_s);
		ask cell where (!each.blocked) {
			neigh_open <- neigh where (!each.blocked);
			nb_open <- length(neigh_open);
			out_v <- list_with(nb_open, 0.0);
			K_face <- [];
			loop nb over: neigh_open {
				K_face << cell_s / (0.5 * (n_man + nb.n_man) * sq_dx);
			}
		}

		// ---- river banks: every land cell touching the estuary is an overbank spill
		// face with its OWN ground as the sill. bank_K folds in the number of estuary
		// faces (each a conduit of width cell_s) and the open-water/land Manning blend.
		ask cell where (!each.blocked) {
			bank_faces <- neigh count each.is_river;
			is_bank <- bank_faces > 0;
			if (is_bank) { bank_K <- bank_faces * cell_s / (0.5 * (n_water + n_man) * sq_dx); }
		}
		bank_cells <- cell where each.is_bank;

		// ---- buildings: bind to a cell, drop those outside the DEM domain
		ask building { my_cell <- first(cell overlapping location); }
		ask building where (each.my_cell = nil) { do die; }

		// ---- baseline (no lakes here, so the land starts dry: water0 = 0)
		ask cell { water0 <- water_h; stirred <- false; idle_hr <- false; idle_subs <- 0; }

		// ---- measured stage record (default forcing)
		matrix sm <- matrix(stage_file);
		loop r over: rows_list(sm) {
			string ds <- string(r[0]);
			if (length(ds) > 0 and ds != "datetime") { stage_series << float(r[1]); }
		}
		n_stage <- length(stage_series);
		if (n_stage = 0) {
			write "WARNING: DongHoiStage2020_hourly.csv empty.";
		} else {
			write "Measured gauge: " + n_stage + " h, " + (min(stage_series) with_precision 2)
				+ ".." + (max(stage_series) with_precision 2) + " m (gauge datum)";
		}

		// ---- observed discharge record (dates are M/D/YYYY H:MM)
		matrix qm <- matrix(discharge_file);
		loop r over: rows_list(qm) {
			string ds <- string(r[0]);
			float qv <- float(r[1]);
			if (length(ds) > 0 and qv > 0.0) {
				list<string> parts <- ds split_with " ";
				list<string> dmy <- first(parts) split_with "/";
				int hh <- 0;
				if (length(parts) > 1) { hh <- int(first(parts[1] split_with ":")); }
				q_dates  << date([int(dmy[2]), int(dmy[0]), int(dmy[1]), hh, 0, 0]);
				q_values << qv;
			}
		}
		if (empty(q_values)) {
			if (use_discharge_csv) { use_discharge_csv <- false; }
			write "WaterDischarge_HamNinh.csv empty -> measured gauge will drive the stage.";
		} else {
			q_min <- min(q_values); q_max <- max(q_values);
			write "Discharge: " + length(q_values) + " values, " + first(q_dates) + ".." + last(q_dates)
				+ ", " + (q_min with_precision 1) + "-" + (q_max with_precision 1) + " m3/s"
				+ (use_discharge_csv ? "  (DRIVING the stage via rating curve)" : "  (loaded for context; gauge is driving)");
		}

		// ---- colours
		float zmin <- valid_cells min_of each.z;
		float zmax <- valid_cells max_of each.z;
		ask cell where (!each.is_nodata) {
			float t <- (z - zmin) / max(0.001, zmax - zmin);
			base_color <- is_river ? rgb(70, 130, 180)
				: rgb(int(110 + 120 * t), int(140 + 40 * t), int(80 + 40 * t));
			color <- base_color;
		}

		// prime the stage and paint the estuary for the opening frame
		river_stage  <- stage_at(starting_date);
		water_level  <- river_stage + datum_offset;
		do refresh_visuals;

		// ---- init report
		write "grid " + dem_name + ": " + length(cell) + " cells of " + int(cell_s) + " m, z "
			+ (zmin with_precision 1) + " .. " + (zmax with_precision 1) + " m (mean " + (z_mean with_precision 2) + ")";
		write "river cells: " + length(river_cells) + " | bank (spill) cells: " + length(bank_cells)
			+ " | urban cells: " + (cell count each.is_urban) + " | buildings in domain: " + length(building);
		write "datum_offset " + datum_offset + " m -> water surface L = stage " + (river_stage with_precision 2)
			+ " m -> " + (water_level with_precision 2) + " m at start (peak gauge "
			+ ((n_stage > 0 ? max(stage_series) : peak_stage) with_precision 2) + " m).";
		write "dyke-free overbank mode: the plain floods wherever L exceeds the bank ground - no breach needed.";
		write "Init done. Simulation: " + starting_date + " -> " + end_date;
	}

	// ================================================================== forcing
	// observed discharge interpolated in time
	float discharge_at (date d) {
		if (empty(q_values)) { return 0.0; }
		if (d <= first(q_dates)) { return first(q_values); }
		if (d >= last(q_dates))  { return last(q_values); }
		loop i from: 1 to: length(q_dates) - 1 {
			if (d <= q_dates[i]) {
				float f <- (d - q_dates[i - 1]) / max(1.0, q_dates[i] - q_dates[i - 1]);
				return q_values[i - 1] + f * (q_values[i] - q_values[i - 1]);
			}
		}
		return last(q_values);
	}

	// measured-gauge stage (default forcing), linearly interpolated by hours since record_start
	float measured_stage_at (date d) {
		if (n_stage = 0) { return base_stage; }
		float hrs <- (d - record_start) / 3600.0;
		if (hrs <= 0.0) { return first(stage_series); }
		int ih <- int(hrs);
		if (ih >= n_stage - 1) { return last(stage_series); }
		return stage_series[ih] + (hrs - ih) * (stage_series[ih + 1] - stage_series[ih]);
	}

	// discharge -> Manning rating curve -> stage (when use_discharge_csv); else the gauge
	float stage_at (date d) {
		if (use_discharge_csv and !empty(q_values)) {
			float q <- discharge_at(d);
			float fq <- (q ^ rating_exponent - q_min ^ rating_exponent)
			          / max(1e-6, q_max ^ rating_exponent - q_min ^ rating_exponent);
			return base_stage + (peak_stage - base_stage) * min(1.0, max(0.0, fq));
		}
		return measured_stage_at(d);
	}

	reflex update_forcing {
		river_discharge <- discharge_at(current_date);
		river_stage     <- stage_at(current_date);
		water_level     <- river_stage + datum_offset;
	}

	// ================================================================== water engine
	// Phase A (parallel): every wet land cell proposes outflow volumes to lower
	// neighbours. Diffusion-wave Manning flux on each face, capped per face at cap_frac
	// of the head difference and scaled so a donor never sends more than it stores.
	// (an action, not a reflex, so the same code could be reused at init)
	action flow_step (float dts) {
		// flow_list holds only land cells (river/no-data are blocked and never enter
		// neigh_open), so no blocked test is needed in the hot loop
		ask (flow_list where (!each.idle_hr)) parallel: true {
			has_out <- false;
			float totV <- 0.0;
			float max_dW <- 0.0;
			// nb_open = 0 -> walled in, can never flow; the guard also matters because
			// GAML's `loop from: 0 to: -1` runs DOWNWARD (would execute at i = 0)
			if (nb_open > 0 and water_h > EPS_WET) {
				float wsA <- z + water_h;
				loop i from: 0 to: nb_open - 1 {
					out_v[i] <- 0.0;
					cell nb <- neigh_open[i];
					float nbz <- nb.z;
					float dW <- wsA - (nbz + nb.water_h);
					if (dW > EPS_HEAD) {
						if (dW > max_dW) { max_dW <- dW; }
						float heff <- wsA - max(z, nbz);
						if (heff > EPS_WET) {
							float v <- min((heff ^ 1.6667) * sqrt(dW) * K_face[i] * dts, dW * cap_vol);
							if (v > 0.0) { out_v[i] <- v; totV <- totV + v; }
						}
					}
				}
			}
			if (totV > 0.0) {
				// donor caps: its storage, AND it may not drop below its lowest neighbour's
				// level (prevents multi-face over-send -> ping-pong -> checkerboard films)
				float availV <- min(water_h * cell_a, max_dW * cap_vol);
				if (totV > availV) {
					float sc <- availV / totV;
					loop i from: 0 to: nb_open - 1 { out_v[i] <- out_v[i] * sc; }
				}
				has_out <- true;
				idle_subs <- 0;
			} else if (idle_sleep) {
				// two consecutive zero-flux sub-steps -> sleep until pushed or next hourly rebuild
				idle_subs <- idle_subs + 1;
				if (idle_subs >= 2) { idle_hr <- true; }
			}
		}

		// Phase B (sequential push, donors only): deliver volumes and wake the receivers.
		ask flow_list {
			if (has_out) {
				float tot <- 0.0;
				loop i from: 0 to: nb_open - 1 {
					float v <- out_v[i];
					if (v > 0.0) {
						cell nb <- neigh_open[i];
						nb.water_h <- nb.water_h + v * inv_cell_a;
						tot <- tot + v;
						nb.stirred <- true;
						nb.idle_hr <- false;
						nb.idle_subs <- 0;
						if (!nb.in_active) { nb.in_active <- true; pending_cells << nb; }
						if (!nb.in_flow)   { nb.in_flow   <- true; pending_flow  << nb; }
					}
				}
				water_h <- max(0.0, water_h - tot * inv_cell_a);
				moved_V_sub <- moved_V_sub + tot;
				stirred <- true;
			}
		}
		if (!empty(pending_cells)) { active_cells <<+ pending_cells; pending_cells <- []; }
		if (!empty(pending_flow))  { flow_list    <<+ pending_flow;  pending_flow  <- []; }
	}

	// River-bank overbank exchange (sequential, few cells): the dyke-free replacement
	// for the Hanoi ABM's weir gates. The estuary is a fixed head L; each bank land cell
	// exchanges with it through the SAME diffusion-wave Manning flux as overland flow,
	// with the bank's own ground as the sill (no structure -> no weir coefficient).
	//   river -> land : spill in   when L > land surface, depth over the sill = L - z
	//   land  -> river: drain back when land surface > L, depth available  = water_h
	// Volumes are booked to river_in_V / river_out_V so the mass ledger closes exactly.
	action bank_exchange (float dts) {
		float L <- water_level;
		loop bc over: bank_cells {
			float wl <- bc.z + bc.water_h;
			if (L > wl + EPS_HEAD) {                                   // spill in
				float heff <- L - bc.z;                               // flow depth over the bank ground
				if (heff > EPS_WET) {
					float dW <- L - wl;
					float v <- min(bc.bank_K * (heff ^ 1.6667) * sqrt(dW) * dts, dW * cap_vol);
					if (v > 0.0) {
						bc.water_h <- bc.water_h + v * inv_cell_a;
						river_in_V  <- river_in_V + v;
						moved_V_sub <- moved_V_sub + v;
						bc.stirred <- true; bc.idle_hr <- false; bc.idle_subs <- 0;
						if (!bc.in_active) { bc.in_active <- true; active_cells << bc; }
						if (!bc.in_flow)   { bc.in_flow   <- true; flow_list    << bc; }
					}
				}
			} else if (wl > L + EPS_HEAD and bc.water_h > EPS_WET) {   // drain back
				float dW <- wl - L;
				float v <- min(min(bc.bank_K * (bc.water_h ^ 1.6667) * sqrt(dW) * dts, bc.water_h * cell_a), dW * cap_vol);
				if (v > 0.0) {
					bc.water_h <- max(0.0, bc.water_h - v * inv_cell_a);
					river_out_V <- river_out_V + v;
					moved_V_sub <- moved_V_sub + v;
					bc.stirred <- true; bc.idle_hr <- false; bc.idle_subs <- 0;
					if (!bc.in_flow) { bc.in_flow <- true; flow_list << bc; }
				}
			}
		}
	}

	// one cycle = one hour: run the engine in sub-steps, stop early once the hour goes quiet
	reflex hydraulics {
		// hourly activity rebuild: only land cells that moved water last hour (+ their ring
		// of unblocked neighbours) are scheduled; still ponds and dry land cost nothing.
		// Bank cells are looped every sub-step in bank_exchange, so they re-arm themselves.
		ask flow_list { in_flow <- false; }
		list<cell> seeds <- active_cells where each.stirred;
		ask active_cells { stirred <- false; }
		flow_list <- [];
		loop c over: seeds {
			if (!c.in_flow) { c.in_flow <- true; c.idle_hr <- false; c.idle_subs <- 0; flow_list << c; }
			loop nb over: c.neigh_open {
				if (!nb.in_flow) { nb.in_flow <- true; nb.idle_hr <- false; nb.idle_subs <- 0; flow_list << nb; }
			}
		}
		cap_vol <- cell_a * cap_frac;                 // keeps a mid-run cap_frac change effective
		int n_sub <- max(1, int(step / sub_dt_s));
		float dts <- step / n_sub;
		loop i from: 1 to: n_sub {
			moved_V_sub <- 0.0;
			do bank_exchange dts: dts;                // inject / drain at the banks ...
			do flow_step dts: dts;                    // ... then spread overland
			if (moved_V_sub < 0.5 * dts) { break; }   // < 0.5 m3/s moving anywhere -> equilibrium
		}
	}

	// ================================================================== bookkeeping (per cycle = per hour)
	reflex bookkeeping {
		// exposure
		ask building {
			depth_w <- my_cell = nil ? 0.0 : max(0.0, my_cell.water_h - my_cell.water0);
			status <- depth_w > flood_thr ? 2 : (depth_w > wet_thr ? 1 : 0);
		}
		n_bldg_wet     <- building count (each.status = 1);
		n_bldg_flooded <- building count (each.status = 2);

		// indicators & exact ledger (the river bank is the only source/sink)
		stored_V <- (active_cells sum_of each.water_h) * cell_a;
		float expected <- river_in_V - river_out_V;
		mass_err_pct <- stored_V = 0.0 ? 0.0 : 100.0 * (stored_V - expected) / max(1.0, expected);
		list<cell> wet <- active_cells where (each.water_h - each.water0 > wet_thr);
		flooded_km2 <- length(wet) * cell_a / 1e6;
		flood_volume_mm3 <- (wet sum_of (each.water_h - each.water0)) * cell_a / 1e6;
		peak_flooded_km2 <- max(peak_flooded_km2, flooded_km2);
		peak_bldg_flooded <- max(peak_bldg_flooded, n_bldg_flooded);
		river_Q_now <- (expected - prev_net_river_V) / 3600.0;        // hourly mean net bank flow
		prev_net_river_V <- expected;

		do refresh_visuals;

		if (current_date.hour mod 3 = 0) {
			write "" + current_date + " | stage " + (river_stage with_precision 2) + " m (L "
				+ (water_level with_precision 2) + ") | flooded " + (flooded_km2 with_precision 2)
				+ " km2 | banks " + int(river_Q_now) + " m3/s | bldg flooded " + n_bldg_flooded
				+ " | mass err " + (mass_err_pct with_precision 3) + " %";
		}
	}

	action refresh_visuals {
		ask river_cells { water_h <- max(0.0, water_level - z); }     // display-only estuary depth
		ask river_cells + active_cells {
			if (water_h > 0.02) {
				float t <- min(1.0, water_h / 3.0);
				color <- rgb(int(70 - 40 * t), int(150 - 90 * t), int(210 + 30 * t));
			} else {
				color <- base_color;
			}
		}
	}

	// every simulated 12 h, drop long-dry cells from the active list
	reflex prune when: every(12 #cycle) {
		list<cell> dried <- active_cells where (each.water_h <= EPS_WET);
		if (!empty(dried)) {
			ask dried { in_active <- false; }
			active_cells <- active_cells - dried;
		}
	}

	reflex stop_simulation when: current_date >= end_date {
		sim_finished <- true;
		write "" + current_date + "  record finished. Flooded area " + (flooded_km2 with_precision 2)
			+ " km2 (peak " + (peak_flooded_km2 with_precision 2) + " km2), peak buildings flooded "
			+ peak_bldg_flooded + ", final mass err " + (mass_err_pct with_precision 3) + " %.";
		if (auto_pause) { do pause; }
	}
}

// ====================================================================== species
grid cell file: dem_file neighbors: 4
	use_regular_agents: false use_individual_shapes: false use_neighbors_cache: true schedules: [] {
	float z <- 0.0;                        // raw DEM elevation (m, EPSG:3857)
	float water_h <- 0.0;                  // water depth (m)
	float water0 <- 0.0;                   // depth right after init (dry land here) - metrics are relative to it
	bool is_river <- false;
	bool is_nodata <- false;
	bool is_urban <- false;
	bool is_bank <- false;                 // land cell touching the estuary (overbank spill face)
	bool blocked <- false;                 // river or no-data: excluded from the land flux
	bool in_active <- false;
	bool has_out <- false;
	bool stirred <- false;                 // moved/received water this hour -> stays scheduled next hour
	bool in_flow <- false;                 // currently in flow_list
	bool idle_hr <- false;                 // asleep for the rest of this hour (zero-flux cell)
	int  idle_subs <- 0;                   // consecutive zero-flux sub-steps
	list<cell> neigh <- [];
	list<cell> neigh_open <- [];           // neighbours that are not blocked (precomputed at init)
	int  nb_open <- 0;
	list<float> K_face <- [];              // per-open-face conductance w/(n_face*sqrt(dx))
	list<float> out_v <- [];
	float n_man <- 0.05;
	float bank_K <- 0.0;                    // overbank conductance to the estuary (folds in the number of river faces)
	int  bank_faces <- 0;
	rgb  base_color <- #grey;
}

species river_poly schedules: [] {
	aspect default {
		draw shape color: rgb(70, 130, 180, 90);     // translucent: real depth is painted on the cells below
		draw shape.contour color: #steelblue width: 2;
	}
}

species building schedules: [] {
	cell my_cell <- nil;
	float depth_w <- 0.0;
	int status <- 0;                       // 0 dry, 1 wet, 2 flooded
	aspect default { draw shape color: status = 2 ? #red : (status = 1 ? #orange : rgb(120, 120, 125)); }
}

// ====================================================================== experiments
experiment donghoi_abm type: gui {
	parameter "Use observed discharge (rating curve)" var: use_discharge_csv category: "Forcing";
	parameter "Rating curve exponent" var: rating_exponent category: "Forcing";
	parameter "Base river stage (m)" var: base_stage category: "Forcing";
	parameter "Peak river stage (m)" var: peak_stage category: "Forcing";
	parameter "Gauge datum offset (m) - CALIBRATE" var: datum_offset min: -15.0 max: 5.0 category: "Forcing";
	parameter "Hydraulic sub-step (s)" var: sub_dt_s min: 15.0 max: 600.0 category: "Engine";
	parameter "Transfer cap per sub-step" var: cap_frac min: 0.1 max: 0.5 category: "Engine";
	parameter "Idle-cell sleep (within hour)" var: idle_sleep category: "Engine";
	parameter "n open water" var: n_water category: "Hydraulics";
	parameter "n floodplain" var: n_field category: "Hydraulics";
	parameter "n urban" var: n_urban category: "Hydraulics";
	parameter "Wet / flooded-extent threshold (m)" var: wet_thr min: 0.02 max: 0.5 category: "Engine";
	parameter "Auto pause at end" var: auto_pause category: "Engine";

	output {
		layout #split;
		display "Flood map" type: 2d antialias: false {
			grid cell;
			species river_poly refresh: false;
			species building;
			graphics "info" {
				draw string(current_date) + "   stage: " + (river_stage with_precision 2)
					+ " m (L " + (water_level with_precision 2) + ")   flooded: " + (flooded_km2 with_precision 1) + " km2"
					at: {world.shape.width * 0.02, world.shape.height * 0.03}
					color: #white font: font("SansSerif", 14, #bold);
			}
		}
		display "Forcing & extent" type: 2d {
			chart "River forcing" type: series size: {1.0, 0.5} position: {0.0, 0.0} {
				data "stage (m)" value: river_stage color: #blue marker: false;
				data "water surface L (m)" value: water_level color: #steelblue marker: false;
				data "discharge (1000 m3/s)" value: river_discharge / 1000.0 color: #grey marker: false;
			}
			chart "Inundation" type: series size: {1.0, 0.5} position: {0.0, 0.5} {
				data "flooded area (km2)" value: flooded_km2 color: #navy marker: false;
				data "bank flow (100 m3/s)" value: river_Q_now / 100.0 color: #red marker: false;
			}
		}
		display "Impact & balance" type: 2d {
			chart "Buildings" type: series size: {1.0, 0.5} position: {0.0, 0.0} {
				data "wet (>0.1 m)" value: n_bldg_wet color: #orange marker: false;
				data "flooded (>1 m)" value: n_bldg_flooded color: #red marker: false;
			}
			chart "Mass balance" type: series size: {1.0, 0.5} position: {0.0, 0.5} {
				data "stored (Mm3)" value: stored_V / 1e6 color: #blue marker: false;
				data "ledger error (%)" value: mass_err_pct color: #black marker: false;
			}
		}
		monitor "Date" value: string(current_date);
		monitor "River discharge (m3/s)" value: river_discharge with_precision 1;
		monitor "Stage (m)" value: river_stage with_precision 2;
		monitor "Water surface L (m)" value: water_level with_precision 2;
		monitor "Bank flow (m3/s)" value: int(river_Q_now);
		monitor "Flooded area (km2)" value: flooded_km2 with_precision 2;
		monitor "Peak flooded (km2)" value: peak_flooded_km2 with_precision 2;
		monitor "Flood volume (Mm3)" value: flood_volume_mm3 with_precision 2;
		monitor "Buildings wet / flooded" value: "" + n_bldg_wet + " / " + n_bldg_flooded;
		monitor "Mass error (%)" value: mass_err_pct with_precision 3;
	}
}

// Display-free production / calibration experiment (numbers only, fastest).
experiment donghoi_abm_fast type: gui {
	parameter "Use observed discharge (rating curve)" var: use_discharge_csv category: "Forcing";
	parameter "Peak river stage (m)" var: peak_stage category: "Forcing";
	parameter "Gauge datum offset (m) - CALIBRATE" var: datum_offset min: -15.0 max: 5.0 category: "Forcing";
	parameter "Hydraulic sub-step (s)" var: sub_dt_s min: 15.0 max: 600.0 category: "Engine";
	parameter "n floodplain" var: n_field category: "Hydraulics";

	output {
		display "Time series" type: 2d {
			chart "Dong Hoi flood event" type: series x_label: "hours since 07-10 00:00" {
				data "River stage (m)" value: river_stage color: #blue marker: false;
				data "Water surface L (m)" value: water_level color: #steelblue marker: false;
				data "Flooded area (km2)" value: flooded_km2 color: #red marker: false;
				data "Flood volume (Mm3)" value: flood_volume_mm3 color: #darkorange marker: false;
			}
		}
		monitor "Date" value: current_date;
		monitor "Stage (m)" value: river_stage with_precision 2;
		monitor "Water surface L (m)" value: water_level with_precision 2;
		monitor "Flooded area (km2)" value: flooded_km2 with_precision 2;
		monitor "Peak flooded (km2)" value: peak_flooded_km2 with_precision 2;
		monitor "Flood volume (Mm3)" value: flood_volume_mm3 with_precision 2;
		monitor "Mass error (%)" value: mass_err_pct with_precision 3;
	}
}

experiment sensitivity type: batch repeat: 1 keep_seed: true until: sim_finished {
	parameter "n floodplain" var: n_field among: [0.035, 0.05, 0.08];
	parameter "Gauge datum offset (m)" var: datum_offset among: [-6.0, -5.185, -4.5];
	parameter "Hydraulic sub-step (s)" var: sub_dt_s among: [30.0, 60.0];
	parameter "Auto pause at end" var: auto_pause among: [false];

	reflex results {
		ask simulations {
			write "RESULT n_field=" + n_field + " datum=" + datum_offset + " dt=" + sub_dt_s
				+ " -> peak flooded " + (peak_flooded_km2 with_precision 2) + " km2, peak buildings flooded "
				+ peak_bldg_flooded + ", mass err " + (mass_err_pct with_precision 3) + " %";
		}
	}
}
