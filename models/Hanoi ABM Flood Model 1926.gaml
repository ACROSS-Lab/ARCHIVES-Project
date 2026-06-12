/**
* Name: Hanoi ABM Flood Model 1926
* Description:
*   Clean-sheet AGENT-BASED model of the July-August 1926 Red River flood in Hanoi.
*
*   Multi-agent architecture (no code inherited from previous models in this project):
*     - grid `cell`     : ~69k terrain agents (mnt-gz40.asc archival MNT). Water moves between
*                         cell agents with DIFFUSION-WAVE Manning fluxes (the default engine of
*                         HEC-RAS 2D), two-phase and strictly mass-conservative. Dyke crests are
*                         only partially captured by the raster, which is why dykes are explicit
*                         barrier agents rather than DEM features.
*     - `dyke`          : 27 polyline agents (Dykes.shp), finite-state machine INTACT -> BREACHED
*                         (-> RESEALED). On its archival break date (BREAK/DATE attributes) a dyke
*                         turns its crest cells into bidirectional BROAD-CRESTED WEIR gates with a
*                         Villemonte submergence correction - the HEC-RAS levee-break pattern
*                         (weir structure + 2D friction routing). The DEM is never carved.
*     - `river`         : 1 polygon agent (RedRiver1925.shp). The observed hourly stage
*                         (RedRiverStage1926_hourly.csv, 1926-07-20 .. 08-15) is the upstream head
*                         of every weir gate; river cells are display-only water.
*     - `building`      : 4806 agents (Buildings1925.shp), DRY / WET / FLOODED exposure states.
*     - `lake`          : 1093 agents (Lakes1925.shp), initial standing water.
*
*   Strict dyke rule: water can ONLY pass from the river to the floodplain through a breached
*   segment. Cells crossed by a dyke line block all flow; the one-cell ring of floodplain in
*   contact with the river polygon is auto-sealed at init (guards against 50 m rasterization
*   gaps), so "no overtopping unless breached" holds by construction and the mass ledger is exact.
*
*   No reseal dates exist in the input data, so breaches stay open by default and casiers
*   re-equilibrate with the falling river (8.95 m at the end of the record). Optional reseal
*   switch (dates from archival research, NOT from the input data) reproduces trapped ponds.
*
*   Time: ONE CYCLE = ONE HOUR, aligned with the hourly stage record so simulated and observed
*   series can be compared row by row. Inside each cycle the engine runs in hydraulic sub-steps
*   (sub_dt_s, default 60 s) and exits early once flows die out, so quiet hours are nearly free.
*   Every face transfer is capped at cap_frac of the head difference per sub-step and by donor
*   storage -> unconditionally stable; halve sub_dt_s for final runs and check convergence.
*
* Inputs (all in ../includes/): mnt-gz40.asc (archival MNT; mnt-gz25.asc = fine mode), Dykes.shp,
*   RedRiver1925.shp, Lakes1925.shp, Buildings1925.shp, RedRiverStage1926_hourly.csv,
*   WaterDischarge.csv (daily event hydrograph, context/sanity only).
*   CAUTION: mnt-gz50.asc and mnt-gz10.asc were overwritten with a modern DEM (z -2.6..43 m) whose
*   datum does not match the 1926 gauge; their .aux.xml statistics describe the old content.
* Author: Nguyen Thanh Do, June 2026.
*/

model HanoiABMFlood1926

global {

	// ------------------------------------------------------------------ parameters
	// DEM: only mnt-gz50.asc / mnt-gz25.asc hold the archival MNT (z 1.9..16.1 m, same datum as the
	// 1926 gauge). mnt-gz50.asc and mnt-gz10.asc were overwritten with a modern DEM (z -2.6..43 m,
	// stale .aux.xml statistics) whose datum is inconsistent with the stage record - do not use them.
	string dem_name <- "mnt-gz40.asc";        // "mnt-gz25.asc" = fine mode (~2.6x more cells)
	string sim_start <- "1926-07-26 00:00:00";// skip quiet pre-breach days (stage record starts 07-20)
	float  sub_dt_s <- 60.0;                  // hydraulic sub-step (s) inside each 1-hour cycle; 30 for final runs, 120-300 for previews
	float  cap_frac <- 0.5;                   // fraction of a head difference equalised per sub-step (0.5 = pair-stable max)

	float  n_water <- 0.025;                  // Manning n: open water (river / lake cells)
	float  n_field <- 0.05;                   // Manning n: floodplain / fields
	float  n_urban <- 0.12;                   // Manning n: cells containing 1925 buildings

	float  weir_coef <- 1.5;                  // broad-crested weir coefficient (SI)
	float  breach_width_fraction <- 1.0;      // share of each breached crest cell acting as gate
	int    breach_hour <- 6;                  // hour of day breaches open (not in data)
	float  breach_invert_offset <- 0.0;       // m, added to per-breach computed invert
	float  datum_offset <- 0.0;               // m, gauge datum minus DEM datum

	float  lake_initial_depth <- 0.5;         // m of water seeded in lake cells
	float  lake_settle_hours <- 12.0;         // init-time settling of the seeded lakes (slab -> still ponds)
	float  wet_thr <- 0.1;                    // m, building WET threshold
	float  flood_thr <- 1.0;                  // m, building FLOODED threshold

	bool   enable_reseal <- false;            // reseal dates are NOT in the input data (archival)
	string reseal_spec <- "Gia-Quat:08-08,Ai-Mo:08-08,Lam-Giu:12-08";
	bool   auto_pause <- true;                // pause the GUI when the stage record ends

	// ------------------------------------------------------------------ input files
	file dykes_shp     <- shape_file("../includes/Dykes.shp");
	file river_shp     <- shape_file("../includes/RedRiver1925.shp");
	file lakes_shp     <- shape_file("../includes/Lakes1925.shp");
	file buildings_shp <- shape_file("../includes/Buildings1925.shp");
	geometry shape <- envelope(grid_file("../includes/" + dem_name));

	// ------------------------------------------------------------------ time
	date record_start  <- date("1926-07-20 00:00:00");   // first row of the stage CSV
	date starting_date <- date(sim_start);
	date end_date      <- date("1926-08-16 00:00:00");   // end of the hourly stage record
	float step <- 3600.0;                                // ONE CYCLE = ONE HOUR (matches the stage record)
	date next_tick <- date(sim_start) + 3600.0;          // next hourly bookkeeping tick
	bool sim_finished <- false;
	date first_break_time <- nil;

	// ------------------------------------------------------------------ forcing series
	list<float> stage_hourly <- [];           // index = hours since 1926-07-20 00:00
	list<float> discharge_daily <- [];        // index = days  since 1926-07-20 (context only)
	float stage_now <- 0.0;                   // gauge stage (m), interpolated
	float discharge_now <- 0.0;               // river discharge (m3/s), daily steps

	// ------------------------------------------------------------------ engine state
	float cell_s <- 50.0;                     // cell size (m), measured at init
	float cell_a <- 2500.0;                   // cell area (m2), measured at init
	float EPS_WET <- 1e-4;                    // m, minimum depth considered wet
	float EPS_HEAD <- 1e-3;                   // m, minimum head difference that drives flow
	list<cell> active_cells <- [];            // wet cells (and cells that ever got water) - metrics
	list<cell> pending_cells <- [];
	list<cell> flow_list <- [];               // cells with MOVING water (+ ring) - the engine only visits these
	list<cell> pending_flow <- [];
	list<cell> river_cells <- [];
	list<cell> breached_cells <- [];

	// ------------------------------------------------------------------ ledger & indicators
	float lake_init_V <- 0.0;                 // m3 seeded in lakes at init
	float river_in_V <- 0.0;                  // m3 delivered by weir gates (river -> land)
	float river_out_V <- 0.0;                 // m3 drained back through gates (land -> river)
	float weir_Q_now <- 0.0;                  // m3/s, net through all gates (+ = into land), hourly average
	float prev_net_weir_V <- 0.0;             // m3, gate net volume at the previous hourly tick
	float moved_V_sub <- 0.0;                 // m3 moved in the current sub-step (early-exit test)
	float stored_V <- 0.0;                    // m3 currently on the floodplain (incl. lakes)
	float mass_err_pct <- 0.0;                // ledger closure error, should stay ~0
	float flooded_km2 <- 0.0;                 // area with > wet_thr of water (excl. river)
	int   n_breached <- 0;
	int   n_bldg_wet <- 0;
	int   n_bldg_flooded <- 0;
	float peak_flooded_km2 <- 0.0;
	int   peak_bldg_flooded <- 0;

	init {
		// ---- terrain
		ask cell {
			z <- grid_value;
			neigh <- neighbors;
			nb_n <- length(neigh);
			out_v <- list_with(nb_n, 0.0);
		}
		cell one <- first(cell);
		cell_a <- one.shape.area;
		cell_s <- sqrt(cell_a);
		float z_mean <- cell mean_of each.z;
		if (z_mean > 9.0) {
			write "WARNING: mean ground = " + (z_mean with_precision 2) + " m -> this is the MODERN DEM vintage "
				+ "(mnt-gz50/gz10 were overwritten; their .aux.xml stats are stale). Use mnt-gz40.asc or "
				+ "mnt-gz25.asc, which hold the archival MNT on the 1926 gauge datum.";
		}

		// ---- vector agents
		create river from: river_shp;
		create lake from: lakes_shp;
		create building from: buildings_shp;
		create dyke from: dykes_shp with: [
			commune::string(read("Commune")),
			break_s::string(read("BREAK")),
			date_s::string(read("DATE"))
		];

		// ---- cell flags (order matters: river, dykes, seal ring, then derived)
		river_cells <- cell overlapping first(river).shape;
		ask river_cells { is_river <- true; }
		ask dyke {
			my_cells <- cell overlapping self;
			ask my_cells { is_dyke <- true; }
		}
		ask river_cells {                      // auto-seal: the floodplain ring touching the river
			loop nb over: neigh {
				if (!nb.is_river and !nb.is_dyke) { nb.is_sealed <- true; }
			}
		}
		ask cell { blocked <- is_river or is_dyke or is_sealed; }

		ask lake { ask cell overlapping self { is_lake <- true; } }
		ask building {
			my_cell <- first(cell overlapping location);
			ask cell overlapping self { is_urban <- true; }
		}
		ask cell {
			n_man <- (is_river or is_lake) ? n_water : (is_urban ? n_urban : n_field);
		}

		// ---- dykes: parse break dates, compute per-breach invert from nearby floodplain ground
		ask dyke {
			breakable <- break_s = "YES";
			if (breakable) {
				list<string> p <- date_s split_with "-";          // "28-07" = DD-MM
				break_time <- date("1926-" + p[1] + "-" + p[0] + " 00:00:00") + breach_hour * 3600.0;
				// invert = LOW land-side ground near the breach (15th percentile of the 60-250 m
				// ring), not the mean: the ring also contains dyke shoulders and berm mounds.
				list<cell> ring <- ((cell overlapping (shape + 250.0)) - (cell overlapping (shape + 60.0)))
					where (!each.is_dyke and !each.is_river and !each.is_sealed);
				if (empty(ring)) {
					invert_lvl <- (my_cells mean_of each.z) - 3.0 + breach_invert_offset;
				} else {
					list<float> zz <- (ring collect each.z) sort_by each;
					invert_lvl <- zz[int(0.15 * length(zz))] + breach_invert_offset;
				}
			}
			if (enable_reseal and breakable) {
				loop spec over: reseal_spec split_with "," {
					list<string> kv <- spec split_with ":";
					if (kv[0] = commune) {
						list<string> rp <- kv[1] split_with "-";
						reseal_time <- date("1926-" + rp[1] + "-" + rp[0] + " 00:00:00");
					}
				}
			}
		}

		// ---- initial water: lakes
		ask cell where (each.is_lake and !each.blocked) {
			water_h <- lake_initial_depth;
			in_active <- true;
		}
		active_cells <- cell where each.in_active;
		lake_init_V <- (active_cells sum_of each.water_h) * cell_a;
		// the flat slab seeded above is not an equilibrium: let it settle into the lake
		// depressions NOW so the run starts from still ponds and the metric baseline is clean
		int settle_iters <- int(lake_settle_hours * 3600.0 / 120.0);
		flow_list <- list(active_cells);
		ask flow_list { in_flow <- true; }
		loop times: settle_iters { do flow_step dts: 120.0; }
		// settled state = metric baseline; clear `stirred` so the ponds start ASLEEP and the
		// engine has zero work until the first weir push wakes cells up
		ask cell { water0 <- water_h; stirred <- false; }

		// ---- forcing series
		matrix ms <- matrix(csv_file("../includes/RedRiverStage1926_hourly.csv", ",", true));
		loop v over: ms column_at 1 { stage_hourly << float(v); }
		matrix md <- matrix(csv_file("../includes/WaterDischarge.csv", ",", true));
		discharge_daily <- list_with(40, 0.0);
		loop r from: 0 to: length(md column_at 0) - 1 {
			list<string> dp <- string(md[0, r]) split_with "/";   // "7/22/1926 0:00"
			int off <- (int(dp[0]) = 7) ? int(dp[1]) - 20 : 11 + int(dp[1]);
			if (off >= 0 and off < 40) { discharge_daily[off] <- float(md[1, r]); }
		}
		stage_now <- stage_hourly[0];

		// ---- colors
		float zmin <- cell min_of each.z;
		float zmax <- cell max_of each.z;
		ask cell {
			float t <- (z - zmin) / max(0.001, zmax - zmin);
			base_color <- is_dyke ? rgb(70, 65, 60) : rgb(int(110 + 120 * t), int(140 + 40 * t), int(80 + 40 * t));
			color <- base_color;
		}
		do refresh_visuals;

		// ---- init report
		write "=== Hanoi ABM Flood 1926 - initialisation ===";
		write "grid " + dem_name + ": " + length(cell) + " cells of " + int(cell_s) + " m, z "
			+ (cell min_of each.z) + " .. " + (cell max_of each.z) + " m";
		write "river cells: " + length(river_cells) + " | dyke cells: " + (cell count each.is_dyke)
			+ " | auto-sealed ring: " + (cell count each.is_sealed)
			+ " | lake cells: " + (cell count each.is_lake) + " | urban cells: " + (cell count each.is_urban);
		write "stage record: " + length(stage_hourly) + " h from " + record_start + ", simulation starts " + starting_date
			+ " (start " + stage_hourly[0] + " m, max " + max(stage_hourly) + " m)";
		ask dyke where each.breakable {
			write "breach scheduled: " + commune + " on " + break_time + ", invert " + (invert_lvl with_precision 2)
				+ " m, gate width " + int(length(my_cells) * cell_s * breach_width_fraction) + " m"
				+ (reseal_time != nil ? ", reseal " + string(reseal_time) : "");
		}
		first_break_time <- (dyke where each.breakable) min_of each.break_time;
		write "lakes settled for " + int(lake_settle_hours) + " h at init ("
			+ int(lake_settle_hours * 30) + " iterations); flood metrics are relative to this state.";
		write "strict dyke mode: floodplain receives water ONLY through breached gates"
			+ " -> everything stays dry until the first breach (" + first_break_time + ").";
	}

	// ================================================================== forcing
	reflex update_forcing {
		float hrs <- (current_date - record_start) / 3600.0;
		int ih <- int(hrs);
		if (ih >= length(stage_hourly) - 1) {
			stage_now <- last(stage_hourly);
		} else {
			stage_now <- stage_hourly[ih] + (hrs - ih) * (stage_hourly[ih + 1] - stage_hourly[ih]);
		}
		int id <- min(int(hrs / 24.0), length(discharge_daily) - 1);
		discharge_now <- discharge_daily[id];
		if (!sim_finished and current_date >= end_date) {
			sim_finished <- true;
			write "" + current_date + "  stage record finished - simulation complete.";
			if (auto_pause) { do pause; }
		}
	}

	// ================================================================== water engine
	// Phase A (parallel): every wet cell proposes outflow volumes to lower neighbours.
	// Diffusion-wave Manning flux on each face, capped per face at cap_frac of the head
	// difference and scaled so a donor never sends more than it stores.
	// (an action rather than a reflex so init can also run it to settle the seeded lakes)
	action flow_step (float dts) {
		ask flow_list parallel: true {
			has_out <- false;
			if (water_h > EPS_WET and !blocked) {
				float wsA <- z + water_h;
				float totV <- 0.0;
				loop i from: 0 to: nb_n - 1 {
					out_v[i] <- 0.0;
					cell nb <- neigh[i];
					if (!nb.blocked) {
						float dW <- wsA - (nb.z + nb.water_h);
						if (dW > EPS_HEAD) {
							float heff <- wsA - max(z, nb.z);
							if (heff > EPS_WET) {
								float q <- (heff ^ 1.6667) * sqrt(dW / cell_s) / (0.5 * (n_man + nb.n_man)) * cell_s;
								float v <- min(q * dts, dW * cell_a * cap_frac);
								if (v > 0.0) { out_v[i] <- v; totV <- totV + v; }
							}
						}
					}
				}
				if (totV > 0.0) {
					float availV <- water_h * cell_a;
					if (totV > availV) {
						float sc <- availV / totV;
						loop i from: 0 to: nb_n - 1 { out_v[i] <- out_v[i] * sc; }
					}
					has_out <- true;
				}
			}
		}

		// Phase B (sequential push): donors deliver their volumes and wake up dry receivers.
		ask flow_list {
			if (has_out) {
				float tot <- 0.0;
				loop i from: 0 to: nb_n - 1 {
					float v <- out_v[i];
					if (v > 0.0) {
						cell nb <- neigh[i];
						nb.water_h <- nb.water_h + v / cell_a;
						tot <- tot + v;
						nb.stirred <- true;
						if (!nb.in_active) { nb.in_active <- true; pending_cells << nb; }
						if (!nb.in_flow) { nb.in_flow <- true; pending_flow << nb; }
					}
				}
				water_h <- max(0.0, water_h - tot / cell_a);
				moved_V_sub <- moved_V_sub + tot;
				stirred <- true;
			}
		}
		if (!empty(pending_cells)) {
			active_cells <- active_cells + pending_cells;
			pending_cells <- [];
		}
		if (!empty(pending_flow)) {
			flow_list <- flow_list + pending_flow;
			pending_flow <- [];
		}
	}

	// one cycle = one hour: run the engine in sub-steps, stop early once the hour goes quiet
	reflex hydraulics {
		// hourly activity rebuild: only cells that moved water last hour (+ their ring of
		// unblocked neighbours) are scheduled; still ponds and dry land cost nothing
		ask flow_list { in_flow <- false; }
		list<cell> seeds <- active_cells where each.stirred;
		ask active_cells { stirred <- false; }
		flow_list <- [];
		loop c over: seeds {
			if (!c.in_flow) { c.in_flow <- true; flow_list << c; }
			loop nb over: c.neigh { if (!nb.in_flow and !nb.blocked) { nb.in_flow <- true; flow_list << nb; } }
		}
		int n_sub <- max(1, int(step / sub_dt_s));
		float dts <- step / n_sub;
		loop i from: 1 to: n_sub {
			moved_V_sub <- 0.0;
			do flow_step dts: dts;
			if (!empty(breached_cells)) { do weir_step dts: dts; }
			if (moved_V_sub < 0.5 * dts) { break; }      // < 0.5 m3/s moving anywhere -> equilibrium
		}
	}

	// Weir gates (sequential, few cells): bidirectional broad-crested weir between the river
	// stage and the floodplain side of each breached crest cell, Villemonte submergence.
	action weir_step (float dts) {
		float se <- stage_now + datum_offset;
		float bw <- cell_s * breach_width_fraction;
		loop bc over: breached_cells {
			list<cell> cands <- bc.rx_cells;
			if (!empty(cands)) {
				cell rc <- cands with_min_of (each.z + each.water_h);
				float wl <- rc.z + rc.water_h;
				if (se > wl + EPS_HEAD) {                            // river -> floodplain
					float hu <- se - bc.invert_z;
					if (hu > 0.0) {
						float hd <- max(0.0, wl - bc.invert_z);
						float f <- (hd >= hu) ? 0.0 : (1.0 - (hd / hu) ^ 1.5) ^ 0.385;
						// a gate may fill its receiver up to 90% of the way to the river level in
						// one step: the source is a Dirichlet boundary, so this cannot oscillate
						float v <- min(weir_coef * bw * (hu ^ 1.5) * f * dts, 0.9 * (se - wl) * cell_a);
						if (v > 0.0) {
							rc.water_h <- rc.water_h + v / cell_a;
							rc.stirred <- true;
							if (!rc.in_active) { rc.in_active <- true; active_cells << rc; }
							if (!rc.in_flow) { rc.in_flow <- true; flow_list << rc; }
							river_in_V <- river_in_V + v;
							moved_V_sub <- moved_V_sub + v;
						}
					}
				} else {                                             // floodplain -> river (drain-back)
					cell sc_ <- cands with_max_of (each.z + each.water_h);
					float wl2 <- sc_.z + sc_.water_h;
					if (wl2 > se + EPS_HEAD and sc_.water_h > EPS_WET) {
						float hu <- wl2 - bc.invert_z;
						if (hu > 0.0) {
							float hd <- max(0.0, se - bc.invert_z);
							float f <- (hd >= hu) ? 0.0 : (1.0 - (hd / hu) ^ 1.5) ^ 0.385;
							float v <- min(min(weir_coef * bw * (hu ^ 1.5) * f * dts, sc_.water_h * cell_a),
								0.9 * (wl2 - se) * cell_a);
							if (v > 0.0) {
								sc_.water_h <- sc_.water_h - v / cell_a;
								sc_.stirred <- true;
								if (!sc_.in_flow) { sc_.in_flow <- true; flow_list << sc_; }
								river_out_V <- river_out_V + v;
								moved_V_sub <- moved_V_sub + v;
							}
						}
					}
				}
			}
		}
	}

	// ================================================================== hourly bookkeeping
	reflex hourly_tick when: current_date >= next_tick {
		next_tick <- next_tick + 3600.0;

		// exposure
		ask building {
			depth_w <- my_cell = nil ? 0.0 : max(0.0, my_cell.water_h - my_cell.water0);
			status <- depth_w > flood_thr ? 2 : (depth_w > wet_thr ? 1 : 0);
			color <- status = 2 ? #red : (status = 1 ? #orange : rgb(120, 120, 125));
		}
		n_bldg_wet <- building count (each.status = 1);
		n_bldg_flooded <- building count (each.status = 2);

		// indicators & exact ledger (only weir gates connect river and floodplain)
		stored_V <- (active_cells sum_of each.water_h) * cell_a;
		float expected <- lake_init_V + river_in_V - river_out_V;
		mass_err_pct <- stored_V = 0.0 ? 0.0 : 100.0 * (stored_V - expected) / max(1.0, expected);
		flooded_km2 <- (active_cells count (each.water_h - each.water0 > wet_thr)) * cell_a / 1e6;
		n_breached <- length(breached_cells) > 0 ? dyke count (each.state_s = "BREACHED") : 0;
		peak_flooded_km2 <- max(peak_flooded_km2, flooded_km2);
		peak_bldg_flooded <- max(peak_bldg_flooded, n_bldg_flooded);
		weir_Q_now <- ((river_in_V - river_out_V) - prev_net_weir_V) / 3600.0;   // hourly mean gate flow
		prev_net_weir_V <- river_in_V - river_out_V;

		do refresh_visuals;

		if (current_date.hour mod 6 = 0) {
			write "" + current_date + " | stage " + (stage_now with_precision 2) + " m | breached " + n_breached
				+ " | flooded " + (flooded_km2 with_precision 2) + " km2 | gates " + int(weir_Q_now)
				+ " m3/s | bldg flooded " + n_bldg_flooded + " | mass err " + (mass_err_pct with_precision 3) + " %";
		}
	}

	action refresh_visuals {
		float se <- stage_now + datum_offset;
		ask river_cells { water_h <- max(0.0, se - z); }            // display-only river water
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
	reflex prune when: every(12 #cycle) {                 // every 12 simulated hours
		list<cell> dried <- active_cells where (each.water_h <= EPS_WET and !each.is_lake);
		if (!empty(dried)) {
			ask dried { in_active <- false; }
			active_cells <- active_cells - dried;
		}
	}
}

// ====================================================================== species
species dyke schedules: dyke where each.breakable {
	string commune;
	string break_s;
	string date_s;
	bool breakable <- false;
	string state_s <- "INTACT";
	date break_time <- nil;
	date reseal_time <- nil;
	float invert_lvl <- 0.0;
	list<cell> my_cells <- [];

	reflex check_break when: state_s = "INTACT" and breakable and current_date >= break_time {
		state_s <- "BREACHED";
		ask my_cells {
			is_breached <- true;
			invert_z <- myself.invert_lvl;
			// receivers within ~2 cells: a 40-50 m dyke rasterisation is often 2 cells thick,
			// so direct neighbours of a crest cell can all be blocked
			rx_cells <- (cell overlapping (shape + 1.6 * cell_s)) where (!each.blocked);
		}
		if (my_cells all_match empty(each.rx_cells)) {
			write "WARNING: breach at " + commune + " found no unblocked receiver cells!";
		}
		breached_cells <- breached_cells + my_cells;
		write "" + current_date + "  *** DYKE BREACH at " + commune + " (invert "
			+ (invert_lvl with_precision 2) + " m) ***";
	}

	reflex check_reseal when: state_s = "BREACHED" and reseal_time != nil and current_date >= reseal_time {
		state_s <- "RESEALED";
		ask my_cells { is_breached <- false; }
		breached_cells <- breached_cells - my_cells;
		write "" + current_date + "  dyke RESEALED at " + commune;
	}

	aspect default {
		draw shape color: state_s = "BREACHED" ? #red : (state_s = "RESEALED" ? #orange : #green) width: 3;
	}
}

species river schedules: [] {
	aspect default {
		draw shape color: rgb(70, 130, 180, 90);     // translucent: actual depth is on the cells below
		draw shape.contour color: #steelblue width: 2;
	}
}

species lake schedules: [] {
	aspect default { draw shape.contour color: rgb(80, 140, 200) width: 1; }
}

species building schedules: [] {
	cell my_cell <- nil;
	float depth_w <- 0.0;
	int status <- 0;                       // 0 dry, 1 wet, 2 flooded
	rgb color <- rgb(120, 120, 125);
	aspect default { draw shape color: color; }
}

grid cell file: grid_file("../includes/" + dem_name) neighbors: 4
	use_regular_agents: false use_individual_shapes: false use_neighbors_cache: true schedules: [] {
	float z <- 0.0;
	float water_h <- 0.0;
	bool is_river <- false;
	bool is_dyke <- false;
	bool is_sealed <- false;
	bool is_lake <- false;
	bool is_urban <- false;
	bool is_breached <- false;
	bool blocked <- false;
	bool in_active <- false;
	bool has_out <- false;
	float invert_z <- 0.0;
	float water0 <- 0.0;                   // depth right after init (lakes): flood metrics are relative to it
	bool stirred <- false;                 // moved/received water this hour -> stays scheduled next hour
	bool in_flow <- false;                 // currently in flow_list
	list<cell> rx_cells <- [];             // weir receivers, filled when this crest cell breaches
	float n_man <- 0.05;
	int nb_n <- 0;
	list<cell> neigh <- [];
	list<float> out_v <- [];
	rgb base_color <- #grey;
}

// ====================================================================== experiments
experiment flood_1926 type: gui {
	parameter "DEM file" var: dem_name among: ["mnt-gz40.asc", "mnt-gz25.asc"] category: "Terrain";
	parameter "Simulation start" var: sim_start among: ["1926-07-20 00:00:00", "1926-07-26 00:00:00"] category: "Time";
	parameter "Hydraulic sub-step (s)" var: sub_dt_s min: 15.0 max: 600.0 category: "Time";
	parameter "Transfer cap per sub-step" var: cap_frac min: 0.1 max: 0.5 category: "Time";
	parameter "n open water" var: n_water category: "Hydraulics";
	parameter "n floodplain" var: n_field category: "Hydraulics";
	parameter "n urban" var: n_urban category: "Hydraulics";
	parameter "Weir coefficient" var: weir_coef min: 0.8 max: 2.2 category: "Breach";
	parameter "Gate width fraction" var: breach_width_fraction min: 0.2 max: 1.0 category: "Breach";
	parameter "Breach hour of day" var: breach_hour min: 0 max: 23 category: "Breach";
	parameter "Invert offset (m)" var: breach_invert_offset min: -2.0 max: 2.0 category: "Breach";
	parameter "Enable reseal (archival dates)" var: enable_reseal category: "Breach";
	parameter "Gauge datum offset (m)" var: datum_offset min: -2.0 max: 2.0 category: "Hydraulics";
	parameter "Initial lake depth (m)" var: lake_initial_depth min: 0.0 max: 2.0 category: "Initial state";
	parameter "Lake settling at init (h)" var: lake_settle_hours min: 0.0 max: 48.0 category: "Initial state";
	parameter "Auto pause at end" var: auto_pause category: "Time";

	output {
		layout #split;
		display "Flood map" type: 2d antialias: false {
			grid cell;
			species river;
			species lake;
			species building;
			species dyke;
		}
		display "Forcing & extent" type: 2d {
			chart "River forcing" type: series size: {1.0, 0.5} position: {0.0, 0.0} {
				data "stage (m)" value: stage_now color: #blue marker: false;
				data "discharge (1000 m3/s)" value: discharge_now / 1000.0 color: #grey marker: false;
			}
			chart "Inundation" type: series size: {1.0, 0.5} position: {0.0, 0.5} {
				data "flooded area (km2)" value: flooded_km2 color: #navy marker: false;
				data "gate flow (100 m3/s)" value: weir_Q_now / 100.0 color: #red marker: false;
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
		monitor "Stage (m)" value: stage_now with_precision 2;
		monitor "Dykes breached" value: n_breached;
		monitor "First breach" value: string(first_break_time);
		monitor "Gate flow (m3/s)" value: int(weir_Q_now);
		monitor "Flooded area (km2)" value: flooded_km2 with_precision 2;
		monitor "Buildings wet / flooded" value: "" + n_bldg_wet + " / " + n_bldg_flooded;
		monitor "Stored volume (Mm3)" value: (stored_V / 1e6) with_precision 2;
		monitor "Mass error (%)" value: mass_err_pct with_precision 3;
	}
}

experiment sensitivity type: batch repeat: 1 keep_seed: true until: sim_finished {
	parameter "n floodplain" var: n_field among: [0.035, 0.05, 0.08];
	parameter "Weir coefficient" var: weir_coef among: [1.0, 1.5, 2.0];
	parameter "Invert offset (m)" var: breach_invert_offset among: [-0.5, 0.0, 0.5];
	parameter "Auto pause at end" var: auto_pause among: [false];
	//method exhaustive;

	reflex results {
		ask simulations {
			write "RESULT n_field=" + n_field + " Cw=" + weir_coef + " invert_off=" + breach_invert_offset
				+ " -> peak flooded " + (peak_flooded_km2 with_precision 2) + " km2, peak buildings flooded "
				+ peak_bldg_flooded;
		}
	}
}
