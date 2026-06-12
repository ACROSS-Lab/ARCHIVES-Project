/**
* Name: Hanoi FastFlood Model LISEM OpenLISEM Optimized V5
* Author: Thành Đô Nguyễn
*
* V5 of the OPTIMIZED variant of "Hanoi FastFlood Model LISEM OpenLISEM.gaml".
* On top of O1-O29 of Optimized / V2 / V3 / V4 (kept unchanged, summarized
* below), V5 adds four exact optimizations:
*
*   O30 [exact*] Incremental active-list extension. When the flood front
*       crossed into a previously-dry block, V4 rebuilt ALL flat active lists
*       from scratch (O(active cells) appends) - and during the expansion days
*       the front enters a new 8x8 block every few iterations, so the rebuild
*       could fire nearly every iteration at 2-3x the cost of a relaxation
*       pass. V5 keeps the from-scratch rebuild only in the hourly wetness
*       recompute (the only place blocks can DEactivate) and extends the lists
*       incrementally mid-hour: a per-block wet_processed flag records that a
*       wet block's neighbours are already active, so a dirty iteration only
*       scans the ~500 block flags and appends the members of newly activated
*       blocks. The active SET is identical (wet blocks + their 4-neighbours,
*       same as V4's mid-hour rebuild, which could also only grow the set).
*       (*) Only the ORDER of cells appended mid-hour differs from a
*       from-scratch rebuild, so on expansion hours the floating-point
*       convergence sum can differ in its last bits (deterministic run-to-run;
*       a convergence decision can only flip for sums within ~1 ulp of
*       relax_tolerance). The hourly rebuild keeps V4's exact bid order.
*   O31 [exact] List diet. River dh_abs is provably always 0.0 in forcing
*       mode (pass 2 never writes river cells since O21, track_dh is false
*       during the static refinement, and adding exact zeros never changes a
*       float sum), so the convergence sum runs over the LAND lists with a
*       bitwise-identical result. The now-redundant active_interior /
*       active_boundary lists are dropped entirely (active_all serves the
*       scans and the monitor): 5 lists instead of 7 in every rebuild.
*   O32 [exact] Pass-1 comparison reuse: q_dir tests dz_dir > 0 instead of
*       re-testing w_dir > 0 (w_dir = sqrt(dz_dir) is > 0 exactly when
*       dz_dir > 0): 4 comparisons per wet cell per iteration removed from
*       the hottest loop.
*   O33 [exact] Boundary-cell asks run sequentially: the domain-edge lists
*       hold at most ~700 cells (2x(219+142)-4 = 718 on this grid), where the
*       parallel fork/join dispatch costs more than the work it distributes.
*       Results are identical: these passes only write the cell's own
*       attributes (plus the idempotent true-only block flags).
*
* Inherited optimizations (full details in the Optimized / V2 / V3 / V4
* headers; all [exact] except O6):
*   O1  flux-split simplification (Manning factors cancel; w = sqrt(dz))
*   O2  relaxation passes 2+3 merged (gather pass commits h directly)
*   O3  hot-loop bodies inlined in the ask blocks
*   O4  minimal grid agents (no regular agents / shapes / neighbour cache)
*   O5  static-hazard display layer frozen (refresh: false)
*   O6  recolor only on wet-state change, terrain colors precomputed [approx]
*   O7  forced river depth precomputed once per hourly step
*   O8  relax_tolerance parameter; peak velocity cached for the monitor
*   O9  adaptive hourly budget on quiet hours
*   O10 cached water-surface level wsl = z_dyn + h
*   O11 dh_abs computed only on convergence-check iterations (track_dh)
*   O12 static vector layers drawn once in a frozen graphics layer
*   O13 nested-max velocity, color refresh inlined in bookkeeping
*   O14 interior/boundary split (nil checks only on the ~700 edge cells)
*   O15 dry-cell flux guard (had_q)
*   O16 block-level activity skipping (active = wet blocks + neighbours)
*   O17 stale-flux leak of O16 closed in the hourly wetness scan
*   O18 bookkeeping and hourly scans restricted to the active lists
*   O19 rebuild_active_lists appends in place (<<+)
*   O20 flux split divides once per cell (q_norm)
*   O21 river cells removed from pass 2 (h/wsl committed once per hour)
*   O22 schedules: [] on the grid and every species (no scheduler stepping)
*   O23 fast-sweeping inner min without list allocation
*   O24 block_size default 8 (active set hugs the flood front tighter)
*   O25 river-interior cells skipped in pass 1 (their flux has no reader)
*   O26 boundary asks skipped when their list is empty
*   O27 hourly wetness / stale-flux scan runs as a parallel ask
*   O28 rebuild marks blocks with an active flag (no remove_duplicates)
*   O29 b-shape fixed-point closure breaks early on convergence (init)
*   ... and the display-free experiment fastflood_1926_fast for production /
*       calibration runs.
*
* Everything else - method, calibration, data - is the original model:
* a GAML port of the FastFlood method of
*   van den Bout, B., Jetten, V.G., van Westen, C.J., Lombardo, L. (2023)
*   "A breakthrough in fast flood simulation", Env. Mod. & Soft. 168, 105787
* and of its reference implementation in LISEM-main/lisem/algorithms/raster/
* rasterfastflow.h (AS_FastFlood, AccuFluxDiffusive, SteadyStateCorrection2),
* applied to the 1926 Red River dyke-breach flood of Hanoi:
*
*  A) STATIC PIPELINE (paper sections 2.1-2.4, computed once at init):
*     fast-sweeping hydrological correction (appendix A, eq. 21-22), D4
*     steady-state flow accumulation (eq. 1), Manning inversion (eq. 2),
*     partial steady-state compensation (eq. 4, 5, 6, 12, 14), and the
*     adaptive pressure-driven refinement of section 2.4 (the final
*     AccuFluxDiffusiveCP pass of AS_FastFlood: 100 iterations, courant 0.1).
*
*  B) DYNAMIC QUASI-STEADY SOLVER (port of AccuFluxDiffusive):
*     artificial-velocity relaxation (a courant fraction of each wet cell's
*     water moves to downslope neighbours of the water surface z+h per
*     iteration, equilibrium-limited by 0.25 (wsl - wsl_nb)), with river
*     cells carrying a forced stage (the do_forced/HForced mechanism) from
*     the observed 1926 discharge record via a Manning rating curve, and
*     dyke breaches opening on the dates recorded in Dykes.shp. Each hourly
*     step relaxes to the steady state of the current forcing
*     (relax_to_convergence), as the method requires.
*
* Input data (includes/): mnt-gz50.asc (50 m DEM, dyke crests included),
* RedRiver1925.shp, Buildings1925.shp, Lakes1925.shp, 5_arrival_time.shp,
* Dykes.shp (BREAK/DATE/Commune), WaterDischarge.csv (observed daily
* discharge 22-07..06-08 1926, peak 30 000 m3/s on 30-07).
*
* Calibration (sources in includes/): observed 1926 stage hydrograph at Hanoi
* (Gourou, Le Tonkin fig. 9; CIA GS 66-5 fig. 5): peak 11.93 m end of July,
* 8-day rise 7.0 -> 11.93 m from ~22-07, recession ~8.3 m in early August
* (reproduced by the rating curve with exponent 0.6). Three breaches on the
* left bank (Ai-Mo / Gia-Quat / Lam-Giu communes); CIA dates them 30-07,
* Dykes.shp records 28/29-07 (kept as the data-driven choice).
*/
model HanoiFastFloodLISEMOpenLISEMOptimizedV5

global {

	// ------------------------------------------------------------------
	// Input data
	// ------------------------------------------------------------------
	file dem_file       <- grid_file("../includes/mnt-gz50.asc");
	file river_file     <- shape_file("../includes/RedRiver1925.shp");
	file buildings_file <- shape_file("../includes/Buildings1925.shp");
	file lakes_file     <- shape_file("../includes/Lakes1925.shp");
	file points_file    <- shape_file("../includes/5_arrival_time.shp");
	file dykes_file     <- shape_file("../includes/Dykes.shp");
	file discharge_file <- csv_file("../includes/WaterDischarge.csv", ",", true);
	geometry shape <- envelope(dem_file);

	// ------------------------------------------------------------------
	// Simulated time : the July-August 1926 Red River flood
	// ------------------------------------------------------------------
	date starting_date <- date("1926-07-20 00:00:00");
	date end_date      <- date("1926-08-08 00:00:00");
	float step <- 1 #h;

	// ------------------------------------------------------------------
	// FastFlood solver parameters (names follow rasterfastflow.h)
	// ------------------------------------------------------------------
	// fraction of stored water a cell exports per relaxation iteration
	// (the "artificial velocity" of AccuFluxDiffusive; decays to courant/8 in-cycle)
	float courant_fastflood <- 0.15 min: 0.01 max: 1.0;
	// base relaxation iterations of the quasi-steady solver per 1 h step
	int sub_iterations <- 60 min: 1;
	// The method assumes the flow field reaches the steady state of the current
	// forcing: the C++ uses iter = max(rows, cols) (~219 here) full-grid
	// iterations per solve, and the paper's Maas levee case is a single solve
	// run to convergence. With relax_to_convergence, each hourly step keeps
	// iterating past sub_iterations (up to max_relax_iterations) until the mean
	// water-depth change per iteration drops below relax_tolerance - so narrow
	// breach openings actually convey their steady-state discharge.
	bool  relax_to_convergence <- true;
	int   max_relax_iterations <- 250 min: 1;
	float relax_tolerance <- 0.0002; // m, mean |dh| per iteration over the grid
	// O9: adaptive budget - on quiet hours convergence checks start early
	bool  adaptive_relaxation <- true;
	int   min_check_iterations <- 10;  // first possible convergence check on quiet hours
	float quiet_stage_delta <- 0.01;   // m, stage change under which an hour counts as quiet
	// O16: block-level activity skipping
	bool block_skipping <- true;
	// O24: 8 hugs the flood front tighter than 16 (fewer dry cells in the
	// active set); skipping is exact for ANY block size (applied at init)
	int  block_size <- 8 min: 4;       // cells per block side
	// equilibrium limiter of the C++ code: q <= stab_factor * (wsl - wsl_neighbour)
	float stab_factor <- 0.25;
	float h_eps <- 0.001;           // m, cells under this depth do not route water
	float flood_threshold <- 0.05;  // m, depth considered "flooded" (arrival, area, color)

	// Manning surface roughness (s m^-1/3)
	float n_land     <- 0.06;
	float n_building <- 0.15;
	float n_lake     <- 0.035;
	float n_river    <- 0.03;

	// ------------------------------------------------------------------
	// River stage forcing.
	// Observed 1926 discharge (WaterDischarge.csv) interpolated in time and
	// converted to stage with a Manning-type rating curve mapped onto
	// [base_stage, peak_stage]; calibrated against the observed 1926
	// hydrograph at Hanoi (see header). Fallback: gaussian hydrograph.
	// ------------------------------------------------------------------
	bool  use_discharge_csv <- true;
	float rating_exponent <- 0.6;  // Manning h ~ Q^(3/5)
	float base_stage <- 7.0;    // m, stage at the lowest recorded discharge (23 ft, ~22-07)
	float peak_stage <- 11.93;  // m, observed 1926 peak stage at Hanoi (Gourou fig. 9)
	date  peak_date  <- date("1926-07-30 00:00:00");
	float sigma_rise_days <- 3.5; // gaussian width of the rising limb (fallback)
	float sigma_fall_days <- 3.0; // gaussian width of the falling limb (fallback)
	float river_stage <- base_stage;
	float river_discharge <- 0.0; // m3/s, current interpolated discharge
	list<date>  q_dates  <- [];
	list<float> q_values <- [];
	float q_min <- 1.0; float q_max <- 2.0;

	// ------------------------------------------------------------------
	// Dyke breaching (historical record: see the original model header;
	// three left-bank breaches, CIA dates 30-07, Dykes.shp 28/29-07).
	// Breach geometry: the embankment is 2-3 cells wide in the 50 m DEM
	// while the polyline only crosses the crest chain, so the invert is
	// taken from the lowest protected-side ground within
	// breach_search_radius, and a corridor of breach_cut_halfwidth is cut
	// through the full embankment width (scenario definition, OpenLISEM
	// FlowBarriers-style; not solver physics).
	// ------------------------------------------------------------------
	int   breach_hour <- 6;        // breaches recorded per day open at this hour
	float breach_floor_min <- 2.0; // m, breach invert never drops below this
	float breach_freeboard <- 0.2; // m, breach invert sits this much above the land side
	float breach_search_radius <- 300.0; // m, search radius for the invert ground level
	float breach_cut_halfwidth <- 60.0;  // m, half-width of the corridor cut through the embankment

	// ------------------------------------------------------------------
	// Static FastFlood hazard map (paper pipeline, computed at init)
	// ------------------------------------------------------------------
	bool  compute_static_hazard <- true;
	float design_rain_mmh <- 50.0;     // mm/h uniform design rainfall
	float design_duration_h <- 6.0;    // h, event duration for the compensation factor
	float dz_min_correction <- 0.001;  // m, minimum elevation increase per cell (eq. 21 delta)
	int   max_correction_sweeps <- 12; // fast-sweeping rounds (4 directional passes each)
	// step iv: pressure-driven refinement pass (flow3 of AS_FastFlood: 100 iterations, courant 0.1)
	int   static_refine_iterations <- 100;
	float static_refine_courant <- 0.1;

	// ------------------------------------------------------------------
	// Bookkeeping
	// ------------------------------------------------------------------
	int grid_cols; int grid_rows;
	float cell_dx <- 50.0;
	float z_min <- 0.0; float z_max <- 1.0;
	list<cell> river_cells;
	list<cell> interior_cells; // O14: all 4 neighbours exist - no nil checks needed
	list<cell> boundary_cells; // O14: domain-edge cells (free-draining virtual neighbours)
	list<cell> interior_land_cells; // O21: interior_cells minus river (pass-2 full-grid fallback)
	list<cell> boundary_land_cells; // O21: boundary_cells minus river
	list<cell> interior_p1_cells;   // O25: interior land + bank-river cells (pass-1 full-grid fallback)
	list<cell> boundary_p1_cells;   // O25: boundary land + bank-river cells
	// O16 state (O31: active_interior / active_boundary dropped - the land
	// lists carry the convergence sum and active_all serves the scans)
	int n_bx; int n_by;                  // blocks per row / column
	list<flow_block> blocks_all <- [];
	list<cell> active_interior_land <- []; // O21: active interior minus river (pass 2 + dh sum)
	list<cell> active_boundary_land <- []; // O21: active boundary minus river (pass 2 + dh sum)
	list<cell> active_interior_p1 <- [];   // O25: active interior land + bank river (pass 1)
	list<cell> active_boundary_p1 <- [];   // O25: active boundary land + bank river (pass 1)
	list<cell> active_all <- [];         // O18: interior + boundary of active blocks
	list<cell> all_cells <- [];          // O18: cached full population (fallback scans)
	bool blocks_dirty <- false;          // a block became wet since the last list rebuild
	// O18: true while the active lists are current (block skipping ran last cycle);
	// false after any full-grid cycle, so stale lists are never used for scans
	bool blocks_in_use <- false;
	float total_cell_count <- 1.0;
	float courant_now <- 0.15;
	// shared relaxation-iteration switches (FlowSource and do_forced of AccuFluxDiffusive)
	float relax_source_m <- 0.0;     // m of water added per cell per iteration
	bool  relax_force_river <- true; // river cells carry the forced stage (committed hourly, O21)
	float flooded_area_km2 <- 0.0;
	float flood_volume_mm3 <- 0.0;  // million m3
	int breaches_open <- 0;
	float static_flooded_km2 <- 0.0;
	float u_peak_max <- 0.0;        // m/s, domain maximum of u_peak (cached for the monitor)
	// O9/O11 relaxation state
	bool  track_dh <- false;            // pass 2 records dh_abs only when true (O11)
	float last_stage <- -999.0;         // stage of the previous hour
	bool  last_hour_converged <- false; // previous hour reached steady state
	bool  breach_this_hour <- false;    // a breach opened in the current cycle
	int   relax_iters_done <- 0;        // iterations used by the last hourly step

	init {
		write "=== Hanoi FastFlood (LISEM / OpenLISEM port, OPTIMIZED V5) ===";
		grid_cols <- 1 + max(cell collect each.grid_x);
		grid_rows <- 1 + max(cell collect each.grid_y);
		cell_dx <- first(cell).shape.width;
		write "Grid: " + grid_cols + " x " + grid_rows + " cells of " + (cell_dx with_precision 2) + " m";

		// --- elevation, neighbour references, inline pit removal -------------
		// (AccuFluxDiffusive raises pit cells to their lowest neighbour: the
		//  "pit" term computed inside every solver loop of the C++ code)
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
			z_fill <- z + pit;
			z_dyn <- z_fill;
			n_man <- n_land;
		}
		z_min <- cell min_of each.z;
		z_max <- cell max_of each.z;
		// O14: split the grid once into interior and boundary cells
		boundary_cells <- cell where (each.nE = nil or each.nW = nil or each.nN = nil or each.nS = nil);
		interior_cells <- cell where (each.nE != nil and each.nW != nil and each.nN != nil and each.nS != nil);
		all_cells <- list(cell); // O18: cached once for full-grid fallback scans
		total_cell_count <- float(length(cell));

		// O16: tile the grid into blocks and wire cells <-> blocks
		n_bx <- 1 + ((grid_cols - 1) div block_size);
		n_by <- 1 + ((grid_rows - 1) div block_size);
		create flow_block number: n_bx * n_by returns: created_blocks;
		blocks_all <- created_blocks;
		loop i from: 0 to: length(blocks_all) - 1 { blocks_all[i].bid <- i; }
		ask cell {
			my_block <- blocks_all[(grid_y div block_size) * n_bx + (grid_x div block_size)];
		}
		loop ce over: interior_cells { ce.my_block.members_interior << ce; }
		loop ce over: boundary_cells { ce.my_block.members_boundary << ce; }
		ask flow_block {
			int bx <- bid mod n_bx;
			int by <- bid div n_bx;
			if bx > 0        { nbrs << blocks_all[bid - 1]; }
			if bx < n_bx - 1 { nbrs << blocks_all[bid + 1]; }
			if by > 0        { nbrs << blocks_all[bid - n_bx]; }
			if by < n_by - 1 { nbrs << blocks_all[bid + n_bx]; }
		}
		write "Blocks: " + n_bx + " x " + n_by + " of " + block_size + "x" + block_size + " cells";

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

		// dykes first: their cells are excluded from river forcing so that the
		// forced water level never sits on a crest cell
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
		write "Dykes: " + length(dyke) + " segments, " + (dyke count (each.will_break)) + " breach during the event";

		ask lake     { ask cell overlapping self { is_lake <- true;  n_man <- n_lake; } }
		ask building { ask cell overlapping self { n_man <- n_building; } }
		ask river_poly {
			ask cell overlapping self {
				if !is_dyke { is_river <- true; n_man <- n_river; }
			}
		}
		river_cells <- cell where each.is_river;
		write "River cells: " + length(river_cells) + ", lake cells: " + (cell count (each.is_lake));

		// O21: land (non-river) pass-2 lists - full-grid and per-block versions
		interior_land_cells <- interior_cells where (!each.is_river);
		boundary_land_cells <- boundary_cells where (!each.is_river);
		// O25: bank classification - a river cell whose flux can actually be
		// consumed (>= 1 non-river neighbour, whose pass-2 gather reads it).
		// River-interior cells feed only river cells, which pass 2 skips
		// (O21), so their pass-1 work has no reader in forcing mode.
		ask river_cells {
			is_bank <- (nE != nil and !nE.is_river) or (nW != nil and !nW.is_river)
			        or (nN != nil and !nN.is_river) or (nS != nil and !nS.is_river);
		}
		interior_p1_cells <- interior_cells where (!each.is_river or each.is_bank);
		boundary_p1_cells <- boundary_cells where (!each.is_river or each.is_bank);
		ask flow_block {
			members_interior_land <- members_interior where (!each.is_river);
			members_boundary_land <- members_boundary where (!each.is_river);
			members_interior_p1 <- members_interior where (!each.is_river or each.is_bank);
			members_boundary_p1 <- members_boundary where (!each.is_river or each.is_bank);
		}
		write "Pass-1 cells: " + (length(interior_p1_cells) + length(boundary_p1_cells)) + " of "
			+ int(total_cell_count) + " (" + (river_cells count (!each.is_bank))
			+ " river-interior cells skipped, O25)";

		ask observation_point { my_cell <- first(cell overlapping self); }

		// --- observed discharge record (WaterDischarge.csv) -------------------
		// rows like "7/22/1926 0:00,10500" : date M/d/yyyy H:mm, discharge m3/s
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
			write "WaterDischarge.csv empty or unreadable -> falling back to the gaussian hydrograph";
		} else {
			q_min <- min(q_values);
			q_max <- max(q_values);
			write "Discharge record: " + length(q_values) + " values, " + first(q_dates) + " .. "
				+ last(q_dates) + ", " + q_min + " - " + q_max + " m3/s";
		}

		// --- static FastFlood hazard map (paper pipeline) ---------------------
		// run before the river is wetted: the refinement pass borrows the cell
		// water-depth field and resets it to zero afterwards
		if compute_static_hazard {
			do static_fastflood;
		}

		// --- precomputed terrain colors (O6) -----------------------------------
		ask cell {
			float shade <- (z - z_min) / max(0.001, z_max - z_min);
			terrain_color <- is_dyke
				? rgb(110, 80, 60)
				: rgb(70 + int(150 * shade), 80 + int(130 * shade), 60 + int(110 * shade));
			color <- terrain_color;
		}

		// --- initial state ----------------------------------------------------
		river_discharge <- discharge_at(starting_date);
		river_stage <- stage_at(starting_date);
		ask river_cells {
			h_forced <- max(0.0, river_stage - z_dyn);
			h <- h_forced;
		}
		ask cell { wsl <- z_dyn + h; do refresh_color; } // O10: prime the wsl cache
		write "Init done. Simulation: " + starting_date + " -> " + end_date;
	}

	// ======================================================================
	//  River forcing: observed discharge -> rating curve -> stage
	// ======================================================================
	// linear interpolation of the observed daily discharge record
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

	// Manning-type rating curve: maps the observed discharge range onto
	// [base_stage, peak_stage] with stage ~ Q^rating_exponent
	float stage_at (date d) {
		if use_discharge_csv and !empty(q_values) {
			float q <- discharge_at(d);
			float fq <- (q ^ rating_exponent - q_min ^ rating_exponent)
			          / max(1e-6, q_max ^ rating_exponent - q_min ^ rating_exponent);
			return base_stage + (peak_stage - base_stage) * min(1.0, max(0.0, fq));
		}
		// fallback: synthetic gaussian hydrograph
		float t_days <- (d - peak_date) / 86400.0;
		float sigma <- t_days < 0.0 ? sigma_rise_days : sigma_fall_days;
		return base_stage + (peak_stage - base_stage) * exp(-0.5 * (t_days / sigma) ^ 2);
	}

	reflex update_stage {
		river_discharge <- discharge_at(current_date);
		river_stage <- stage_at(current_date);
	}

	// ======================================================================
	//  Dyke breaching on the dates recorded in Dykes.shp
	// ======================================================================
	reflex open_breaches {
		int before <- breaches_open;
		ask dyke where (each.will_break and !each.opened and current_date >= each.breach_date) {
			do open_breach;
		}
		breaches_open <- dyke count (each.opened);
		if breaches_open > before { breach_this_hour <- true; } // O9: full budget this hour
	}

	// ======================================================================
	//  Quasi-steady FastFlood relaxation (port of AccuFluxDiffusive)
	//  One iteration = two passes (O2, O3):
	//   pass 1: each wet cell exports a courant fraction of its water to the
	//           downslope neighbours of the water surface, split by
	//           w = sqrt(dz) (the Manning factors cancel in the normalized
	//           split, O1), limited by stab_factor * dz. In forcing mode it
	//           runs on land + bank-river cells only (O25: river-interior
	//           flux has no reader); the full grid in source mode.
	//   pass 2: gather incoming fluxes and commit h directly. In forcing
	//           mode it runs on LAND cells only (O21: river h/wsl are
	//           committed once per hour); in source mode (static refinement)
	//           it runs on all cells with the rain source.
	// ======================================================================
	// O30: append one block's members to every active list and mark it
	// listed (active). Used by the mid-hour incremental extension; each
	// block is appended at most once per hour.
	action activate_block (flow_block b) {
		b.active <- true;
		// O19: '<<+' appends in place; O31: 5 lists (land, p1, all)
		active_interior_land <<+ b.members_interior_land; // O21
		active_boundary_land <<+ b.members_boundary_land; // O21
		active_interior_p1 <<+ b.members_interior_p1;     // O25
		active_boundary_p1 <<+ b.members_boundary_p1;     // O25
		active_all <<+ b.members_interior;
		active_all <<+ b.members_boundary;
	}

	// O16: flat cell lists of the currently active blocks (wet + 4-neighbours).
	// Full from-scratch rebuild - hourly only (the one place blocks can
	// DEactivate); collects in bid order, exactly as V4 did.
	action rebuild_active_lists {
		// O28: mark active blocks with a flag and collect in one pass over
		// the blocks (no remove_duplicates, no temporary concatenated lists).
		// O30: wet blocks are marked processed here (their neighbours are
		// activated below), so mid-hour extensions skip them.
		ask flow_block { active <- false; wet_processed <- false; }
		loop b over: blocks_all {
			if b.wet {
				b.active <- true;
				b.wet_processed <- true;
				loop nb over: b.nbrs { nb.active <- true; }
			}
		}
		active_interior_land <- [];
		active_boundary_land <- [];
		active_interior_p1 <- [];
		active_boundary_p1 <- [];
		active_all <- [];
		loop b over: blocks_all {
			if b.active {
				// O19: '<<+' appends in place ('+' would copy the accumulated list per block)
				active_interior_land <<+ b.members_interior_land; // O21
				active_boundary_land <<+ b.members_boundary_land; // O21
				active_interior_p1 <<+ b.members_interior_p1;     // O25
				active_boundary_p1 <<+ b.members_boundary_p1;     // O25
				active_all <<+ b.members_interior;
				active_all <<+ b.members_boundary;
			}
		}
	}

	// O30: mid-hour incremental extension - the front entered new block(s).
	// Within the hour the wet set only grows, so the active set only grows:
	// extending is equivalent to V4's from-scratch rebuild (same set,
	// appended at the end of the lists instead of in bid order). A newly wet
	// block was itself already active (water only reaches cells the gather
	// pass visits, i.e. cells of active blocks), so what is genuinely new
	// are its not-yet-active neighbours; the active flag guards against
	// double-appending either way.
	action extend_active_lists {
		loop b over: blocks_all {
			if (b.wet and !b.wet_processed) {
				b.wet_processed <- true;
				if !b.active { do activate_block(b); } // defensive; see invariant above
				loop nb over: b.nbrs {
					if !nb.active { do activate_block(nb); }
				}
			}
		}
	}

	action relax_iteration (float courant_val) {
		courant_now <- courant_val;
		// O16: skipping is exact only without a distributed source (rain wets
		// every cell, so the static refinement always runs on the full grid)
		bool use_blocks <- block_skipping and relax_source_m <= 0.0;
		// O25: pass 1 skips river-interior cells in forcing mode (their flux
		// has no reader, see init); source mode keeps the full grid because
		// pass 2 then visits every cell
		list<cell> p_int <- relax_force_river
			? (use_blocks ? active_interior_p1 : interior_p1_cells)
			: interior_cells;
		list<cell> p_bnd <- relax_force_river
			? (use_blocks ? active_boundary_p1 : boundary_p1_cells)
			: boundary_cells;
		// ---- pass 1, interior cells (O14: no nil checks) ----
		ask p_int parallel: true {
			if h > h_eps {
				// O10: wsl of every cell is maintained in pass 2, so a single
				// attribute read per neighbour suffices here
				float dz_e <- wsl - nE.wsl;
				float dz_w <- wsl - nW.wsl;
				float dz_n <- wsl - nN.wsl;
				float dz_s <- wsl - nS.wsl;
				float w_e <- dz_e > 0.0 ? sqrt(dz_e) : 0.0;
				float w_w <- dz_w > 0.0 ? sqrt(dz_w) : 0.0;
				float w_n <- dz_n > 0.0 ? sqrt(dz_n) : 0.0;
				float w_s <- dz_s > 0.0 ? sqrt(dz_s) : 0.0;
				float w_tot <- w_e + w_w + w_n + w_s;
				if w_tot > 0.0 {
					float q_norm <- courant_now * h / w_tot; // O20: one division per cell
					// equilibrium limiter of the C++ code (anti-oscillation):
					// no neighbour may receive more than what levels both surfaces
					// (O32: dz > 0 exactly when w > 0 - no second test on w)
					q_e <- dz_e > 0.0 ? min(q_norm * w_e, stab_factor * dz_e) : 0.0;
					q_w <- dz_w > 0.0 ? min(q_norm * w_w, stab_factor * dz_w) : 0.0;
					q_n <- dz_n > 0.0 ? min(q_norm * w_n, stab_factor * dz_n) : 0.0;
					q_s <- dz_s > 0.0 ? min(q_norm * w_s, stab_factor * dz_s) : 0.0;
					had_q <- true;
				} else if had_q { // O15
					q_e <- 0.0; q_w <- 0.0; q_n <- 0.0; q_s <- 0.0; had_q <- false;
				}
			} else if had_q { // O15: reset once when the cell dries out
				q_e <- 0.0; q_w <- 0.0; q_n <- 0.0; q_s <- 0.0; had_q <- false;
			}
		}
		// ---- pass 1, boundary cells (nil neighbours drain freely: dz = h;
		//      O33: sequential - the edge lists are too small to parallelize) ----
		if !empty(p_bnd) { // O26
			ask p_bnd {
				if h > h_eps {
					float dz_e <- nE = nil ? h : wsl - nE.wsl;
					float dz_w <- nW = nil ? h : wsl - nW.wsl;
					float dz_n <- nN = nil ? h : wsl - nN.wsl;
					float dz_s <- nS = nil ? h : wsl - nS.wsl;
					float w_e <- dz_e > 0.0 ? sqrt(dz_e) : 0.0;
					float w_w <- dz_w > 0.0 ? sqrt(dz_w) : 0.0;
					float w_n <- dz_n > 0.0 ? sqrt(dz_n) : 0.0;
					float w_s <- dz_s > 0.0 ? sqrt(dz_s) : 0.0;
					float w_tot <- w_e + w_w + w_n + w_s;
					if w_tot > 0.0 {
						float q_norm <- courant_now * h / w_tot; // O20
						// O32: dz tests reused
						q_e <- dz_e > 0.0 ? min(q_norm * w_e, stab_factor * dz_e) : 0.0;
						q_w <- dz_w > 0.0 ? min(q_norm * w_w, stab_factor * dz_w) : 0.0;
						q_n <- dz_n > 0.0 ? min(q_norm * w_n, stab_factor * dz_n) : 0.0;
						q_s <- dz_s > 0.0 ? min(q_norm * w_s, stab_factor * dz_s) : 0.0;
						had_q <- true;
					} else if had_q {
						q_e <- 0.0; q_w <- 0.0; q_n <- 0.0; q_s <- 0.0; had_q <- false;
					}
				} else if had_q {
					q_e <- 0.0; q_w <- 0.0; q_n <- 0.0; q_s <- 0.0; had_q <- false;
				}
			}
		}
		// ---- pass 2 (gather + commit) ----
		if relax_force_river {
			// O21: forcing mode - river h/wsl were committed at the start of the
			// hourly step (as the C++ does: H is set from HForced before the
			// iteration loop of AccuFluxDiffusive) and nothing changes them
			// within the hour, so pass 2 runs on land cells only, without the
			// per-cell river test and without the source term (relax_source_m
			// is always 0 in forcing mode).
			list<cell> g_int <- use_blocks ? active_interior_land : interior_land_cells;
			list<cell> g_bnd <- use_blocks ? active_boundary_land : boundary_land_cells;
			// interior land cells (O14: no nil checks)
			ask g_int parallel: true {
				float hv <- max(0.0, h - (q_e + q_w + q_n + q_s)
					+ nE.q_w + nW.q_e + nN.q_s + nS.q_n);
				if track_dh { dh_abs <- abs(hv - h); }          // O11
				h <- hv;
				wsl <- z_dyn + hv;                              // O10
				// O16: a newly wetted cell activates its block (true-only write)
				if hv > h_eps and !my_block.wet { my_block.wet <- true; blocks_dirty <- true; }
			}
			// boundary land cells (O33: sequential)
			if !empty(g_bnd) { // O26
				ask g_bnd {
					float hv <- max(0.0, h - (q_e + q_w + q_n + q_s)
						+ (nE = nil ? 0.0 : nE.q_w) + (nW = nil ? 0.0 : nW.q_e)
						+ (nN = nil ? 0.0 : nN.q_s) + (nS = nil ? 0.0 : nS.q_n));
					if track_dh { dh_abs <- abs(hv - h); }
					h <- hv;
					wsl <- z_dyn + hv;
					if hv > h_eps and !my_block.wet { my_block.wet <- true; blocks_dirty <- true; }
				}
			}
		} else {
			// source mode (static refinement): all cells gather + rain source.
			// use_blocks is false whenever a source is active, so p_int/p_bnd
			// are the full interior/boundary lists here.
			ask p_int parallel: true {
				float hv <- max(0.0, h - (q_e + q_w + q_n + q_s)
					+ nE.q_w + nW.q_e + nN.q_s + nS.q_n + relax_source_m);
				if track_dh { dh_abs <- abs(hv - h); }
				h <- hv;
				wsl <- z_dyn + hv;
				if hv > h_eps and !my_block.wet { my_block.wet <- true; blocks_dirty <- true; }
			}
			// (O33: sequential boundary pass)
			ask p_bnd {
				float hv <- max(0.0, h - (q_e + q_w + q_n + q_s)
					+ (nE = nil ? 0.0 : nE.q_w) + (nW = nil ? 0.0 : nW.q_e)
					+ (nN = nil ? 0.0 : nN.q_s) + (nS = nil ? 0.0 : nS.q_n)
					+ relax_source_m);
				if track_dh { dh_abs <- abs(hv - h); }
				h <- hv;
				wsl <- z_dyn + hv;
				if hv > h_eps and !my_block.wet { my_block.wet <- true; blocks_dirty <- true; }
			}
		}
	}

	// Each hour the water field is relaxed towards the steady state belonging
	// to the current stage / breach configuration; with relax_to_convergence
	// the iterations continue until that steady state is actually reached.
	reflex fastflood_relax {
		relax_force_river <- true;
		relax_source_m <- 0.0;
		// O7 + O21: the forced river depth is constant within the hourly step,
		// and pass 2 no longer touches river cells - commit their depth and
		// water surface here, once (the C++ ordering: H from HForced before
		// the iteration loop)
		ask river_cells parallel: true {
			h_forced <- max(0.0, river_stage - z_dyn);
			h <- h_forced;
			wsl <- z_dyn + h_forced;
		}
		// O16: recompute block wetness from scratch once per hour (lazy
		// deactivation of blocks that dried out), then build the active lists
		if block_skipping {
			ask flow_block { wet <- false; }
			// O18: wet cells are provably confined to the previous active set
			// (inactive cells cannot gain water), so when the lists are current
			// only they need scanning; full grid otherwise (first hour, or
			// block skipping was off last cycle).
			// O17: dry cells that still carry outgoing flux (dried in the very
			// last iteration of the previous hour) are cleared HERE, before
			// their block can be deactivated - otherwise an active neighbour
			// across the block boundary would keep gathering that stale flux
			// (phantom water). Restores the invariant: inactive => q = 0.
			// O27: parallel ask - the wet <- true writes are idempotent (same
			// benign-race pattern as the pass-2 activation writes) and the
			// q/had_q resets touch only the cell's own attributes.
			ask (blocks_in_use ? active_all : all_cells) parallel: true {
				if h > h_eps {
					if !my_block.wet { my_block.wet <- true; }
				} else if had_q {
					q_e <- 0.0; q_w <- 0.0; q_n <- 0.0; q_s <- 0.0;
					had_q <- false;
				}
			}
			// river cells that the rising stage wets THIS hour must also
			// activate their block (their block may have been inactive, so the
			// scan above did not visit them)
			ask river_cells where (each.h_forced > h_eps) { my_block.wet <- true; }
			do rebuild_active_lists;
			blocks_dirty <- false;
			blocks_in_use <- true;
		} else {
			blocks_in_use <- false; // O18: active lists go stale on full-grid cycles
		}
		// O9: on quiet hours (previous hour converged, no breach, stage nearly
		// unchanged) convergence checks may start almost immediately
		bool quiet <- adaptive_relaxation and relax_to_convergence and last_hour_converged
			and !breach_this_hour and abs(river_stage - last_stage) < quiet_stage_delta;
		int check_from <- quiet ? min_check_iterations : sub_iterations;
		int n_max <- relax_to_convergence ? max_relax_iterations : sub_iterations;
		bool converged <- false;
		loop it from: 1 to: n_max {
			// O16/O30: the front crossed into a new block - extend the active
			// lists incrementally (the from-scratch rebuild is hourly only)
			if block_skipping and blocks_dirty {
				do extend_active_lists;
				blocks_dirty <- false;
			}
			// in-cycle decay of the artificial velocity, as in the C++ code:
			// courant_here = courant (1-progress) + progress courant/8
			float progress <- it / n_max;
			// O11: record dh only on iterations whose convergence check runs
			track_dh <- relax_to_convergence and it >= check_from and (it mod 10 = 0);
			do relax_iteration(courant_fastflood * (1.0 - progress) + progress * courant_fastflood / 8.0);
			if track_dh {
				// O16: inactive cells have exactly dh = 0, so summing the active
				// lists over the full cell count equals the full-grid mean.
				// O31: river dh_abs is identically 0.0 in forcing mode (never
				// written since O21), so the land lists carry the whole sum -
				// bitwise the same value as summing land + river
				float dh_mean <- block_skipping
					? ((active_interior_land sum_of each.dh_abs) + (active_boundary_land sum_of each.dh_abs)) / total_cell_count
					: (cell mean_of each.dh_abs);
				if dh_mean < relax_tolerance {
					converged <- true;
					relax_iters_done <- it;
					break;
				}
			}
		}
		if !converged { relax_iters_done <- n_max; }
		last_hour_converged <- converged;
		last_stage <- river_stage;
		breach_this_hour <- false;
		track_dh <- false;
	}

	// ======================================================================
	//  Bookkeeping: arrival times, statistics, colors, stop condition
	// ======================================================================
	reflex bookkeeping {
		// O18: when block skipping ran this cycle (fastflood_relax precedes this
		// reflex, so blocks_in_use is current), every cell outside the active
		// lists has h <= h_eps: it cannot be flooded or change color, and its
		// arrival_h / h_peak / u_peak are frozen - scanning only the active
		// lists is exact. Full grid otherwise.
		list<cell> scan <- blocks_in_use ? active_all : all_cells;
		ask scan parallel: true {
			if h > flood_threshold and !is_river {
				if arrival_h < 0.0 { arrival_h <- (current_date - starting_date) / 3600.0; }
				if h > h_peak { h_peak <- h; }
			}
			// peak flow velocity (paper fig. 6D; C++ Vel map; OpenLISEM Vmax.map):
			// diffusive-wave Manning velocity on the steepest water-surface slope
			// (O13: nested max - no list allocation; O10: cached wsl)
			if h > h_eps {
				float dzmax <- max(max(
					max(nE = nil ? h : wsl - nE.wsl, nW = nil ? h : wsl - nW.wsl),
					max(nN = nil ? h : wsl - nN.wsl, nS = nil ? h : wsl - nS.wsl)), 0.0);
				float u <- (h ^ (2.0 / 3.0)) * sqrt(dzmax / cell_dx) / n_man;
				if u > u_peak { u_peak <- u; }
			}
			// O13: color refresh inlined (recolor only on wet-state change, O6)
			if h > flood_threshold {
				float fc <- min(1.0, h / 4.0);
				color <- rgb(int(150 * (1 - fc)), int(190 * (1 - fc) + 30), int(180 + 75 * fc));
				was_wet_color <- true;
			} else if was_wet_color {
				color <- terrain_color;
				was_wet_color <- false;
			}
		}
		ask observation_point where (each.arrival_h < 0.0) {
			if my_cell != nil and my_cell.h > flood_threshold {
				arrival_h <- (current_date - starting_date) / 3600.0;
				write "Observation point " + pid + " reached by the flood on " + current_date
					+ " (h = " + (my_cell.h with_precision 2) + " m)";
			}
		}
		list<cell> wet <- scan where (each.h > flood_threshold and !each.is_river);
		flooded_area_km2 <- length(wet) * cell_dx * cell_dx / 1e6;
		flood_volume_mm3 <- (wet sum_of each.h) * cell_dx * cell_dx / 1e6;
		// O18: u_peak is frozen outside the active set, so the domain maximum
		// is the previous maximum folded with the maximum over active cells
		u_peak_max <- blocks_in_use
			? (empty(scan) ? u_peak_max : max(u_peak_max, scan max_of each.u_peak))
			: (cell max_of each.u_peak);
	}

	reflex stop_simulation when: current_date >= end_date {
		write "End of event. Flooded area: " + (flooded_area_km2 with_precision 2) + " km2";
		ask observation_point {
			write "Point " + pid + " arrival: " + (arrival_h < 0.0 ? "never" : string(arrival_h with_precision 1) + " h after 20-07 00:00");
		}
		do pause;
	}

	// ======================================================================
	//  STATIC FASTFLOOD PIPELINE  (paper sections 2.1 - 2.4 + appendix A)
	// ======================================================================
	action static_fastflood {
		write "Static FastFlood: hydrological correction (fast sweeping)...";
		// 1. fast-sweeping hydro-correction: zcorr monotonically increasing away
		//    from the domain boundary, slope at least dz_min_correction per cell
		ask cell {
			bool edge <- nE = nil or nW = nil or nN = nil or nS = nil;
			zcorr <- edge ? z_fill : z_fill + 1e6;
		}
		// the four directional visiting orders of the Fast Sweeping Method (fig. 11)
		list<list<cell>> sweep_orders <- [
			cell sort_by (float(each.grid_y * grid_cols + each.grid_x)),
			cell sort_by (float(each.grid_y * grid_cols - each.grid_x)),
			cell sort_by (float(-(each.grid_y * grid_cols) + each.grid_x)),
			cell sort_by (float(-(each.grid_y * grid_cols) - each.grid_x))
		];
		int sweeps <- 0;
		bool changed <- true;
		loop while: (changed and sweeps < max_correction_sweeps) {
			changed <- false;
			sweeps <- sweeps + 1;
			loop ord over: sweep_orders {
				loop ce over: ord {
					// O23: nested binary min - no per-visit list allocation
					float zmin_nb <- min(
						min(ce.nE = nil ? ce.zcorr : ce.nE.zcorr, ce.nW = nil ? ce.zcorr : ce.nW.zcorr),
						min(ce.nN = nil ? ce.zcorr : ce.nN.zcorr, ce.nS = nil ? ce.zcorr : ce.nS.zcorr));
					float znew <- max(ce.z_fill, zmin_nb + dz_min_correction);
					if znew < ce.zcorr - 1e-6 {
						ce.zcorr <- znew;
						changed <- true;
					}
				}
			}
		}
		write "Static FastFlood: corrected in " + sweeps + " sweep rounds.";

		// 2. D4 drainage network on the corrected DEM (appendix A: the
		//    multi-directional network converted to steepest descent), and
		//    steady-state flow accumulation in a single elevation-ordered pass
		ask cell {
			cell best <- nil;
			float zb <- zcorr;
			loop nb over: [nE, nW, nN, nS] {
				if nb != nil and nb.zcorr < zb { zb <- nb.zcorr; best <- nb; }
			}
			downstream <- best;
			slope_ss <- max(0.001, (zcorr - (best = nil ? zcorr : best.zcorr)) / cell_dx);
			af1 <- 1.0;
		}
		list<cell> ordered <- cell sort_by (-each.zcorr);
		loop ce over: ordered {
			if ce.downstream != nil { ce.downstream.af1 <- ce.downstream.af1 + ce.af1; }
		}
		// AF(AF(1)) for the mean upstream travel distance (paper eq. 5)
		ask cell { af2 <- af1; }
		loop ce over: ordered {
			if ce.downstream != nil { ce.downstream.af2 <- ce.downstream.af2 + ce.af2; }
		}

		// 3. invert accumulation to steady-state flow height (paper eq. 2)
		float rain_ms <- design_rain_mmh / 1000.0 / 3600.0; // m/s
		ask cell {
			float q_ss <- af1 * rain_ms * cell_dx * cell_dx;  // m3/s through this cell
			h_static <- (q_ss * n_man / (cell_dx * sqrt(slope_ss))) ^ 0.6;
		}

		// 4. partial steady-state compensation (paper eq. 4, 5, 6, 12)
		list<float> u_vals <- (cell where (each.h_static > 0.001)) collect
			((each.h_static ^ (2.0 / 3.0)) * sqrt(each.slope_ss) / each.n_man);
		float u_mean <- max(0.05, empty(u_vals) ? 0.1 : mean(u_vals));
		float duration_s <- design_duration_h * 3600.0;
		ask cell {
			if af1 > 1.5 {
				float mean_s <- cell_dx * af2 / af1;             // <s>, eq. 5
				float b <- 1.0;
				loop times: 15 {                                  // closure of eq. 4 + eq. 6
					float smax_i <- cell_dx * (af1 ^ (1.0 / (1.0 + b)));
					float ratio <- min(0.95, max(0.36, mean_s / smax_i));
					float bnew <- min(10.0, max(-0.9, (2.0 * ratio - 1.0) / (1.0 - ratio)));
					// O29: fixed-point converged - stop early, as the C++
					// do_break does (b moves < 0.001 from here on)
					bool done <- abs(bnew - b) < 0.001;
					b <- bnew;
					if done { break; }
				}
				b_shape <- b;
				smax <- cell_dx * (af1 ^ (1.0 / (1.0 + b)));
				float s_ss <- min(1.0, duration_s * u_mean / max(cell_dx, smax));
				f_ss <- s_ss ^ (1.0 + b);                         // eq. 12
			} else {
				b_shape <- 0.0; smax <- cell_dx; f_ss <- 1.0;
			}
			float q_c <- f_ss * af1 * rain_ms * cell_dx * cell_dx; // eq. 14
			h_static <- (q_c * n_man / (cell_dx * sqrt(slope_ss))) ^ 0.6;
		}

		// 5. step iv (paper section 2.4): adaptive pressure-driven refinement.
		// Seed the diffusive relaxation solver with the compensated inverted
		// heights and the design rain as flow source, and let it spread the
		// water beyond the D4 network. This is the GAML equivalent of the
		// final pass of AS_FastFlood:
		//   flow2 = AccuFluxDiffusiveCP(DEM, Rain, flowinv, Zero, SS, 100, 0.1, ...)
		write "Static FastFlood: pressure-driven refinement (" + static_refine_iterations + " iterations)...";
		ask cell {
			h <- h_static;
			wsl <- z_dyn + h; // O10
			q_e <- 0.0; q_w <- 0.0; q_n <- 0.0; q_s <- 0.0; had_q <- false;
		}
		relax_force_river <- false;       // pluvial map: no forced river stage
		relax_source_m <- rain_ms;        // FlowSource analogue, m per iteration
		loop it from: 1 to: static_refine_iterations {
			float progress <- it / static_refine_iterations;
			do relax_iteration(static_refine_courant * (1.0 - progress) + progress * static_refine_courant / 8.0);
		}
		// harvest the refined map and restore the dynamic state
		ask cell {
			h_static <- h;
			h <- 0.0; wsl <- z_dyn; q_e <- 0.0; q_w <- 0.0; q_n <- 0.0; q_s <- 0.0; had_q <- false;
		}
		relax_source_m <- 0.0;
		relax_force_river <- true;

		static_flooded_km2 <- (cell count (each.h_static > flood_threshold)) * cell_dx * cell_dx / 1e6;
		write "Static FastFlood: hazard map done (R = " + design_rain_mmh + " mm/h, t = "
			+ design_duration_h + " h, mean u = " + (u_mean with_precision 2)
			+ " m/s, flooded " + (static_flooded_km2 with_precision 2) + " km2).";
	}
}

// ==========================================================================
//  Raster domain (O4: minimal grid agents - we use explicit neighbour refs;
//  O22: cells have no reflexes, so they are removed from the scheduler)
// ==========================================================================
grid cell file: dem_file neighbors: 4 use_regular_agents: false use_individual_shapes: false use_neighbors_cache: false schedules: [] {
	// terrain
	float z;       // raw DEM elevation (contains the dyke crests)
	float z_fill;  // pit-adjusted elevation (inline pit term of AccuFluxDiffusive)
	float z_dyn;   // elevation used by the dynamic solver (lowered when a dyke breaches)
	float n_man;   // Manning roughness
	bool is_river <- false;
	bool is_lake  <- false;
	bool is_dyke  <- false;
	bool is_bank  <- false;  // O25: river cell with >= 1 non-river neighbour

	// dynamic solver state
	float h <- 0.0;          // water depth (m)
	float wsl <- 0.0;        // m, cached water-surface level z_dyn + h (O10)
	float h_forced <- 0.0;   // m, forced river depth for the current hour (O7)
	bool had_q <- false;     // O15: cell carried outgoing flux in the last iteration
	flow_block my_block;     // O16: the activity block this cell belongs to
	float q_e <- 0.0; float q_w <- 0.0; float q_n <- 0.0; float q_s <- 0.0;
	float h_peak <- 0.0;
	float u_peak <- 0.0;     // m/s, peak diffusive-wave Manning velocity (Vmax)
	float dh_abs <- 0.0;     // m, |dh| of the last relaxation iteration (convergence check)
	float arrival_h <- -1.0; // hours after simulation start when first flooded

	// neighbour references (E/W = +x/-x, S/N = +y/-y in grid coordinates)
	cell nE; cell nW; cell nN; cell nS;

	// display state (O6)
	rgb terrain_color <- #gray;
	bool was_wet_color <- false;

	// static pipeline state
	float zcorr;
	cell downstream;
	float slope_ss <- 0.001;
	float af1 <- 1.0;     // flow accumulation AF(1) (contributing cells)
	float af2 <- 1.0;     // AF(AF(1)), for the mean travel distance
	float b_shape <- 0.0; // catchment shape parameter b
	float smax <- 50.0;   // maximum travel distance (m)
	float f_ss <- 1.0;    // partial steady-state compensation factor
	float h_static <- 0.0;

	// O6: recolor only when the wet state changes; terrain colors precomputed
	action refresh_color {
		if h > flood_threshold {
			float f <- min(1.0, h / 4.0);
			color <- rgb(int(150 * (1 - f)), int(190 * (1 - f) + 30), int(180 + 75 * f));
			was_wet_color <- true;
		} else if was_wet_color {
			color <- terrain_color;
			was_wet_color <- false;
		}
	}
}

// ==========================================================================
//  O16: activity blocks (block_size x block_size tiles of the grid)
// ==========================================================================
species flow_block schedules: [] {
	int bid;
	bool wet <- false;
	bool active <- false;        // O28: in the active lists (wet or wet-adjacent)
	bool wet_processed <- false; // O30: this wet block's neighbours are already active
	list<cell> members_interior <- [];
	list<cell> members_boundary <- [];
	list<cell> members_interior_land <- []; // O21: non-river members (pass 2)
	list<cell> members_boundary_land <- []; // O21
	list<cell> members_interior_p1 <- [];   // O25: land + bank-river members (pass 1)
	list<cell> members_boundary_p1 <- [];   // O25
	list<flow_block> nbrs <- [];
}

// ==========================================================================
//  Vector species (O22: no reflexes anywhere - all behavior is driven by
//  global asks, so none of these species needs scheduling)
// ==========================================================================
species dyke schedules: [] {
	string break_s; string date_s; string commune;
	bool will_break <- false;
	bool opened <- false;
	date breach_date;
	list<cell> my_cells;

	// open the breach: the invert is the lowest protected-side ground within
	// breach_search_radius (the floodplain behind the dyke, not the embankment
	// shoulder), and the cut is a corridor of breach_cut_halfwidth around the
	// polyline so it crosses the full width of the embankment in the DEM
	action open_breach {
		opened <- true;
		list<cell> search_zone <- cell overlapping (shape + breach_search_radius);
		list<cell> ground <- search_zone where (!each.is_dyke and !each.is_river);
		float target <- empty(ground)
			? (my_cells min_of each.z_fill) - 5.0
			: (ground min_of each.z_fill) + breach_freeboard;
		target <- max(breach_floor_min, target);
		list<cell> corridor <- (cell overlapping (shape + breach_cut_halfwidth)) where (!each.is_river);
		ask corridor { z_dyn <- min(z_dyn, target); wsl <- z_dyn + h; } // keep the wsl cache valid (O10)
		write "BREACH at " + commune + " on " + current_date + " (invert lowered to "
			+ (target with_precision 2) + " m, " + length(corridor) + " cells cut)";
	}

	aspect default {
		draw shape color: opened ? #red : (will_break ? #orange : #darkgreen) width: 3;
	}
}

species river_poly schedules: [] {
	aspect default { draw shape color: rgb(70, 130, 180, 120) border: #steelblue; }
}

species lake schedules: [] {
	aspect default { draw shape color: rgb(150, 200, 230, 150); }
}

species building schedules: [] {
	aspect default { draw shape color: rgb(90, 90, 90); }
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
experiment fastflood_1926 type: gui {
	parameter "Courant fraction (artificial velocity)" var: courant_fastflood category: "Solver";
	parameter "Relaxation iterations per hour" var: sub_iterations category: "Solver";
	parameter "Relax to convergence (steady state)" var: relax_to_convergence category: "Solver";
	parameter "Max relaxation iterations" var: max_relax_iterations category: "Solver";
	parameter "Convergence tolerance (m)" var: relax_tolerance category: "Solver";
	parameter "Adaptive budget on quiet hours" var: adaptive_relaxation category: "Solver";
	parameter "Block-level activity skipping" var: block_skipping category: "Solver";
	parameter "Block size (cells, applied at init)" var: block_size category: "Solver";
	parameter "Manning n land" var: n_land category: "Roughness";
	parameter "Manning n buildings" var: n_building category: "Roughness";
	parameter "Use observed discharge (WaterDischarge.csv)" var: use_discharge_csv category: "Forcing";
	parameter "Rating curve exponent" var: rating_exponent category: "Forcing";
	parameter "Base river stage (m)" var: base_stage category: "Forcing";
	parameter "Peak river stage (m)" var: peak_stage category: "Forcing";
	parameter "Breach hour of day" var: breach_hour category: "Breaching";
	parameter "Breach invert search radius (m)" var: breach_search_radius category: "Breaching";
	parameter "Breach cut half-width (m)" var: breach_cut_halfwidth category: "Breaching";
	parameter "Compute static FastFlood hazard map" var: compute_static_hazard category: "Static pipeline";
	parameter "Design rainfall (mm/h)" var: design_rain_mmh category: "Static pipeline";
	parameter "Design event duration (h)" var: design_duration_h category: "Static pipeline";
	parameter "Refinement iterations (step iv)" var: static_refine_iterations category: "Static pipeline";

	output {
		layout #split;

		display "Flood simulation" type: 2d background: #black {
			grid cell;
			// O12: the landscape never changes - draw it once and freeze it
			graphics "static landscape" refresh: false {
				loop rp over: river_poly { draw rp.shape color: rgb(70, 130, 180, 120) border: #steelblue; }
				loop lk over: lake { draw lk.shape color: rgb(150, 200, 230, 110); }
				loop b over: building { draw b.shape color: rgb(90, 90, 90); }
			}
			species dyke;
			species observation_point;
			graphics "info" {
				draw string(current_date) + "   stage: " + (river_stage with_precision 2)
					+ " m   flooded: " + (flooded_area_km2 with_precision 1) + " km2"
					at: {world.shape.width * 0.02, world.shape.height * 0.03}
					color: #white font: font("SansSerif", 16, #bold);
			}
		}

		// O5: the hazard layer is computed once at init and never changes -
		// freeze it so it is not rescanned and redrawn every frame
		display "Static FastFlood hazard" type: 2d background: #black {
			grid cell transparency: 0.55;
			graphics "static hazard" refresh: false {
				loop ce over: cell where (each.h_static > flood_threshold) {
					float f <- min(1.0, ce.h_static / 3.0);
					draw ce.shape color: rgb(int(255 * f), int(180 * (1 - f)), int(255 * (1 - f) * 0.4 + 120));
				}
			}
			species dyke;
			species river_poly transparency: 0.6;
		}

		display "Time series" type: 2d {
			chart "1926 flood event" type: series x_label: "hours since 20-07 00:00" {
				data "River stage (m)" value: river_stage color: #blue;
				data "Discharge (1000 m3/s)" value: river_discharge / 1000.0 color: #darkblue;
				data "Flooded area (km2)" value: flooded_area_km2 color: #red;
				data "Flood volume (10^6 m3)" value: flood_volume_mm3 color: #darkorange;
			}
		}

		monitor "Date" value: current_date;
		monitor "Discharge (m3/s)" value: river_discharge with_precision 0;
		monitor "River stage (m)" value: river_stage with_precision 2;
		monitor "Breaches open" value: breaches_open;
		monitor "Relax iterations (last hour)" value: relax_iters_done;
		// O31: active_all holds interior + boundary members of active blocks
		monitor "Active cells" value: block_skipping ? length(active_all) : int(total_cell_count);
		monitor "Flooded area (km2)" value: flooded_area_km2 with_precision 2;
		monitor "Peak velocity (m/s)" value: u_peak_max with_precision 2;
		monitor "Flood volume (10^6 m3)" value: flood_volume_mm3 with_precision 2;
		monitor "Arrival times (h)" value: observation_point collect (string(each.pid) + ": "
			+ (each.arrival_h < 0.0 ? "-" : string(each.arrival_h with_precision 1)));
	}
}

// Production / calibration experiment: identical simulation, but no grid
// raster display - redrawing the full cell grid every cycle is typically the
// largest fixed cost of a GUI run once the solver is tight. Only the
// time-series chart and the monitors are kept. Use fastflood_1926 when you
// want to watch the flood; use this one when you want the numbers fast.
experiment fastflood_1926_fast type: gui {
	parameter "Courant fraction (artificial velocity)" var: courant_fastflood category: "Solver";
	parameter "Relaxation iterations per hour" var: sub_iterations category: "Solver";
	parameter "Relax to convergence (steady state)" var: relax_to_convergence category: "Solver";
	parameter "Max relaxation iterations" var: max_relax_iterations category: "Solver";
	parameter "Convergence tolerance (m)" var: relax_tolerance category: "Solver";
	parameter "Adaptive budget on quiet hours" var: adaptive_relaxation category: "Solver";
	parameter "Block-level activity skipping" var: block_skipping category: "Solver";
	parameter "Block size (cells, applied at init)" var: block_size category: "Solver";
	parameter "Manning n land" var: n_land category: "Roughness";
	parameter "Manning n buildings" var: n_building category: "Roughness";
	parameter "Use observed discharge (WaterDischarge.csv)" var: use_discharge_csv category: "Forcing";
	parameter "Rating curve exponent" var: rating_exponent category: "Forcing";
	parameter "Base river stage (m)" var: base_stage category: "Forcing";
	parameter "Peak river stage (m)" var: peak_stage category: "Forcing";
	parameter "Breach hour of day" var: breach_hour category: "Breaching";
	parameter "Breach invert search radius (m)" var: breach_search_radius category: "Breaching";
	parameter "Breach cut half-width (m)" var: breach_cut_halfwidth category: "Breaching";
	parameter "Compute static FastFlood hazard map" var: compute_static_hazard category: "Static pipeline";
	parameter "Design rainfall (mm/h)" var: design_rain_mmh category: "Static pipeline";
	parameter "Design event duration (h)" var: design_duration_h category: "Static pipeline";
	parameter "Refinement iterations (step iv)" var: static_refine_iterations category: "Static pipeline";

	output {
		display "Time series" type: 2d {
			chart "1926 flood event" type: series x_label: "hours since 20-07 00:00" {
				data "River stage (m)" value: river_stage color: #blue;
				data "Discharge (1000 m3/s)" value: river_discharge / 1000.0 color: #darkblue;
				data "Flooded area (km2)" value: flooded_area_km2 color: #red;
				data "Flood volume (10^6 m3)" value: flood_volume_mm3 color: #darkorange;
			}
		}

		monitor "Date" value: current_date;
		monitor "Discharge (m3/s)" value: river_discharge with_precision 0;
		monitor "River stage (m)" value: river_stage with_precision 2;
		monitor "Breaches open" value: breaches_open;
		monitor "Relax iterations (last hour)" value: relax_iters_done;
		// O31: active_all holds interior + boundary members of active blocks
		monitor "Active cells" value: block_skipping ? length(active_all) : int(total_cell_count);
		monitor "Flooded area (km2)" value: flooded_area_km2 with_precision 2;
		monitor "Peak velocity (m/s)" value: u_peak_max with_precision 2;
		monitor "Flood volume (10^6 m3)" value: flood_volume_mm3 with_precision 2;
		monitor "Arrival times (h)" value: observation_point collect (string(each.pid) + ": "
			+ (each.arrival_h < 0.0 ? "-" : string(each.arrival_h with_precision 1)));
	}
}
