/**
* Name: Hanoi FastFlood Model LISEM OpenLISEM
* Author: Thành Đô Nguyễn
*
* A GAML port of the FastFlood method of:
*   van den Bout, B., Jetten, V.G., van Westen, C.J., Lombardo, L. (2023)
*   "A breakthrough in fast flood simulation", Environmental Modelling and Software 168, 105787.
* and of its reference implementation in the LISEM source tree
*   (LISEM-main/lisem/algorithms/raster/rasterfastflow.h, functions AS_FastFlood,
*    AccuFluxDiffusive, AS_Accumulate2DT, SteadyStateCorrection2),
* applied to the 1926 Red River dyke-breach flood of Hanoi.
*
* The model contains the two pillars of the FastFlood method:
*
*  A) STATIC PIPELINE (paper sections 2.1-2.4, computed once at init):
*     1. Hydrological correction of the DEM with a fast-sweeping scheme that makes the
*        elevation field monotonically increasing away from the domain boundary
*        (paper Appendix A, eq. 21-22). The result is depression-free by construction.
*     2. Steady-state flow accumulation over the resulting flow network (paper eq. 1).
*     3. Inversion of accumulated discharge to flow height with Manning's law:
*        h = (q n / (dx sqrt(S)))^(3/5)                                  (paper eq. 2)
*     4. Compensation for partial steady state from the catchment shape parameter b:
*        <s> = dx AF(AF(1))/AF(1)  (eq. 5),  smax = dx AF(1)^(1/(1+b))   (eq. 6)
*        f_ss = CDF(s_ss) = s_ss^(1+b)                                   (eq. 12)
*     5. Adaptive pressure-driven refinement (paper section 2.4, step iv): the
*        compensated inverted heights seed the iterative diffusive solver (with the
*        design rain as flow source), which spreads them beyond the D4 network --
*        the GAML equivalent of the final AccuFluxDiffusiveCP pass ("flow3", 100
*        iterations at courant 0.1) of AS_FastFlood.
*
*  B) DYNAMIC QUASI-STEADY SOLVER (port of AccuFluxDiffusive in rasterfastflow.h):
*     an iterative scheme with artificial velocity: every relaxation iteration each wet
*     cell exports a fraction (courant) of its stored water to its downslope neighbours
*     (slope of the water surface z+h, as in the adaptive pressure-driven refinement of
*     paper section 2.4). Outflow is split between the 4 neighbours proportionally to the
*     diffusive-wave Manning velocity u = h^(2/3) sqrt(S_ws) / n (paper eq. 17), and is
*     limited near equilibrium by q <= 0.25 (wsl - wsl_neighbour), exactly like the
*     stabilisation step of the C++ code. River cells carry a forced water level
*     (the do_forced/HForced mechanism of AccuFluxDiffusive), driven by a stage
*     hydrograph of the July-August 1926 Red River flood. Dykes recorded in Dykes.shp
*     with BREAK = YES open on their DATE (28-07 / 29-07), lowering the crest cells and
*     letting the river pour into the protected plain.
*
* Input data (includes/):
*   mnt-gz50.asc        50 m DEM (dyke crests are present in the elevation data)
*   RedRiver1925.shp    Red River water surface polygon
*   Buildings1925.shp   buildings of 1925 Hanoi  (used to raise Manning roughness)
*   Lakes1925.shp       lakes/ponds of 1925      (used to lower Manning roughness)
*   5_arrival_time.shp  5 observation points where flood arrival time is recorded
*   Dykes.shp           dyke segments, attributes BREAK (YES/NO), DATE (dd-MM), Commune
*   WaterDischarge.csv  observed daily Red River discharge (m3/s) 22-07 .. 06-08 1926,
*                       peak 30 000 m3/s on 30-07; converted to stage with a
*                       Manning-type rating curve (stage forcing of the river cells)
*/
model HanoiFastFloodLISEMOpenLISEM

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
	// relaxation iterations of the quasi-steady solver per 1 h simulation step
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
	// equilibrium limiter of the C++ code: q <= stab_factor * (wsl - wsl_neighbour)
	float stab_factor <- 0.25;
	float h_eps <- 0.001;           // m, cells under this depth do not route water
	float flood_threshold <- 0.05;  // m, depth considered "flooded" (arrival, area)

	// Manning surface roughness (s m^-1/3)
	float n_land     <- 0.06;
	float n_building <- 0.15;
	float n_lake     <- 0.035;
	float n_river    <- 0.03;

	// ------------------------------------------------------------------
	// River stage forcing.
	// Primary mode: the observed 1926 discharge record (WaterDischarge.csv,
	// daily values, peak 30 000 m3/s on 30-07) interpolated in time and
	// converted to stage with a Manning-type rating curve that maps the
	// observed discharge range onto [base_stage, peak_stage]:
	//   stage(Q) = base + (peak - base) (Q^e - Qmin^e)/(Qmax^e - Qmin^e)
	// (e = 0.6, the Manning depth-discharge exponent).
	// Calibration against the observed 1926 stage hydrograph at Hanoi
	// (Gourou, Le Tonkin, fig. 9; reproduced in CIA GS 66-5 "Levees in the
	// Red River Delta", fig. 5):
	//   - peak 11.93 m at end of July, then the highest level ever recorded;
	//   - the river rose 23 ft -> 39 ft (7.0 m -> 11.93 m) in the 8 days
	//     before the peak, so ~7.0 m around 22-07 = the first CSV record;
	//   - recession ~8.3 m in the first days of August: the rating curve
	//     with e = 0.6 gives 8.34 m for the last CSV record (15 000 m3/s
	//     on 06-08), matching the observed limb.
	// Fallback mode (use_discharge_csv = false): gaussian stage hydrograph.
	// ------------------------------------------------------------------
	bool  use_discharge_csv <- true;
	float rating_exponent <- 0.6;  // Manning h ~ Q^(3/5)
	float base_stage <- 7.0;    // m, stage at the lowest recorded discharge (23 ft, ~22-07)
	float peak_stage <- 15.00;  // m, observed 1926 peak stage at Hanoi (Gourou fig. 9)
	date  peak_date  <- date("1926-07-30 00:00:00");
	float sigma_rise_days <- 3.5; // gaussian width of the rising limb (fallback)
	float sigma_fall_days <- 3.0; // gaussian width of the falling limb (fallback)
	float river_stage <- base_stage;
	float river_discharge <- 0.0; // m3/s, current interpolated discharge
	list<date>  q_dates  <- [];
	list<float> q_values <- [];
	float q_min <- 1.0; float q_max <- 2.0;

	// ------------------------------------------------------------------
	// Dyke breaching.
	// Historical record (CIA GS 66-5 and CIA/BI GB 66-20, quoting French
	// sources): three breaches occurred in the levees not far from Hanoi,
	// on the LEFT bank (the Gia Lam side - matching the Ai-Mo, Gia-Quat and
	// Lam-Giu communes of Dykes.shp); two were closed on 8 and 12 August,
	// one remained open; over 250 000 acres were flooded. The CIA documents
	// date the breaches to 30 July; Dykes.shp records 28-07/29-07, which is
	// kept here as the data-driven choice (the river stood within ~1 m of
	// its peak on those days).
	// Breach geometry: in the 50 m DEM the embankment is 2-3 cells wide while
	// the Dykes.shp polyline only crosses the 1-cell crest chain. The breach
	// therefore (a) takes its invert from the lowest protected-side ground
	// within breach_search_radius of the segment (the floodplain at 5-8 m,
	// not the embankment shoulder), and (b) cuts a corridor of
	// breach_cut_halfwidth around the polyline so the opening goes through
	// the full embankment width. This is scenario definition (OpenLISEM
	// FlowBarriers-style user-specified barrier change), not solver physics.
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
	float courant_now <- 0.15;
	// shared relaxation-iteration switches (FlowSource and do_forced of AccuFluxDiffusive)
	float relax_source_m <- 0.0;     // m of water added per cell per iteration
	bool  relax_force_river <- true; // reset river cells to the forced stage each iteration
	float flooded_area_km2 <- 0.0;
	float flood_volume_mm3 <- 0.0;  // million m3
	int breaches_open <- 0;
	float static_flooded_km2 <- 0.0;

	init {
		write "=== Hanoi FastFlood (LISEM / OpenLISEM port) ===";
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

		// --- initial state ----------------------------------------------------
		river_discharge <- discharge_at(starting_date);
		river_stage <- stage_at(starting_date);
		ask river_cells { h <- max(0.0, river_stage - z_dyn); }
		ask cell { do update_color; }
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
		ask dyke where (each.will_break and !each.opened and current_date >= each.breach_date) {
			do open_breach;
		}
		breaches_open <- dyke count (each.opened);
	}

	// ======================================================================
	//  Quasi-steady FastFlood relaxation (port of AccuFluxDiffusive)
	// ======================================================================
	// one relaxation iteration of the artificial-velocity scheme; shared by
	// the dynamic simulation and the static refinement pass (step iv)
	action relax_iteration (float courant_val) {
		courant_now <- courant_val;
		ask cell parallel: true { do compute_outflux; }
		ask cell parallel: true { do apply_flux; }
		ask cell parallel: true {
			float hv <- (relax_force_river and is_river) ? max(0.0, river_stage - z_dyn) : h_new;
			dh_abs <- abs(hv - h);
			h <- hv;
		}
	}

	// Each hour the water field is relaxed towards the steady state belonging
	// to the current stage / breach configuration; with relax_to_convergence
	// the iterations continue until that steady state is actually reached.
	reflex fastflood_relax {
		relax_force_river <- true;
		relax_source_m <- 0.0;
		int n_max <- relax_to_convergence ? max_relax_iterations : sub_iterations;
		loop it from: 1 to: n_max {
			// in-cycle decay of the artificial velocity, as in the C++ code:
			// courant_here = courant (1-progress) + progress courant/8
			float progress <- it / n_max;
			do relax_iteration(courant_fastflood * (1.0 - progress) + progress * courant_fastflood / 8.0);
			// steady-state check every 10 iterations once the base budget is spent
			if relax_to_convergence and it >= sub_iterations and (it mod 10 = 0) {
				if (cell mean_of each.dh_abs) < relax_tolerance { break; }
			}
		}
	}

	// ======================================================================
	//  Bookkeeping: arrival times, statistics, colors, stop condition
	// ======================================================================
	reflex bookkeeping {
		ask cell parallel: true {
			if h > flood_threshold and !is_river {
				if arrival_h < 0.0 { arrival_h <- (current_date - starting_date) / 3600.0; }
				if h > h_peak { h_peak <- h; }
			}
			// peak flow velocity (paper fig. 6D; C++ Vel map; OpenLISEM Vmax.map):
			// diffusive-wave Manning velocity on the steepest water-surface slope
			if h > h_eps {
				float wsl <- z_dyn + h;
				float dzmax <- max([
					nE = nil ? h : wsl - (nE.z_dyn + nE.h),
					nW = nil ? h : wsl - (nW.z_dyn + nW.h),
					nN = nil ? h : wsl - (nN.z_dyn + nN.h),
					nS = nil ? h : wsl - (nS.z_dyn + nS.h), 0.0]);
				float u <- (h ^ (2.0 / 3.0)) * sqrt(dzmax / cell_dx) / n_man;
				if u > u_peak { u_peak <- u; }
			}
			do update_color;
		}
		ask observation_point where (each.arrival_h < 0.0) {
			if my_cell != nil and my_cell.h > flood_threshold {
				arrival_h <- (current_date - starting_date) / 3600.0;
				write "Observation point " + pid + " reached by the flood on " + current_date
					+ " (h = " + (my_cell.h with_precision 2) + " m)";
			}
		}
		list<cell> wet <- cell where (each.h > flood_threshold and !each.is_river);
		flooded_area_km2 <- length(wet) * cell_dx * cell_dx / 1e6;
		flood_volume_mm3 <- (wet sum_of each.h) * cell_dx * cell_dx / 1e6;
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
					float zmin_nb <- min([
						ce.nE = nil ? ce.zcorr : ce.nE.zcorr,
						ce.nW = nil ? ce.zcorr : ce.nW.zcorr,
						ce.nN = nil ? ce.zcorr : ce.nN.zcorr,
						ce.nS = nil ? ce.zcorr : ce.nS.zcorr]);
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
					b <- min(10.0, max(-0.9, (2.0 * ratio - 1.0) / (1.0 - ratio)));
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
			h_new <- 0.0; q_e <- 0.0; q_w <- 0.0; q_n <- 0.0; q_s <- 0.0;
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
			h <- 0.0; h_new <- 0.0; q_e <- 0.0; q_w <- 0.0; q_n <- 0.0; q_s <- 0.0;
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
//  Raster domain
// ==========================================================================
grid cell file: dem_file neighbors: 4 {
	// terrain
	float z;       // raw DEM elevation (contains the dyke crests)
	float z_fill;  // pit-adjusted elevation (inline pit term of AccuFluxDiffusive)
	float z_dyn;   // elevation used by the dynamic solver (lowered when a dyke breaches)
	float n_man;   // Manning roughness
	bool is_river <- false;
	bool is_lake  <- false;
	bool is_dyke  <- false;

	// dynamic solver state
	float h <- 0.0;      // water depth (m)
	float h_new <- 0.0;
	float q_e <- 0.0; float q_w <- 0.0; float q_n <- 0.0; float q_s <- 0.0;
	float h_peak <- 0.0;
	float u_peak <- 0.0;     // m/s, peak diffusive-wave Manning velocity (Vmax)
	float dh_abs <- 0.0;     // m, |dh| of the last relaxation iteration (convergence check)
	float arrival_h <- -1.0; // hours after simulation start when first flooded

	// neighbour references (E/W = +x/-x, S/N = +y/-y in grid coordinates)
	cell nE; cell nW; cell nN; cell nS;

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

	// ------------------------------------------------------------------
	// Pass 1: export a courant fraction of the stored water to downslope
	// neighbours, weighted by the diffusive-wave Manning velocity on the
	// water-surface slope (AccuFluxDiffusive: GetVelocity + weight split).
	// Nil neighbours behave like the C++ OUTORMV case: ground at own z with
	// no water, so the domain boundary drains freely.
	// ------------------------------------------------------------------
	action compute_outflux {
		q_e <- 0.0; q_w <- 0.0; q_n <- 0.0; q_s <- 0.0;
		if h > h_eps {
			float wsl <- z_dyn + h;
			float dz_e <- nE = nil ? h : wsl - (nE.z_dyn + nE.h);
			float dz_w <- nW = nil ? h : wsl - (nW.z_dyn + nW.h);
			float dz_n <- nN = nil ? h : wsl - (nN.z_dyn + nN.h);
			float dz_s <- nS = nil ? h : wsl - (nS.z_dyn + nS.h);
			float h23 <- h ^ (2.0 / 3.0);
			float w_e <- dz_e > 0.0 ? h23 * sqrt(dz_e / cell_dx) / n_man : 0.0;
			float w_w <- dz_w > 0.0 ? h23 * sqrt(dz_w / cell_dx) / n_man : 0.0;
			float w_n <- dz_n > 0.0 ? h23 * sqrt(dz_n / cell_dx) / n_man : 0.0;
			float w_s <- dz_s > 0.0 ? h23 * sqrt(dz_s / cell_dx) / n_man : 0.0;
			float w_tot <- w_e + w_w + w_n + w_s;
			if w_tot > 0.0 {
				float q_out <- courant_now * h;
				// equilibrium limiter of the C++ code (anti-oscillation):
				// no neighbour may receive more than what levels both surfaces
				q_e <- w_e > 0.0 ? min(q_out * w_e / w_tot, stab_factor * dz_e) : 0.0;
				q_w <- w_w > 0.0 ? min(q_out * w_w / w_tot, stab_factor * dz_w) : 0.0;
				q_n <- w_n > 0.0 ? min(q_out * w_n / w_tot, stab_factor * dz_n) : 0.0;
				q_s <- w_s > 0.0 ? min(q_out * w_s / w_tot, stab_factor * dz_s) : 0.0;
			}
		}
	}

	// ------------------------------------------------------------------
	// Pass 2: gather incoming fluxes (HN = H + QIN - QOUT of the C++ code)
	// ------------------------------------------------------------------
	action apply_flux {
		float q_in <- (nE = nil ? 0.0 : nE.q_w) + (nW = nil ? 0.0 : nW.q_e)
		            + (nN = nil ? 0.0 : nN.q_s) + (nS = nil ? 0.0 : nS.q_n);
		h_new <- max(0.0, h - (q_e + q_w + q_n + q_s) + q_in + relax_source_m);
	}

	action update_color {
		float shade <- (z - z_min) / max(0.001, z_max - z_min);
		rgb terrain <- rgb(70 + int(150 * shade), 80 + int(130 * shade), 60 + int(110 * shade));
		if is_dyke and h <= flood_threshold {
			color <- rgb(110, 80, 60);
		} else if h > 0.02 {
			float f <- min(1.0, h / 4.0);
			color <- rgb(int(150 * (1 - f)), int(190 * (1 - f) + 30), int(180 + 75 * f));
		} else {
			color <- terrain;
		}
	}
}

// ==========================================================================
//  Vector species
// ==========================================================================
species dyke {
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
		ask corridor { z_dyn <- min(z_dyn, target); }
		write "BREACH at " + commune + " on " + current_date + " (invert lowered to "
			+ (target with_precision 2) + " m, " + length(corridor) + " cells cut)";
	}

	aspect default {
		draw shape color: opened ? #red : (will_break ? #orange : #darkgreen) width: 3;
	}
}

species river_poly {
	aspect default { draw shape color: rgb(70, 130, 180, 120) border: #steelblue; }
}

species lake {
	aspect default { draw shape color: rgb(150, 200, 230, 150); }
}

species building {
	aspect default { draw shape color: rgb(90, 90, 90); }
}

species observation_point {
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
//  Experiment
// ==========================================================================
experiment fastflood_1926 type: gui {
	parameter "Courant fraction (artificial velocity)" var: courant_fastflood category: "Solver";
	parameter "Relaxation iterations per hour" var: sub_iterations category: "Solver";
	parameter "Relax to convergence (steady state)" var: relax_to_convergence category: "Solver";
	parameter "Max relaxation iterations" var: max_relax_iterations category: "Solver";
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
			species river_poly;
			species lake transparency: 0.4;
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

		display "Static FastFlood hazard" type: 2d background: #black {
			grid cell transparency: 0.55;
			graphics "static hazard" {
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
		monitor "Flooded area (km2)" value: flooded_area_km2 with_precision 2;
		monitor "Peak velocity (m/s)" value: (cell max_of each.u_peak) with_precision 2;
		monitor "Flood volume (10^6 m3)" value: flood_volume_mm3 with_precision 2;
		monitor "Arrival times (h)" value: observation_point collect (string(each.pid) + ": "
			+ (each.arrival_h < 0.0 ? "-" : string(each.arrival_h with_precision 1)));
	}
}
