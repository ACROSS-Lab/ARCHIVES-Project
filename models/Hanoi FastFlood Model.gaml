/**
* Hanoi 1926 Red River flood — FastFlood steady-state model
* =========================================================================
* Implements: van den Bout, Jetten, van Westen & Lombardo (2023),
*   "A breakthrough in fast flood simulation", Environ. Model. Softw. 168, 105787.
*   (includes/1-s2.0-S1364815223001731-main.pdf)
*
* This is NOT a time-stepping dynamic model. It estimates the PEAK flood
* (depth, velocity, arrival order) as a COMPENSATED STEADY STATE, which is why
* it is fast and free of the CFL instabilities of the dynamic models. We run the
* paper's discharge-boundary mode — the same setup as their dijkring-41 levee
* breach case — driven by the documented 1926 dyke breaches:
*
*   (1) Fast-Sweeping hydrological DEM correction  -> depression-free zc + flow net
*       (Appendix A, Eq 21-23; Planchon-Darboux fill accelerated by 4-way sweeps)
*   (2) Steady-state flow accumulation of the breach inflow  -> discharge AF  (Eq 1)
*       plus AF(1) and AF(AF(1)) for the catchment-shape parameter
*   (3) Manning inversion (Eq 2) + partial-steady-state compensation (Eq 3-14)
*       -> per-cell shape b, s_max, s_ss, factor f_ss, compensated discharge q_c,
*          and the first-guess steady-state flow height h_af
*   (4) Adaptive diffusive-wave inundation refinement (Eq 15-19) on the REAL DEM,
*       fed by a weir source at the breaches, draining at the OPEN domain edges,
*       relaxed to equilibrium (inflow = outflow) with an artificial-velocity
*       acceleration that ramps down to 1.
*   (5) Arrival-time field along the flow network (seeded with the breach dates)
*       -> arrival ORDER at the 5 observation points, vs the documented sequence.
*
* SOURCE  : dyke breaches (Dykes.shp BREAK/DATE). Weir inflow Q = C*L*(stage-crest)^1.5
*           at the 1926 peak river stage. Gia-Quat breaks 28-07, the rest 29-07.
* DATA    : mnt-gz50.asc (219x142, ~50 m), RedRiver1925.shp, Dykes.shp,
*           5_arrival_time.shp, WaterDischarge.csv.
* OUTPUT  : peak depth map, velocity, arrival time + a written arrival-order report.
*
* All heavy work (steps 1-3) is precomputed ONCE in init; the displayed reflex
* runs the step-4 relaxation so you watch the basin fill to its peak.
*/
model HanoiFastFlood

global {
	// ------------------------------------------------------------- inputs
	file dem_file    <- grid_file("../includes/mnt-gz50.asc");
	file river_file  <- file("../includes/RedRiver1925.shp");
	file dykes_file  <- file("../includes/Dykes.shp");
	file points_file <- file("../includes/5_arrival_time.shp");
	csv_file q_csv   <- csv_file("../includes/WaterDischarge.csv", ",", true);

	geometry shape <- envelope(dem_file);

	// ------------------------------------------------- grid / geometry
	float dx <- 49.9736;                 // cell size (m), from the DEM header
	int   nb_cols;
	int   nb_rows;

	// --------------------------------------- (1) hydrological correction
	float delta_slope <- 0.001;          // min elevation rise per cell (m) — Eq 21 delta
	int   sweep_rounds <- 40;            // 4-direction rounds for every sweeping solve

	// --------------------------------------------- physics / calibration
	float manning   <- 0.06;             // floodplain Manning n (calibration knob)
	float gravity   <- 9.81;
	float min_depth <- 0.01;             // wet/dry threshold (m)

	// -------------------------------------------- (4) breach weir source
	float weir_C  <- 1.7;                // broad-crested weir coeff (SI: Q = C*L*H^1.5)
	float h_peak  <- 13.3;               // 1926 peak river stage (m) — calibration knob
	float h_base  <- 11.5;               // baseline stage (m) — context only
	int   day0_july <- 22;               // discharge series starts 22 July 1926
	float river_stage;                   // steady peak stage that drives the weir
	float Q_peak;                        // peak river discharge (context, from CSV)

	// ------------------------------ (3) partial steady-state parameters
	float event_days <- 8.0;             // representative event duration (days) for s_ss
	float v_mean     <- 0.3;             // a-priori mean velocity (m/s): s_ss = v_mean * t

	// --------------------------------------- (4) diffusive relaxation
	// Gauss-Seidel, flux-limited redistribution: each transfer is capped to a
	// quarter of the head difference and to the water available, so it stays
	// stable on flat ground (an explicit diffusive update blows up there).
	int   iters_per_cycle <- 3;          // redistribution passes per displayed cycle
	float relax_dt  <- 2.0;              // pseudo-step (s); the caps keep it stable for any value
	int   max_cycles <- 800;
	int   cycle <- 0;                    // RELAXATION ITERATION counter — NOT real time
	float conv_tol  <- 0.005;            // converged when max |dh| per cycle < this (m)
	float max_dh <- 1.0;                 // convergence monitor (max |dh| this cycle)
	bool  converged <- false;

	// --------------------------------------- precomputed cell sets
	list<cell> breach_cells <- [];
	list<cell> sample_cells <- [];

	init {
		// --- discharge series (context: report the 1926 peak) ---
		matrix data <- matrix(q_csv);
		list<float> Q_series <- [];
		loop i from: 0 to: data.rows - 1 {
			string s <- string(data[1, i]);
			if (s != nil and s != "" and s != "m3/s") { add float(s) to: Q_series; }
		}
		Q_peak <- empty(Q_series) ? 0.0 : max(Q_series);
		river_stage <- h_peak;

		// --- grid dimensions + neighbours + per-cell init ---
		nb_cols <- (cell max_of each.grid_x) + 1;
		nb_rows <- (cell max_of each.grid_y) + 1;
		ask cell {
			z <- grid_value;
			active <- z > -1000.0;                       // guard NODATA (-9999); this DEM has none
			zc <- z; h <- 0.0;
			is_edge <- (grid_x = 0 or grid_y = 0 or grid_x = nb_cols - 1 or grid_y = nb_rows - 1);
			if (grid_x < nb_cols - 1) { nE <- cell grid_at {grid_x + 1, grid_y}; }
			if (grid_x > 0)           { nW <- cell grid_at {grid_x - 1, grid_y}; }
			if (grid_y < nb_rows - 1) { nS <- cell grid_at {grid_x, grid_y + 1}; }
			if (grid_y > 0)           { nN <- cell grid_at {grid_x, grid_y - 1}; }
		}

		// --- dyke breaches -> datable openings (the water source) ---
		create dyke_seg from: dykes_file with: [
			brk::string(read("BREAK")), dnum_s::string(read("DATE")), commune::string(read("Commune"))
		];
		ask dyke_seg where (each.brk = "YES") {
			int dday <- int(first(dnum_s split_with "-"));   // "28-07" -> 28
			open_time <- (dday - myself.day0_july) * 86400.0;
			ask cell overlapping self {
				is_breach <- true;
				breach_open <- (breach_open < 0.0) ? myself.open_time : min(breach_open, myself.open_time);
			}
		}
		breach_cells <- cell where each.is_breach;
		write "Breach cells: " + length(breach_cells)
			+ " (earliest opens day " + ((empty(breach_cells)) ? -1.0 : (breach_cells min_of each.breach_open) / 86400.0) + ")";

		// --- river polygon -> fixed-stage reservoir cells (the water SOURCE side) ---
		// Held at the river stage and EXCLUDED from floodplain redistribution, so the
		// river feeds the protected area ONLY through the breach weir (dykes hold).
		create river_area from: river_file;
		ask river_area { ask cell overlapping self { is_river <- true; } }
		write "River cells: " + length(cell where each.is_river);

		// --- observation points (record order = ids 1..5) ---
		list<arrival_pt> pts <- [];
		create arrival_pt from: points_file returns: pts;
		loop i from: 0 to: length(pts) - 1 {
			ask (pts at i) {
				sid <- i + 1;
				cell c <- first(cell overlapping self);
				if (c != nil) { c.is_sample <- true; c.sample_id <- sid; }
			}
		}
		sample_cells <- cell where each.is_sample;

		// ============================ FASTFLOOD PIPELINE (steps 1-3) ============================
		write "[1/4] Fast-Sweeping DEM hydro-correction ...";
		do correct_dem;            // (1) depression-free zc
		do velocity_field;         // (2a) flow directions + accumulation weights (Eq 23)
		do compute_weir_source;    // breach inflow discharge per cell
		write "[2/4] Steady-state flow accumulation ...";
		do accumulate_field(1);    // AF(1)  : upstream cell count
		do accumulate_field(2);    // AF(AF(1))
		do accumulate_field(0);    // AF(q)  : steady-state discharge from the breach (Eq 1)
		write "[3/4] Manning inversion + partial-steady-state compensation ...";
		do compensate;             // (3) b, s_max, s_ss, f_ss, q_c, h_af (Eq 2-14)
		do seed_inundation;        // warm-start the relaxation from h_af
		write "[4/4] Diffusive-wave inundation relaxation (running) ...";

		// --- colour the static terrain once ---
		ask cell where each.active {
			int gv <- int(min(235, max(60, (z - 4) * 11)));
			base_color <- rgb(gv, gv, gv);
			color <- base_color;
		}
		do recolor;
	}

	// =====================================================================
	// (1) FAST-SWEEPING DEM CORRECTION  (Appendix A, Eq 20-22)
	// Planchon-Darboux depression fill: zc = max(z, min_neighbour(zc) + delta),
	// boundary cells fixed at z (outlets). Four directional sweeps per round use
	// the natural distance ordering so it converges in a handful of rounds.
	// =====================================================================
	action correct_dem {
		float BIG <- 1e9;
		ask cell where each.active { zc <- is_edge ? z : BIG; }
		loop r from: 0 to: sweep_rounds - 1 {
			loop dir from: 0 to: 3 {
				bool rx <- (dir = 1 or dir = 3);
				bool ry <- (dir = 2 or dir = 3);
				loop ix from: 0 to: nb_cols - 1 {
					int gx <- rx ? (nb_cols - 1 - ix) : ix;
					loop iy from: 0 to: nb_rows - 1 {
						int gy <- ry ? (nb_rows - 1 - iy) : iy;
						cell c <- cell grid_at {gx, gy};
						if (c.active and not c.is_edge) {
							float m <- BIG;
							loop nb over: [c.nE, c.nW, c.nS, c.nN] {
								if (nb != nil and nb.active and nb.zc < m) { m <- nb.zc; }
							}
							float cand <- max(c.z, m + delta_slope);
							if (cand < c.zc) { c.zc <- cand; }
						}
					}
				}
			}
		}
	}

	// =====================================================================
	// (2a) VELOCITY FIELD + ACCUMULATION WEIGHTS  (Eq 23)
	// Downslope unit vector of the corrected DEM. Outflow is split between the two
	// downslope neighbours in proportion to |fx|,|fy| (a D-infinity-style network).
	// =====================================================================
	action velocity_field {
		ask cell where each.active {
			float zE <- (nE != nil and nE.active) ? nE.zc : zc;
			float zW <- (nW != nil and nW.active) ? nW.zc : zc;
			float zS <- (nS != nil and nS.active) ? nS.zc : zc;
			float zN <- (nN != nil and nN.active) ? nN.zc : zc;
			float gxg <- (zE - zW) / (2 * dx);           // d(zc)/dx
			float gyg <- (zS - zN) / (2 * dx);           // d(zc)/dy  (S is +y)
			float gmag <- sqrt(gxg * gxg + gyg * gyg);
			if (gmag < 1e-9) {
				fx <- 0.0; fy <- 0.0; slope <- delta_slope / dx;
				wE <- 0.0; wW <- 0.0; wS <- 0.0; wN <- 0.0;
			} else {
				fx <- -gxg / gmag;                       // downslope unit x
				fy <- -gyg / gmag;                       // downslope unit y
				slope <- gmag;                           // |grad zc| (m/m), downslope slope
				float wsum <- abs(fx) + abs(fy);
				wE <- max(0.0,  fx) / wsum; wW <- max(0.0, -fx) / wsum;
				wS <- max(0.0,  fy) / wsum; wN <- max(0.0, -fy) / wsum;
			}
		}
	}

	// breach inflow discharge per cell: weir over the dyke crest (= cell elevation)
	action compute_weir_source {
		ask cell { q_in <- 0.0; }
		ask breach_cells {
			float head <- max(0.0, myself.river_stage - z);
			q_in <- myself.weir_C * myself.dx * (head ^ 1.5);   // m3/s entering this cell
		}
		write "Total breach inflow Q = " + (breach_cells sum_of each.q_in) with_precision 1
			+ " m3/s  (peak river Q on record = " + Q_peak + " m3/s)";
	}

	// =====================================================================
	// (2b) FLOW ACCUMULATION via sweeping  (Eq 1)
	// Pull formulation: each cell = its own source + the share of each upslope
	// neighbour's value that drains into it. Repeated 4-way sweeps converge.
	//   mode 0 -> af  (source q_in)   : steady-state discharge
	//   mode 1 -> af1 (source 1)      : upstream cell count
	//   mode 2 -> afa (source af1)    : AF(AF(1))
	// =====================================================================
	action accumulate_field (int mode) {
		ask cell where each.active {
			if      (mode = 0) { af  <- q_in; }
			else if (mode = 1) { af1 <- 1.0; }
			else               { afa <- af1; }
		}
		loop r from: 0 to: sweep_rounds - 1 {
			loop dir from: 0 to: 3 {
				bool rx <- (dir = 1 or dir = 3);
				bool ry <- (dir = 2 or dir = 3);
				loop ix from: 0 to: nb_cols - 1 {
					int gx <- rx ? (nb_cols - 1 - ix) : ix;
					loop iy from: 0 to: nb_rows - 1 {
						int gy <- ry ? (nb_rows - 1 - iy) : iy;
						cell c <- cell grid_at {gx, gy};
						if (c.active) {
							float inflow <- 0.0;
							if (c.nW != nil and c.nW.active) { inflow <- inflow + (mode = 0 ? c.nW.af : (mode = 1 ? c.nW.af1 : c.nW.afa)) * c.nW.wE; }
							if (c.nE != nil and c.nE.active) { inflow <- inflow + (mode = 0 ? c.nE.af : (mode = 1 ? c.nE.af1 : c.nE.afa)) * c.nE.wW; }
							if (c.nN != nil and c.nN.active) { inflow <- inflow + (mode = 0 ? c.nN.af : (mode = 1 ? c.nN.af1 : c.nN.afa)) * c.nN.wS; }
							if (c.nS != nil and c.nS.active) { inflow <- inflow + (mode = 0 ? c.nS.af : (mode = 1 ? c.nS.af1 : c.nS.afa)) * c.nS.wN; }
							float src <- (mode = 0 ? c.q_in : (mode = 1 ? 1.0 : c.af1));
							float val <- src + inflow;
							if      (mode = 0) { c.af  <- val; }
							else if (mode = 1) { c.af1 <- val; }
							else               { c.afa <- val; }
						}
					}
				}
			}
		}
	}

	// =====================================================================
	// (3) MANNING INVERSION + PARTIAL-STEADY-STATE COMPENSATION  (Eq 2-14)
	// Catchment shape b solved per cell from Eq 4 & 6 combined:
	//     afa/af1 = (1+b)/(2+b) * af1^(1/(1+b))
	// Then s_max (Eq 6), s_ss = v_mean*t (event distance), f_ss (Eq 12-13),
	// compensated discharge q_c (Eq 14), and first-guess height h_af (Eq 2).
	// =====================================================================
	action compensate {
		ask cell where each.active {
			float bsol <- 0.5;
			if (af1 > 1.0) {
				float ratio <- afa / af1;
				float blo <- 0.001; float bhi <- 0.999;
				float flo <- ((1 + blo) / (2 + blo)) * (af1 ^ (1.0 / (1 + blo))) - ratio;
				loop k from: 0 to: 40 {
					float bm <- 0.5 * (blo + bhi);
					float fm <- ((1 + bm) / (2 + bm)) * (af1 ^ (1.0 / (1 + bm))) - ratio;
					if ((fm > 0.0) = (flo > 0.0)) { blo <- bm; flo <- fm; } else { bhi <- bm; }
				}
				bsol <- 0.5 * (blo + bhi);
			}
			bshape <- bsol;
			smax <- dx * (af1 ^ (1.0 / (1.0 + bshape)));              // Eq 6 rearranged
			sss  <- v_mean * (event_days * 86400.0);                 // event travel distance
			float sn <- (smax > 0.0) ? min(1.0, sss / smax) : 1.0;
			fss  <- (smax > sss) ? (sn ^ (1.0 + bshape)) : 1.0;      // Eq 12-13
			qc   <- fss * af;                                        // Eq 14 (af is already vol/time)
			float S <- max(slope, 1e-4);
			h_af <- (qc > 0.0) ? (((qc * manning) / (dx * sqrt(S))) ^ 0.6) : 0.0;   // Eq 2
		}
	}

	// start dry; the breach weir source fills the basin during relaxation.
	// (h_af from step 3 is kept as the FastFlood steady-state estimate / output,
	//  but is NOT used as the seed — on flat ground its slope->0 makes it ill-posed.)
	action seed_inundation {
		ask cell where each.active { h <- 0.0; hmax <- 0.0; }
		// fill the river channel to the stage (display + the head the weir sees)
		ask cell where (each.is_river and not each.is_breach) {
			h <- max(0.0, myself.river_stage - z);
			hmax <- h;
		}
	}

	// =====================================================================
	// (4) ADAPTIVE DIFFUSIVE-WAVE INUNDATION (Eq 15-19) — displayed relaxation
	// Two-phase per iteration: compute face discharges (analytic Manning velocity,
	// Eq 17, water-surface form), then update depth by mass balance with the breach
	// source and free outflow at the open domain edges. Acceleration ramps to 1.
	// =====================================================================
	reflex relax when: not converged and cycle < max_cycles {
		ask cell where each.active { h_prev <- h; }
		// hold the river channel at its stage (Dirichlet reservoir)
		ask cell where (each.is_river and not each.is_breach) { h <- max(0.0, river_stage - z); }
		loop times: iters_per_cycle {
			// inject the breach weir inflow (discharge -> depth) ...
			ask breach_cells {
				h <- h + relax_dt * q_in / (dx * dx);
				if (h > hmax) { hmax <- h; }
			}
			// ... then push water downhill, highest water-surface first, so it
			// cascades in one Gauss-Seidel pass (sequential -> race-free, stable).
			// River cells are excluded (reservoir feeds only via the breach weir).
			list<cell> wet <- (cell where (each.active and each.h > min_depth
				and (each.is_breach or not each.is_river))) sort_by (-(each.z + each.h));
			ask wet { do redistribute; }
		}
		max_dh <- (cell where each.active) max_of (abs(each.h - each.h_prev));
		cycle <- cycle + 1;
		if (cycle > 2 and max_dh < conv_tol) { do finalize; }
		do recolor;
	}

	action recolor {
		ask cell where each.active parallel: true {
			if (h > min_depth) {
				color <- rgb(0, int(max(40, 170 - h * 18)), 255);            // blue, darker = deeper
			} else {
				color <- base_color;
			}
		}
	}

	// =====================================================================
	// (5) FINALIZE — velocity, arrival-time field, arrival-order report
	// =====================================================================
	action finalize {
		converged <- true;
		// peak velocity from Manning on the converged water surface: v = (1/n) h^(2/3) sqrt(|grad eta|)
		ask cell where each.active {
			if (h > min_depth) {
				float etaC <- z + h;
				float eE <- (nE != nil and nE.active) ? (nE.z + nE.h) : etaC;
				float eW <- (nW != nil and nW.active) ? (nW.z + nW.h) : etaC;
				float eS <- (nS != nil and nS.active) ? (nS.z + nS.h) : etaC;
				float eN <- (nN != nil and nN.active) ? (nN.z + nN.h) : etaC;
				float gmag <- sqrt(((eE - eW) / (2 * dx)) ^ 2 + ((eS - eN) / (2 * dx)) ^ 2);
				vmag <- (1.0 / manning) * (h ^ 0.6667) * sqrt(gmag);
			} else { vmag <- 0.0; }
		}
		do arrival_time;
		do report;
		write "Converged after " + cycle + " cycles (max dh = " + max_dh + " m). Paused.";
		do pause;
	}

	// arrival time: multi-source Eikonal from the breaches (seeded with the breach
	// dates) propagated through wet cells at the local flow velocity.
	action arrival_time {
		float BIG <- 1e12;
		ask cell where each.active { t_arr <- is_breach ? breach_open : BIG; }
		loop r from: 0 to: sweep_rounds - 1 {
			loop dir from: 0 to: 3 {
				bool rx <- (dir = 1 or dir = 3);
				bool ry <- (dir = 2 or dir = 3);
				loop ix from: 0 to: nb_cols - 1 {
					int gx <- rx ? (nb_cols - 1 - ix) : ix;
					loop iy from: 0 to: nb_rows - 1 {
						int gy <- ry ? (nb_rows - 1 - iy) : iy;
						cell c <- cell grid_at {gx, gy};
						if (c.active and not c.is_breach) {
							float best <- c.t_arr;
							loop nb over: [c.nW, c.nE, c.nN, c.nS] {
								if (nb != nil and nb.active and nb.t_arr < BIG and nb.h > min_depth and nb.vmag > 1e-3) {
									float tt <- nb.t_arr + dx / max(nb.vmag, 1e-3);
									if (tt < best) { best <- tt; }
								}
							}
							c.t_arr <- best;
						}
					}
				}
			}
		}
	}

	action report {
		write "==================== FASTFLOOD ARRIVAL ORDER ====================";
		list<cell> samples <- sample_cells sort_by each.sample_id;
		ask samples {
			if (t_arr >= 1e11) {
				write "  pt" + sample_id + " : NOT reached (z=" + z with_precision 1 + ")";
			} else {
				write "  pt" + sample_id + " : " + (t_arr / 86400.0) with_precision 2
					+ " days  (peak depth=" + hmax with_precision 2 + " m, z=" + z with_precision 1 + ")";
			}
		}
		list<cell> reached <- (samples where (each.t_arr < 1e11)) sort_by each.t_arr;
		write "  ORDER reached: " + (reached collect ("pt" + each.sample_id));
		write "  Documented breaches: Gia-Quat 28-07 (earliest), then Ai-Mo / Lam-Giu / Gia-Quat 29-07";
		write "================================================================";
	}
}

// ====================================================================== grid
// Terrain + all FastFlood state. frequency:0 -> driven by the global reflexes
// (no per-cell scheduling), like the project's other grid models.
grid cell file: dem_file neighbors: 4 frequency: 0
	use_regular_agents: false use_individual_shapes: false {
	float z;                 // original bed elevation (m)
	float zc;                // hydrologically-corrected elevation (m)
	float h;                 // flow depth (m)
	float hmax;              // peak flow depth over the relaxation (m)
	float h_af;              // inverted steady-state flow height (first guess, Eq 2)

	// flow network (from corrected DEM)
	float fx; float fy;      // downslope unit vector
	float slope;             // |grad zc| (m/m)
	float wE; float wW; float wS; float wN;   // normalised outflow weights

	// accumulation fields
	float af;                // steady-state discharge (Eq 1)
	float af1;               // AF(1) upstream cell count
	float afa;               // AF(AF(1))

	// compensation
	float bshape;            // catchment shape parameter b
	float smax;              // max catchment distance (m)
	float sss;               // steady-state travel distance (m)
	float fss;               // partial-steady-state factor
	float qc;                // compensated peak discharge

	// source + diffusive solver
	float q_in;              // breach inflow discharge (m3/s) — 0 except breach cells
	float vmag;              // velocity magnitude at peak (m/s)
	float t_arr <- 1e12;     // flood arrival time (s)
	float h_prev;            // depth at start of cycle (convergence check)

	bool active   <- true;
	bool is_edge  <- false;
	bool is_breach <- false;
	bool is_river  <- false;     // river channel: fixed-stage reservoir (the source side)
	bool is_sample <- false;
	int  sample_id <- 0;
	float breach_open <- -1.0;   // breach opening time (s since day0)

	rgb base_color <- #gray;
	rgb color <- #gray;
	cell nE; cell nW; cell nS; cell nN;

	// Gauss-Seidel diffusive redistribution (Eq 15-19, flux-limited form).
	// Push water to every lower-water-surface neighbour with the Manning flux
	//   q = (1/n) hflow^(5/3) sqrt(Sw),   capped to:
	//   (a) a quarter of the head difference  -> no overshoot, stable on flat ground;
	//   (b) the water actually available      -> depth never goes negative.
	// A missing/inactive neighbour is free outfall (open boundary): outside eta = bed.
	action redistribute {
		if (h > min_depth) {
			loop nb over: [nE, nW, nS, nN] {
				float etaC <- z + h;
				float zN <- z; float etaN <- z; bool outside <- true;
				if (nb != nil and nb.active) { zN <- nb.z; etaN <- nb.z + nb.h; outside <- false; }
				if (etaC > etaN) {
					float hf <- etaC - max(z, zN);
					if (hf > min_depth) {
						float dhead <- etaC - etaN;
						float q  <- (hf ^ 1.6667) / manning * sqrt(dhead / dx);   // m2/s
						float dV <- relax_dt * q / dx;                            // depth (m)
						dV <- min(dV, 0.25 * dhead);                              // no overshoot
						dV <- min(dV, h - min_depth);                             // keep depth >= 0
						if (dV > 0.0) {
							h <- h - dV;
							if (not outside) {
								nb.h <- nb.h + dV;
								if (nb.h > nb.hmax) { nb.hmax <- nb.h; }
							}
						}
					}
				}
			}
			if (h > hmax) { hmax <- h; }
		}
	}
}

// ====================================================================== other
species dyke_seg {
	string brk;
	string dnum_s;
	string commune;
	float  open_time <- -1.0;
	aspect default { draw shape color: (brk = "YES") ? #red : rgb(120, 80, 40) width: 2; }
}

species arrival_pt {
	int sid;
	aspect default {
		draw circle(90) color: #yellow border: #black;
		draw string(sid) color: #black size: 14 at: location + {0, -120};
	}
}

species river_area {
	aspect default { draw shape color: rgb(40, 90, 160, 60) border: rgb(30, 70, 140); }
}

experiment HanoiFastFlood type: gui {
	parameter "Manning n"              var: manning      min: 0.02 max: 0.12 step: 0.005;
	parameter "Peak river stage (m)"   var: h_peak       min: 11.5 max: 16.0 step: 0.1;
	parameter "Weir coefficient C"     var: weir_C       min: 1.0  max: 2.2  step: 0.05;
	parameter "Hydro-correct slope (m/cell)" var: delta_slope min: 0.0001 max: 0.01 step: 0.0005;
	parameter "Event duration (days)"  var: event_days   min: 1.0  max: 16.0 step: 1.0;
	parameter "Relaxation dt (s)"      var: relax_dt     min: 0.5  max: 20.0 step: 0.5;

	output {
		display "Flood depth" type: 2d {
			grid cell;
			species river_area aspect: default;
			species dyke_seg aspect: default;
			species arrival_pt aspect: default;
			overlay position: {10, 10} size: {340 #px, 80 #px} background: #black transparency: 0.5 {
				draw "relaxation iter " + world.cycle + (world.converged ? "  (converged)" : "  (solving)")
					at: {15 #px, 24 #px} color: #white font: font("Helvetica", 13, #bold);
				draw "max change " + (world.max_dh with_precision 4) + " m   (steady-state solver, not a clock)"
					at: {15 #px, 48 #px} color: rgb(185, 185, 185) font: font("Helvetica", 11, #plain);
			}
		}
		monitor "Relaxation iter"  value: cycle;
		monitor "Converged?"       value: converged;
		monitor "Max |dh| (m)"     value: max_dh with_precision 4;
		monitor "Wet cells"        value: cell count (each.active and each.h > min_depth);
		monitor "Max depth (m)"    value: (cell where each.active max_of each.hmax) with_precision 2;
		monitor "Breach Q (m3/s)"  value: (breach_cells sum_of each.q_in) with_precision 0;
	}
}
