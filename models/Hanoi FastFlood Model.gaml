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
*       -> per-cell shape b, s_max, s_ss (s_ss from the model's OWN mean velocity),
*          factor f_ss, compensated discharge q_c, and the compensated steady-state
*          flow height h_af. THIS field IS the flood estimate (steps 1-3 drive it).
*   (4) Diffusive-wave inundation on the REAL DEM = paper's Eq 17 velocity + Eq 19
*       h^(5/3) discharge-conserving update, flux-limited for stability. The floodplain
*       starts DRY; the breach (held at the routed river stage) feeds it, so the flood
*       SPREADS from the river and FILLS the polder — closed domain edges by default
*       (the dijkring-41 storage case; set drain_edges for open through-flow). It self-
*       limits when the basin level reaches the river stage. ('exact_eq19' = raw Jacobi.)
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

	// --------------------- discharge boundary: WaterDischarge.csv -> river stage
	// Q_peak from the CSV is routed down the channel (Manning normal depth) to a
	// north-inlet -> south-outlet water-surface profile; that stage drives the breach
	// weir AND the foreland reservoir. (Replaces the hand-set peak stage.)
	float weir_C  <- 1.7;                // broad-crested weir coeff (SI: Q = C*L*H^1.5)
	float river_n <- 0.03;               // CHANNEL Manning n (for the Q -> stage routing)
	float stage_cal <- 0.78;             // calibration multiplier (Q-routed ~17 m -> ~13.3 m known peak)
	bool  allow_overtopping <- false;    // also exchange where the routed stage tops a (non-breach) dyke
	float h_peak  <- 13.3;               // fallback inlet stage (m) if the CSV carries no Q
	float h_base  <- 11.5;               // baseline stage (m) — context only
	int   day0_july <- 22;               // discharge series starts 22 July 1926
	float river_stage;                   // representative (inlet) stage — set by compute_river_stage
	float stage_north;                   // Q-driven stage at the north inlet (m)
	float stage_south;                   // Q-driven stage at the south outlet (m)
	float Q_peak;                        // peak river discharge from WaterDischarge.csv (drives the stage)
	int   foreland_reach <- 50;          // flood-stage river fills this many cells inland (to the dykes)

	// ------------------------------ (3) partial steady-state parameters
	float event_days <- 8.0;             // representative event duration (days) for s_ss
	float v_mean     <- 0.3;             // mean flow velocity (m/s) — recomputed from the field
	float fss_event  <- 1.0;             // discharge-weighted event compensation (scales the source)

	// --------------------------------------- (4) diffusive relaxation
	// Gauss-Seidel, flux-limited redistribution: each transfer is capped to a
	// quarter of the head difference and to the water available, so it stays
	// stable on flat ground (an explicit diffusive update blows up there).
	// NB: 'cycle' is a GAMA BUILT-IN (the sim step), so the relaxation counter
	// must be named something else -> relax_iter.
	int   iters_per_cycle <- 6;          // redistribution passes per displayed cycle
	float relax_dt  <- 8.0;              // pseudo-step (s); the caps keep it stable for any value
	// false = POLDER / basin-fill (closed edges -> the flood spreads in and fills, like the
	//   paper's dijkring-41 levee case); true = open downstream boundary (through-flow/conveyance).
	bool  drain_edges <- false;
	bool  exact_eq19 <- false;           // true -> paper-literal Eq 17/19 update (fragile); false -> stable GS
	float eq19_dt   <- 0.5;              // explicit step (s) for the exact_eq19 mode (keep small)
	int   max_cycles <- 3000;            // hard cap — we still finalize (report) when hit
	int   relax_iter <- 0;               // relaxation counter — NOT real time
	float conv_tol  <- 0.02;             // depth-converged when max |dh| per cycle < this (m)
	float max_dh <- 1.0;                 // convergence monitor (max |dh| this cycle)
	// flood-EXTENT convergence: stop when the inundated footprint stops growing,
	// because depth keeps creeping up long after the extent (what we care about) settles.
	int   wet_prev <- -1;
	int   stable_count <- 0;
	int   stable_needed <- 15;           // consecutive extent-stable cycles to call it settled
	bool  converged <- false;
	bool  flag_disconnect <- false;      // when true, recolor paints truly-disconnected wet cells red

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
		ask dyke_seg { ask cell overlapping self { is_dyke <- true; } }   // whole dyke line = barrier
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
		ask river_area { ask cell overlapping self { is_river <- true; is_channel <- true; } }
		do route_channel;          // inject Q at the inlet, route N->S with continuity (sets rstage + q_in spill)
		do build_reservoir;        // grow the channel out to the dykes (using the routed stage)
		write "River + foreland reservoir cells: " + length(cell where each.is_river)
			+ "  (channel grown " + foreland_reach + " cells max to the dyke line)";

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
		// (breach source q_in already set by route_channel = the river spill, with continuity)
		write "[2/4] Steady-state flow accumulation of the river spill ...";
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

	// grow the river reservoir from the channel out to the dyke line: take every low
	// cell (z < stage) 4-connected to the channel, stopping at dyke cells and high
	// ground. Bounded to 'foreland_reach' rings so a gap in the dyke ring can't run away.
	action build_reservoir {
		loop k from: 1 to: foreland_reach {
			list<cell> grow <- [];
			ask cell where (each.active and not each.is_river and not each.is_dyke and each.z < each.rstage) {
				if ((nE != nil and nE.is_river) or (nW != nil and nW.is_river)
				  or (nS != nil and nS.is_river) or (nN != nil and nN.is_river)) {
					add self to: grow;
				}
			}
			if (empty(grow)) { break; }
			ask grow { is_river <- true; }
		}
	}

	// --------------------- DISCHARGE-BOUNDARY river-corridor routing ---------------------
	// Strict version of the paper's discharge BC: inject Q_total at the north inlet and
	// sweep DOWNSTREAM (north -> south) station by station (one grid_y row of channel cells
	// = one cross-section). At each station:
	//   stage = bed + Manning normal depth for the CURRENT channel discharge & local width;
	//   breach/overtopping spill = weir(stage - crest)  -> becomes the floodplain source q_in;
	//   CONTINUITY: the channel discharge is reduced by that spill before the next station.
	// The remaining discharge leaves at the south outlet. So the stage is locally varying and
	// the river loses water as it spills (true 1D-channel-with-lateral-outflow + 2D floodplain).
	action route_channel {
		ask cell { q_in <- 0.0; rstage <- 0.0; }
		list<cell> chan <- cell where each.is_channel;
		if (empty(chan)) {
			// no channel cells: fall back to a flat stage and a weir at the mapped breaches
			ask cell { rstage <- h_peak; }
			ask breach_cells { q_in <- myself.weir_C * myself.dx * (max(0.0, h_peak - z) ^ 1.5); }
			stage_north <- h_peak; stage_south <- h_peak; river_stage <- h_peak;
		} else {
			int ymin <- chan min_of each.grid_y;                 // north inlet (top row)
			int ymax <- chan max_of each.grid_y;                 // south outlet (bottom row)
			float len <- max(dx, (ymax - ymin) * dx);
			float bed_n <- (chan where (each.grid_y <= ymin + 2)) mean_of each.z;
			float bed_s <- (chan where (each.grid_y >= ymax - 2)) mean_of each.z;
			float Sc <- max(1e-4, abs(bed_n - bed_s) / len);     // channel longitudinal slope
			float Qch <- (Q_peak > 0.0) ? Q_peak : (weir_C * 1000.0);
			float Qin0 <- Qch;
			float st <- bed_n;                                   // running stage (carried through gaps)
			stage_north <- 0.0; stage_south <- 0.0;
			loop gy from: ymin to: ymax {
				list<cell> row <- chan where (each.grid_y = gy);
				if (not empty(row)) {
					float Wrow <- max(dx, length(row) * dx);     // channel width at this station (m)
					float bedrow <- row mean_of each.z;
					float d <- stage_cal * ((Qch * river_n / (Wrow * sqrt(Sc))) ^ 0.6);   // normal depth
					st <- bedrow + d;
					if (stage_north = 0.0) { stage_north <- st; }
					stage_south <- st;
				}
				// local river stage along this station (carried st where the channel is absent)
				ask cell where (each.grid_y = gy) { rstage <- st; }
				// overtopping: non-breach dyke cells topped by the stage become spill points
				if (allow_overtopping) {
					ask cell where (each.is_dyke and not each.is_breach and each.grid_y = gy and (st > each.z + 0.1)) {
						is_breach <- true; breach_open <- 0.0;
					}
				}
				// spill through breach/overtop cells at this station (weir) -> floodplain source
				ask cell where (each.is_breach and each.grid_y = gy) {
					q_in <- myself.weir_C * myself.dx * (max(0.0, st - z) ^ 1.5);
				}
				float spill_raw <- (cell where (each.is_breach and each.grid_y = gy)) sum_of each.q_in;
				// CONTINUITY: cannot spill more than is left in the channel -> cap & rescale
				float spill <- min(spill_raw, Qch);
				if (spill_raw > spill and spill_raw > 0.0) {
					float fac <- spill / spill_raw;
					ask cell where (each.is_breach and each.grid_y = gy) { q_in <- q_in * fac; }
				}
				Qch <- Qch - spill;                              // river loses exactly the (capped) spill
			}
			river_stage <- max(stage_north, stage_south);
			breach_cells <- cell where each.is_breach;
			write "Channel routing: inlet stage=" + (stage_north with_precision 1) + " m, outlet="
				+ (stage_south with_precision 1) + " m;  Q in=" + (Qin0 with_precision 0)
				+ " -> out=" + (Qch with_precision 0) + " m3/s;  total spill to floodplain="
				+ ((breach_cells sum_of each.q_in) with_precision 0) + " m3/s";
		}
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
		// --- catchment shape b (Eq 4 & 6) and s_max, per cell (independent of v_mean) ---
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
		}
		// --- f_ss, q_c, h_af; iterate twice so s_ss uses the model's OWN mean velocity ---
		//     (the paper: "s_ss estimated using the duration of the event and the average
		//      flow velocities" — not an a-priori constant).
		loop pass from: 1 to: 2 {
			ask cell where each.active {
				sss  <- v_mean * (event_days * 86400.0);                 // event travel distance
				float sn <- (smax > 0.0) ? min(1.0, sss / smax) : 1.0;
				fss  <- (smax > sss) ? (sn ^ (1.0 + bshape)) : 1.0;      // Eq 12-13
				qc   <- fss * af;                                        // Eq 14 (af is already vol/time)
				float S <- max(slope, 1e-4);
				h_af <- (qc > 0.0) ? (((qc * manning) / (dx * sqrt(S))) ^ 0.6) : 0.0;   // Eq 2
			}
			// discharge-weighted mean Manning velocity  v = (1/n) h_af^(2/3) sqrt(S)
			float wsum <- (cell where each.active) sum_of (each.af);
			if (wsum > 0.0) {
				float vw <- (cell where each.active) sum_of
					(each.af * (1.0 / manning) * (each.h_af ^ 0.6667) * sqrt(max(each.slope, 1e-4)));
				v_mean <- max(0.05, min(5.0, vw / wsum));
			}
		}
		// representative event compensation (discharge-weighted) -> scales the sustained source
		float wsum2 <- (cell where each.active) sum_of (each.af);
		fss_event <- (wsum2 > 0.0) ? ((cell where each.active) sum_of (each.fss * each.af)) / wsum2 : 1.0;
		write "Compensation: mean v = " + (v_mean with_precision 2) + " m/s, s_ss = "
			+ ((v_mean * event_days * 86400.0) / 1000.0) with_precision 1
			+ " km, event f_ss = " + (fss_event with_precision 3);
	}

	// Floodplain starts DRY so the flood visibly SPREADS from the breaches and FILLS the
	// polder (basin-fill mode); the diffusive solver does the spreading. The river channel
	// + foreland are filled to the routed stage = the source reservoir behind the dykes.
	// (h_af from steps 1-3 is still computed as the steady-state estimate / for the velocity
	//  field, but is not used as the seed — the fill is what the user wants to watch.)
	action seed_inundation {
		ask cell where each.active { h <- 0.0; hmax <- 0.0; }
		ask cell where (each.is_river and not each.is_breach) {
			h <- max(0.0, rstage - z);
			hmax <- h;
		}
	}

	// =====================================================================
	// (4) DIFFUSIVE-WAVE REFINEMENT (Eq 15-19) — relax seeded h_af to steady state
	// Each cycle: hold the fixed-head BOUNDARIES (river/foreland reservoir + the breach
	// at the f_ss-scaled stage), then sweep wet cells (highest surface first) and
	// redistribute to lower neighbours with the flux-limited Manning flux. The breach
	// is a boundary, re-held each pass; from a DRY floodplain it spreads in and FILLS the
	// polder (domain edges are walls unless drain_edges) until the basin reaches the river
	// stage. relax_dt sets only the speed. Converges on FLOODPLAIN depth/extent stability.
	// =====================================================================
	reflex relax when: not converged and relax_iter < max_cycles {
		ask cell where each.active { h_prev <- h; }
		// fixed-head boundary: river channel + foreland at the stage (Dirichlet reservoir)
		ask cell where (each.is_river and not each.is_breach) { h <- max(0.0, rstage - z); }
		if (exact_eq19) {
			// ---- paper-literal mode: Eq 17 velocity (potential z + 1/2 h^2) + Eq 19
			//      (h^(5/3) discharge-conserving) two-phase Jacobi update. No flux limiter. ----
			loop times: iters_per_cycle {
				ask breach_cells { h <- fss_event * max(0.0, rstage - z); }
				ask cell where (each.active and not each.is_river and not (each.is_dyke and not each.is_breach))
					parallel: true { do flux_eq17; }
				ask cell where (each.active and not each.is_river and not each.is_dyke)
					parallel: true { do update_eq19; }
			}
		} else {
			// ---- stable mode: flux-limited Gauss-Seidel, highest water-surface first.
			//      Recompute the wet set EACH pass so the flood FRONT advances every iteration
			//      (filling a polder needs the front to travel far; once-per-cycle made it
			//      crawl). Non-breach dyke cells are walls (excluded here AND in redistribute). ----
			loop times: iters_per_cycle {
				// breach = fixed-head boundary (f_ss-scaled river stage), re-held each pass.
				ask breach_cells { h <- fss_event * max(0.0, rstage - z); if (h > hmax) { hmax <- h; } }
				list<cell> wet <- (cell where (each.active and each.h > min_depth
					and (each.is_breach or (not each.is_river and not each.is_dyke)))) sort_by (-(each.z + each.h));
				ask wet { do redistribute; }
			}
		}
		// convergence on the FLOODPLAIN depth (the held boundaries are excluded)
		max_dh <- (cell where (each.active and not each.is_breach and not each.is_river))
			max_of (abs(each.h - each.h_prev));
		relax_iter <- relax_iter + 1;
		// flood-extent convergence
		int wet_now <- cell count (each.active and each.h > min_depth and not each.is_river);
		if (wet_prev >= 0 and abs(wet_now - wet_prev) <= max(3, int(0.001 * wet_now))) {
			stable_count <- stable_count + 1;
		} else { stable_count <- 0; }
		wet_prev <- wet_now;
		// POLDER fill: converge when DEPTH settles (basin actually full). Conveyance: when
		// EXTENT settles (depth there creeps forever). The hard cap always finalizes too.
		bool done <- drain_edges ? (stable_count >= stable_needed) : (max_dh < conv_tol);
		if (relax_iter > 5 and (done or relax_iter >= max_cycles - 1)) { do finalize; }
		do recolor;
	}

	action recolor {
		ask cell where each.active parallel: true {
			if (h > min_depth) {
				if (flag_disconnect and not connected and not is_river) {
					color <- #red;                                          // wet but NO path back to a breach/river
				} else {
					color <- rgb(0, int(max(40, 170 - h * 18)), 255);       // blue, darker = deeper
				}
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
		do check_connectivity;       // diagnostic: are there any truly-floating wet cells?
		flag_disconnect <- true;
		do recolor;                  // paint disconnected wet cells red, if any
		do report;
		write "Settled after " + relax_iter + " relaxation iters (max depth change "
			+ (max_dh with_precision 4) + " m/iter). Paused.";
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
						// only a cell that is itself FLOODED can have an arrival time
						if (c.active and not c.is_breach and c.h > min_depth) {
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

	// DIAGNOSTIC: flood-fill connectivity from the breaches/river through wet cells.
	// Any wet floodplain cell with no such path is genuinely floating (a bug) -> red.
	action check_connectivity {
		ask cell where each.active { connected <- (is_river or is_breach); }
		loop r from: 0 to: sweep_rounds - 1 {
			loop dir from: 0 to: 3 {
				bool rx <- (dir = 1 or dir = 3);
				bool ry <- (dir = 2 or dir = 3);
				loop ix from: 0 to: nb_cols - 1 {
					int gx <- rx ? (nb_cols - 1 - ix) : ix;
					loop iy from: 0 to: nb_rows - 1 {
						int gy <- ry ? (nb_rows - 1 - iy) : iy;
						cell c <- cell grid_at {gx, gy};
						if (c.active and c.h > min_depth and not c.connected) {
							loop nb over: [c.nW, c.nE, c.nN, c.nS] {
								if (nb != nil and nb.active and nb.connected and nb.h > min_depth) { c.connected <- true; }
							}
						}
					}
				}
			}
		}
		int discon <- cell count (each.active and each.h > min_depth and not each.is_river and not each.connected);
		write "Connectivity check: " + discon + " wet cell(s) have NO path back to a breach/river"
			+ (discon = 0 ? "  -> ALL flood water is connected (thin links can hide it visually)." : "  -> shown in RED.");
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
	float FE; float FS;      // h^(5/3)-flux on the E and S faces (exact_eq19 mode)
	float rstage;            // local river water-surface elevation (Q-driven, N->S profile)

	bool active   <- true;
	bool is_edge  <- false;
	bool is_breach <- false;
	bool is_dyke   <- false;     // any dyke cell — a barrier for the flood-stage foreland fill
	bool is_channel <- false;    // the river CHANNEL (RedRiver polygon) — used to route Q -> stage
	bool is_river  <- false;     // channel + flood-stage foreland: fixed-stage reservoir (source side)
	bool is_sample <- false;
	bool connected <- false;     // diagnostic: has a wet path back to a breach / the river
	int  sample_id <- 0;
	float breach_open <- -1.0;   // breach opening time (s since day0)

	rgb base_color <- #gray;
	rgb color <- #gray;
	cell nE; cell nW; cell nS; cell nN;

	// Gauss-Seidel DIFFUSIVE-WAVE redistribution = the paper's Eq 17 + Eq 19.
	// Eq 17 velocity (Manning, driven by the WATER-SURFACE slope d(z+h)/dx — the
	//   dimensionally-correct reading of the pressure term; the printed d(1/2 h^2)/dx
	//   is a mis-render, since z is a slope and 1/2 h^2 is metres):
	//     u = (hflow^(2/3)/n) * sqrt(Sw),   Sw = (etaC - etaN)/dx
	// Eq 19 conserves the DISCHARGE quantity psi = h^(5/3) (not depth): we transfer psi,
	//   so h is NOT volume-conserved (the paper's explicit choice). Flux-limited by:
	//   (a) a quarter of the head difference (no overshoot) and (b) the water available.
	// A missing/inactive neighbour is free outfall (open boundary): outside eta = bed.
	action redistribute {
		if (h > min_depth) {
			loop nb over: [nE, nW, nS, nN] {
				// non-breach dyke cells are WALLS (dykes hold except at breaches); the domain
				// edge is a wall too unless drain_edges (POLDER fill vs open-conveyance).
				bool wall <- (nb != nil and nb.is_dyke and not nb.is_breach);
				bool isedge <- (nb = nil or not nb.active);
				if (not wall and (not isedge or drain_edges)) {
					float etaC <- z + h;
					float zN   <- isedge ? z : nb.z;                          // open edge: outside surface = bed
					float etaN <- isedge ? z : (nb.z + nb.h);
					if (etaC > etaN) {
						float hf <- etaC - max(z, zN);
						if (hf > min_depth) {
							float dhead <- etaC - etaN;                               // Eq-17 water-surface slope
							float u  <- (hf ^ 0.6667) / manning * sqrt(dhead / dx);   // Eq-17 velocity (m/s)
							float dpsi <- relax_dt * (h ^ 1.6667) * u / dx;           // Eq-19 flux of psi = h^(5/3)
							// flux limiter (the paper's flux-to-volume cap), via the equivalent depth:
							float dVcap <- min(0.25 * dhead, h - min_depth);          // no overshoot + available
							float dpsicap <- (h ^ 1.6667) - ((h - dVcap) ^ 1.6667);
							dpsi <- min(dpsi, dpsicap);
							if (dpsi > 0.0) {
								h <- (max(0.0, (h ^ 1.6667) - dpsi)) ^ 0.6;           // Eq-19: C loses psi
								if (not isedge) {
									nb.h <- ((nb.h ^ 1.6667) + dpsi) ^ 0.6;           // N gains the same psi
									if (nb.h > nb.hmax) { nb.hmax <- nb.h; }
								}
							}
						}
					}
				}
			}
			if (h > hmax) { hmax <- h; }
		}
	}

	// -------- paper-literal Eq 17 / Eq 19 update (optional 'exact_eq19' mode) --------
	// Eq 17 analytic velocity with the depth-integrated pressure potential phi = z + 1/2 h^2:
	//   u = sign(G) * (hflow^(2/3)/n) * sqrt(|G|),   G = -d(phi)/dx
	// Phase 1 stores the h^(5/3)-flux on the E and S faces (upwind). No flux limiter, so it
	// can go unstable on flat ground -> h is clamped to 50 m so instability is visible, not NaN.
	action flux_eq17 {
		// EAST face
		bool edgeE <- (nE = nil or not nE.active);
		bool wallE <- edgeE ? (not drain_edges) : (nE.is_river or (nE.is_dyke and not nE.is_breach));
		if (wallE) { FE <- 0.0; }
		else {
			float zN <- edgeE ? z : nE.z;          // open edge: outside surface = bed (free outfall)
			float hN <- edgeE ? 0.0 : nE.h;
			float G  <- -(((zN + 0.5 * hN * hN) - (z + 0.5 * h * h)) / dx);
			float hface <- max(z + h, zN + hN) - max(z, zN);
			if (hface > min_depth and abs(G) > 1e-12) {
				float u <- ((G > 0.0) ? 1.0 : -1.0) * (hface ^ 0.6667) / manning * sqrt(abs(G));
				FE <- ((u >= 0.0) ? (h ^ 1.6667) : (hN ^ 1.6667)) * u;
			} else { FE <- 0.0; }
		}
		// SOUTH face
		bool edgeS <- (nS = nil or not nS.active);
		bool wallS <- edgeS ? (not drain_edges) : (nS.is_river or (nS.is_dyke and not nS.is_breach));
		if (wallS) { FS <- 0.0; }
		else {
			float zN2 <- edgeS ? z : nS.z;
			float hN2 <- edgeS ? 0.0 : nS.h;
			float G2  <- -(((zN2 + 0.5 * hN2 * hN2) - (z + 0.5 * h * h)) / dx);
			float hface2 <- max(z + h, zN2 + hN2) - max(z, zN2);
			if (hface2 > min_depth and abs(G2) > 1e-12) {
				float u2 <- ((G2 > 0.0) ? 1.0 : -1.0) * (hface2 ^ 0.6667) / manning * sqrt(abs(G2));
				FS <- ((u2 >= 0.0) ? (h ^ 1.6667) : (hN2 ^ 1.6667)) * u2;
			} else { FS <- 0.0; }
		}
	}

	// Phase 2: Eq 19 conservative update of psi = h^(5/3), then h = psi^(3/5).
	action update_eq19 {
		float psi <- h ^ 1.6667;
		float Fw <- (nW != nil and nW.active) ? nW.FE : 0.0;   // flux from west face into me
		float Fn <- (nN != nil and nN.active) ? nN.FS : 0.0;   // flux from north face into me
		float net <- (Fw - FE) + (Fn - FS);
		float psinew <- max(0.0, psi + eq19_dt * net / dx);
		h <- min(50.0, psinew ^ 0.6);                          // clamp so instability is visible, not NaN
		if (h > hmax) { hmax <- h; }
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
	parameter "Manning n (floodplain)" var: manning      min: 0.02 max: 0.12 step: 0.005;
	parameter "Weir coefficient C"     var: weir_C       min: 1.0  max: 2.2  step: 0.05;
	parameter "Foreland reach (cells)" var: foreland_reach min: 0 max: 120 step: 2;
	parameter "Hydro-correct slope (m/cell)" var: delta_slope min: 0.0001 max: 0.01 step: 0.0005;
	parameter "Event duration (days)"  var: event_days   min: 1.0  max: 16.0 step: 1.0;
	parameter "Relaxation dt (s)"      var: relax_dt     min: 0.5  max: 20.0 step: 0.5;
	parameter "Drain at edges (conveyance vs polder fill)" var: drain_edges;
	parameter "Channel Manning n"        var: river_n     min: 0.02 max: 0.06 step: 0.005 category: "Discharge boundary";
	parameter "Stage calibration x"      var: stage_cal   min: 0.5  max: 1.5  step: 0.05  category: "Discharge boundary";
	parameter "Allow overtopping"        var: allow_overtopping category: "Discharge boundary";
	parameter "Stage fallback (no Q, m)" var: h_peak      min: 11.5 max: 16.0 step: 0.1   category: "Discharge boundary";
	parameter "Exact Eq-17/19 solver"  var: exact_eq19   category: "Compare: paper-literal (fragile)";
	parameter "  Eq-19 dt (s)"         var: eq19_dt      min: 0.05 max: 2.0  step: 0.05 category: "Compare: paper-literal (fragile)";

	output {
		display "Flood depth" type: 2d {
			grid cell;
			species river_area aspect: default;
			species dyke_seg aspect: default;
			species arrival_pt aspect: default;
			overlay position: {10, 10} size: {340 #px, 80 #px} background: #black transparency: 0.5 {
				draw "relaxation iter " + world.relax_iter + (world.converged ? "  (settled)" : "  (solving)")
					at: {15 #px, 24 #px} color: #white font: font("Helvetica", 13, #bold);
				draw "max change " + (world.max_dh with_precision 4) + " m   (steady-state solver, not a clock)"
					at: {15 #px, 48 #px} color: rgb(185, 185, 185) font: font("Helvetica", 11, #plain);
			}
		}
		monitor "Relaxation iter"  value: relax_iter;
		monitor "Extent-stable cycles" value: stable_count;
		monitor "Converged?"       value: converged;
		monitor "Max |dh| (m)"     value: max_dh with_precision 4;
		monitor "Wet cells"        value: cell count (each.active and each.h > min_depth);
		monitor "Max depth (m)"    value: (cell where each.active max_of each.hmax) with_precision 2;
		monitor "Breach Q (m3/s)"  value: (breach_cells sum_of each.q_in) with_precision 0;
	}
}
