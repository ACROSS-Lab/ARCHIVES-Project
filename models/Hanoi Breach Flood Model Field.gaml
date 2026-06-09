/**
* Hanoi 1926 breach flood — FIELD port of "Hanoi Breach Flood Model.gaml".
*
* Same physics (local-inertial, Bates 2010), same breach mechanism and Q->stage
* forcing, same dynamic frontier. The difference is purely computational: the state
* (z, h, qx, qy) lives in `field`s accessed by [col,row], and the hot loop iterates
* the frontier as parallel int lists (f_col / f_row) — no agents, no action dispatch
* per cell, no point allocation. The grid `cell` is kept ONLY for one-time geometric
* tagging in init (overlapping the river polygon / breach lines / points).
*
* Honest expectation: this is a constant-factor speedup (~2-4x) over the grid model,
* not order-of-magnitude — GAMA has no vectorised field-stencil arithmetic, so the
* timestep count (CFL-bound) is unchanged. For calibration, also run experiment
* `Fast` (no display) and/or a coarser DEM.
*/
model HanoiBreachFloodField

global {
	// ---------------------------------------------------------------- inputs
	file dem_file     <- grid_file("../includes/mnt-gz50.asc");
	file river_file   <- file("../includes/RedRiver1925.shp");
	file dykes_file   <- file("../includes/Dykes.shp");
	file points_file  <- file("../includes/5_arrival_time.shp");
	csv_file q_csv    <- csv_file("../includes/WaterDischarge.csv", ",", true);

	geometry shape <- envelope(dem_file);

	// ----------------------------------------------------------- state fields
	field Zf;                            // bed elevation (breaches lower it)
	field Hf;                            // water depth
	field QXf;                           // unit-width flux to EAST  (m2/s)
	field QYf;                           // unit-width flux to SOUTH (m2/s)
	field WSf;                           // water-surface elevation for the display mesh
	int   ncols;
	int   nrows;

	// ------------------------------------------------- hydraulic parameters
	float dx        <- 49.9736;
	float gravity   <- 9.81;
	float dt        <- 4.0;
	int   substeps  <- 120;
	float step      <- (dt * substeps) #s;
	float manning   <- 0.06;
	float min_depth <- 0.01;

	// -------------------------------------------------------- stage forcing
	list<float> Q_series;
	float Q_peak;
	float h_peak <- 13.3;
	float h_base <- 11.5;
	float river_stage <- h_base;
	float breach_floor <- 9.0;

	// ------------------------------------------------------------- time
	float sim_seconds  <- 0.0;
	float sim_days     <- 12.0;
	float total_seconds <- sim_days * 86400.0;
	int   day0_july <- 22;

	// ----------------------------------------- frontier + sets (flat indices)
	list<int>  f_col;  list<int>  f_row;        // the dynamic active set
	list<bool> inFront;                         // membership, linear index r*ncols+c
	list<bool> isRiverL;                        // river mask, linear index
	list<int>  src_col;  list<int>  src_row;    // river-edge (stage source) cells
	list<int>  brk_col;  list<int>  brk_row;    // breach cells
	list<float> brk_time; list<bool> brk_open;
	list<int>  smp_col;  list<int>  smp_row;  list<int> smp_id;  list<float> smp_arr;

	// ----------------------------------------- hoisted constants + monitors
	float manning_sq;  float gdt;  float gdt_nsq;  float pow73 <- 7.0 / 3.0;
	float max_depth <- 0.0;  int wet_count <- 0;

	init {
		// --- discharge ---
		matrix data <- matrix(q_csv);
		loop i from: 0 to: data.rows - 1 {
			string s <- string(data[1, i]);
			if (s != nil and s != "" and s != "m3/s") { add float(s) to: Q_series; }
		}
		Q_peak <- max(Q_series);

		// --- state fields ---
		Zf <- field(dem_file);
		ncols <- Zf.columns;
		nrows <- Zf.rows;
		Hf  <- field(ncols, nrows);
		QXf <- field(ncols, nrows);
		QYf <- field(ncols, nrows);
		WSf <- field(ncols, nrows);
		inFront  <- list_with(ncols * nrows, false);
		isRiverL <- list_with(ncols * nrows, false);
		write "Grid " + ncols + " x " + nrows + ", peak Q = " + Q_peak + " m3/s";

		// --- river polygon -> river mask (grid used only for geometry) ---
		create river_area from: river_file;
		ask river_area { ask cell overlapping self { is_river <- true; } }
		ask cell where each.is_river { isRiverL[grid_y * ncols + grid_x] <- true; }

		// --- breaches: cells under a BREAK=YES segment, with their open time ---
		create dyke_seg from: dykes_file with: [
			brk::string(read("BREAK")), dnum_s::string(read("DATE")), commune::string(read("Commune"))
		];
		ask dyke_seg where (each.brk = "YES") {
			int dday <- int(first(dnum_s split_with "-"));
			open_time <- (dday - myself.day0_july) * 86400.0;
			ask cell overlapping self {
				if (not is_breach) { is_breach <- true; breach_time <- myself.open_time; }
				else { breach_time <- min(breach_time, myself.open_time); }
			}
		}
		ask cell where each.is_breach {
			add grid_x to: brk_col;  add grid_y to: brk_row;
			add breach_time to: brk_time;  add false to: brk_open;
		}

		// --- observation points (record order = ids 1..5) ---
		list<arrival_pt> pts <- [];
		create arrival_pt from: points_file returns: pts;
		loop i from: 0 to: length(pts) - 1 {
			ask (pts at i) {
				sid <- i + 1;
				cell c <- first(cell overlapping self);
				if (c != nil) {
					add c.grid_x to: smp_col;  add c.grid_y to: smp_row;
					add sid to: smp_id;  add -1.0 to: smp_arr;
				}
			}
		}

		// --- source = river cell with at least one non-river neighbour ---
		ask cell where each.is_river {
			int c <- grid_x;  int r <- grid_y;  bool edge <- false;
			if (c < ncols - 1 and not isRiverL[r * ncols + c + 1])     { edge <- true; }
			if (c > 0         and not isRiverL[r * ncols + c - 1])     { edge <- true; }
			if (r < nrows - 1 and not isRiverL[(r + 1) * ncols + c])   { edge <- true; }
			if (r > 0         and not isRiverL[(r - 1) * ncols + c])   { edge <- true; }
			if (edge) { add c to: src_col;  add r to: src_row; }
		}
		write "Source cells: " + length(src_col) + ", breach cells: " + length(brk_col);

		// --- hoist constants ---
		manning_sq <- manning ^ 2;
		gdt <- gravity * dt;
		gdt_nsq <- gdt * manning_sq;

		// --- seed the frontier with the sources and their neighbours ---
		loop k from: 0 to: length(src_col) - 1 {
			int c <- src_col[k];  int r <- src_row[k];
			do try_add(c: c, r: r);
			if (c < ncols - 1) { do try_add(c: c + 1, r: r); }
			if (c > 0)         { do try_add(c: c - 1, r: r); }
			if (r < nrows - 1) { do try_add(c: c, r: r + 1); }
			if (r > 0)         { do try_add(c: c, r: r - 1); }
		}
	}

	// add a cell to the frontier if not already in it (O(1) via the membership flag)
	action try_add(int c, int r) {
		int lin <- r * ncols + c;
		if (not inFront[lin]) { inFront[lin] <- true; add c to: f_col; add r to: f_row; }
	}

	// ------------------------------------------------- per-cycle simulation
	reflex simulate when: sim_seconds < total_seconds {
		do update_stage;
		loop times: substeps {
			if (sim_seconds >= total_seconds) { break; }
			do open_due_breaches;
			do source_bc;          // impose stage before the flux pass
			do flux_pass;
			do depth_pass;
			do record_arrivals;
			sim_seconds <- sim_seconds + dt;
		}
		do update_display_field;
	}

	action update_stage {
		float dnum <- sim_seconds / 86400.0;
		int d0 <- min(int(dnum), length(Q_series) - 1);
		int d1 <- min(d0 + 1, length(Q_series) - 1);
		float frac <- dnum - d0;
		float Q <- Q_series[d0] * (1 - frac) + Q_series[d1] * frac;
		river_stage <- h_base + (h_peak - h_base) * (Q / Q_peak);
	}

	action open_due_breaches {
		loop k from: 0 to: length(brk_col) - 1 {
			if (not brk_open[k] and sim_seconds >= brk_time[k]) {
				int c <- brk_col[k];  int r <- brk_row[k];
				Zf[c, r] <- min(Zf[c, r], breach_floor);
				brk_open[k] <- true;
				do try_add(c: c, r: r);
				if (c < ncols - 1) { do try_add(c: c + 1, r: r); }
				if (c > 0)         { do try_add(c: c - 1, r: r); }
				if (r < nrows - 1) { do try_add(c: c, r: r + 1); }
				if (r > 0)         { do try_add(c: c, r: r - 1); }
			}
		}
	}

	action source_bc {
		loop k from: 0 to: length(src_col) - 1 {
			Hf[src_col[k], src_row[k]] <- max(0.0, river_stage - Zf[src_col[k], src_row[k]]);
		}
	}

	// PHASE 1: local-inertial momentum on each frontier cell's E and S faces
	action flux_pass {
		int n <- length(f_col);
		loop k from: 0 to: n - 1 {
			int c <- f_col[k];  int r <- f_row[k];  int lin <- r * ncols + c;
			float zc <- Zf[c, r];  float etaC <- zc + Hf[c, r];
			// east face
			if (c < ncols - 1 and not (isRiverL[lin] and isRiverL[lin + 1])) {
				float ze <- Zf[c + 1, r];
				float hflow <- max(etaC, ze + Hf[c + 1, r]) - max(zc, ze);
				if (hflow > min_depth) {
					float slope <- ((ze + Hf[c + 1, r]) - etaC) / dx;
					float q <- QXf[c, r];
					QXf[c, r] <- (q - gdt * hflow * slope) / (1 + gdt_nsq * abs(q) / (hflow ^ pow73));
				} else { QXf[c, r] <- 0.0; }
			} else { QXf[c, r] <- 0.0; }
			// south face
			if (r < nrows - 1 and not (isRiverL[lin] and isRiverL[lin + ncols])) {
				float zs <- Zf[c, r + 1];
				float hflow2 <- max(etaC, zs + Hf[c, r + 1]) - max(zc, zs);
				if (hflow2 > min_depth) {
					float slope2 <- ((zs + Hf[c, r + 1]) - etaC) / dx;
					float q2 <- QYf[c, r];
					QYf[c, r] <- (q2 - gdt * hflow2 * slope2) / (1 + gdt_nsq * abs(q2) / (hflow2 ^ pow73));
				} else { QYf[c, r] <- 0.0; }
			} else { QYf[c, r] <- 0.0; }
		}
	}

	// PHASE 2: mass balance; collect newly-wet cells, then grow the frontier
	action depth_pass {
		int n <- length(f_col);
		list<int> nc <- [];  list<int> nr <- [];
		loop k from: 0 to: n - 1 {
			int c <- f_col[k];  int r <- f_row[k];
			float inx <- (c > 0 ? QXf[c - 1, r] : 0.0) - QXf[c, r];
			float iny <- (r > 0 ? QYf[c, r - 1] : 0.0) - QYf[c, r];
			float h0 <- Hf[c, r];
			float hn <- max(0.0, h0 + dt * (inx + iny) / dx);
			Hf[c, r] <- hn;
			if (h0 <= min_depth and hn > min_depth) { add c to: nc;  add r to: nr; }
		}
		loop k from: 0 to: length(nc) - 1 {
			int c <- nc[k];  int r <- nr[k];
			if (c < ncols - 1) { do try_add(c: c + 1, r: r); }
			if (c > 0)         { do try_add(c: c - 1, r: r); }
			if (r < nrows - 1) { do try_add(c: c, r: r + 1); }
			if (r > 0)         { do try_add(c: c, r: r - 1); }
		}
	}

	action record_arrivals {
		loop k from: 0 to: length(smp_col) - 1 {
			if (smp_arr[k] < 0 and Hf[smp_col[k], smp_row[k]] > min_depth) {
				smp_arr[k] <- sim_seconds;
			}
		}
	}

	// build the water-surface field for the mesh; dry cells park just under terrain
	action update_display_field {
		do source_bc;
		max_depth <- 0.0;  wet_count <- 0;
		loop k from: 0 to: length(f_col) - 1 {
			int c <- f_col[k];  int r <- f_row[k];
			float hh <- Hf[c, r];
			if (hh > min_depth) {
				WSf[c, r] <- Zf[c, r] + hh;
				wet_count <- wet_count + 1;
				if (hh > max_depth) { max_depth <- hh; }
			} else {
				WSf[c, r] <- Zf[c, r] - 1.0;
			}
		}
	}

	// stop at the horizon, or early if every point is already wet
	reflex finish when: sim_seconds >= total_seconds or min(smp_arr) >= 0 {
		write "==================== ARRIVAL RESULTS ====================";
		loop k from: 0 to: length(smp_id) - 1 {
			if (smp_arr[k] < 0) {
				write "  pt" + smp_id[k] + " : NOT reached";
			} else {
				write "  pt" + smp_id[k] + " : " + (smp_arr[k] / 86400.0) with_precision 2 + " days";
			}
		}
		list<int> idxs <- [];
		loop k from: 0 to: length(smp_arr) - 1 { add k to: idxs; }
		list<int> reached <- (idxs where (smp_arr[each] >= 0)) sort_by (smp_arr[each]);
		write "  ORDER reached: " + (reached collect ("pt" + smp_id[each]));
		write "=========================================================";
		do pause;
	}
}

// grid kept ONLY for one-time geometric tagging in init (overlapping). Not used in
// the hot loop, not displayed.
grid cell file: dem_file neighbors: 4 frequency: 0
	use_regular_agents: false use_individual_shapes: false {
	bool  is_river  <- false;
	bool  is_breach <- false;
	float breach_time <- -1.0;
}

species river_area {
	aspect default { draw shape color: rgb(60, 110, 200) border: #blue; }
}

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

experiment HanoiBreachField type: gui {
	parameter "Manning n"             var: manning      min: 0.02 max: 0.12 step: 0.005;
	parameter "Peak stage h_peak (m)" var: h_peak       min: 11.5 max: 15.0 step: 0.1;
	parameter "Base stage h_base (m)" var: h_base       min: 9.0  max: 12.5 step: 0.1;
	parameter "Breach floor z (m)"    var: breach_floor min: 6.0  max: 12.0 step: 0.25;
	parameter "Hydraulic dt (s)"      var: dt           min: 1.0  max: 8.0  step: 0.5;
	parameter "Simulate days"         var: sim_days     min: 8.0  max: 16.0 step: 1.0;

	output {
		display "Flood" type: opengl {
			mesh Zf  scale: 1 grayscale: true smooth: false triangulation: true refresh: true;
			mesh WSf scale: 1 color: rgb(70, 130, 220, 200) smooth: false triangulation: true refresh: true;
			species dyke_seg aspect: default;
			species arrival_pt aspect: default;
		}
		monitor "Sim day"        value: sim_seconds / 86400.0 with_precision 2;
		monitor "River stage (m)" value: river_stage with_precision 2;
		monitor "Frontier size"  value: length(f_col);
		monitor "Wet cells"      value: wet_count;
		monitor "Max depth (m)"  value: max_depth with_precision 2;
	}
}

// No-display run for fast calibration; arrival results still print to the console.
experiment Fast type: gui {
	parameter "Manning n"             var: manning      min: 0.02 max: 0.12 step: 0.005;
	parameter "Peak stage h_peak (m)" var: h_peak       min: 11.5 max: 15.0 step: 0.1;
	parameter "Base stage h_base (m)" var: h_base       min: 9.0  max: 12.5 step: 0.1;
	parameter "Breach floor z (m)"    var: breach_floor min: 6.0  max: 12.0 step: 0.25;
	parameter "Simulate days"         var: sim_days     min: 8.0  max: 16.0 step: 1.0;
	output { }
}
