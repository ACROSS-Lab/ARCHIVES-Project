/**
* Hanoi 1926 Red River flood — breach-driven 2D local-inertial (LISFLOOD-FP) model
*
* Mechanism: the river dyke BREACHES on documented dates (Dykes.shp BREAK/DATE:
*   Gia-Quat 28-07, then Ai-Mo / Lam-Giu / Gia-Quat 29-07). Water at the river
*   stage pours through the opened breach cells into the protected floodplain and
*   spreads. We record the first-wet time at the 5 observation points
*   (5_arrival_time.shp) to compare arrival ORDER against the breach dates.
*
* Forcing: WaterDischarge.csv (daily Q, 1926) -> river stage via
*   h(t) = h_base + (h_peak - h_base) * Q(t)/Q_peak   (the DEM has no channel
*   bathymetry, so a Manning rating curve is not derivable; this is the agreed
*   data-driven mapping. h_peak / h_base are calibration knobs.)
*
* Numerics: explicit local-inertial scheme (Bates, Horritt & Fewtrell 2010) on the
*   DEM grid. Cells do NOT self-schedule (frequency: 0); the global reflex runs the
*   two-phase (flux, then depth) update so all fluxes use a consistent state.
*/
model HanoiBreachFlood

global {
	// ---------------------------------------------------------------- inputs
	file dem_file     <- grid_file("../includes/mnt-gz50.asc");
	file river_file   <- file("../includes/RedRiver1925.shp");
	file dykes_file   <- file("../includes/Dykes.shp");
	file points_file  <- file("../includes/5_arrival_time.shp");
	csv_file q_csv    <- csv_file("../includes/WaterDischarge.csv", ",", true);

	geometry shape <- envelope(dem_file);

	// ------------------------------------------------- hydraulic parameters
	float dx        <- 49.9736;          // cell size (m), from the DEM header
	float gravity   <- 9.81;
	float dt        <- 4.0;              // hydraulic timestep (s) — CFL ~ dx/sqrt(g*h)
	int   substeps  <- 60;               // hydraulic steps per displayed cycle
	float step      <- (dt * substeps) #s;
	float manning   <- 0.06;             // floodplain roughness (calibration knob)
	float min_depth <- 0.01;             // wetting / dry threshold (m)

	// -------------------------------------------------------- stage forcing
	list<float> Q_series;
	float Q_peak;
	float h_peak <- 13.3;                // 1926 Hanoi peak stage (m) — calibration knob
	float h_base <- 11.5;                // danger-level baseline (m) — calibration knob
	float river_stage <- h_base;

	// --------------------------------------------------------- breach setup
	float breach_floor <- 9.0;           // crest elevation a breach is cut down to (m)

	// ------------------------------------------------------------- bookkeeping
	int   nb_cols;
	int   nb_rows;
	float sim_seconds  <- 0.0;           // elapsed hydraulic time (sub-step granularity)
	float total_seconds <- 16.0 * 86400.0;
	int   day0_july <- 22;               // discharge series starts 22 July 1926

	init {
		// --- discharge series (skip header, guard non-numeric) ---
		matrix data <- matrix(q_csv);
		loop i from: 0 to: data.rows - 1 {
			string s <- string(data[1, i]);
			if (s != nil and s != "" and s != "m3/s") {
				Q_series <+ float(s);
			}
		}
		Q_peak <- max(Q_series);
		write "Loaded " + length(Q_series) + " daily discharges, peak = " + Q_peak + " m3/s";

		// --- grid dimensions + per-cell init ---
		nb_cols <- (cell max_of each.grid_x) + 1;
		nb_rows <- (cell max_of each.grid_y) + 1;
		ask cell {
			z <- grid_value;
			active <- z > -1000.0;       // guard NODATA (-9999); this DEM has none
			h <- 0.0; qx <- 0.0; qy <- 0.0;
			if (grid_x < nb_cols - 1) { nE <- cell grid_at {grid_x + 1, grid_y}; }
			if (grid_x > 0)           { nW <- cell grid_at {grid_x - 1, grid_y}; }
			if (grid_y < nb_rows - 1) { nS <- cell grid_at {grid_x, grid_y + 1}; }
			if (grid_y > 0)           { nN <- cell grid_at {grid_x, grid_y - 1}; }
		}

		// --- river polygon -> stage boundary cells ---
		create river_area from: river_file;
		ask river_area { ask cell overlapping self { is_river <- true; } }
		write "River boundary cells: " + length(cell where each.is_river);

		// --- dyke segments; breached ones become datable openings ---
		create dyke_seg from: dykes_file with: [
			brk::string(read("BREAK")), dnum_s::string(read("DATE")), commune::string(read("Commune"))
		];
		ask dyke_seg where (each.brk = "YES") {
			int dday <- int(first(dnum_s split_with "-"));   // "28-07" -> 28
			open_time <- (dday - myself.day0_july) * 86400.0;
			ask cell overlapping self {
				is_breach <- true;
				orig_z <- z;
				breach_time <- (breach_time < 0) ? myself.open_time : min(breach_time, myself.open_time);
			}
		}
		write "Breach cells: " + length(cell where each.is_breach)
			+ " (earliest opens day " + ((cell where each.is_breach) min_of each.breach_time) / 86400.0 + ")";

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
	}

	// ------------------------------------------------- per-cycle simulation
	reflex simulate when: sim_seconds < total_seconds {
		loop times: substeps {
			if (sim_seconds >= total_seconds) { break; }
			do update_stage;
			do open_due_breaches;
			ask cell where each.active { do compute_flux; }
			ask cell where each.active { do update_depth; }
			ask cell where each.is_river { h <- max(0.0, river_stage - z); }   // Dirichlet reservoir
			ask cell where each.is_sample {
				if (h > min_depth and arrival < 0) { arrival <- myself.sim_seconds; }
			}
			sim_seconds <- sim_seconds + dt;
		}
		do recolor;
	}

	action update_stage {
		float dnum <- sim_seconds / 86400.0;                 // 0-based day index
		int d0 <- min(int(dnum), length(Q_series) - 1);
		int d1 <- min(d0 + 1, length(Q_series) - 1);
		float frac <- dnum - d0;
		float Q <- Q_series[d0] * (1 - frac) + Q_series[d1] * frac;
		river_stage <- h_base + (h_peak - h_base) * (Q / Q_peak);
	}

	action open_due_breaches {
		ask cell where (each.is_breach and not each.opened and sim_seconds >= each.breach_time) {
			z <- min(z, myself.breach_floor);
			opened <- true;
		}
	}

	action recolor {
		ask cell where each.active {
			if (h > min_depth) {
				color <- rgb(0, int(max(40, 170 - h * 18)), 255);   // blue, darker = deeper
			} else {
				int gv <- int(min(235, max(60, (z - 4) * 11)));     // grey terrain by elevation
				color <- rgb(gv, gv, gv);
			}
		}
	}

	reflex finish when: sim_seconds >= total_seconds {
		write "==================== ARRIVAL RESULTS ====================";
		list<cell> samples <- (cell where each.is_sample) sort_by each.sample_id;
		ask samples {
			if (arrival < 0) {
				write "  pt" + sample_id + " : NOT reached (z=" + z with_precision 1 + ")";
			} else {
				write "  pt" + sample_id + " : " + (arrival / 86400.0) with_precision 2
					+ " days  (z=" + z with_precision 1 + ")";
			}
		}
		list<cell> reached <- (samples where (each.arrival >= 0)) sort_by each.arrival;
		write "  ORDER reached: " + (reached collect ("pt" + each.sample_id));
		write "=========================================================";
		do pause;
	}
}

// ====================================================================== grid
// Terrain + flood state. Elevation from the DEM via grid_value. frequency:0 so
// the cells are driven explicitly by the global reflex (no per-cell scheduling).
grid cell file: dem_file neighbors: 4 frequency: 0
	use_regular_agents: false use_individual_shapes: false {
	float z;                 // bed elevation (m)
	float h;                 // water depth (m)
	float qx;                // unit-width flux to EAST neighbour (m2/s)
	float qy;                // unit-width flux to SOUTH neighbour (m2/s)
	bool  active   <- true;
	bool  is_river <- false;
	bool  is_breach <- false;
	bool  opened   <- false;
	bool  is_sample <- false;
	int   sample_id <- 0;
	float orig_z;
	float breach_time <- -1.0;
	float arrival <- -1.0;
	rgb   color <- #gray;
	cell  nE; cell nW; cell nS; cell nN;

	// local-inertial momentum update on the two outgoing faces (E and S)
	action compute_flux {
		if (nE != nil and nE.active) {
			float etaC <- z + h;  float etaE <- nE.z + nE.h;
			float hflow <- max(etaC, etaE) - max(z, nE.z);
			if (hflow > min_depth) {
				float slope <- (etaE - etaC) / dx;
				qx <- (qx - gravity * hflow * dt * slope)
					/ (1 + gravity * dt * manning ^ 2 * abs(qx) / (hflow ^ (7.0 / 3.0)));
			} else { qx <- 0.0; }
		} else { qx <- 0.0; }

		if (nS != nil and nS.active) {
			float etaC <- z + h;  float etaS <- nS.z + nS.h;
			float hflow <- max(etaC, etaS) - max(z, nS.z);
			if (hflow > min_depth) {
				float slope <- (etaS - etaC) / dx;
				qy <- (qy - gravity * hflow * dt * slope)
					/ (1 + gravity * dt * manning ^ 2 * abs(qy) / (hflow ^ (7.0 / 3.0)));
			} else { qy <- 0.0; }
		} else { qy <- 0.0; }
	}

	// mass balance: net of (west-in - east-out) + (north-in - south-out)
	action update_depth {
		float inx <- (nW != nil ? nW.qx : 0.0) - qx;
		float iny <- (nN != nil ? nN.qy : 0.0) - qy;
		h <- max(0.0, h + dt * (inx + iny) / dx);
	}
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

experiment HanoiBreach type: gui {
	parameter "Manning n"            var: manning      min: 0.02 max: 0.12 step: 0.005;
	parameter "Peak stage h_peak (m)" var: h_peak       min: 11.5 max: 15.0 step: 0.1;
	parameter "Base stage h_base (m)" var: h_base       min: 9.0  max: 12.5 step: 0.1;
	parameter "Breach floor z (m)"    var: breach_floor min: 6.0  max: 12.0 step: 0.25;
	parameter "Hydraulic dt (s)"      var: dt           min: 1.0  max: 8.0  step: 0.5;

	output {
		display "Flood" type: 2d {
			grid cell;
			species river_area aspect: default transparency: 0.4;
			species dyke_seg aspect: default;
			species arrival_pt aspect: default;
			overlay position: {10, 10} size: {220 #px, 70 #px} background: #black transparency: 0.4 {
				draw "day " + (world.sim_seconds / 86400.0) with_precision 2
					+ "   stage " + world.river_stage with_precision 2 + " m"
					at: {20 #px, 30 #px} color: #white font: font("Helvetica", 14, #bold);
			}
		}
		monitor "Sim day"        value: sim_seconds / 86400.0 with_precision 2;
		monitor "River stage (m)" value: river_stage with_precision 2;
		monitor "Wet cells"      value: length(cell where (each.h > min_depth));
		monitor "Max depth (m)"  value: (cell max_of each.h) with_precision 2;
		monitor "Breaches open"  value: length(cell where each.opened);
	}
}
