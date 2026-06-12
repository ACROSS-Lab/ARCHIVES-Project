/**
* Hanoi 1926 Red River flood — hourly stage-driven spreading model
*
* One simulation step = ONE HOUR, 20 July 00:00 to 15 August 23:00 1926
* (648 steps). Each hour the Red River is set to the reconstructed gauge
* stage at Hanoi (RedRiverStage1926_hourly.csv — daily readings digitized
* from Gourou, "Le Tonkin" fig. 9, reproduced in CIA-RDP79T01018A000900010001-4
* and CIA-RDP79T01019A000400120001-3, interpolated to hourly with a
* monotone cubic; record peak 11.93 m on 30 July).
*
* Spreading rule (hourly fill-and-spill, volume-conserving):
*   - the POOL = the river plus contiguous water already equalized with it;
*     the pool follows the stage instantly (it is one body of water);
*   - NOTHING dry floods instantly: every dry cell on the pool margin below
*     the stage is a spill GATE into its basin, and basins fill lowest-first
*     with one hour of broad-crested-weir flow (Q = C * b * head^1.5) per
*     gate — over a ground sill (b = 50 m), through an open breach
*     (b = breach_cell_width), or over an overtopped crest (b = 50 m);
*   - a dyke segment with BREAK=YES opens at its DATE (Gia-Quat 28-07,
*     Ai-Mo / Lam-Giu / Gia-Quat 29-07), its cells cut to the breach floor;
*     per the Jan-1966 CIA report two of the three breaches were resealed
*     (8 and 12 August) — they become walls again; Lam-Giu stays open;
*   - unbroken dykes block until the stage exceeds the segment crest
*     (= highest DEM value along the segment);
*   - a basin still below the river level keeps being weir-fed until it
*     catches up, then merges with the pool; when its connection is lost
*     (stage below the sill, breach resealed) it stays ponded — the trapped
*     casier water described in the reports.
*/
model HanoiHourlyStageFlood1926

global {
	// ---------------------------------------------------------------- inputs
	file dem_file       <- grid_file("../includes/mnt-gz25.asc");
	file river_file     <- file("../includes/RedRiver1925.shp");
	file lakes_file     <- file("../includes/Lakes1925.shp");
	file buildings_file <- file("../includes/Buildings1925.shp");
	file dykes_file     <- file("../includes/Dykes.shp");
	csv_file stage_csv  <- csv_file("../includes/RedRiverStage1926_hourly.csv", ",", true);

	geometry shape <- envelope(dem_file);
	date starting_date <- date([1926, 7, 20, 0, 0, 0]);
	float step <- 1 #h;

	// -------------------------------------------------------- stage forcing
	list<float> stage_series <- [];      // one value per hour, cycle-indexed
	int   n_stage;
	float river_stage <- 0.0;            // raw gauge stage (Hanoi gauge datum)
	float stage_offset <- 0.0;           // gauge-to-DEM datum correction (m): the water
	                                     // surface compared to the DEM everywhere is
	                                     // river_stage + stage_offset

	// --------------------------------------------------------- breach setup
	float breach_floor <- 8.5;           // crest cut-down once breached (m)
	int   breach_hour  <- 6;             // hour of the DATE day a breach opens
	int   day0         <- 20;            // simulation starts 20 July
	// communes whose breach was resealed, with the August day of closure
	map<string, int> reseal_day <- ["Gia-Quat"::8, "Ai-Mo"::12];   // Lam-Giu stays open

	// ------------------------------------------------------- dyke overtopping
	// A dyke segment's crest = highest DEM value along it (single cells under
	// the line read low — rasterization artifacts, not real gaps). Water flows
	// over an unbroken dyke when the stage exceeds that crest.
	bool  allow_overtop <- true;
	float crest_trim    <- 0.0;          // lower every crest by this much (what-if runs)

	// ------------------------------------------------------ breach hydraulics
	// A breach passes water like a broad-crested weir, Q = C * b * head^1.5,
	// so a casier fills lowest-first over hours-days (CIA Jan-66: 250,000
	// acres through three breaches) instead of equalizing in one step.
	// CALIBRATED 2026-06-12 against the historical anchors (CIA reports, Ross
	// 2023, Verheij & Rohde 2006): offset 0 (1926 peak 11.93 m sits just under
	// the crests, city saved, breaches do the damage); casier inundation
	// builds over ~3 days from the 28 Jul breach (peak area 31 Jul-1 Aug,
	// lagging the stage peak by ~17 h); ponding persists for weeks; total
	// breach width 24 cells x 3 m = 72 m across the three openings.
	float weir_coef <- 1.7;              // broad-crested weir coefficient (metric)
	float breach_cell_width <- 3.0;      // effective flow width per breach cell (m)
	// the weir caps VOLUME per hour; these two cap how fast that volume may
	// turn into AREA (friction-limited sheet flow over the flat plain):
	// stored volume meeting a flat terrace must spread gradually, not at once
	float front_speed <- 300.0;          // m/h lateral advance into dry ground
	float max_spread_km2h <- 0.3;        // new land claimed per hour, ALL basins together (km2)
	float spread_budget <- 0.0;          // cells of new land left this hour (internal)
	bool  instant_fill <- false;         // init only: seed the equilibrium state

	// ------------------------------------------------------------ thresholds
	float river_bed_z <- 6.5;            // channel burned into the DEM (m): the DEM has
	                                     // no bathymetry, so river cells are lowered to
	                                     // this bed so the whole 1925 channel is wet at
	                                     // the 7.02 m starting stage (display only —
	                                     // river cells are stage sources regardless)
	float flood_depth_min <- 0.05;       // depth (m) below which a cell shows dry
	float bldg_flood_depth <- 0.10;      // depth (m) flagging a building as flooded
	float lake_init_depth <- 0.5;        // initial ponded depth in the 1925 lakes (m)

	// ------------------------------------------------------------ bookkeeping
	list<cell> all_cells   <- [];
	list<cell> river_cells <- [];
	list<cell> breach_cells <- [];
	float dx;                            // cell size (m), read from the DEM
	float cell_area;                     // m2
	float flooded_km2 <- 0.0;
	float max_flooded_km2 <- 0.0;
	int   flooded_bldgs <- 0;

	init {
		// --- hourly stage series ---
		matrix data <- matrix(stage_csv);
		loop i from: 0 to: data.rows - 1 {
			add float(data[1, i]) to: stage_series;
		}
		n_stage <- length(stage_series);
		write "Stage series: " + n_stage + " hourly values, peak "
			+ max(stage_series) + " m (30 July 1926)";

		// --- grid wiring (cell size from the loaded DEM, never hard-coded) ---
		int nb_cols <- (cell max_of each.grid_x) + 1;
		int nb_rows <- (cell max_of each.grid_y) + 1;
		dx <- shape.width / nb_cols;
		cell_area <- dx * dx;
		write "DEM: " + nb_cols + "x" + nb_rows + " cells of " + (dx with_precision 2) + " m";
		ask cell {
			z <- grid_value;
			orig_z <- z;
			active <- z > -1000.0;       // guard NODATA
			if (grid_x < nb_cols - 1) { add (cell grid_at {grid_x + 1, grid_y}) to: neigh; }
			if (grid_x > 0)           { add (cell grid_at {grid_x - 1, grid_y}) to: neigh; }
			if (grid_y < nb_rows - 1) { add (cell grid_at {grid_x, grid_y + 1}) to: neigh; }
			if (grid_y > 0)           { add (cell grid_at {grid_x, grid_y - 1}) to: neigh; }
		}
		all_cells <- cell where each.active;

		// --- dykes FIRST so river marking can exclude dyke cells ---
		create dyke_seg from: dykes_file with: [
			brk::string(read("BREAK")), date_s::string(read("DATE")), commune::string(read("Commune"))
		];
		ask dyke_seg {
			if (brk = "YES") {
				int dday <- int(first(date_s split_with "-"));            // "28-07" -> 28
				open_cycle <- (dday - day0) * 24 + breach_hour;
				if (reseal_day contains_key commune) {
					// 11 days of July remain after the 20th; reseal at noon
					close_cycle <- (11 + reseal_day[commune]) * 24 + 12;
				}
			}
			// crest from the untouched DEM (grid_value), then raise the whole
			// segment to it so low artifact cells cannot leak.
			// The tiny buffer makes the rasterized line 4-CONNECTED: without
			// it, the chain touches diagonally and 4-neighbour flow slips
			// between two corner-touching dyke cells (a sieve at 25 m).
			list<cell> dcells <- cell overlapping (shape + (dx * 0.1));
			crest <- (dcells max_of each.grid_value) - crest_trim;
			ask dcells {
				is_dyke <- true;
				z <- max(z, myself.crest);
				orig_z <- z;
				if (myself.brk = "YES") {
					is_breach <- true;
					open_cycle <- (open_cycle < 0) ? myself.open_cycle : min(open_cycle, myself.open_cycle);
					if (myself.close_cycle < 0) { never_close <- true; }
					else { close_cycle <- max(close_cycle, myself.close_cycle); }
				}
			}
		}
		ask cell where each.never_close { close_cycle <- -1; }
		breach_cells <- all_cells where each.is_breach;
		write "Dyke cells: " + (all_cells count each.is_dyke)
			+ "  breach cells: " + length(breach_cells)
			+ " (open cycles " + ((breach_cells collect each.open_cycle)) + ")";

		// --- river: the stage reservoir (dyke cells excluded so walls hold) ---
		create river_area from: river_file;
		ask river_area {
			ask (cell overlapping self) where (each.active and not each.is_dyke) {
				is_river <- true;
			}
		}
		river_cells <- all_cells where each.is_river;
		ask river_cells {                // burn the channel (no bathymetry in the DEM)
			z <- min(z, river_bed_z);
			orig_z <- z;
		}
		write "River cells: " + length(river_cells);

		// --- lakes: initialized AT REST — one flat surface per lake holding
		//     the nominal volume, capped at the spill rim so an undisturbed
		//     lake never moves or overflows on its own ---
		create lake from: lakes_file;
		loop lk over: list(lake) {
			list<cell> lcells <- (cell overlapping lk)
				where (each.active and not each.is_dyke and not each.is_river);
			if (not empty(lcells)) {
				ask lcells { is_lake <- true; }
				// spill rim: lowest ground on the lake's outer border
				list<cell> rim <- remove_duplicates(lcells accumulate each.neigh)
					where (each.active and not each.is_lake);
				float spill <- empty(rim) ? 1e9 : (rim min_of each.z);
				// flat surface storing ~lake_init_depth per cell, never above the rim
				list<cell> srt <- lcells sort_by each.z;
				float lvl <- min(fill_level(srt, length(lcells) * lake_init_depth * cell_area),
					spill);
				ask lcells where (each.z < lvl) {
					wet <- true;
					level <- lvl;
					h <- lvl - z;
				}
			}
		}

		// --- buildings: remember their cells once ---
		create building from: buildings_file;
		ask building {
			my_cells <- cell overlapping self;
		}
		write "Buildings: " + length(building);

		// --- seed the 20 Jul 00:00 equilibrium (instant: it is the *initial
		//     condition*, not a propagating flood wave) ---
		river_stage <- first(stage_series);
		instant_fill <- true;
		do spread_water;
		instant_fill <- false;

		// --- static terrain colouring (water recoloured per cycle) ---
		ask all_cells { do recolor_self; }
	}

	// ------------------------------------------------- one cycle = one hour
	reflex hourly when: cycle < n_stage {
		river_stage <- stage_series[cycle];
		do update_breaches;
		do spread_water;
		do update_buildings;
		ask all_cells { do recolor_self; }

		flooded_km2 <- (all_cells count (each.h > flood_depth_min
			and not each.is_river and not each.is_lake)) * cell_area / 1e6;
		max_flooded_km2 <- max(max_flooded_km2, flooded_km2);
	}

	// open breaches on their DATE; reseal the two repaired ones in August
	action update_breaches {
		ask breach_cells {
			bool should_be_open <- (cycle >= open_cycle)
				and (close_cycle < 0 or cycle < close_cycle);
			if (should_be_open and not is_open) {
				is_open <- true;
				z <- min(orig_z, breach_floor);
				write string(current_date, "dd MMM HH:mm") + "  BREACH OPEN  (crest "
					+ orig_z with_precision 1 + " m -> " + z with_precision 1 + " m)";
			}
			if (not should_be_open and is_open) {
				is_open <- false;
				z <- orig_z;                 // repaired: wall again at full crest
				wet <- false; h <- 0.0; level <- 0.0;
				write string(current_date, "dd MMM HH:mm") + "  BREACH RESEALED";
			}
		}
	}

	// Volume-conserving fill-and-spill. The RIVER is the only instantly-driven
	// water (it is the gauge boundary condition). Every body of water on land
	// is a persistent REGION holding its own volume: it exchanges with the
	// river through weir gates (bank, breach, overtopped crest) in BOTH
	// directions, and each hour it redistributes its stored volume lowest-
	// first — so water delivered during the rise KEEPS SPREADING after the
	// river peak, and drains only as fast as its gates allow.
	action spread_water {
		float W <- river_stage + stage_offset;   // datum-corrected water surface
		spread_budget <- max_spread_km2h * 1e6 / cell_area;   // global hourly area budget
		ask all_cells { connected <- false; rid <- -1; is_gate <- false; adv <- 0; }

		ask river_cells {
			connected <- true;
			wet <- true;
			level <- W;
			h <- max(0.0, W - z);
		}

		if (instant_fill) {
			// initialization only: classic connected bathtub = the equilibrium
			list<cell> frontier <- list(river_cells);
			loop while: not empty(frontier) {
				list<cell> nxt <- [];
				ask frontier {
					loop nb over: neigh {
						if (nb.active and not nb.connected and not nb.is_dyke and nb.z < W) {
							nb.connected <- true;
							add nb to: nxt;
						}
					}
				}
				frontier <- nxt;
			}
			ask all_cells where (each.connected and not each.is_river) {
				wet <- true;
				level <- W;
				h <- W - z;
			}
		} else {
			// ---- gates: every interface where water can pass this hour
			list<cell> gates <- [];
			// river-bank gates: land beside the river that should fill
			// (below stage) or drain (standing above stage)
			ask river_cells {
				loop nb over: neigh {
					if (nb.active and not nb.connected and not nb.is_gate and not nb.is_dyke
						and ((nb.z < W and (not nb.wet or nb.level < W - 0.01))
							or (nb.wet and nb.level > W + 0.01))) {
						nb.is_gate <- true;
						add nb to: gates;
					}
				}
			}
			// dyke gates: open breach / overtoppable crest with water against
			// it on either side (river at stage, or a region standing higher)
			ask all_cells where (each.is_dyke and not each.is_gate
				and (each.is_open or allow_overtop)) {
				bool active_gate <- false;
				loop nb over: neigh {
					if ((nb.connected and z < W) or (nb.wet and nb.level > z)) {
						active_gate <- true;
					}
				}
				if (active_gate) {
					is_gate <- true;
					add self to: gates;
				}
			}
			do fill_regions(W, gates);
		}

		// ---- ponded leftovers (wet, in no region, no gate) keep their level
		ask all_cells where (each.wet and not each.connected) {
			if (level <= z or is_dyke) { wet <- false; h <- 0.0; level <- 0.0; }
			else { h <- level - z; }
		}
	}

	// Build the regions: first those reached by gates (they receive/lose one
	// hour of weir flow), then every remaining wet body (no gate this hour —
	// volume unchanged, but it still redistributes and keeps creeping).
	// A dyke gate NEVER merges its two sides: the side with the higher water
	// surface is the source (it can only LOSE water through the weir), the
	// lower side is fed — each is its own region with its own volume.
	action fill_regions (float w, list<cell> gates) {
		int cur <- 0;
		list<cell> todo <- copy(gates);
		loop while: not empty(todo) {
			cell g0 <- first(todo);
			todo >- g0;

			if (g0.is_dyke) {
				// source level = highest water surface standing at the gate
				float src <- 0.0;
				loop nb over: g0.neigh {
					if (nb.connected) { src <- max(src, w); }
					else if (nb.wet) { src <- max(src, nb.level); }
				}
				list<cell> lo_seeds <- [];
				list<cell> hi_seeds <- [];
				loop nb over: g0.neigh {
					if (nb.active and not nb.is_dyke and not nb.connected and nb.rid < 0
						and (nb.z < w or (nb.wet and nb.level > nb.z))) {
						float surf <- nb.wet ? nb.level : nb.z;
						if (surf < src - 0.01) {
							nb.rid <- cur;
							nb.adv <- nb.wet ? 0 : 1;
							add nb to: lo_seeds;
						} else if (nb.wet) {
							nb.rid <- cur + 1;
							nb.adv <- 0;
							add nb to: hi_seeds;
						}
					}
				}
				if (not empty(lo_seeds)) {       // receiving side: gains via weir
					list<cell> region <- grow_region(w, lo_seeds, cur);
					list<cell> feeders <- [g0];
					loop c over: copy(todo) {
						if (c.rid = cur or (c.neigh one_matches (each.rid = cur))) {
							add c to: feeders;
							todo >- c;
						}
					}
					do apply_fill(w, region, feeders);
				}
				if (not empty(hi_seeds)) {       // source side: may drain via weir
					list<cell> region <- grow_region(w, hi_seeds, cur + 1);
					do apply_fill(w, region, [g0]);
				}
				cur <- cur + 2;
			} else {
				g0.rid <- cur;
				g0.adv <- g0.wet ? 0 : 1;
				list<cell> region <- grow_region(w, [g0], cur);
				list<cell> feeders <- [g0];
				loop c over: copy(todo) {
					if (c.rid = cur or (c.neigh one_matches (each.rid = cur))) {
						add c to: feeders;
						todo >- c;
					}
				}
				do apply_fill(w, region, feeders);
				cur <- cur + 1;
			}
		}

		// gateless wet bodies: redistribute in place (stored volume keeps
		// spreading outward after the river has dropped — no inflow needed).
		// No river drive here, and LAKES NEVER SEED a region: they start in
		// equilibrium and only external floodwater may move them (adjacent
		// ponds are separated by sub-cell bunds, so they must not equalize
		// with each other on their own).
		list<cell> strays <- all_cells where (each.wet and not each.connected
			and each.rid < 0 and not each.is_dyke and not each.is_lake);
		loop while: not empty(strays) {
			cell s0 <- first(strays);
			s0.rid <- cur;
			s0.adv <- 0;
			list<cell> region <- grow_region(-1e9, [s0], cur);
			do apply_fill(w, region, []);
			cur <- cur + 1;
			strays <- strays where (each.rid < 0);
		}
	}

	// BFS growth of one region: through wet cells freely; into dry cells if
	// they are below the DRIVING head (river stage for gate-fed regions;
	// nothing for isolated ponds/lakes) OR below the water standing right
	// next to them, within the hourly front budget. Never across dykes,
	// never into the river, never into another region (rid).
	list<cell> grow_region (float drive, list<cell> seeds, int cur) {
		int max_adv <- max(1, int(front_speed / dx));
		list<cell> region <- copy(seeds);
		list<cell> frontier <- copy(seeds);
		loop while: not empty(frontier) {
			list<cell> nxt <- [];
			ask frontier {
				loop nb over: neigh {
					if (nb.active and not nb.is_dyke and not nb.connected and nb.rid < 0
						and (nb.z < drive
							or (nb.wet and nb.level > nb.z)
							or (wet and nb.z < level))) {
						// flood-water cells conduct freely (one body), but LAKE
						// water consumes front budget: a lake is absorbed
						// slice-by-slice over hours, so its deep storage cannot
						// swallow the region's level (and dry its margins) in
						// a single step
						int nadv <- (nb.wet and not nb.is_lake) ? 0 : adv + 1;
						if (nadv <= max_adv) {
							nb.rid <- cur;
							nb.adv <- nadv;
							add nb to: nxt;
						}
					}
				}
			}
			region <- region + nxt;
			frontier <- nxt;
		}
		return region;
	}

	// One hour of weir exchange for a region, then lowest-first
	// redistribution of its (conserved) volume, with the area limiter.
	action apply_fill (float w, list<cell> region, list<cell> feeders) {
		// present level of the water that is ACTUALLY there (wet cells only:
		// dry candidate cells must not distort it)
		list<cell> wets <- region where each.wet;
		float v <- sum(wets collect max(0.0, each.level - each.z)) * cell_area;
		float l_old <- empty(wets) ? (region min_of each.z)
			: fill_level(wets sort_by each.z, v);

		float vin <- 0.0;
		loop c over: feeders {
			// breach width is calibrated per 50 m of dyke (total ~72 m over the
			// three 1926 openings), so it scales with the DEM resolution
			float wdt <- (c.is_dyke and c.is_open) ? (breach_cell_width * dx / 50.0) : dx;
			if (w > l_old) {
				float head <- w - max(c.z, l_old);
				if (head > 0) { vin <- vin + weir_coef * wdt * (head ^ 1.5) * 3600.0; }
			} else if (l_old > w) {
				float head <- l_old - max(c.z, w);
				if (head > 0) { vin <- vin - weir_coef * wdt * (head ^ 1.5) * 3600.0; }
			}
		}

		// STORAGE admission: newly reached dry ground BELOW the present level
		// (deep pockets the front just arrived at) is admitted only as fast
		// as the inflow can feed it AND within the hourly area budget —
		// otherwise redistribution would drain the flooded margins into the
		// pocket and the area would collapse while the river is still
		// rising. A draining region admits none; a gateless pond may creep.
		int allow <- max(0, int(spread_budget));
		if (not empty(feeders)) {
			float spare <- max(0.0, vin);
			list<cell> deep <- (region where (not each.wet and each.z < l_old))
				sort_by (- each.z);
			loop c over: deep {
				float need <- (l_old - c.z) * cell_area;
				if (spare >= need and allow > 0) {
					spare <- spare - need;
					allow <- allow - 1;
				} else {
					region >- c;                     // waits for a later hour
				}
			}
		}

		list<cell> sorted_cells <- region sort_by each.z;
		float capacity <- sum(region collect max(0.0, w - each.z)) * cell_area;
		float vt <- (vin >= 0.0) ? min(v + vin, max(capacity, v)) : max(capacity, v + vin);
		float lvl <- fill_level(sorted_cells, vt);

		// AREA limiter: the hourly budget of new land is GLOBAL across all
		// basins; held-back volume is conserved and spreads in later hours.
		// While FILLING, the clamp may defer new area but must never push
		// the level below the water already present. ('allow' already
		// accounts for the deep cells admitted above.)
		list<cell> drys_above <- sorted_cells where (not each.wet and each.z >= l_old);
		if (length(drys_above) > allow) {
			float zcap <- (allow = 0) ? l_old : (drys_above at allow).z;
			lvl <- min(lvl, (vin >= 0.0) ? max(zcap, min(lvl, l_old)) : zcap);
		}
		spread_budget <- spread_budget
			- (sorted_cells count (not each.wet and each.z < lvl));

		ask region {
			connected <- true;               // managed: skip the ponded pass
			if (z < lvl) { wet <- true; level <- lvl; h <- lvl - z; }
			else { wet <- false; h <- 0.0; level <- 0.0; }
		}
		// dyke gates show the flow passing over them
		ask feeders where each.is_dyke {
			connected <- true;
			wet <- true;
			level <- max(w, l_old);
			h <- max(0.0, max(w, l_old) - z);
		}
	}

	// water-surface level that stores volume vol in cells sorted by ground z
	float fill_level (list<cell> srt, float vol) {
		int n <- length(srt);
		float lvl <- first(srt).z;
		float rem <- vol / cell_area;        // metres of depth x cells
		int k <- 0;
		loop while: (rem > 0.0 and k < n) {
			float znext <- (k < n - 1) ? (srt at (k + 1)).z : 1e9;
			float cap_layer <- (znext - lvl) * (k + 1);
			if (rem >= cap_layer and k < n - 1) {
				rem <- rem - cap_layer;
				lvl <- znext;
				k <- k + 1;
			} else {
				lvl <- lvl + rem / (k + 1);
				rem <- 0.0;
			}
		}
		return lvl;
	}

	action update_buildings {
		ask building {
			flooded <- my_cells one_matches (each.h > bldg_flood_depth and not each.is_lake);
		}
		flooded_bldgs <- building count each.flooded;
	}

	reflex finish when: cycle >= n_stage {
		write "==================== 1926 FLOOD SUMMARY ====================";
		write "  Peak stage:            " + max(stage_series) + " m (30 July 12:00)";
		write "  Max flooded area:      " + max_flooded_km2 with_precision 2 + " km2 (outside river/lakes)";
		write "  Flooded buildings now: " + flooded_bldgs + " / " + length(building);
		write "  Still-wet cells:       " + (all_cells count (each.wet and each.h > flood_depth_min));
		write "============================================================";
		do pause;
	}
}

// ====================================================================== grid
grid cell file: dem_file neighbors: 4 frequency: 0
	use_regular_agents: false use_individual_shapes: false {
	float z;                             // ground / crest elevation (m)
	float orig_z;
	float h <- 0.0;                      // water depth (m)
	float level <- 0.0;                  // stored water-surface elevation (m)
	bool  active <- true;
	bool  is_river <- false;
	bool  is_lake <- false;
	bool  is_dyke <- false;
	bool  is_breach <- false;
	bool  is_open <- false;              // breach currently open
	bool  never_close <- false;
	bool  wet <- false;
	bool  connected <- false;
	int   rid <- -1;                     // id of the region owning this cell this step
	bool  is_gate <- false;              // water-margin gate this step
	int   adv <- 0;                      // dry-cell distance from water this step
	int   open_cycle <- -1;
	int   close_cycle <- -1;
	list<cell> neigh <- [];
	rgb   color <- #gray;

	action recolor_self {
		if (h > flood_depth_min) {                            // overtopped dykes show water
			color <- rgb(0, int(max(40, 170 - h * 25)), 255); // darker = deeper
		} else if (is_dyke and not is_open) {
			color <- rgb(115, 75, 35);                       // intact dyke wall
		} else {
			int gv <- int(min(235, max(60, (z - 4) * 11)));   // grey terrain
			color <- rgb(gv, gv, gv);
		}
	}
}

species river_area {
	aspect default { draw shape color: rgb(60, 110, 200, 120) border: #blue; }
}

species lake {
	aspect default { draw shape color: rgb(90, 160, 220, 150) border: rgb(50, 100, 180); }
}

species building {
	list<cell> my_cells <- [];
	bool flooded <- false;
	aspect default { draw shape color: flooded ? #red : rgb(40, 40, 40); }
}

species dyke_seg {
	string brk;
	string date_s;
	string commune;
	float crest <- 0.0;
	int open_cycle <- -1;
	int close_cycle <- -1;
	aspect default {
		bool now_open <- (brk = "YES") and cycle >= open_cycle
			and (close_cycle < 0 or cycle < close_cycle);
		draw shape color: now_open ? #red : ((brk = "YES") ? #orange : rgb(115, 75, 35)) width: 3;
	}
}

experiment Hourly1926 type: gui {
	parameter "Breach floor z (m)"        var: breach_floor    min: 6.0  max: 11.0 step: 0.25;
	parameter "River bed z (m)"           var: river_bed_z     min: 3.0  max: 7.0  step: 0.25;
	parameter "Allow crest overtopping"   var: allow_overtop;
	parameter "Crest trim (m)"            var: crest_trim      min: 0.0  max: 3.0  step: 0.1;
	parameter "Stage datum offset (m)"    var: stage_offset    min: -3.0 max: 3.0  step: 0.1;
	parameter "Weir coefficient"          var: weir_coef       min: 0.5  max: 3.0  step: 0.1;
	parameter "Breach width per cell (m)" var: breach_cell_width min: 1.0 max: 50.0 step: 1.0;
	parameter "Flood front speed (m/h)"   var: front_speed     min: 100.0 max: 5000.0 step: 100.0;
	parameter "Max spread (km2/h, total)" var: max_spread_km2h min: 0.1  max: 5.0  step: 0.1;
	parameter "Breach opens at hour"      var: breach_hour     min: 0    max: 23   step: 1;
	parameter "Min display depth (m)"     var: flood_depth_min min: 0.01 max: 0.5  step: 0.01;
	parameter "Building flood depth (m)"  var: bldg_flood_depth min: 0.05 max: 1.0 step: 0.05;

	output {
		display "Flood 1926" type: 2d {
			grid cell;
			species lake aspect: default transparency: 0.3;
			species river_area aspect: default transparency: 0.4;
			species building aspect: default;
			species dyke_seg aspect: default;
			overlay position: {10, 10} size: {320 #px, 90 #px} background: #black transparency: 0.4 {
				draw string(current_date, "dd MMM yyyy  HH:mm")
					at: {20 #px, 30 #px} color: #white font: font("Helvetica", 14, #bold);
				draw "stage " + (river_stage + stage_offset) with_precision 2 + " m"
					+ (stage_offset != 0.0 ? " (gauge " + river_stage with_precision 2 + ")" : "")
					+ "    flooded " + flooded_km2 with_precision 1 + " km2"
					at: {20 #px, 55 #px} color: #white font: font("Helvetica", 13, #plain);
			}
		}
		display "Stage hydrograph" type: 2d {
			chart "Red River stage at Hanoi (m)" type: series {
				data "stage (m)" value: river_stage color: #blue;
				data "flooded km2" value: flooded_km2 color: #red;
			}
		}
		monitor "Date"               value: string(current_date, "dd MMM HH:mm");
		monitor "Stage gauge (m)"    value: river_stage with_precision 2;
		monitor "Stage effective (m)" value: (river_stage + stage_offset) with_precision 2;
		monitor "Flooded area (km2)" value: flooded_km2 with_precision 2;
		monitor "Max flooded (km2)"  value: max_flooded_km2 with_precision 2;
		monitor "Flooded buildings"  value: flooded_bldgs;
		monitor "Breaches open"      value: breach_cells count each.is_open;
	}
}
