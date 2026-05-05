/**
* Name: HanoiFloodingModel3DBFS
* BFS / edge-cell variant of the Hanoi 3D flooding model.
*
* Spread mechanism is ported from the Quang Binh Three Dimensional Flooding
* Model: only edge cells push water to dry neighbours, with an incremental
* edge-update tracker that avoids the full O(N) scan. There is NO flow1/flow2
* selector and NO diffusion rate. Water moves via two complementary rules:
*   - Spread (BFS to dry): edge wet cells fill dry neighbours up to just below
*     the donor's water surface, draining the donor by the same volume.
*   - Equalise (volume-conserving): wet cells pass volume to their already-
*     processed wet neighbours of lower surface, halving the head difference
*     each pass — replaces Quang Binh's artificial rise with real flow.
*
* Time stepping: 1 cycle = 1 hour. Total simulation = 15 days (360 cycles).
* Author: Thanh-Do Nguyen (Nguyen Thanh Do)
* Tags: flooding, hydrology, 3D, mesh, BFS, volumetric
*/
model HanoiFloodingModel3DBFS

global {
	file river_shapefile <- file("../includes/RedRiver1925.shp");
	file lakes_shapefile <- file("../includes/Lakes1925.shp");
	file buildings_shapefile <- file("../includes/Buildings1925.shp");
	int resolution_grille <- 10 among: [10, 25, 40, 50];
	file mnt_csv <- file("../includes/mnt-gz" + resolution_grille + ".csv");
	geometry shape <- envelope("../includes/mnt-gz" + resolution_grille + ".asc");
	int nb_cols <- resolution_grille = 10 ? 1283 : (resolution_grille = 25 ? 513 : (resolution_grille = 40 ? 321 : 257));
	int nb_rows <- resolution_grille = 10 ? 854 : (resolution_grille = 25 ? 342 : (resolution_grille = 40 ? 214 : 171));

	float cell_area <- float(resolution_grille * resolution_grille);

	// Time
	float step <- 1 #h;
	int total_days <- 15;
	int total_cycles <- total_days * 24;

	// Hydrology state
	list<cell> drains;
	list<cell> active_cells;
	list<cell> river_cells;
	list<cell> active_water_cells <- [];
	list<cell> edge_water_cells <- [];
	float max_altitude;
	float min_altitude;
	int neighbours_type <- 8;

	// Parameters
	bool is_raining <- false;
	float rain <- 0.01 min: 0.0 max: 1.0 step: 0.01;                   // m / day
	float water_inp <- 2.0 min: 0.0 max: 10.0 step: 1.0;               // m, initial river depth
	float flow_rate <- 1000000.0 min: 0.0 max: 50000000.0 step: 100000.0; // m^3 / day
	float evaporation <- 0.001 min: 0.0 max: 1.0 step: 0.001;          // m / day
	bool river_input <- true;

	// BFS-specific parameters (from the Quang Binh model)
	float flow_threshold <- 0.01 min: 0.0 max: 1.0 step: 0.001;        // m, min depth needed to spread
	float min_flow_diff <- 0.001 min: 0.0 max: 0.1 step: 0.0001;       // m, min surface-vs-terrain head to push
	float min_visible_depth <- 0.01;                                   // m, render threshold

	string export_dir <- "../exported_results/";

	// 3D mesh fields
	field elevation_map;
	field water_field;

	init {
		matrix mnt_val <- matrix(mnt_csv);
		elevation_map <- field(nb_cols, nb_rows);
		water_field <- field(nb_cols, nb_rows);

		ask cell {
			altitude2 <- mnt_val at {grid_x, grid_y};
			is_inactive <- altitude2 = -9999.0;
			if (is_inactive) {
				altitude2 <- 0;
			}
			elevation_map[grid_x, grid_y] <- altitude2;
		}

		active_cells <- cell where !(each.is_inactive);

		ask active_cells {
			active_neighbours <- (self neighbors_at 1) where !(each.is_inactive);
			if (length(active_neighbours) < neighbours_type) {
				is_drain <- true;
			}
		}

		max_altitude <- active_cells max_of (each.altitude2);
		min_altitude <- active_cells min_of (each.altitude2);
		write "max_altitude : " + max_altitude + " min_altitude : " + min_altitude;
		write "cell_area : " + cell_area + " m^2 (resolution " + resolution_grille + " m)";

		ask cell {
			alt_norm <- 100 * (((altitude2 - min_altitude) / (max_altitude - min_altitude)) ^ 2);
			water_field[grid_x, grid_y] <- min_altitude - 1.0;
		}

		create building from: buildings_shapefile with: [height::rnd(10) + 5.0];

		create river from: river_shapefile {
			cells_concerned <- cell overlapping self;
			if (river_input) {
				ask cells_concerned {
					water_volume <- world.water_inp * world.cell_area;
					is_water <- true;
				}
			}
		}
		river_cells <- remove_duplicates(river accumulate (each.cells_concerned)) where !(each.is_inactive);

		// Lake is loaded for display only; no water seeded.
		create lake from: lakes_shapefile {
			cells_concerned <- cell overlapping self;
		}

		drains <- active_cells where (each.is_drain);

		// Build initial active and edge sets from the river seed.
		active_water_cells <- active_cells where (each.is_water);
		ask active_water_cells {
			water_field[grid_x, grid_y] <- altitude2 + water_height;
			if (!empty(active_neighbours where (!each.is_water))) {
				is_edge_cell <- true;
			}
		}
		edge_water_cells <- active_water_cells where (each.is_edge_cell);

		write "Active cells: " + length(active_cells)
			+ " | River cells: " + length(river_cells)
			+ " | Initial wet cells: " + length(active_water_cells)
			+ " | Initial edge cells: " + length(edge_water_cells)
			+ " | Drain cells: " + length(drains);
		write "Run length: " + total_cycles + " cycles (" + total_days + " days at 1 h/step)";
	}

	reflex pause_at_end when: cycle >= total_cycles {
		write "Reached " + total_days + " days. Pausing.";
		do pause;
	}

	reflex raining when: is_raining {
		float vol_per_cell <- (rain / 24.0) * cell_area;
		list<cell> newly_wet <- [];
		ask active_cells {
			water_volume <- water_volume + vol_per_cell;
			if (!is_water and water_height > flow_threshold) {
				is_water <- true;
				newly_wet <- newly_wet + self;
			}
		}
		if (!empty(newly_wet)) {
			active_water_cells <- active_water_cells + newly_wet;
			do update_edge_cells(newly_wet, []);
		}
	}

	reflex river_inflow when: river_input and !empty(river_cells) {
		float vol_per_cell <- (flow_rate / 24.0) / length(river_cells);
		ask river_cells {
			water_volume <- water_volume + vol_per_cell;
		}
	}

	reflex adding_input_water {
		loop so over: source {
			int n <- length(so.cells_concerned);
			if (n > 0) {
				float per_cell <- (so.water_input / 24.0) / n;
				ask so.cells_concerned {
					water_volume <- water_volume + per_cell;
				}
			}
		}
	}

	// BFS spread (edge -> dry) followed by volume-conserving equalisation among wet cells.
	reflex simulate_water_flow {
		list<cell> new_water_cells <- [];
		list<cell> affected_neighbors <- [];

		// Phase 1: SPREAD from edge cells to dry neighbours.
		ask edge_water_cells {
			if (water_volume > 0 and water_height > flow_threshold) {
				float donor_surface <- altitude2 + water_height;
				list<cell> dry_targets <- active_neighbours where (!each.is_water);
				ask dry_targets {
					float my_top <- altitude2 + (the_digue = nil ? 0.0 : the_digue.height);
					float head_above_terrain <- donor_surface - my_top;
					if (head_above_terrain > world.min_flow_diff) {
						// Quang Binh rule: receiver lands just below donor's surface.
						float target_height <- max(world.flow_threshold, donor_surface - 0.05 - altitude2);
						float fill_volume <- target_height * world.cell_area;
						float available <- min(fill_volume, myself.water_volume);
						if (available > 0) {
							myself.water_volume <- myself.water_volume - available;
							water_volume <- water_volume + available;
							is_water <- true;
							new_water_cells <- new_water_cells + self;
							// Record wet neighbours of the new cell for edge re-evaluation.
							ask active_neighbours where (each.is_water) {
								if (!(affected_neighbors contains self)) {
									affected_neighbors <- affected_neighbors + self;
								}
							}
						}
					}
				}
			}
		}

		// Phase 2: incremental edge-set update (only changed cells, per Quang Binh).
		if (!empty(new_water_cells)) {
			active_water_cells <- active_water_cells + new_water_cells;
			do update_edge_cells(new_water_cells, affected_neighbors);
		}

		// Phase 3: EQUALISE among wet cells (volume-conserving, replaces artificial rise).
		ask active_water_cells {
			already <- false;
		}
		ask shuffle(active_water_cells) {
			do equalise_with_wet;
		}
	}

	// Mark new water cells as edges if they border any dry cell, and re-check
	// affected wet neighbours — those that just lost a dry neighbour stop being
	// edge cells. Mirrors update_edge_cells in the Quang Binh model.
	action update_edge_cells (list<cell> new_water_cells, list<cell> affected_neighbors) {
		ask new_water_cells {
			if (!is_edge_cell and !empty(active_neighbours where (!each.is_water))) {
				is_edge_cell <- true;
				edge_water_cells <- edge_water_cells + self;
			}
		}
		ask affected_neighbors {
			if (is_water and is_edge_cell) {
				if (empty(active_neighbours where (!each.is_water))) {
					is_edge_cell <- false;
					edge_water_cells <- edge_water_cells - self;
				}
			}
		}
	}

	reflex evaporation {
		float vol_per_cell <- (evaporation / 24.0) * cell_area;
		ask active_cells where (each.water_volume > 0) {
			water_volume <- max(0.0, water_volume - vol_per_cell);
		}
	}

	reflex draining {
		ask drains {
			water_volume <- 0.0;
		}
	}

	reflex update_water_field {
		ask active_cells {
			water_field[grid_x, grid_y] <- water_height > min_visible_depth
				? (altitude2 + water_height)
				: (altitude2 - 1.0);
		}
	}

//	reflex export_water_height when: cycle mod 24 = 0 and cycle <= total_cycles {
//		int day <- int(cycle / 24);
//		string day_str <- (day < 10 ? "0" : "") + string(day);
//		string base_name <- export_dir + "water_height_3D_BFS_res" + resolution_grille + "_day" + day_str;
//		ask cell {
//			grid_value <- water_height;
//		}
//		save cell to: base_name + ".tif" format: "geotiff";
//		write "Exported day " + day + " (cycle " + cycle + ") -> " + base_name;
//	}

	species river {
		list<cell> cells_concerned;

		aspect geometry {
			draw shape color: rgb("blue") depth: world.water_inp;
		}
	}

	species lake {
		list<cell> cells_concerned;

		aspect geometry {
			draw shape color: rgb("green") depth: 0.5;
		}
	}

	species building parent: digue frequency: 0 {
		aspect geometry {
			draw shape color: rgb("pink") depth: height;
		}
	}

	species digue {
		float height;
		list<cell> cells_concerned;

		init {
			cells_concerned <- cell overlapping self;
			ask cells_concerned {
				the_digue <- myself;
			}
		}

		aspect geometry {
			draw shape color: rgb("red") depth: height;
		}

		action destroy_digue {
			ask cells_concerned {
				the_digue <- nil;
			}
			do die;
		}

		user_command "Destroy digue" action: destroy_digue;
	}

	species source {
		list<cell> cells_concerned <- cell overlapping self;
		float water_input;  // m^3 / day

		aspect geometry {
			draw shape + 2.0 color: rgb("green") depth: 1.0;
		}
	}

	grid cell width: nb_cols height: nb_rows neighbors: neighbours_type frequency: 0 use_regular_agents: false use_individual_shapes: false {
		bool is_inactive <- false;
		float altitude2;

		// Cached active neighbours - computed once during init
		list<cell> active_neighbours;

		float water_volume <- 0.0 min: 0.0;
		float water_height -> water_volume / world.cell_area;

		bool is_drain <- false;
		bool is_water <- false;
		bool is_edge_cell <- false;
		float alt_norm <- 0.0;
		digue the_digue;
		bool already <- false;

		// Volume-conserving pair equalisation: transfer half the head
		// difference (in volume) to each already-processed lower wet neighbour.
		// One pass per cycle; converges over multiple cycles.
		action equalise_with_wet {
			if (water_volume > 0) {
				float my_surface <- altitude2 + water_height;
				list<cell> wet_lower <- active_neighbours where
					(each.is_water and each.already and (each.altitude2 + each.water_height) < my_surface);
				loop n over: shuffle(wet_lower) {
					float surface_diff <- (altitude2 + water_height) - (n.altitude2 + n.water_height);
					if (surface_diff > 0) {
						float volume_transfer <- min(water_volume, surface_diff * world.cell_area / 2.0);
						water_volume <- water_volume - volume_transfer;
						n.water_volume <- n.water_volume + volume_transfer;
					}
				}
			}
			already <- true;
		}
	}
}

experiment main_gui type: gui {
	parameter "Is raining" var: is_raining;
	parameter "Rain (m/day)" var: rain;
	parameter "Initial river depth (m)" var: water_inp;
	parameter "River flow rate (m^3/day)" var: flow_rate;
	parameter "Evaporation (m/day)" var: evaporation;
	parameter "River input enabled" var: river_input;
	parameter "Flow threshold (m)" var: flow_threshold;
	parameter "Min flow diff (m)" var: min_flow_diff;
	parameter "Min visible depth (m)" var: min_visible_depth;

	action _init_ {
		int resolution <- 0;
		loop while: (!(resolution in [10, 25, 40, 50])) {
			map resolution_input <- user_input_dialog("Choose the resolution of the grid among 10, 25, 40, 50 meters", [choose("Choose a value", int, 10, [10, 25, 40, 50])]);
			resolution <- int(resolution_input["Choose a value"]);
			if (resolution in [10, 25, 40, 50]) {
				create simulation with: [resolution_grille::resolution];
			}
		}
	}

	output {
		display map type: opengl {
			mesh elevation_map
				scale: 1
				grayscale: true
				smooth: false
				triangulation: true;

			mesh water_field
				scale: 1
				color: rgb(100, 150, 255, 180)
				smooth: false
				triangulation: true;

			species river aspect: geometry;
			species lake aspect: geometry;
			species building aspect: geometry;
			species digue aspect: geometry;
		}
	}
}
