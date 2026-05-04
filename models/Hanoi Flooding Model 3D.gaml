/**
* Name: HanoiFloodingModel3D
* 3D variant of the Hanoi 1925 flooding model. Same UI and parameters as the
* original Hanoi Flooding Model, but the flow algorithm is volumetric (m^3 per
* cell) and terrain + water are rendered as OpenGL meshes inspired by the
* Quang Binh Three Dimensional Flooding Model.
*
* Time stepping: 1 cycle = 1 hour. Total simulation = 15 days (360 cycles).
* Author: Thanh-Do Nguyen (Nguyen Thanh Do)
* Tags: flooding, hydrology, 3D, mesh, volumetric
*/
model HanoiFloodingModel3D

global {
	file river_shapefile <- file("../includes/RedRiver1925.shp");
	file lakes_shapefile <- file("../includes/Lakes1925.shp");
	file buildings_shapefile <- file("../includes/Buildings1925.shp");
	int resolution_grille <- 10 among: [10, 25, 40, 50];
	file mnt_csv <- file("../includes/mnt-gz" + resolution_grille + ".csv");
	geometry shape <- envelope("../includes/mnt-gz" + resolution_grille + ".asc");
	int nb_cols <- resolution_grille = 10 ? 1283 : (resolution_grille = 25 ? 513 : (resolution_grille = 40 ? 321 : 257));
	int nb_rows <- resolution_grille = 10 ? 854 : (resolution_grille = 25 ? 342 : (resolution_grille = 40 ? 214 : 171));

	// Each grid cell is a square with side = resolution_grille metres.
	float cell_area <- float(resolution_grille * resolution_grille);

	// Time
	float step <- 1 #h;
	int total_days <- 15;
	int total_cycles <- total_days * 24;

	// Hydrology state
	list<cell> drains;
	list<cell> active_cells;
	list<cell> river_cells;
	float max_altitude;
	float min_altitude;
	int neighbours_type <- 8;

	// Parameters (carry the same names/units as the 2D model where possible)
	bool is_raining <- false;
	float rain <- 0.01 min: 0.0 max: 1.0 step: 0.01;                   // m / day
	float water_inp <- 2.0 min: 0.0 max: 10.0 step: 1.0;               // m, initial water depth on river cells
	float flow_rate <- 1000000.0 min: 0.0 max: 50000000.0 step: 100000.0; // m^3 / day, continuous river inflow
	float evaporation <- 0.001 min: 0.0 max: 1.0 step: 0.001;          // m / day
	bool river_input <- true;
	float diffusion_rate <- 0.8;
	int algo_flowing <- 1 among: [1, 2];
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

		// Cache active neighbours once (is_inactive never changes after init)
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

		// Park the water surface clearly below the terrain everywhere by default
		// so the translucent water mesh is hidden behind the terrain mesh on dry
		// cells. Wet cells get raised above terrain in the loop below.
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
				}
			}
		}
		river_cells <- remove_duplicates(river accumulate (each.cells_concerned)) where !(each.is_inactive);

		// Lake is loaded for display only; no water is seeded here per spec.
		create lake from: lakes_shapefile {
			cells_concerned <- cell overlapping self;
		}

		drains <- active_cells where (each.is_drain);

		ask active_cells where (each.water_volume > 0) {
			water_field[grid_x, grid_y] <- altitude2 + water_height;
		}

		write "Active cells: " + length(active_cells)
			+ " | River cells: " + length(river_cells)
			+ " | Drain cells: " + length(drains);
		write "Run length: " + total_cycles + " cycles (" + total_days + " days at 1 h/step)";
	}

	reflex pause_at_end when: cycle >= total_cycles {
		write "Reached " + total_days + " days. Pausing.";
		do pause;
	}

	reflex raining when: is_raining {
		float vol_per_cell <- (rain / 24.0) * cell_area;
		ask active_cells {
			water_volume <- water_volume + vol_per_cell;
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

	reflex flowing {
		ask active_cells {
			already <- false;
		}

		switch algo_flowing {
			match 1 {
				ask shuffle(active_cells) {
					do flow1;
				}
			}

			match 2 {
				ask (active_cells sort_by ((each.altitude2 + each.water_height + (each.the_digue = nil ? 0.0 : each.the_digue.height)))) {
					do flow2;
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

	reflex export_water_height when: cycle mod 24 = 0 and cycle <= total_cycles {
		int day <- int(cycle / 24);
		string day_str <- (day < 10 ? "0" : "") + string(day);
		string base_name <- export_dir + "water_height_3D_res" + resolution_grille + "_day" + day_str;
		ask cell {
			grid_value <- water_height;
		}
		save cell to: base_name + ".tif" format: "geotiff";
		write "Exported day " + day + " (cycle " + cycle + ") -> " + base_name;
	}

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

		// Volumetric state. water_height is derived so the rest of the model
		// can reason in metres while bookkeeping is conserved in m^3.
		float water_volume <- 0.0 min: 0.0;
		float water_height -> water_volume / world.cell_area;

		bool is_drain <- false;
		float alt_norm <- 0.0;
		float height2;
		digue the_digue;
		bool already <- false;

		action flow1 {
			if (water_volume > 0) {
				list<cell> neighbour_cells_al <- active_neighbours where (each.already);
				if (!empty(neighbour_cells_al)) {
					ask neighbour_cells_al {
						height2 <- altitude2 + water_height + (the_digue = nil ? 0.0 : the_digue.height);
					}

					height2 <- altitude2 + water_height;
					list<cell> flow_cells <- (neighbour_cells_al where (height2 > each.height2));
					if (!empty(flow_cells)) {
						loop flow_cell over: shuffle(flow_cells) {
							float head_diff <- height2 - flow_cell.height2;
							float volume_flowing <- max(0.0, min(head_diff * world.cell_area, water_volume * world.diffusion_rate));
							water_volume <- water_volume - volume_flowing;
							flow_cell.water_volume <- flow_cell.water_volume + volume_flowing;
							height2 <- altitude2 + water_height;
						}
					}
				}
			}

			already <- true;
		}

		action flow2 {
			if (water_volume > 0) {
				list<cell> neighbour_cells_al <- active_neighbours where (each.already);
				if (!empty(neighbour_cells_al)) {
					ask neighbour_cells_al {
						height2 <- altitude2 + water_height + (the_digue = nil ? 0.0 : the_digue.height);
					}

					height2 <- altitude2 + water_height;
					list<cell> flow_cells <- (neighbour_cells_al where (height2 > each.height2));
					if (!empty(flow_cells)) {
						float mean_height <- mean([height2] + flow_cells collect (each.height2));
						loop flow_cell over: flow_cells sort_by (each.height2) {
							float head_to_target <- mean_height - flow_cell.height2;
							float volume_flowing <- max(0.0, min(head_to_target * world.cell_area, water_volume));
							water_volume <- water_volume - volume_flowing;
							flow_cell.water_volume <- flow_cell.water_volume + volume_flowing;
						}
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
	parameter "Diffusion rate" var: diffusion_rate;
	parameter "Flowing algorithm" var: algo_flowing;
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
				scale: 100
				grayscale: true
				smooth: false
				triangulation: true;

			mesh water_field
				scale: 100
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
