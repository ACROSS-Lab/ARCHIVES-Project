/**
* Name: HanoiFloodingModel
* This model is based on the hydrological model developed by Patrick Taillandier in the ARCHIVES Project. The model incorporates modifications to ensure compatibility with GAMA 2025-06.
* Author: Thanh-Do Nguyen (Nguyen Thanh Do)
* Tags:
*/
model HanoiFloodingModel

global {
	file river_shapefile <- file("../includes/RedRiver1925.shp");
	file lakes_shapefile <- file("../includes/Lakes1925.shp");
	file buildings_shapefile <- file("../includes/Buildings1925.shp");
	int resolution_grille <- 10 among: [10, 25, 40, 50];
	file mnt_csv <- file("../includes/mnt-gz" + resolution_grille + ".csv");
	geometry shape <- envelope("../includes/mnt-gz" + resolution_grille + ".asc");
	int nb_cols <- resolution_grille = 10 ? 1283 : (resolution_grille = 25 ? 513 : (resolution_grille = 40 ? 321 : 257));
	int nb_rows <- resolution_grille = 10 ? 854 : (resolution_grille = 25 ? 342 : (resolution_grille = 40 ? 214 : 171));
	list<cell> drains;
	float max_altitude;
	float min_altitude;
	int neighbours_type <- 8;
	bool is_raining <- false;
	float rain <- 0.01 min: 0.0 max: 1.0 step: 0.01;
	float water_inp <- 2.0 min: 0.0 max: 10.0 step: 1.0;
	float evaporation <- 0.001 min: 0.0 max: 1.0 step: 0.001;
	bool river_input <- true;
	float diffusion_rate <- 0.8;
	int algo_flowing <- 1 among: [1, 2];
	list<cell> active_cells;

	init {
		matrix mnt_val <- matrix(mnt_csv);
		ask cell {
			altitude2 <- mnt_val at {grid_x, grid_y};
			is_inactive <- altitude2 = -9999.0;
			if (is_inactive) {
				altitude2 <- 0;
			}
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
		ask cell {
			alt_norm <- 100 * (((altitude2 - min_altitude) / (max_altitude - min_altitude)) ^ 2);
			int val_c <- int(255 * (1.0 - ((altitude2 - min_altitude) / (max_altitude - min_altitude))));
			color_mnt <- rgb([val_c, val_c, val_c]);
		}

		create building from: buildings_shapefile with: [height::rnd(10) + 5.0];
		create river from: river_shapefile {
			if (river_input) {
				ask (cell overlapping self) {
					water_height <- water_inp;
					do update_color;
				}
			}
		}

		create lake from: lakes_shapefile {
			if (river_input) {
				ask (cell overlapping self) {
					water_height <- water_inp;
					do update_color;
				}
			}
		}

		drains <- active_cells where (each.is_drain);
	}

	reflex raining when: is_raining {
		ask active_cells {
			water_height <- water_height + rain;
		}
	}

	reflex adding_input_water {
		loop so over: source {
			ask so.cells_concerned {
				water_height <- water_height + (so.water_input * water_inp);
			}
		}
	}

	reflex flowing {
		ask active_cells {
			already <- false;
		}

		switch algo_flowing {
			match 1 {
				ask (shuffle(active_cells)) {
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
		// Only evaporate and update color for cells that have water
		ask active_cells where (each.water_height > 0) {
			water_height <- water_height - evaporation;
			do update_color;
		}
	}

	reflex draining {
		ask drains {
			water_height <- 0;
		}
	}

	species river {

		aspect geometry {
			draw shape color: rgb("blue") depth: world.water_inp;
		}
	}

	species lake {

		aspect geometry {
			draw shape color: rgb("green") depth: world.water_inp;
		}
	}

	species building parent: digue frequency: 0 {

		aspect geometry {
			draw shape color: rgb("pink") depth: 100;
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
		float water_input;

		aspect geometry {
			draw shape + 2.0 color: rgb("green") depth: water_input;
		}
	}

	grid cell width: nb_cols height: nb_rows neighbors: neighbours_type frequency: 0 use_regular_agents: false use_individual_shapes: false {
		bool is_inactive <- false;
		float altitude2;

		// Cached active neighbours - computed once during init
		list<cell> active_neighbours;

		float water_height <- 0.0 min: 0.0;
		rgb color_mnt;
		bool is_drain <- false;
		float alt_norm <- 0.0;
		float height2;
		digue the_digue;
		bool already <- false;

		action update_color {
			int val_water <- max([0, min([255, int(255 * (1 - (water_height / 1.0)))])]);
			color <- water_height > 0 ? rgb([val_water, val_water, 255]) : color_mnt;
		}

		action flow2 {
			if (water_height > 0) {
				// Use cached neighbours instead of recomputing each call
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
							float water_flowing <- max([0.0, min([(mean_height - flow_cell.height2), water_height])]);
							water_height <- water_height - water_flowing;
							flow_cell.water_height <- flow_cell.water_height + water_flowing;
						}
					}
				}
			}

			already <- true;
		}

		action flow1 {
			if (water_height > 0) {
				// Use cached neighbours instead of recomputing each call
				list<cell> neighbour_cells_al <- active_neighbours where (each.already);
				if (!empty(neighbour_cells_al)) {
					ask neighbour_cells_al {
						height2 <- altitude2 + water_height + (the_digue = nil ? 0.0 : the_digue.height);
					}

					height2 <- altitude2 + water_height;
					list<cell> flow_cells <- (neighbour_cells_al where (height2 > each.height2));
					if (!empty(flow_cells)) {
						loop flow_cell over: shuffle(flow_cells) {
							float water_flowing <- max([0.0, min([(height2 - flow_cell.height2), water_height * diffusion_rate])]);
							water_height <- water_height - water_flowing;
							flow_cell.water_height <- flow_cell.water_height + water_flowing;
							height2 <- altitude2 + water_height;
						}
					}
				}
			}

			already <- true;
		}

		aspect mnt {
			draw shape color: color_mnt depth: alt_norm;
		}

		aspect water {
			rgb c <- water_height > 0 ? color : color_mnt;
			draw shape color: c border: c depth: water_height;
		}
	}
}

experiment main_gui type: gui {
	parameter "Is raining" var: is_raining;
	parameter "Rain quantity" var: rain;
	parameter "Water input" var: water_inp;
	parameter "Evaporation" var: evaporation;
	parameter "River input" var: river_input;
	parameter "Diffusion rate" var: diffusion_rate;
	parameter "Flowing algorithm" var: algo_flowing;

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
			species river aspect: geometry;
			species lake aspect: geometry;
			species building aspect: geometry;
			species digue aspect: geometry;
			grid cell transparency: 0.5 elevation: -0.1;
		}
	}
}
