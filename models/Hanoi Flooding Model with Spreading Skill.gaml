/**
* Name: HanoiFloodingModelwithSpreadingSkill
* Based on the internal empty template. 
* Author: thanhdonguyen
* Tags: 
*/


model HanoiFloodingModelwithSpreadingSkill

/* Insert your model definition here */


global skills: [spreading] { 
    // === DATA FILES ===
    int resolution_grille <- 10 among: [10, 25, 40, 50];
    file dem_file <- file("../includes/mnt-gz" + resolution_grille + ".asc");
    field elevation_map <- field(dem_file);
    geometry shape <- envelope(dem_file);
    
    // Water sources
    file river_shapefile <- file("../includes/RedRiver1925.shp");
    file lakes_shapefile <- file("../includes/Lakes1925.shp");
    list<geometry> water_geometries <- [];
    
    // Obstacles (buildings/dykes)
    file buildings_shapefile <- file("../includes/Buildings1925.shp");
    bool load_buildings_as_obstacles <- true;
    
    // === DISPLAY FIELDS ===
    field display_water_field;
    field display_obstacle_field;
    
    // === SIMULATION PARAMETERS ===
    float initial_water_depth <- 2.0 min: 0.0 max: 10.0;
    float flow_threshold_param <- 0.01;
    float rising_rate_param <- 0.3;
    float building_height <- 5.0 min: 0.0 max: 20.0;
    bool buildings_destroyable <- false;
    
    // === SIMULATION STATE ===
    // bool simulation_active <- false;
    int total_cycles <- 0;
    
    init {
        write "=== HANOI FLOODING MODEL WITH SPREADINGSKILL ===";
        write "Resolution: " + resolution_grille + "m";
        write "DEM dimensions: " + elevation_map.columns + "x" + elevation_map.rows;
        write "World bounds: " + shape;
        
        // Load water geometries from shapefiles
        list<geometry> rivers <- river_shapefile.contents;
        list<geometry> lakes <- lakes_shapefile.contents;
        water_geometries <- rivers + lakes;
        create water from: water_geometries;
        create building from: buildings_shapefile;
        write "Loaded " + length(rivers) + " rivers and " + length(lakes) + " lakes";
        
        // Initialize display fields
        display_water_field <- field(elevation_map.columns, elevation_map.rows);
        display_obstacle_field <- field(elevation_map.columns, elevation_map.rows);
        
        loop i from: 0 to: elevation_map.columns - 1 {
            loop j from: 0 to: elevation_map.rows - 1 {
                display_water_field[i, j] <- 0.0;
                display_obstacle_field[i, j] <- 0.0;
            }
        }
        
        // Initialize grid cells
        ask water_grid {
            water_level <- 0.0;
        }
        
        write "=== SPREADING SKILL INITIALIZATION ===";
        
        // Create obstacle field with proper dimensions
        field obstacle_field_init <- field(elevation_map.columns, elevation_map.rows);
        loop i from: 0 to: elevation_map.columns - 1 {
            loop j from: 0 to: elevation_map.rows - 1 {
                obstacle_field_init[i, j] <- 0.0;
            }
        }
        
        // Initialize spreading grid (calling skill action directly on global)
        do initialize_spreading_grid_with_obstacle_field(
            dem_field: elevation_map,
            obstacle_field: obstacle_field_init,
            water_geometries: water_geometries,
            initial_water_depth: initial_water_depth,
            flow_threshold: flow_threshold_param,
            rising_rate: rising_rate_param,
            min_flow_diff: 0.001,
            equalization_threshold: 0.1
        );
        
        // Sync display fields
        if (water_field != nil) {
            display_water_field <- water_field;
        }
        if (obstacle_field != nil) {
            display_obstacle_field <- obstacle_field;
        }
        
        write "Grid initialized: " + grid_width + "x" + grid_height;
        write "Initial water cells: " + get_active_water_count();
        write "Initial edge cells: " + get_edge_cell_count();
        
        // Load buildings as obstacles if enabled
        if (load_buildings_as_obstacles and buildings_shapefile != nil) {
            write "Loading buildings as obstacles...";
            list<geometry> building_geoms <- buildings_shapefile.contents;
            
            if (length(building_geoms) > 0) {
                bool success <- add_obstacles_from_shapefile(
                    obstacle_shapefile: building_geoms,
                    obstacle_type: "building",
                    height: building_height,
                    uniform_height: true,
                    destroyable: buildings_destroyable,
                    destruction_time: 20.0
                );
                
                if (success) {
                    write "Buildings loaded successfully as obstacles";
                    write "Total obstacle cells: " + get_active_obstacle_count();
                    
                    // Update display field
                    if (obstacle_field != nil) {
                        display_obstacle_field <- obstacle_field;
                    }
                } else {
                    write "WARNING: Failed to load buildings";
                }
            }
        }
        
        write "Initialization complete";
        write "========================================";
    }
    
    // === SIMULATION STEP ===
    reflex simulation_step when: simulation_active {
        do simulate_spreading_step();
        
        // Update display fields
        if (water_field != nil) {
            display_water_field <- water_field;
        }
        
        // Update obstacle field periodically
        if (obstacle_field != nil and (cycle mod 10 = 0)) {
            display_obstacle_field <- obstacle_field;
        }
        
        // Report progress
        if (cycle mod 50 = 0) {
            write "Step " + get_current_step() + ": " + 
                  get_active_water_count() + " water cells, " +
                  get_edge_cell_count() + " edge cells";
            
            if (get_obstacles_under_attack() > 0) {
                write "WARNING: " + get_obstacles_under_attack() + " obstacles under water attack";
            }
        }
    }
    
    reflex count_cycles {
        total_cycles <- total_cycles + 1;
    }
    
    // === CONTROL ACTIONS ===
    action start_sim {
        do start_spreading_simulation();
        simulation_active <- true;
        write "Simulation STARTED";
    }
    
    action stop_sim {
        do stop_spreading_simulation();
        simulation_active <- false;
        write "Simulation STOPPED";
    }
    
    action reset_sim {
        do reset_spreading_simulation(water_geometries, initial_water_depth);
        simulation_active <- false;
        
        // Reinitialize
        field obstacle_field_reset <- field(elevation_map.columns, elevation_map.rows);
        loop i from: 0 to: elevation_map.columns - 1 {
            loop j from: 0 to: elevation_map.rows - 1 {
                obstacle_field_reset[i, j] <- 0.0;
            }
        }
        
        do initialize_spreading_grid_with_obstacle_field(
            dem_field: elevation_map,
            obstacle_field: obstacle_field_reset,
            water_geometries: water_geometries,
            initial_water_depth: initial_water_depth,
            flow_threshold: flow_threshold_param,
            rising_rate: rising_rate_param,
            min_flow_diff: 0.001,
            equalization_threshold: 0.1
        );
        
        // Reload buildings if enabled
        if (load_buildings_as_obstacles) {
            list<geometry> building_geoms <- buildings_shapefile.contents;
            if (length(building_geoms) > 0) {
                do add_obstacles_from_shapefile(
                    obstacle_shapefile: building_geoms,
                    obstacle_type: "building",
                    height: building_height,
                    uniform_height: true,
                    destroyable: buildings_destroyable,
                    destruction_time: 20.0
                );
            }
        }
        
        // Sync display fields
        if (water_field != nil) {
            display_water_field <- water_field;
        }
        if (obstacle_field != nil) {
            display_obstacle_field <- obstacle_field;
        }
        
        write "Simulation RESET";
    }
    
    // === RAIN ACTIONS ===
    action start_monsoon_rain {
        do start_rain(0.01, 1.0);
        write "Monsoon rain started (0.01m/step)";
    }
    
    action start_heavy_rain {
        do start_rain(0.3, 1.5);
        write "Heavy rain started (0.3m/step)";
    }
    
    action stop_rain_action {
        do stop_rain();
        write "Rain stopped";
    }
    
    // === OBSTACLE ACTIONS ===
    action clear_all_buildings {
        do clear_all_obstacles();
        if (obstacle_field != nil) {
            display_obstacle_field <- obstacle_field;
        }
        write "All buildings cleared";
    }
    
    action report_status {
        write "=== HANOI FLOOD SIMULATION STATUS ===";
        write "Simulation: " + (simulation_active ? "ACTIVE" : "STOPPED");
        write "Cycle: " + total_cycles;
        write "Water cells: " + get_active_water_count();
        write "Edge cells: " + get_edge_cell_count();
        write "Obstacle cells: " + get_active_obstacle_count();
        write "Rain active: " + is_rain_active();
        if (is_rain_active()) {
            write "Rain rate: " + get_rain_rate() + "m/step";
        }
        write "====================================";
    }
    
    action report_detailed_status {
        write "=== DETAILED STATUS ===";
        write "SIMULATION:";
        write "  Active: " + is_simulation_active();
        write "  Step: " + get_current_step();
        write "WATER:";
        write "  Water cells: " + get_active_water_count();
        write "  Edge cells: " + get_edge_cell_count();
        write "  Coverage: " + ((get_active_water_count() / (grid_width * grid_height)) * 100.0 with_precision 1) + "%";
        write "RAIN:";
        write "  Active: " + is_rain_active();
        if (is_rain_active()) {
            write "  Rate: " + (get_rain_rate() with_precision 3) + "m/step";
            write "  Intensity: " + (get_rain_intensity() with_precision 1) + "x";
        }
        write "OBSTACLES:";
        write "  Total cells: " + get_active_obstacle_count();
        write "  Under attack: " + get_obstacles_under_attack();
        write "=====================";
    }
}

// Water species
species water {
    rgb color <- rgb("blue");
    float depth <- 0.5;
    
    aspect geometry {
        draw shape color: color depth: depth;
    }
}

// Building species
species building {
    rgb color <- rgb("pink");
    float height <- 0.5;
    
    aspect geometry {
        draw shape color: color depth: height;
    }
}

// Dike species
species dike {
    rgb color <- rgb("red");
    float height <- 5.0;
    
    aspect geometry {
        draw shape color: color depth: height;
    }
}

// Water grid
grid water_grid width:elevation_map.columns height: elevation_map.rows frequency: 1 use_regular_agents: false use_individual_shapes: false {
    float water_level;
    
    init {
    	water_level <- 0.0;
    }
    
    reflex update_water when: simulation_active and (cycle mod 5 = 0) {
        if (display_water_field != nil) {
            water_level <- float(display_water_field[grid_x, grid_y]);
            write self.name + " " + water_level;
        }
        
        if (water_level > 0.01) {
        	color <- rgb("blue");
            draw shape color: color depth: water_level * 100;
        }
    }
    
    aspect default {
        if (display_water_field != nil) {
            water_level <- float(display_water_field[grid_x, grid_y]);
            write self.name + " " + water_level;
        }
        
        if (water_level > 0.01) {
        	color <- rgb("blue");
            draw shape color: color depth: water_level;
        }
    }
}

experiment HanoiFloodSimulation type: gui {
    
    // === PARAMETERS ===
    parameter "Resolution (m)" var: resolution_grille category: "Grid";
    parameter "Initial water depth (m)" var: initial_water_depth category: "Water";
    parameter "Flow threshold (m)" var: flow_threshold_param category: "Water";
    parameter "Rising rate (m/step)" var: rising_rate_param category: "Water";
    parameter "Load buildings as obstacles" var: load_buildings_as_obstacles category: "Obstacles";
    parameter "Building height (m)" var: building_height category: "Obstacles";
    parameter "Buildings destroyable" var: buildings_destroyable category: "Obstacles";
    
    // === INITIALIZATION ===
    action _init_ {
        int resolution <- 0;
        loop while: (!(resolution in [10, 25, 40, 50])) {
            map resolution_input <- user_input_dialog(
                "Choose the resolution of the grid among 10, 25, 40, 50 meters",
                [choose("Choose a value", int, 10, [10, 25, 40, 50])]
            );
            resolution <- int(resolution_input["Choose a value"]);
            if (resolution in [10, 25, 40, 50]) {
                create simulation with: [resolution_grille::resolution];
            }
        }
    }
    
    output {
        display "3D View" type: opengl refresh: true {
            grid water_grid ;
            species water aspect: geometry;
            species building aspect: geometry;
            species dike aspect: geometry;
        }
        
        display "Hanoi Flood Simulation 3D" type: opengl {
            mesh elevation_map scale: 20 triangulation: true grayscale: true transparency: 0.1 refresh: false;
            mesh display_water_field scale: 20 triangulation: true color: rgb(0, 100, 255, 180) refresh: true;
            mesh display_obstacle_field scale: 20 triangulation: true color: rgb(139, 69, 19, 220) refresh: true;
            
            light #ambient intensity: 100;
            
            overlay position: {5, 5} size: {350, 150} background: #black transparency: 0.3 border: #white {
                draw "HANOI FLOODING SIMULATION" at: {10, 20} color: #yellow font: font("Arial", 14, #bold);
                draw ("Status: " + (simulation_active ? "RUNNING" : "STOPPED")) at: {10, 40} 
                    color: (simulation_active ? #green : #red);
                draw ("Water cells: " + world.get_active_water_count()) at: {10, 55} color: #cyan;
                draw ("Obstacles: " + world.get_active_obstacle_count()) at: {10, 70} color: #brown;
                draw ("Rain: " + (world.is_rain_active() ? "Active" : "None")) at: {10, 85} 
                    color: (world.is_rain_active() ? #blue : #gray);
                draw ("Cycle: " + total_cycles) at: {10, 100} color: #white;
            }
        }
        
        display "Simulation Metrics" type: java2D {
            chart "Water Spread and System Status" type: series {
                data "Water Cells" value: world.get_active_water_count() color: #blue;
                data "Edge Cells" value: world.get_edge_cell_count() color: #red;
                data "Rain Rate x100" value: world.get_rain_rate() * 100 color: #purple;
                data "Obstacles Under Attack x10" value: world.get_obstacles_under_attack() * 10 color: #orange;
            }
        }
        
        monitor "Simulation Active" value: simulation_active color: simulation_active ? #green : #red;
        monitor "Cycle" value: total_cycles;
        monitor "Water Cells" value: world.get_active_water_count();
        monitor "Edge Cells" value: world.get_edge_cell_count();
        monitor "Obstacle Cells" value: world.get_active_obstacle_count();
        monitor "Rain Active" value: world.is_rain_active() 
            color: world.is_rain_active() ? #blue : #gray;
        monitor "Water Coverage %" value: 
            ((world.get_active_water_count() / 
            (world.grid_width * world.grid_height)) * 100.0) with_precision 1;
    }
    
    // === USER COMMANDS ===
    user_command "Start Simulation" action: {ask world {do start_sim();}} category: "Simulation";
    user_command "Stop Simulation" action: {ask world {do stop_sim();}} category: "Simulation";
    user_command "Reset Simulation" action: {ask world {do reset_sim();}} category: "Simulation";
    
    user_command "Start Monsoon Rain" action: {ask world {do start_monsoon_rain();}} category: "Rain";
    user_command "Start Heavy Rain" action: {ask world {do start_heavy_rain();}} category: "Rain";
    user_command "Stop Rain" action: {ask world {do stop_rain_action();}} category: "Rain";
    
    user_command "Clear All Buildings" action: {ask world {do clear_all_buildings();}} category: "Obstacles";
    
    user_command "Quick Status" action: {ask world {do report_status();}} category: "Status";
    user_command "Detailed Status" action: {ask world {do report_detailed_status();}} category: "Status";
}