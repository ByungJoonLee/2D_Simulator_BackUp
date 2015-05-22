#pragma once

#include "EULERIAN_FLUID_SOLVER_2D.h"
#include "WORLD_DISCRETIZATION_2D.h"
#include <ctype.h>
#include <boost/chrono.hpp>

class SIMULATION_WORLD
{
public: // Primitive Solvers which are updated
	EULERIAN_FLUID_SOLVER_2D		eulerian_solver;

public: 
	WORLD_DISCRETIZATION_2D			world_discretization;

public: // Essential for script update
	std::string						script_abs_dir;
	std::string						app_abs_dir;
	std::string						script_file_path;
	std::string						script_base_name;

public: // Simulation Properties
	T								dt;
	
public: // Options for Simulation
	bool							air_water_simulation;
	bool							oil_water_simulation;
	bool							vortex_sheet_problem;

public: // Options for Oil Water Simulation
	bool							is_vertical, is_parallel;

public: // Wave length and Wave number for Oil Water Simulation
	T								A_0, alpha;

public: // Multithreading
	MULTITHREADING*					multithreading;

public: // Constructor and Destructor
	SIMULATION_WORLD(void)
		: multithreading(0), air_water_simulation(false), oil_water_simulation(false), is_vertical(false), is_parallel(false)
	{}

	~SIMULATION_WORLD(void)
	{}

public: // Initialization Functions
	void InitializeFromScript(const char* script)
	{
		DELETE_POINTER(multithreading);
		multithreading = new MULTITHREADING;

		SetCurrentDirectory(script_abs_dir.c_str());

		SCRIPT_READER script_reader(script);

		SCRIPT_BLOCK script_block_for_this = script_reader.FindBlock("SIMULATION_WORLD");

		dt = script_block_for_this.GetFloat("dt", (T)0.01);
		
		// Simulation Options
		air_water_simulation = script_block_for_this.GetBoolean("air_water_simulation", false);
		oil_water_simulation = script_block_for_this.GetBoolean("oil_water_simulation", false);
		vortex_sheet_problem = script_block_for_this.GetBoolean("vortex_sheet_problem", false);

		if (oil_water_simulation)
		{
			A_0 = script_block_for_this.GetFloat("A_0", (T)1);
			alpha = script_block_for_this.GetFloat("alpha", (T)1);
		}

		// Multithreading
		multithreading->Initialize(script_block_for_this.GetInteger("number_of_threads"));

		// Display
		cout << "---------------SIMULATION WORLD VARIABLES---------------" << endl; 
		cout << "Dt = " << dt << endl;

		if (air_water_simulation)
		{
			cout << "Air-Water Simulation is activated!" << endl;
		}
		if (oil_water_simulation)
		{
			cout << "Oil-Water Simulation is activated!" << endl;
		}

		if (oil_water_simulation)
		{
			A_0 = script_block_for_this.GetFloat("A_0", (T)1);
			alpha = script_block_for_this.GetFloat("alpha", (T)1);
		}

		cout << "Number of threads: " << multithreading->num_threads << endl;

		// Initialize world discretization
		if (air_water_simulation)
		{
			world_discretization.air_water_simulation = true;
			world_discretization.oil_water_simulation = false;
			world_discretization.Initialize(script_block_for_this.FindBlock("WORLD_DISCRETIZATION_AIR_WATER"));
		}
		if (oil_water_simulation)
		{
			world_discretization.air_water_simulation = false;
			world_discretization.oil_water_simulation = true;
			world_discretization.Initialize(script_block_for_this.FindBlock("WORLD_DISCRETIZATION_OIL_WATER"));
			is_vertical = world_discretization.is_vertical;
			is_parallel = world_discretization.is_parallel;
		}
		if (vortex_sheet_problem)
		{
			world_discretization.Initialize(script_block_for_this.FindBlock("WORLD_DISCRETIZATION"));
		}
		InitializeDynamicsSolversFromScript(script_reader);
		
		SetCurrentDirectory(app_abs_dir.c_str());
	}

	void InitializeDynamicsSolversFromScript(SCRIPT_READER& script_reader)
	{
		eulerian_solver.world_discretization = &world_discretization;
		if (air_water_simulation)
		{
			eulerian_solver.air_water_simulation = air_water_simulation;
			eulerian_solver.oil_water_simulation = oil_water_simulation;
			eulerian_solver.InitializeFromScriptBlock(world_discretization.world_grid, script_reader.FindBlock("FLUID_SOLVER_UNIFORM_AIR_WATER"));
		}
		if (oil_water_simulation)
		{
			eulerian_solver.air_water_simulation = air_water_simulation;
			eulerian_solver.oil_water_simulation = oil_water_simulation;
			eulerian_solver.InitializeFromScriptBlock(world_discretization.world_grid, script_reader.FindBlock("FLUID_SOLVER_UNIFORM_OIL_WATER"));
		}
		if (vortex_sheet_problem)
		{
			eulerian_solver.vortex_sheet_problem = vortex_sheet_problem;
			eulerian_solver.InitializeFromScriptBlock(world_discretization.world_grid, script_reader.FindBlock("FLUID_SOLVER_UNIFORM_AIR_WATER"));
		}
	}

	void AdvanceOneTimeStep(const T& dt_input)
	{
		eulerian_solver.AdvanceOneTimeStep(dt_input);
	}
};
		

