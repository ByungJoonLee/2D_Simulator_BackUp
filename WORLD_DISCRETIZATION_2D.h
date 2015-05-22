#pragma once

#include "COMMON_DEFINITION.h"
#include "ARRAY_2D.h"
#include "SCRIPT_READER.h"
#include "GRID_STRUCTURE_2D.h"

class WORLD_DISCRETIZATION_2D
{
public: // Essential Data
	GRID_STRUCTURE_2D		world_grid;
	GRID_STRUCTURE_2D		world_grid_ghost;

	int						ghost_width;

	bool					use_grid_uniform;

	// Options For Simulation
	bool					air_water_simulation;
	bool					oil_water_simulation;

	bool					large_bubble, small_bubble;
	bool					is_vertical, is_parallel;

public: // Speedup constants
	T						dx, dy;
	T						dx_over_two;
	T						dy_over_two;


public: // Constructor and Destructor
	WORLD_DISCRETIZATION_2D(void)
		: use_grid_uniform(true), air_water_simulation(false), oil_water_simulation(false), large_bubble(false), small_bubble(false), is_vertical(false), is_parallel(false)
	{}

	~WORLD_DISCRETIZATION_2D(void)
	{}

public: // Initialization Function
	void Initialize(const SCRIPT_BLOCK& script_block)
	{
		ghost_width = script_block.GetInteger("ghost_width", 1);
		large_bubble = script_block.GetBoolean("large_bubble", (bool)false);
		small_bubble = script_block.GetBoolean("small_bubble", (bool)false);

		if (air_water_simulation)
		{
			if (large_bubble)
			{
				world_grid.InitializeFromBlock(script_block.FindBlock("GRID_STRUCTURE_2D_LARGE"));
				world_grid_ghost.Initialize(world_grid.Enlarged(ghost_width));
				use_grid_uniform = script_block.GetBoolean("use_grid_uniform", (bool)true);
			}
			if (small_bubble)
			{
				world_grid.InitializeFromBlock(script_block.FindBlock("GRID_STRUCTURE_2D_SMALL"));
				world_grid_ghost.Initialize(world_grid.Enlarged(ghost_width));
				use_grid_uniform = script_block.GetBoolean("use_grid_uniform", (bool)true);
			} 
		}
		else if (oil_water_simulation)
		{
			is_vertical = script_block.FindBlock("PIPE_OPTIONS").GetBoolean("is_vertical", (bool)false);
			is_parallel = script_block.FindBlock("PIPE_OPTIONS").GetBoolean("is_parallel", (bool)false);

			if (is_vertical)
			{
				world_grid.InitializeFromBlock(script_block.FindBlock("GRID_STRUCTURE_2D_VERTICAL"));
				cout << "Vertical Pipe Simulation is activated!" << endl;
				world_grid_ghost.Initialize(world_grid.Enlarged(ghost_width));
				use_grid_uniform = script_block.GetBoolean("use_grid_uniform", (bool)true);
			}
			if (is_parallel)
			{
				world_grid.InitializeFromBlock(script_block.FindBlock("GRID_STRUCTURE_2D_PARALLEL"));
				cout << "Parallel Pipe Simulation is activated!" << endl;
				world_grid_ghost.Initialize(world_grid.Enlarged(ghost_width));
				use_grid_uniform = script_block.GetBoolean("use_grid_uniform", (bool)true);
			} 
		}
		else
		{
			world_grid.InitializeFromBlock(script_block.FindBlock("GRID_STRUCTURE_2D"));
			world_grid_ghost.Initialize(world_grid.Enlarged(ghost_width));
			use_grid_uniform = script_block.GetBoolean("use_grid_uniform", (bool)true);
		}			

		// Display
		cout << "---------------DISCRETIZATION VARIABLES---------------" << endl;
		cout << "Start Indices = " << " (" << world_grid.i_start << " ," << world_grid.j_start << ") " << endl;
		cout << "Grid Resolution = " << " (" << world_grid.i_res << " ," << world_grid.j_res << ") " << endl;

		InitializeSpeedupConstants();
	}

	void InitializeSpeedupConstants()
	{
		dx = world_grid.dx;
		dy = world_grid.dy;

		dx_over_two = dx*(T)0.5;
		dy_over_two = dy*(T)0.5;
	}
};

