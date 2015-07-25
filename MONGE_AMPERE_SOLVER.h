#pragma once

#include "COMMON_DEFINITION.h"
#include "WORLD_DISCRETIZATION_2D.h"
#include "LEVELSET_2D.h"

class MONGE_AMPERE_SOLVER
{
public: // World Discretization
	WORLD_DISCRETIZATION_2D*			world_discretization;

public: // Main variables
	FIELD_STRUCTURE_2D<T>				solution;

public: // Convenient variables and references
	GRID_STRUCTURE_2D					base_grid;

public: // For multithreading
	MULTITHREADING*						multithreading;

public: // Constructor and Destructor
	MONGE_AMPERE_SOLVER(void)
		: multithreading(0)
	{}

	~MONGE_AMPERE_SOLVER(void)
	{}

public: // Initialization Function
	void InitializeFromBlock(const SCRIPT_BLOCK* monge_ampere_eqn_block, MULTITHREADING* multithreading_input)
	{
		multithreading = multithreading_input;

		// Grid Initialization
		base_grid.Initialize(world_discretization->world_grid);

		// Initialize Fields
		solution.Initialize(base_grid, 1, multithreading);
	}	
};