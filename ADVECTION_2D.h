#pragma once

#include "LEVELSET_2D.h"
#include "ADVECTION_METHOD_2D.h"

class ADVECTION_2D
{
public: // Essential Data
	// Water Levelset
	LEVELSET_2D&					water_levelset;
	FIELD_STRUCTURE_2D<T>&			water_signed_distance_field;

	FIELD_STRUCTURE_2D<T>*			pressure_field;

	GRID_STRUCTURE_2D&				base_grid;

	FIELD_STRUCTURE_2D<VT>&			water_velocity_field;

	// Vortex Levelset
	LEVELSET_2D&					vortex_levelset;
	FIELD_STRUCTURE_2D<T>&			vortex_signed_distance_field;

	// For MAC GRID
	FIELD_STRUCTURE_2D<T>&			water_velocity_field_mac_x;
	FIELD_STRUCTURE_2D<T>&			water_velocity_field_mac_y;

	// Ghost Value
	FIELD_STRUCTURE_2D<T>&			scalar_field_ghost;
	FIELD_STRUCTURE_2D<VT>&			vector_field_ghost;
	FIELD_STRUCTURE_2D<T>&			velocity_field_mac_ghost_x;
	FIELD_STRUCTURE_2D<T>&			velocity_field_mac_ghost_y;

public: // Options
	// For Levelset Advection
	bool							use_1st_sl, use_5th_weno, use_3rd_eno, gradient_augmented_method;

	// For Velocity Advection
	bool							use_5th_weno_v, use_3rd_eno_v;

	// MAC grid
	bool							use_mac_grid;

	// For Pipe
	bool							is_vertical, is_parallel;

public: // Conditional Variable
	T								epsilon, epsilon_v;

public: // Multithreading
	MULTITHREADING&					multithreading;

public: // Constructor and Destructor
	ADVECTION_2D(LEVELSET_2D& water_levelset_input, FIELD_STRUCTURE_2D<T>& scalar_field_ghost_input, LEVELSET_2D& vortex_levelset_input, FIELD_STRUCTURE_2D<VT>& velocity_field_input, FIELD_STRUCTURE_2D<VT>& vector_field_ghost_input, FIELD_STRUCTURE_2D<T>& velocity_field_mac_x_input, FIELD_STRUCTURE_2D<T>& velocity_field_mac_y_input, FIELD_STRUCTURE_2D<T>& velocity_field_ghost_mac_x_input, FIELD_STRUCTURE_2D<T>& velocity_field_ghost_mac_y_input, MULTITHREADING& multithreading_input)
		:water_levelset(water_levelset_input), scalar_field_ghost(scalar_field_ghost_input), vortex_levelset(vortex_levelset_input), water_velocity_field(velocity_field_input), vector_field_ghost(vector_field_ghost_input), water_velocity_field_mac_x(velocity_field_mac_x_input), water_velocity_field_mac_y(velocity_field_mac_y_input), velocity_field_mac_ghost_x(velocity_field_ghost_mac_x_input), velocity_field_mac_ghost_y(velocity_field_ghost_mac_y_input),
		base_grid(water_levelset_input.grid), pressure_field(0), water_signed_distance_field(water_levelset_input.signed_distance_field), vortex_signed_distance_field(vortex_levelset.signed_distance_field), multithreading(multithreading_input)
	{
		use_1st_sl = false;
		use_5th_weno = false;
		use_3rd_eno = false;
		gradient_augmented_method = false;
		use_5th_weno_v = false;
		use_3rd_eno_v = false;
	}

	~ADVECTION_2D(void)
	{}

public: // Initialization Function
	void InitializeFromBlock(const SCRIPT_BLOCK& outer_block);

public: // Solver
	void Solve_Velocity(const T& dt, const int& thread_id);
	void Solve_Velocity(const T& dt);
	
	void Solve_Levelset(const T& dt, const int& thread_id);
	void Solve_Levelset(const T& dt);
	
	void Solve_Vortex(const T& dt);
	void Solve_Vortex(const T& dt, const int& thread_id);
	
	void ReinitializationBySussman(const T& dt, FIELD_STRUCTURE_2D<T>& sign_function);
	void ReinitializationBySussman(const T& dt, FIELD_STRUCTURE_2D<T>& sign_function, const int& thread_id);
};





