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
	void InitializeFromBlock(const SCRIPT_BLOCK& outer_block)
	{
		const SCRIPT_BLOCK& advection_block(outer_block.FindBlock("LEVELSET_ADVECTION"));
		use_1st_sl = advection_block.GetBoolean("use_1st_sl", false);
		use_5th_weno = advection_block.GetBoolean("use_5th_weno", false);
		use_3rd_eno = advection_block.GetBoolean("use_3rd_eno", false);
		gradient_augmented_method = advection_block.GetBoolean("gradient_augmented_method");

		epsilon = advection_block.GetFloat("epsilon", (float)1e-6);

		const SCRIPT_BLOCK& advection_block_v(outer_block.FindBlock("VELOCITY_ADVECTION"));
		use_5th_weno_v = advection_block_v.GetBoolean("use_5th_weno_v", false);
		use_3rd_eno_v = advection_block_v.GetBoolean("use_3rd_eno_v", false);

		epsilon_v = advection_block_v.GetFloat("epsilon_v", (float)1e-6);

		// Display
		cout << "---------------LEVELSET ADVECTION---------------" << endl;
		if (use_1st_sl)
		{
			cout << "use 1st Semi-Lagrangian: " << "True" << endl;
		}
		else
		{
			cout << "use 1st Semi-Lagrangian: " << "False" << endl;
		}

		if (use_5th_weno)
		{
			cout << "use 5th WENO: " << "True" << endl;
		}
		else
		{
			cout << "use 5th WENO: " << "False" << endl;
		}

		if (use_3rd_eno)
		{
			cout << "use 3rd ENO: " << "True" << endl;
		}
		else
		{
			cout << "use 3rd ENO: " << "False" << endl;
		}

		if (gradient_augmented_method)
		{
			cout << "Gradient Augmented Levelset Method: True" << endl;
		}
		else
		{
			cout << "Gradient Augmented Levelset Method: False" << endl;
		}
		
		cout << "Epsilon: " << epsilon << endl;

		cout << "---------------VELOCITY ADVECTION---------------" << endl;
		if (use_5th_weno_v)
		{
			cout << "use 5th WENO: " << "True" << endl;
		}
		else
		{
			cout << "use 5th WENO: " << "False" << endl;
		}

		if (use_3rd_eno_v)
		{
			cout << "use 3rd ENO: " << "True" << endl;
		}
		else
		{
			cout << "use 3rd ENO: " << "False" << endl;
		}

		cout << "Epsilon: " << epsilon_v << endl;
	}

public: // 
	void Solve_Velocity(const T& dt, const int& thread_id)
	{
		if (use_mac_grid)
		{
			// Prepare For an Advection Step
			velocity_field_mac_ghost_x.FillGhostCellsFrom(water_velocity_field_mac_x.array_for_this, true, thread_id);
			velocity_field_mac_ghost_y.FillGhostCellsFrom(water_velocity_field_mac_y.array_for_this, true, thread_id);

			/*FIELD_STRUCTURE_2D<T> velocity_field_mac_ghost_x_ini, velocity_field_mac_ghost_y_ini;
			velocity_field_mac_ghost_x_ini.Initialize(velocity_field_mac_ghost_x.grid, 2);
			velocity_field_mac_ghost_y_ini.Initialize(velocity_field_mac_ghost_y.grid, 2);

			GRID_ITERATION_2D(velocity_field_mac_ghost_x_ini.grid)
			{
				velocity_field_mac_ghost_x_ini(i, j) = velocity_field_mac_ghost_x(i, j);
			}

			GRID_ITERATION_2D(velocity_field_mac_ghost_y_ini.grid)
			{
				velocity_field_mac_ghost_y_ini(i, j) = velocity_field_mac_ghost_y(i, j);
			}

			velocity_field_mac_ghost_x_ini.FillGhostCellsFrom(water_velocity_field_mac_x.array_for_this, true);
			velocity_field_mac_ghost_y_ini.FillGhostCellsFrom(water_velocity_field_mac_y.array_for_this, true);*/

			// Velocity Advection with MAC grid
			if (use_5th_weno_v)
			{
				ADVECTION_METHOD_2D<T>::WENO5th(water_velocity_field_mac_x, velocity_field_mac_ghost_x, velocity_field_mac_ghost_x, velocity_field_mac_ghost_y, dt, epsilon_v);
				ADVECTION_METHOD_2D<T>::WENO5th(water_velocity_field_mac_y, velocity_field_mac_ghost_y, velocity_field_mac_ghost_x, velocity_field_mac_ghost_y, dt, epsilon_v);
			}
			else if (use_3rd_eno_v)
			{
				ADVECTION_METHOD_2D<T>::ENO3rd(water_velocity_field_mac_x, velocity_field_mac_ghost_x, velocity_field_mac_ghost_x, velocity_field_mac_ghost_y, dt, multithreading, thread_id);
				ADVECTION_METHOD_2D<T>::ENO3rd(water_velocity_field_mac_y, velocity_field_mac_ghost_y, velocity_field_mac_ghost_x, velocity_field_mac_ghost_y, dt, multithreading, thread_id);
				/*ADVECTION_METHOD_2D<T>::ENO3rd(water_velocity_field_mac_x, velocity_field_mac_ghost_x, velocity_field_mac_ghost_x_ini, velocity_field_mac_ghost_y_ini, dt);
				ADVECTION_METHOD_2D<T>::ENO3rd(water_velocity_field_mac_y, velocity_field_mac_ghost_y, velocity_field_mac_ghost_x_ini, velocity_field_mac_ghost_y_ini, dt);*/
			}
		}
	}

	void Solve_Velocity(const T& dt)
	{
		// Prepare For an Advection Step
		velocity_field_mac_ghost_x.FillGhostCellsFrom(water_velocity_field_mac_x.array_for_this, true);
		velocity_field_mac_ghost_y.FillGhostCellsFrom(water_velocity_field_mac_y.array_for_this, true);
		//velocity_field_mac_ghost_x.FillGhostCellsContinuousDerivativesFrom(water_velocity_field_mac_x.array_for_this, true);
		//velocity_field_mac_ghost_y.FillGhostCellsContinuousDerivativesFrom(water_velocity_field_mac_y.array_for_this, true);

		FIELD_STRUCTURE_2D<T> velocity_field_mac_ghost_x_ini, velocity_field_mac_ghost_y_ini;
		velocity_field_mac_ghost_x_ini.Initialize(velocity_field_mac_ghost_x.grid, 2);
		velocity_field_mac_ghost_y_ini.Initialize(velocity_field_mac_ghost_y.grid, 2);

		GRID_ITERATION_2D(velocity_field_mac_ghost_x_ini.grid)
		{
			velocity_field_mac_ghost_x_ini(i, j) = velocity_field_mac_ghost_x(i, j);
		}

		GRID_ITERATION_2D(velocity_field_mac_ghost_y_ini.grid)
		{
			velocity_field_mac_ghost_y_ini(i, j) = velocity_field_mac_ghost_y(i, j);
		}

		velocity_field_mac_ghost_x_ini.FillGhostCellsFrom(water_velocity_field_mac_x.array_for_this, true);
		velocity_field_mac_ghost_y_ini.FillGhostCellsFrom(water_velocity_field_mac_y.array_for_this, true);
		
		// Velocity Advection with MAC grid
		if (use_5th_weno_v)
		{
			ADVECTION_METHOD_2D<T>::WENO5th(water_velocity_field_mac_x, velocity_field_mac_ghost_x, velocity_field_mac_ghost_x, velocity_field_mac_ghost_y, dt, epsilon_v);
			ADVECTION_METHOD_2D<T>::WENO5th(water_velocity_field_mac_y, velocity_field_mac_ghost_y, velocity_field_mac_ghost_x, velocity_field_mac_ghost_y, dt, epsilon_v);
		}
		else if (use_3rd_eno_v)
		{
			/*ADVECTION_METHOD_2D<T>::ENO3rd(water_velocity_field_mac_x, velocity_field_mac_ghost_x, velocity_field_mac_ghost_x, velocity_field_mac_ghost_y, dt);
			ADVECTION_METHOD_2D<T>::ENO3rd(water_velocity_field_mac_y, velocity_field_mac_ghost_y, velocity_field_mac_ghost_x, velocity_field_mac_ghost_y, dt);*/
			ADVECTION_METHOD_2D<T>::ENO3rd(water_velocity_field_mac_x, velocity_field_mac_ghost_x, velocity_field_mac_ghost_x_ini, velocity_field_mac_ghost_y_ini, dt);
			ADVECTION_METHOD_2D<T>::ENO3rd(water_velocity_field_mac_y, velocity_field_mac_ghost_y, velocity_field_mac_ghost_x_ini, velocity_field_mac_ghost_y_ini, dt);
		}
        water_velocity_field_mac_x.FillGhostCellsFrom(water_velocity_field_mac_x.array_for_this, true);
        water_velocity_field_mac_y.FillGhostCellsFrom(water_velocity_field_mac_y.array_for_this, true);
	}

	void Solve_Levelset(const T& dt, const int& thread_id)
	{
		// Water Levelset Advection with MAC Grid
		if (use_5th_weno)
		{
			ADVECTION_METHOD_2D<T>::WENO5th(water_signed_distance_field, scalar_field_ghost, velocity_field_mac_ghost_x, velocity_field_mac_ghost_y, dt, epsilon, multithreading, thread_id);
		}
		else if (use_3rd_eno)
		{
			ADVECTION_METHOD_2D<T>::ENO3rd(water_signed_distance_field, scalar_field_ghost, velocity_field_mac_ghost_x, velocity_field_mac_ghost_y, dt, multithreading, thread_id);
		}
		else if (use_1st_sl)
		{
		
		}
	}

	void Solve_Levelset(const T& dt)
	{
		// Water Levelset Advection with MAC Grid
		if (use_5th_weno)
		{
			if (use_mac_grid)
			{
				// Periodic Boundary Condition
				scalar_field_ghost.FillGhostCellsFrom(water_signed_distance_field.array_for_this, true);

				if (is_vertical)
				{
					for (int i = scalar_field_ghost.i_start; i <= scalar_field_ghost.i_end; i++)
					{
						for (int ghost = 1; ghost <= scalar_field_ghost.ghost_width; ghost++)
						{
							scalar_field_ghost(i, scalar_field_ghost.j_start - ghost) = scalar_field_ghost(i, scalar_field_ghost.j_end - ghost);
							scalar_field_ghost(i, scalar_field_ghost.j_end + ghost) = scalar_field_ghost(i, scalar_field_ghost.j_start + ghost);
						}
					}
				}
				
				ADVECTION_METHOD_2D<T>::WENO5th(water_signed_distance_field, scalar_field_ghost, velocity_field_mac_ghost_x, velocity_field_mac_ghost_y, dt, epsilon);
			}
			else
			{
				// Boundary Condition
				scalar_field_ghost.FillGhostCellsFrom(water_signed_distance_field.array_for_this, true);

				for (int k = scalar_field_ghost.j_start; k <= scalar_field_ghost.j_end; k++)
				{
					for (int ghost = 1; ghost <= scalar_field_ghost.ghost_width; ghost++)
					{
						scalar_field_ghost(scalar_field_ghost.i_start - ghost, k) = scalar_field_ghost(scalar_field_ghost.i_end - ghost, k);
						scalar_field_ghost(scalar_field_ghost.i_end + ghost, k) = scalar_field_ghost(scalar_field_ghost.i_start + ghost, k);
					}
				}
				
				for (int k = scalar_field_ghost.i_start; k <= scalar_field_ghost.i_end; k++)
				{
					for (int ghost = 1; ghost <= scalar_field_ghost.ghost_width; ghost++)
					{
						//scalar_field_ghost(k, scalar_field_ghost.j_start - ghost) = scalar_field_ghost(k, scalar_field_ghost.j_end - ghost) + 2;
						scalar_field_ghost(k, scalar_field_ghost.j_end + ghost) = scalar_field_ghost(k, scalar_field_ghost.j_start + ghost) + 2;
					}
				}

				ADVECTION_METHOD_2D<T>::WENO5th(water_signed_distance_field, scalar_field_ghost, vector_field_ghost, dt, epsilon);
			}
			
		}
		if (use_3rd_eno)
		{
			ADVECTION_METHOD_2D<T>::ENO3rd(water_signed_distance_field, scalar_field_ghost, velocity_field_mac_ghost_x, velocity_field_mac_ghost_y, dt);
		}
		if (use_1st_sl)
		{
		
		}
		if (gradient_augmented_method)
		{
			ADVECTION_METHOD_2D<T>::Gradient_Augmented_Method(water_levelset, scalar_field_ghost, vector_field_ghost, dt);
		}
	}

	void Solve_Vortex(const T& dt)
	{
		// Water Levelset Advection with MAC Grid
		if (use_5th_weno)
		{
			if (use_mac_grid)
			{
				ADVECTION_METHOD_2D<T>::WENO5th(vortex_signed_distance_field, scalar_field_ghost, velocity_field_mac_ghost_x, velocity_field_mac_ghost_y, dt, epsilon);
			}
			else
			{
				// Boundary Condition
				scalar_field_ghost.FillGhostCellsFrom(vortex_signed_distance_field.array_for_this, true);

				for (int k = scalar_field_ghost.j_start; k <= scalar_field_ghost.j_end; k++)
				{
					for (int ghost = 1; ghost <= scalar_field_ghost.ghost_width; ghost++)
					{
						scalar_field_ghost(scalar_field_ghost.i_start - ghost, k) = scalar_field_ghost(scalar_field_ghost.i_end - ghost, k);
						scalar_field_ghost(scalar_field_ghost.i_end + ghost, k) = scalar_field_ghost(scalar_field_ghost.i_start + ghost, k);
					}
				}

				for (int k = scalar_field_ghost.i_start; k <= scalar_field_ghost.i_end; k++)
				{
					for (int ghost = 1; ghost <= scalar_field_ghost.ghost_width; ghost++)
					{
						scalar_field_ghost(k, scalar_field_ghost.j_start - ghost) = scalar_field_ghost(k, scalar_field_ghost.j_end - ghost);
						scalar_field_ghost(k, scalar_field_ghost.j_end + ghost) = scalar_field_ghost(k, scalar_field_ghost.j_start + ghost);
					}
				}

				ADVECTION_METHOD_2D<T>::WENO5th(vortex_signed_distance_field, scalar_field_ghost, vector_field_ghost, dt, epsilon);
			}
			
		}
		if (use_3rd_eno)
		{
			ADVECTION_METHOD_2D<T>::ENO3rd(vortex_signed_distance_field, scalar_field_ghost, velocity_field_mac_ghost_x, velocity_field_mac_ghost_y, dt);
		}
	}

	void ReinitializationBySussman(const T& dt, const int& thread_id)
	{
		ADVECTION_METHOD_2D<T>::WENO5thReinitialization(water_signed_distance_field, scalar_field_ghost, dt, epsilon, base_grid.dx, multithreading, thread_id); 
		water_levelset.FillGhostCellsContinuousDerivativesFrom(water_levelset.arr, true, thread_id);
	}

	void ReinitializationBySussman(const T& dt, FIELD_STRUCTURE_2D<T>& sign_function, FIELD_STRUCTURE_2D<T>& phi_0)
	{
		ADVECTION_METHOD_2D<T>::WENO5thReinitialization(water_signed_distance_field, scalar_field_ghost, dt, epsilon, sign_function); 
		//ADVECTION_METHOD_2D<T>::GodunovReinitialization(water_signed_distance_field, scalar_field_ghost, dt, epsilon, sign_function, phi_0); 
		//ADVECTION_METHOD_2D<T>::SubcellFixedReinitialization(water_signed_distance_field, scalar_field_ghost, dt, epsilon, sign_function, phi_0); 
		water_levelset.FillGhostCellsContinuousDerivativesFrom(water_levelset.arr, true);
	}
};





