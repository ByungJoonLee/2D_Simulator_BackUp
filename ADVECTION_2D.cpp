#include "ADVECTION_2D.h"

// Initialization Function
void ADVECTION_2D::InitializeFromBlock(const SCRIPT_BLOCK& outer_block)
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

// Solver
void ADVECTION_2D::Solve_Velocity(const T& dt, const int& thread_id)
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

void ADVECTION_2D::Solve_Velocity(const T& dt)
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

void ADVECTION_2D::Solve_Levelset(const T& dt, const int& thread_id)
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

void ADVECTION_2D::Solve_Levelset(const T& dt)
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

void ADVECTION_2D::Solve_Vortex(const T& dt)
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

void ADVECTION_2D::Solve_Vortex(const T& dt, const int& thread_id)
{
	// Water Levelset Advection with MAC Grid
	if (use_5th_weno)
	{
		if (use_mac_grid)
		{
			ADVECTION_METHOD_2D<T>::WENO5th(vortex_signed_distance_field, scalar_field_ghost, velocity_field_mac_ghost_x, velocity_field_mac_ghost_y, dt, epsilon, multithreading, thread_id);
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
		ADVECTION_METHOD_2D<T>::ENO3rd(vortex_signed_distance_field, scalar_field_ghost, velocity_field_mac_ghost_x, velocity_field_mac_ghost_y, dt, multithreading, thread_id);
	}
}

void ADVECTION_2D::ReinitializationBySussman(const T& dt, FIELD_STRUCTURE_2D<T>& sign_function)
{
	ADVECTION_METHOD_2D<T>::WENO5thReinitialization(water_signed_distance_field, scalar_field_ghost, dt, epsilon, sign_function); 
	//ADVECTION_METHOD_2D<T>::GodunovReinitialization(water_signed_distance_field, scalar_field_ghost, dt, epsilon, sign_function, phi_0); 
	//ADVECTION_METHOD_2D<T>::SubcellFixedReinitialization(water_signed_distance_field, scalar_field_ghost, dt, epsilon, sign_function, phi_0); 
	water_levelset.FillGhostCellsContinuousDerivativesFrom(water_levelset.arr, true);
}
void ADVECTION_2D::ReinitializationBySussman(const T& dt, FIELD_STRUCTURE_2D<T>& sign_function, const int& thread_id)
{
	ADVECTION_METHOD_2D<T>::WENO5thReinitialization(water_signed_distance_field, scalar_field_ghost, dt, epsilon, sign_function, multithreading, thread_id); 
	//ADVECTION_METHOD_2D<T>::GodunovReinitialization(water_signed_distance_field, scalar_field_ghost, dt, epsilon, sign_function, phi_0); 
	//ADVECTION_METHOD_2D<T>::SubcellFixedReinitialization(water_signed_distance_field, scalar_field_ghost, dt, epsilon, sign_function, phi_0); 
	water_levelset.FillGhostCellsContinuousDerivativesFrom(water_levelset.arr, true);
}