#pragma once

#include "COMMON_DEFINITION.h"
#include "FIELD_STRUCTURE_2D.h"
#include "DYNAMIC_ARRAY.h"

class LEVELSET_2D
{
public: // Essential Data
	FIELD_STRUCTURE_2D<T>			signed_distance_field;
	FIELD_STRUCTURE_2D<T>			scalar_field_ghost;
	
	// For test
	FIELD_STRUCTURE_2D<T>			exact_solution;

	GRID_STRUCTURE_2D&				grid;
	ARRAY<GRID_STRUCTURE_2D>&		partial_grids;
	ARRAY<GRID_STRUCTURE_2D>&		partial_grids_ghost;
	int&							ghost_width;
	ARRAY_2D<T>						&phi, &arr, phi_true;

	FIELD_STRUCTURE_2D<VT>			normal;
	FIELD_STRUCTURE_2D<VT>			tangential;
	FIELD_STRUCTURE_2D<T>			curvature;

	// For Gradient Agmented Method
	FIELD_STRUCTURE_2D<VT>			gradient;

	int								sweep_direction;

public: // Speedup Variables
	ARRAY_2D<bool>					fixed;

public: // The way of calculating curvature vector
	bool							curvature_by_normal_vector;
	bool							curvature_by_levelset;

public: // Options for Axisymmetric object
	bool							is_axisymmetric;

public: // Multithreading
	MULTITHREADING*					multithreading;

public: // Constructors and Destructor
	LEVELSET_2D(void)
		: grid(signed_distance_field.grid), partial_grids(signed_distance_field.partial_grids), partial_grids_ghost(signed_distance_field.partial_grids_ghost), phi(signed_distance_field.array_for_this), arr(signed_distance_field.array_for_this), phi_true(signed_distance_field.array_for_this), sweep_direction(0), ghost_width(signed_distance_field.ghost_width), is_axisymmetric(false)
	{}
	
	~LEVELSET_2D(void)
	{}

public: // Initialization Functions
	void Initialize(const GRID_STRUCTURE_2D& grid_input, const int& ghost_width_input, MULTITHREADING* multithreading_input = 0)
	{
		assert(ghost_width_input >= 2);

		multithreading = multithreading_input;
		signed_distance_field.Initialize(grid_input, ghost_width_input, true, false, multithreading_input);
		scalar_field_ghost.Initialize(grid_input, ghost_width_input, true, false, multithreading_input);
		exact_solution.Initialize(grid_input, ghost_width_input, true, false, multithreading_input);

		fixed.Initialize(grid.i_start - ghost_width, grid.j_start - ghost_width, grid.i_res + 2*ghost_width, grid.j_res + 2*ghost_width, false);

		normal.Initialize(grid_input, 2, false, true, multithreading_input);
		tangential.Initialize(grid_input, 1, false, true, multithreading_input);
		curvature.Initialize(grid_input, 2, true, false, multithreading_input);

		phi.AssignAllValues(grid.dx);
		phi_true.Initialize(phi.i_start, phi.j_start, phi.i_res, phi.j_res, true);

		gradient.Initialize(grid_input, 2, false, true, multithreading_input);

		sweep_direction = 0;

		curvature_by_normal_vector = false;
		curvature_by_levelset = false;
	}

	void Initialize(const int& i_start_input, const int& j_start_input, const int& i_res_input, const int& j_res_input, const T& x_min_input, const T& y_min_input, const T& x_max_input, const T& y_max_input, const int& ghost_width_input, MULTITHREADING* multithreading_input = 0)
	{
		assert(ghost_width_input >= 2);

		multithreading = multithreading_input;
		
		signed_distance_field.Initialize(i_res_input, j_res_input, i_start_input, j_start_input, x_min_input, y_min_input, x_max_input, y_max_input, ghost_width_input, true, false, multithreading_input);
		scalar_field_ghost.Initialize(i_res_input, j_res_input, i_start_input, j_start_input, x_min_input, y_min_input, x_max_input, y_max_input, ghost_width_input, true, false, multithreading_input);
		exact_solution.Initialize(i_res_input, j_res_input, i_start_input, j_start_input, x_min_input, y_min_input, x_max_input, y_max_input, ghost_width_input, true, false, multithreading_input);

		fixed.Initialize(grid.i_start - ghost_width, grid.j_start - ghost_width, grid.i_res + 2*ghost_width, grid.j_res + 2*ghost_width, false);

		normal.Initialize(i_res_input + 2, j_res_input + 2, i_start_input - 1, j_start_input - 1, x_min_input - grid.dx, y_min_input - grid.dy, x_max_input + grid.dx, y_max_input + grid.dy, 0, false, true, multithreading_input);
		tangential.Initialize(i_res_input + 2, j_res_input + 2, i_start_input - 1, j_start_input - 1, x_min_input - grid.dx, y_min_input - grid.dy, x_max_input + grid.dx, y_max_input + grid.dy, 0, false, true, multithreading_input);
		curvature.Initialize(i_res_input, j_res_input, i_start_input, j_start_input, x_min_input, y_min_input, x_max_input, y_max_input, 0, true, false, multithreading_input);

		phi.AssignAllValues(grid.dx);
		phi_true.Initialize(phi.i_start, phi.j_start, phi.i_res, phi.j_res, true);

		gradient.Initialize(i_res_input + 2, j_res_input + 2, i_start_input - 1, j_start_input - 1, x_min_input - grid.dx, y_min_input - grid.dy, x_max_input + grid.dx, y_max_input + grid.dy, 0, true, false, multithreading_input);

		sweep_direction = 0;

		curvature_by_normal_vector = false;
		curvature_by_levelset = false;
	}

public: // Operator Overloading
	inline T& operator ()(const VI& ix) const
	{
		return phi(ix);
	}

	inline T& operator ()(const int& i, const int& j) const
	{
		return phi(i, j);
	}

	inline T operator ()(const VT& position) const
	{
		return BiLinearInterpolation(position);
	}

public: // Member Functions
	inline VT CellCenter(const int& i, const int& j) const
	{
		return grid.CellCenter(i, j);
	}

	inline VT GridPoint(const int& i, const int& j) const
	{
		return grid.GridPoint(i, j);
	}

	inline T BiLinearInterpolation(const VT& position) const
	{
		return signed_distance_field.BilinearInterpolation(position);
	}

	void Translate(const VT& deviation)
	{
		signed_distance_field.Translate(deviation);
		normal.Translate(deviation);
		tangential.Translate(deviation);
		curvature.Translate(deviation);
	}

	void SweepFill(const VI& start_idx, const VI& end_idx, const T& target_value, const T& assign_value)
	{
		int i = start_idx.i;
		int j = start_idx.j;
		
		if (start_idx.j < end_idx.j)
		{
			for (; j <= end_idx.j ; j++)
			{
				T& phi_ = phi(i, j);
				if (phi_ == target_value)
				{
					phi_ = assign_value;
				}
				else
				{
					break;
				}
			}
		}
		else
		{
			for (; j >= end_idx.j ; j--)
			{
				T& phi_ = phi(i, j);
				if (phi_ == target_value)
				{
					phi_ = assign_value;
				}
				else
				{
					break;
				}
			}
		}
	}

	inline bool InsideOBB(const VT& position) const
	{
		return grid.Inside(position);
	}

	inline const T SignedDistance(const VT& position) const
	{
		return BiLinearInterpolation(position);
	}

	void AssignAllValuesLevelset(const T& value)
	{
		for (int j = grid.j_start; j <= grid.j_end; j++)
		{
			for (int i = grid.i_start; i <= grid.i_end; i++)
			{
				phi(i, j) = value;
			}
		}
	}

	void CopyAllValuesFrom(const LEVELSET_2D& levelset_input)
	{
		for (int j = grid.j_start; j <= grid.j_end; j++)
		{
			for (int i = grid.i_start; i <= grid.i_end; i++)
			{
				phi(i, j) = levelset_input.phi(i, j);
			}
		}
	}

	void FillGhostCellsFrom(ARRAY_2D<T>& phi_real, const bool& copy_real_data)
	{
		signed_distance_field.FillGhostCellsFrom(phi_real, copy_real_data);
	}

	void FillGhostCellsFromPointer(ARRAY_2D<T>* phi_real, const bool& copy_real_data)
	{
		signed_distance_field.FillGhostCellsFrom(*phi_real, copy_real_data);
	}

	void FillGhostCellsContinuousDerivativesFrom(ARRAY_2D<T>& phi_real, const bool& copy_real_data)
	{
		signed_distance_field.FillGhostCellsContinuousDerivativesFrom(phi_real, copy_real_data);
	}

	void FillGhostCellsContinuousDerivativesFrom(ARRAY_2D<T>& phi_real, const bool& copy_real_data, const int& thread_id)
	{
		signed_distance_field.FillGhostCellsContinuousDerivativesFrom(phi_real, copy_real_data, thread_id);
	}

	void FillGhostCellsContinuousDerivativesFromPointer(ARRAY_2D<T>* phi_real, const bool& copy_real_data)
	{
		signed_distance_field.FillGhostCellsContinuousDerivativesFrom(*phi_real, copy_real_data);
	}

	void FastSweepingMethodOriginal(const int& sweep_number)
	{
		// Initialize fixed field as false
		signed_distance_field.AssignAllValue<bool>(fixed, false);

		const T nb_width = (T)2*grid.dx;

		// Fix the levelset value
		GRID_ITERATION_2D(grid)
		{
			/*T& phi_center(phi(i, j));
			bool& fix_center(fixed(i, j));

			if (phi_center > (T)0)
			{
				if (phi(i+1, j) <= (T)0)
				{
					fix_center = true;
					fixed(i+1, j) = true;
				}
				else if (phi(i-1, j) <= (T)0)
				{
					fix_center = true;
					fixed(i-1, j) = true;
				}
				else if (phi(i, j+1) <= (T)0)
				{
					fix_center = true;
					fixed(i, j+1) = true;
				}
				else if (phi(i, j-1) <= (T)0)
				{
					fix_center = true;
					fixed(i, j-1) = true;
				}
			}
			else
			{
				if (phi(i+1, j) > (T)0)
				{
					fix_center = true;
					fixed(i+1, j) = true;
				}
				else if (phi(i-1, j) > (T)0)
				{
					fix_center = true;
					fixed(i-1, j) = true;
				}
				else if (phi(i, j+1) > (T)0)
				{
					fix_center = true;
					fixed(i, j+1) = true;
				}
				else if (phi(i, j-1) > (T)0)
				{
					fix_center = true;
					fixed(i, j-1) = true;
				}
			}*/
			
			if (abs(arr(i, j)) < nb_width)
			{
				fixed(i, j) = true;
			}
			else
			{
				fixed(i, j) = false;
			}
		}
		
		int i_start(grid.i_start), j_start(grid.j_start);
		int i_end(grid.i_end), j_end(grid.j_end);

		FillGhostCellsContinuousDerivativesFrom(phi, true);
		//FillGhostCellsFrom(phi, false);

		sweep_direction = 0;

		while (sweep_direction != sweep_number)
		{
			switch (sweep_direction % 4)
			{
			case 0:
				Sweep(i_start, j_start, i_end, j_end);
				break;
			
			case 1:
				Sweep(i_end, j_start, i_start, j_end);
				break;

			case 2:
				Sweep(i_end, j_end, i_start, j_start);
				break;

			case 3:
				Sweep(i_start, j_end, i_end, j_start);
				break;
			}

			sweep_direction++;
		}
		
		FillGhostCellsContinuousDerivativesFrom(phi, false);
	}

	void FastSweepingMethod(void)
	{
		const T nb_width = (T)20*grid.dx;

		signed_distance_field.AssignAllValue<bool>(fixed, false);

		for (int j = grid.j_start; j <= grid.j_end ; j++)
		{
			for (int i = grid.i_start; i <= grid.i_end; i++)
			{
				bool& fix_center(fixed(i, j));

				if (fix_center == true)
				{
					continue;
				}

				T& phi_center(phi(i, j));

				// Fix values of the cells too far from the interface
				
				if (abs(phi_center) > nb_width)
				{
					fix_center = true;
					continue;
				}
				
				// Fix (do not update) interfacial cells
				if (phi_center > (T)0)
				{
					if (phi(i+1, j) <= (T)0)
					{
						fix_center = true;
						fixed(i+1, j) = true;
					}
					else if (phi(i-1, j) <= (T)0)
					{
						fix_center = true;
						fixed(i-1, j) = true;
					}
					else if (phi(i, j+1) <= (T)0)
					{
						fix_center = true;
						fixed(i, j+1) = true;
					}
					else if (phi(i, j-1) <= (T)0)
					{
						fix_center = true;
						fixed(i, j-1) = true;
					}
				}
				else
				{
					if (phi(i+1, j) > (T)0)
					{
						fix_center = true;
						fixed(i+1, j) = true;
					}
					else if (phi(i-1, j) > (T)0)
					{
						fix_center = true;
						fixed(i-1, j) = true;
					}
					else if (phi(i, j+1) > (T)0)
					{
						fix_center = true;
						fixed(i, j+1) = true;
					}
					else if (phi(i, j-1) > (T)0)
					{
						fix_center = true;
						fixed(i, j-1) = true;
					}
				}

				// Assign large magnitude values to non-interfacial cells
				if (fix_center == false)
				{
					if (phi_center > (T)0)
					{
						phi_center = nb_width;
					}
					else
					{
						phi_center = -nb_width;	
					}
				}
			}
		}

		FillGhostCellsFrom(phi, false);

		int i_start(grid.i_start), j_start(grid.j_start);
		int i_end(grid.i_end), j_end(grid.j_end);

		while (sweep_direction != 4)
		{
			switch (sweep_direction % 4)
			{
			case 0:
				Sweep(i_start, j_start, i_end, j_end);
				break;
			
			case 1:
				Sweep(i_end, j_start, i_start, j_end);
				break;

			case 2:
				Sweep(i_end, j_end, i_start, j_start);
				break;

			case 3:
				Sweep(i_start, j_end, i_end, j_start);
				break;
			}

			sweep_direction++;
		}
		
	}
		
	void Sweep(const int& i_start, const int& j_start, const int& i_end, const int& j_end)
	{
		const T one_over_two((T)1/(T)2);
		const T h(grid.dx), hh(h*h), hh2(hh*(T)2);

		T a, b;
		T a1, a2;
		T update;

		// Traverse orders
		int ip(1), jp(1);						// Increasing order
		
		if (i_start > i_end)
		{
			ip = -1;							// Decreasing order
		}
		
		if (j_start > j_end)
		{
			jp = -1;							// Decreasing order
		}

		int i, j;
		j = j_start - jp;
		while (true)
		{
			j += jp;
			i = i_start - ip;
			while (true)
			{
				i += ip;

				T& phi_center(phi(i, j));
				if (fixed(i, j) == false)
				{
					if (phi_center > (T)0)
					{
						// From the paper of Hongkai Zhao
						a = MIN(phi(i-1, j), phi(i+1, j));
						b = MIN(phi(i, j-1), phi(i, j+1));
						
						INCREASING_SORT2(a, b, a1, a2);
						
						if (abs(a1 - a2) >= h)
						{
							update = a1 + h;
						}
						else
						{
							update = one_over_two*(a1 + a2 + sqrt(hh2 - POW2(a1- a2)));
						}

						phi_center = MIN(update, phi_center);
					}
					else if (phi_center <= (T)0)
					{
						// From the paper of Hongkai Zhao
						a = MAX(phi(i-1, j), phi(i+1, j));
						b = MAX(phi(i, j-1), phi(i, j+1));
					
						INCREASING_SORT2(a, b, a1, a2);

						if (abs(a1 - a2) >= h)
						{
							update = a2 - h;
						}
						else
						{
							update = one_over_two*(a2 + a1 - sqrt(hh2 - POW2(a2 - a1)));
						}

						phi_center = MAX(update, phi_center);
					}
				}
				if (i == i_end)
				{
					break;
				}
			}
			if (j == j_end)
			{
				break;
			}
		}
	}

	void ComputeNormals()
	{
		for (int j = grid.j_start; j <= grid.j_end; j++)
		{
			for (int i = grid.i_start; i <= grid.i_end; i++)
			{
				VT& nor(normal(i, j));

				T nor_x = (phi(i+1, j) - phi(i-1, j))*signed_distance_field.one_over_2dx, nor_y = (phi(i, j+1) - phi(i, j-1))*signed_distance_field.one_over_2dy;
				T mag_of_normal = sqrt(POW2(nor_x) + POW2(nor_y));

				if (mag_of_normal != 0)
				{
					nor.x = nor_x;
					nor.y = nor_y;
				}
				else
				{
					nor.x = (phi(i+1, j) - phi(i, j))*signed_distance_field.one_over_dx;
					nor.y = (phi(i, j+1) - phi(i, j))*signed_distance_field.one_over_dy;
				}
				nor.MakeThisUnit();			
			}	
		}

		/*GRID_ITERATION_2D(grid)
		{
			VT& nor(normal(i, j));

			T nor_x_r = (phi(i+1, j) - phi(i, j))*signed_distance_field.one_over_dx, nor_x_l = (phi(i, j) - phi(i - 1, j))*signed_distance_field.one_over_dx, nor_y_u = (phi(i, j+1) - phi(i, j))*signed_distance_field.one_over_dy, nor_y_d = (phi(i, j) - phi(i, j-1))*signed_distance_field.one_over_dy;
			T nor_x = (T)0.5*(nor_x_r + nor_x_l), nor_y = (T)0.5*(nor_y_u + nor_y_d);
			T mag_of_normal = sqrt(POW2(nor_x) + POW2(nor_y));

			if (mag_of_normal != 0)
			{
				nor.x = nor_x;
				nor.y = nor_y;
			}
			else
			{
				nor.x = (phi(i+1, j) - phi(i, j))*signed_distance_field.one_over_dx;
				nor.y = (phi(i, j+1) - phi(i, j))*signed_distance_field.one_over_dy;
			}
			nor.MakeThisUnit();				
		}*/
	}

	void ComputeGradient()
	{
		for (int j = grid.j_start; j <= grid.j_end; j++)
		{
			for (int i = grid.i_start; i <= grid.i_end; i++)
			{
				gradient(i, j).x = (phi(i+1, j) - phi(i-1, j))*signed_distance_field.one_over_2dx;
				gradient(i, j).y = (phi(i, j+1) - phi(i, j-1))*signed_distance_field.one_over_2dy;
			}	
		}

		//gradient.FillGhostCellsContinuousDerivativesFrom(gradient.array_for_this, true);

		/*GRID_ITERATION_2D(grid)
		{
			VT& nor(normal(i, j));

			T nor_x_r = (phi(i+1, j) - phi(i, j))*signed_distance_field.one_over_dx, nor_x_l = (phi(i, j) - phi(i - 1, j))*signed_distance_field.one_over_dx, nor_y_u = (phi(i, j+1) - phi(i, j))*signed_distance_field.one_over_dy, nor_y_d = (phi(i, j) - phi(i, j-1))*signed_distance_field.one_over_dy;
			T nor_x = (T)0.5*(nor_x_r + nor_x_l), nor_y = (T)0.5*(nor_y_u + nor_y_d);
			T mag_of_normal = sqrt(POW2(nor_x) + POW2(nor_y));

			if (mag_of_normal != 0)
			{
				nor.x = nor_x;
				nor.y = nor_y;
			}
			else
			{
				nor.x = (phi(i+1, j) - phi(i, j))*signed_distance_field.one_over_dx;
				nor.y = (phi(i, j+1) - phi(i, j))*signed_distance_field.one_over_dy;
			}
			nor.MakeThisUnit();				
		}*/
	}

	// For 2D case
	void ComputeTangential(void)
	{
		ComputeNormals();
		
		int i(0), j(0);
		int i_start(tangential.i_start), i_end(tangential.i_end), j_start(tangential.j_start), j_end(tangential.j_end);

		LOOPS_2D(i, j, i_start, j_start, i_end, j_end)
		{
			tangential.array_for_this(i, j).x = normal.array_for_this(i, j).y;
			tangential.array_for_this(i, j).y = -normal.array_for_this(i, j).x;
		}
	}

	// Need to consider 1. Boundary Condition when curvature is defined as the 1st one. 
	void ComputeCurvatures(void)
	{
		if (is_axisymmetric)
		{
			FIELD_STRUCTURE_2D<T> phi_r, phi_x, deno;
			phi_r.Initialize(signed_distance_field.grid, 2);
			phi_x.Initialize(signed_distance_field.grid, 2);
			deno.Initialize(signed_distance_field.grid, 2);

			GRID_ITERATION_2D(signed_distance_field.grid)
			{
				phi_r(i, j) = (phi(i + 1, j) - phi(i - 1, j))*signed_distance_field.one_over_2dx;
				phi_x(i, j) = (phi(i, j + 1) - phi(i, j - 1))*signed_distance_field.one_over_2dy;
				deno(i, j) = sqrt(POW2(phi_r(i, j)) + POW2(phi_x(i, j)));
			}

			GRID_ITERATION_2D(signed_distance_field.grid)
			{
				T r_coor = signed_distance_field.grid.x_min + i*signed_distance_field.grid.dx, r_coor_r = signed_distance_field.x_min + (i + 1)*signed_distance_field.grid.dx, r_coor_l = signed_distance_field.x_min + (i - 1)*signed_distance_field.grid.dx;
				T one_over_r_coor = (T)1/r_coor;

				T curv = one_over_r_coor*signed_distance_field.one_over_2dx*(r_coor_r*phi_r(i + 1, j)/deno(i + 1, j) - r_coor_l*phi_r(i - 1, j)/deno(i - 1, j));
				curv += signed_distance_field.one_over_2dx*(phi_x(i, j + 1)/deno(i, j + 1) - phi_x(i, j - 1)/deno(i, j - 1));

				curvature(i, j) = curv;
			}
		}
		else
		{
			T tolerance = (T)1/min(grid.dx, grid.dy);
		
			// Curvature by the definition of divergence of normal
			if (curvature_by_normal_vector == true)
			{
				normal.FillGhostCellsFrom(normal.array_for_this, true);

				GRID_ITERATION_2D(grid)
				{
					T curv = (normal(i+1, j).x - normal(i-1, j).x)*signed_distance_field.one_over_2dx;
					curv += (normal(i, j+1).y - normal(i, j+1).y)*signed_distance_field.one_over_2dy;
				
					if (abs(curv) <= tolerance)
					{
						curvature(i, j) = -curv;
					}
					else
					{
						curvature(i, j) = -tolerance;
					}
				}
			}
			// Curvature calculated by levelset
			if (curvature_by_levelset == true)
			{
				FIELD_STRUCTURE_2D<T> phi_x, phi_y, phi_xx, phi_yy, phi_xy;
				phi_x.Initialize(signed_distance_field.grid, 2);
				phi_y.Initialize(signed_distance_field.grid, 2);
				phi_xx.Initialize(signed_distance_field.grid, 2);
				phi_yy.Initialize(signed_distance_field.grid, 2);
				phi_xy.Initialize(signed_distance_field.grid, 2);
	
				GRID_ITERATION_2D(signed_distance_field.grid)
				{
					phi_x(i, j) = (phi(i + 1, j) - phi(i - 1, j))*signed_distance_field.one_over_2dx;
					phi_y(i, j) = (phi(i, j + 1) - phi(i, j - 1))*signed_distance_field.one_over_2dy;
					phi_xx(i, j) = (phi(i + 1, j) - 2*phi(i, j) + phi(i - 1, j))*signed_distance_field.one_over_dx2;
					phi_yy(i, j) = (phi(i, j + 1) - 2*phi(i, j) + phi(i, j - 1))*signed_distance_field.one_over_dy2;
					phi_xy(i, j) = (phi(i + 1, j + 1) - phi(i + 1, j - 1) - phi(i - 1, j + 1) + phi(i - 1, j - 1))*signed_distance_field.one_over_2dx*signed_distance_field.one_over_2dy;		
				}
			
				T px, py, pxx, pyy, pxy;
			
				for (int j = grid.j_start; j <= grid.j_end; j++)
				{
					for (int i = grid.i_start; i <= grid.i_end; i++)
					{
						px = phi_x(i, j);
						py = phi_y(i, j);
						pxx = phi_xx(i, j);
						pyy = phi_yy(i, j);
						pxy = phi_xy(i, j);
						
						T magnitude = sqrt(POW2(px) + POW2(py));
						T deno = POW3(magnitude);
						
						T curv(0);
	
						if (deno != 0)
						{
							T one_over_deno = (T)1/deno;
							curv = -(POW2(px)*pyy - (T)2*px*py*pxy + POW2(py)*pxx)*one_over_deno;
						}
						else
						{
							px = (phi(i + 1, j) - phi(i, j))*grid.one_over_dx, py = (phi(i, j + 1) - phi(i, j))*grid.one_over_dy;
							magnitude = sqrt(POW2(px) + POW2(py));
							deno = POW3(magnitude);
							if (deno == 0)
							{
								cout << "Denominator cannot be zero!!" << endl;
								exit(0);
							}
							T one_over_deno = (T)1/deno;
							curv = -(POW2(px)*pyy - (T)2*px*py*pxy + POW2(py)*pxx)*one_over_deno;
						}
		
						if (abs(curv) <= tolerance)
						{
							curvature(i, j) = curv;
						}
						else if (curv > tolerance)
						{
							curvature(i, j) = tolerance;
						}
						else if (curv < -tolerance)
						{
							curvature(i, j) = -tolerance;
						}
	
						//curvature(i, j) = curv;
					}
				}
			}
		}
	}
	
	inline const VT Normal(const VT& position) const
	{
		return normal.BilinearInterpolation(position);
	}

	inline void UnitNormal(const VT& position, VT& normal_output) const
	{
		normal_output = normal.BilinearInterpolation(position);
		normal_output.MakeThisUnit();
	}

	inline T Curvature(const VT& position) const
	{
		return curvature.BilinearInterpolation(position); 
	}

public: // Error Estimate Functions (For the Test)
	T ErrorInf(void)
	{
		T err_inf(0);

		for (int j = grid.j_start; j <= grid.j_end; j++)
		{
			for (int i = grid.i_start; i <= grid.i_end; i++)
			{
				T temp = abs(phi(i,j) - phi_true(i,j));
				err_inf = MAX(err_inf, temp);
			}
		}

		return err_inf;
	}

	T Error2Norm(void)
	{
		T err_2norm(0);

		for (int j = grid.j_start; j <= grid.j_end; j++)
		{
			for (int i = grid.i_start; i <= grid.i_end; i++)
			{
				err_2norm += POW2(abs(phi(i,j) - phi_true(i,j)));
			}
		}

		err_2norm = grid.dx*sqrt(err_2norm);

		return err_2norm;
	}

	T Error1Norm(void)
	{
		T err_1norm(0);

		for (int j = grid.j_start; j <= grid.j_end; j++)
		{
			for (int i = grid.i_start; i <= grid.i_end; i++)
			{
				err_1norm += abs(phi(i,j) - phi_true(i,j));
			}
		}

		err_1norm = grid.dx*grid.dx*err_1norm;

		return err_1norm;
	}
};

template<class TT>
static void ComputeCurvature(const FIELD_STRUCTURE_2D<TT>& levelset_phi, FIELD_STRUCTURE_2D<TT>& curvature)
{
	FIELD_STRUCTURE_2D<T> phi_x, phi_y, phi_xx, phi_yy, phi_xy;
	phi_x.Initialize(levelset_phi.grid, 2);
	phi_y.Initialize(levelset_phi.grid, 2);
	phi_xx.Initialize(levelset_phi.grid, 2);
	phi_yy.Initialize(levelset_phi.grid, 2);
	phi_xy.Initialize(levelset_phi.grid, 2);

	GRID_ITERATION_2D(levelset_phi.grid)
	{
		phi_x(i, j) = (levelset_phi(i + 1, j) - levelset_phi(i - 1, j))*levelset_phi.one_over_2dx;
		phi_y(i, j) = (levelset_phi(i, j + 1) - levelset_phi(i, j - 1))*levelset_phi.one_over_2dy;
		phi_xx(i, j) = (levelset_phi(i + 1, j) - 2*levelset_phi(i, j) + levelset_phi(i - 1, j))*levelset_phi.one_over_dx2;
		phi_yy(i, j) = (levelset_phi(i, j + 1) - 2*levelset_phi(i, j) + levelset_phi(i, j - 1))*levelset_phi.one_over_dy2;
		phi_xy(i, j) = (levelset_phi(i + 1, j + 1) - levelset_phi(i + 1, j - 1) - levelset_phi(i - 1, j + 1) + levelset_phi(i - 1, j - 1))*levelset_phi.one_over_2dx*levelset_phi.one_over_2dy;		
	}
		
	T px, py, pxx, pyy, pxy;
	
	for (int j = levelset_phi.grid.j_start; j <= levelset_phi.grid.j_end; j++)
	{
		for (int i = levelset_phi.grid.i_start; i <= levelset_phi.grid.i_end; i++)
		{
			px = phi_x(i, j);
			py = phi_y(i, j);
			pxx = phi_xx(i, j);
			pyy = phi_yy(i, j);
			pxy = phi_xy(i, j);
				
			T magnitude = sqrt(POW2(px) + POW2(py));
			T deno = POW3(magnitude);
					
			T curv(0);

			if (deno != 0)
			{
				T one_over_deno = (T)1/deno;
				curv = -(POW2(px)*pyy - 2*px*py*pxy + POW2(py)*pxx)*one_over_deno;
			}
			else
			{
				px = (levelset_phi(i + 1, j) - levelset_phi(i, j))*levelset_phi.grid.one_over_dx, py = (levelset_phi(i, j + 1) - levelset_phi(i, j))*levelset_phi.grid.one_over_dy;
				magnitude = sqrt(POW2(px) + POW2(py));
				deno = POW3(magnitude);
				if (deno == 0)
				{
					cout << "Denominator cannot be zero!!" << endl;
					exit(0);
				}
				T one_over_deno = (T)1/deno;
				curv = -(POW2(px)*pyy - 2*px*py*pxy + POW2(py)*pxx)*one_over_deno;
			}
			
			/*T tolerance = min(levelset_phi.grid.one_over_dx, levelset_phi.grid.one_over_dy);

			if (abs(curv) <= tolerance)
			{
				curvature(i, j) = curv;
			}
			else if (curv > tolerance)
			{
				curvature(i, j) = tolerance;
			}
			else if (curv < -tolerance)
			{
				curvature(i, j) = -tolerance;
			}*/
			curvature(i, j) = curv;
		}
	}
}


