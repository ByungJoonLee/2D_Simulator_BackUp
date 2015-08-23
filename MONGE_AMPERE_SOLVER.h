#pragma once

#include "COMMON_DEFINITION.h"
#include "WORLD_DISCRETIZATION_2D.h"
#include "LEVELSET_2D.h"
#include "CSR_MATRIX.h"
#include "LINEAR_SOLVER.h"
#include "GAUSS_SEIDEL_METHOD.h"
#include "POISSON_SOLVER.h"
#include "BICGSTAB_METHOD.h"

class MONGE_AMPERE_SOLVER
{
public: // World Discretization
	WORLD_DISCRETIZATION_2D*			world_discretization;

public: // Main variables
	FIELD_STRUCTURE_2D<T>				u;
	FIELD_STRUCTURE_2D<VT>				grad_u;

	FIELD_STRUCTURE_2D<VT>				true_solution;

	FIELD_STRUCTURE_2D<T>				density_x;
	FIELD_STRUCTURE_2D<T>				density_y;
	FIELD_STRUCTURE_2D<T>				density_h_x;
	FIELD_STRUCTURE_2D<T>				density_h_y;

	ARRAY<FIELD_STRUCTURE_2D<T>>		MA_array;
	FIELD_STRUCTURE_2D<T>				MA;

	FIELD_STRUCTURE_2D<int>				bc;
	FIELD_STRUCTURE_2D<T>				RHS;
	FIELD_STRUCTURE_2D<T>				RHS_h;

public: // Convenient variables and references
	GRID_STRUCTURE_2D					base_grid;

public: // For multithreading
	MULTITHREADING*						multithreading;

public: // For Newton Method
	CSR_MATRIX<T>						jacobian;
	VECTOR_ND<T>						x;
	VECTOR_ND<T>						b;

	VECTOR_ND<T>						Jb;
	int									num_iteration, max_iteration_for_sub_linear_solver, max_iteration_for_Newton_Method;
	T									delta;
	
	T									alpha;

	enum POISSON_SOLVER_TYPE            sub_linear_solver_type;
	LINEAR_SOLVER*						sub_linear_solver;
	
	enum POISSON_SOLVER_TYPE			poisson_solver_type;
	POISSON_SOLVER						poisson_solver;

	T									tolerance;

public: // Test
	int									test_number;

public: // Constructor and Destructor
	MONGE_AMPERE_SOLVER(void)
		: multithreading(0), sub_linear_solver(0)
	{}

	~MONGE_AMPERE_SOLVER(void)
	{}

public: // Initialization Function
	void InitializeFromBlock(const SCRIPT_BLOCK& monge_ampere_eqn_block, MULTITHREADING* multithreading_input)
	{
		multithreading = multithreading_input;

		// Grid Initialization
		base_grid.Initialize(world_discretization->world_grid);

		// Initialize Fields
		u.Initialize(base_grid, 1, multithreading);
		grad_u.Initialize(base_grid, 1, multithreading);

		true_solution.Initialize(base_grid, 1, multithreading);

		density_x.Initialize(base_grid, 1, multithreading);
		density_y.Initialize(base_grid, 1, multithreading);
		density_h_x.Initialize(base_grid, 1, multithreading);
		density_h_y.Initialize(base_grid, 1, multithreading);
		
		bc.Initialize(base_grid, 1, multithreading);
		RHS.Initialize(base_grid, 1, multithreading);
		RHS_h.Initialize(base_grid, 1, multithreading);

		// We limited here as compact scheme
		MA_array.Initialize(2);
		MA.Initialize(base_grid, 1, multithreading);

		MA_array[0].Initialize(base_grid, 1, multithreading);
		MA_array[1].Initialize(base_grid, 1, multithreading);

		// For practical reason
		delta = base_grid.dx2;

		tolerance = monge_ampere_eqn_block.GetFloat("tolerance", (T)0.001);
		max_iteration_for_sub_linear_solver = monge_ampere_eqn_block.GetInteger("max_iteration_for_sub_linear_solver", (int)100);
		max_iteration_for_Newton_Method = monge_ampere_eqn_block.GetInteger("max_iteration_for_Newton_Method", (int)100);
		
		// For sub solver
		const char* sub_linear_solver_type_input = monge_ampere_eqn_block.GetString("sub_linear_solver_type", "Null");

		if (!strcmp(sub_linear_solver_type_input, "GS"))
		{
			sub_linear_solver_type = GS;
		}
		else
		{
			sub_linear_solver_type = GS;
		}
		if (!strcmp(sub_linear_solver_type_input, "BICG"))
		{
			sub_linear_solver_type = BICG;
		}
		else
		{
			sub_linear_solver_type = GS;
		}

		InitializeLinearSolver(sub_linear_solver_type);

		cout << "-------------------For Sublinear Solver-------------------" << endl;
		cout << "tolerance: " << tolerance << endl;
		cout << "max_iteration: " << max_iteration_for_sub_linear_solver << endl;
				
		switch (sub_linear_solver_type)
		{
		case NO_SOLVER:
			cout << "sublinear solver type: " << "NO SOLVER" << endl;
			break;

		case GS:
			cout << "sublinear solver type: " << "GS" << endl;
			break;

		case BICG:
			cout << "sublinear solver type: " << "BICG" << endl;
			break;

		default:
			break;
		}

		// For having initial value of Newton Method
		const char* poisson_solver_type_input = monge_ampere_eqn_block.GetString("poisson_solver_type", "Null");
		if (!strcmp(poisson_solver_type_input, "CG"))
		{
			poisson_solver_type = CG;
		}
		else
		{
			poisson_solver_type = CG;
		}
		if (!strcmp(poisson_solver_type_input, "PCG"))
		{
			poisson_solver_type = PCG;
		}
		else
		{
			poisson_solver_type = CG;
		}

		alpha = monge_ampere_eqn_block.GetFloat("alpha", (T)1);

		cout << "-------------------For Newton's Method-------------------" << endl;
		cout << "tolerance: " << tolerance << endl;
		cout << "max_iteration: " << max_iteration_for_Newton_Method << endl;
				
		switch (poisson_solver_type)
		{
		case NO_SOLVER:
			cout << "poisson solver type: " << "NO SOLVER" << endl;
			break;

		case CG:
			cout << "poisson solver type: " << "CG" << endl;
			break;

		case PCG:
			cout << "poisson solver type: " << "PCG" << endl;
			break;

		default:
			break;
		}

		// Initialize Poisson solvers
		poisson_solver.Initialize(tolerance, max_iteration_for_sub_linear_solver, 0, multithreading);
		poisson_solver.InitializeLinearSolver(poisson_solver_type);
	}

	void InitializeLinearSolver(const POISSON_SOLVER_TYPE linear_solver_type);

public: // Member Function
	void SolveThreaded(const int& thread_id)
	{
		SetupDensity(thread_id);
		SetupInitialForNewtonMethod(thread_id);
		SetupBoundaryCondition(bc, thread_id);
		NewtonMethod(jacobian, x, b, thread_id);
		ComputeGradient();
	}
	
	void AssignStencil(const FIELD_STRUCTURE_2D<T>& rho_x, const FIELD_STRUCTURE_2D<T>& rho_y_1, const FIELD_STRUCTURE_2D<T>& rho_y_2, const int& thread_id)
	{
		// Speedup Variable
		T one_over_dx2 = u.grid.one_over_dx2, one_over_dy2 = u.grid.one_over_dy2;
		T one_over_2dx = u.grid.one_over_2dx, one_over_2dy = u.grid.one_over_2dy;

		BEGIN_GRID_ITERATION_2D(u.partial_grids[thread_id])
		{
			T dx, dy, dxx, dyy, dv, dvp, dvv, dvpvp;
			
			// x, y derivatives
			dx = one_over_2dx*(u(i + 1, j) - u(i - 1, j));
			dy = one_over_2dy*(u(i, j + 1) - u(i, j - 1));
			dxx = one_over_dx2*(u(i + 1, j) - 2*u(i, j) + u(i - 1, j));
			dyy = one_over_dy2*(u(i, j + 1) - 2*u(i, j) + u(i, j - 1));

			// v, vp derivatives
			dv    = (T)1/sqrt(2)*one_over_2dx*(u(i + 1, j + 1) - u(i - 1, j - 1));
			dvp   = (T)1/sqrt(2)*one_over_2dy*(u(i + 1, j - 1) - u(i - 1, j + 1));
			dvv   = (T)1/2*one_over_dx2*(u(i + 1, j + 1) + u(i - 1, j - 1) - 2*u(i, j));
			dvpvp = (T)1/2*one_over_dx2*(u(i + 1, j - 1) + u(i - 1, j + 1) - 2*u(i, j));

			// MA_1 and MA_2 - Need to add the constant term u_0
			(MA_array.values[0])(i, j) = max(dxx, delta)*max(dyy, delta) - min(dxx, delta) - min(dyy, delta) - rho_x(i, j)/rho_y_1(i, j);
			(MA_array.values[1])(i, j) = max(dvv, delta)*max(dvpvp, delta) - min(dvv, delta) - min(dvpvp, delta) - rho_x(i, j)/rho_y_2(i, j);
			
			MA(i, j) = min((MA_array.values[0])(i, j), (MA_array.values[1])(i, j));
		}
		END_GRID_ITERATION_2D;
	}

	void CalculateJacobian(CSR_MATRIX<T>& J, VECTOR_ND<T>& b_vector, const int& thread_id)
	{
		AssignStencil(density_x, density_y, density_y, thread_id);
		
		const int num_all_full_cells = AssignSequentialindicesToFullCells(bc, thread_id);
		const int nnz = CountNonZeroElements(bc, thread_id);
		
		HEAD_THREAD_WORK(J.Initialize(num_all_full_cells, nnz, multithreading));
		HEAD_THREAD_WORK(b_vector.Initialize(num_all_full_cells));
		HEAD_THREAD_WORK(multithreading->SplitDomainIndex1D(0, num_all_full_cells));

		BEGIN_HEAD_THREAD_WORK
		{
			J.start_ix[0] = 0;
			J.end_ix[0] = multithreading->sync_value_int[0] - 1;
			J.prev_row_array[0] = -1;
			J.values_ix_array[0] = 0;
			for (int thread_id = 1; thread_id < multithreading->num_threads; thread_id++)
			{
				J.start_ix[thread_id] = J.end_ix[thread_id - 1] + 1;
				J.end_ix[thread_id] = J.end_ix[thread_id - 1] + multithreading->sync_value_int[thread_id];
				J.prev_row_array[thread_id] = -1;
				J.values_ix_array[thread_id] = J.start_ix[thread_id];
			}
		}
		END_HEAD_THREAD_WORK

		// Speedup Variable
		T one_over_dx  = u.grid.one_over_dx, one_over_dy = u.grid.one_over_dy;
		T one_over_dx2 = u.grid.one_over_dx2, one_over_dy2 = u.grid.one_over_dy2;
		T one_over_2dx = u.grid.one_over_2dx, one_over_2dy = u.grid.one_over_2dy;
		
		int i_start = u.grid.i_start, i_end = u.grid.i_end;
		int j_start = u.grid.j_start, j_end = u.grid.j_end;
		
		// For fixed point
		int mid_i(0), mid_j(0);
			
		mid_i = (T)0.5*(i_start + i_end);
		mid_j = (T)0.5*(j_start + j_end);

		BEGIN_GRID_ITERATION_2D(u.partial_grids[thread_id])
		{
			T dxx, dyy, dvv, dvpvp;

			dxx = one_over_dx2*(u(i + 1, j) - 2*u(i, j) + u(i - 1, j));
			dyy = one_over_dy2*(u(i, j + 1) - 2*u(i, j) + u(i, j - 1));
			dvv   = (T)1/2*one_over_dx2*(u(i + 1, j + 1) + u(i - 1, j - 1) - 2*u(i, j));
			dvpvp = (T)1/2*one_over_dx2*(u(i + 1, j - 1) + u(i - 1, j + 1) - 2*u(i, j));

			if (bc(i, j) < 0)
			{
				continue;
			}

			// Need to add the derivative of source term
			//if (i == i_start)
			//{
			//	if (j == j_start)
			//	{
			//		// Define b_vector
			//		T H_star_n;

			//		H_star_n = max(-u.grid.x_min, -u.grid.y_min);
			//		
			//		T dxp = one_over_dx*(u(i + 1, j) - u(i, j)), dyp = one_over_dy*(u(i, j + 1) - u(i, j));
			//		
			//		T b_vector_value = -(T)1/sqrt(2)*dxp - (T)1/sqrt(2)*dyp - H_star_n;

			//		b_vector[bc(i, j)] = b_vector_value;

			//		// Assign Jacobian
			//		T u_r, u_u, u_c;

			//		u_r = -(T)1/sqrt(2)*one_over_dx;
			//		u_u = -(T)1/sqrt(2)*one_over_dx;
			//		u_c = (T)2/sqrt(2)*one_over_dx;

			//		if (bc(i + 1, j) > -1)
			//		{
			//			J.AssignValue(bc(i, j), bc(i + 1, j), u_r, thread_id);
			//		}
			//		if (bc(i, j + 1) > -1)
			//		{
			//			J.AssignValue(bc(i, j), bc(i, j + 1), u_u, thread_id);
			//		}

			//		J.AssignValue(bc(i, j), bc(i, j), u_c, thread_id);
			//	}
			//	else if (j == j_end)
			//	{
			//		// Define b_vector
			//		T H_star_n;

			//		H_star_n = max(-u.grid.x_min, u.grid.y_max);
			//		
			//		T dxp = one_over_dx*(u(i + 1, j) - u(i, j)), dym = one_over_dy*(u(i, j) - u(i, j - 1));
			//		
			//		T b_vector_value = -(T)1/sqrt(2)*dxp + (T)1/sqrt(2)*dym - H_star_n;

			//		b_vector[bc(i, j)] = b_vector_value;

			//		// Assign Jacobian
			//		T u_r, u_d, u_c;

			//		u_r = -(T)1/sqrt(2)*one_over_dx;
			//		u_d = -(T)1/sqrt(2)*one_over_dx;
			//		u_c = (T)2/sqrt(2)*one_over_dx;

			//		if (bc(i + 1, j) > -1)
			//		{
			//			J.AssignValue(bc(i, j), bc(i + 1, j), u_r, thread_id);
			//		}
			//		if (bc(i, j - 1) > -1)
			//		{
			//			J.AssignValue(bc(i, j), bc(i, j - 1), u_d, thread_id);
			//		}

			//		J.AssignValue(bc(i, j), bc(i, j), u_c, thread_id);
			//	}
			//	else
			//	{
			//		// Define b_vector
			//		T H_star_n;

			//		H_star_n = -u.grid.x_min;
			//		
			//		T dxp = one_over_dx*(u(i + 1, j) - u(i, j));
			//		
			//		T b_vector_value = -(T)1/sqrt(2)*dxp - H_star_n;

			//		b_vector[bc(i, j)] = b_vector_value;

			//		// Assign Jacobian
			//		T u_r, u_c;

			//		u_r = -one_over_dx;
			//		u_c = one_over_dx;

			//		if (bc(i + 1, j) > -1)
			//		{
			//			J.AssignValue(bc(i, j), bc(i + 1, j), u_r, thread_id);
			//		}

			//		J.AssignValue(bc(i, j), bc(i, j), u_c, thread_id);
			//	}
			//}
			//else if (i == i_end)
			//{
			//	if (j == j_start)
			//	{
			//		// Define b_vector
			//		T H_star_n;

			//		H_star_n = max(u.grid.x_max, -u.grid.y_min);
			//		
			//		T dxm = one_over_dx*(u(i, j) - u(i - 1, j)), dyp = one_over_dy*(u(i, j + 1) - u(i, j));
			//		
			//		T b_vector_value = (T)1/sqrt(2)*dxm - (T)1/sqrt(2)*dyp - H_star_n;

			//		b_vector[bc(i, j)] = b_vector_value;

			//		// Assign Jacobian
			//		T u_l, u_u, u_c;

			//		u_l = -(T)1/sqrt(2)*one_over_dx;
			//		u_u = -(T)1/sqrt(2)*one_over_dx;
			//		u_c = (T)2/sqrt(2)*one_over_dx;

			//		if (bc(i - 1, j) > -1)
			//		{
			//			J.AssignValue(bc(i, j), bc(i - 1, j), u_l, thread_id);
			//		}
			//		if (bc(i, j + 1) > -1)
			//		{
			//			J.AssignValue(bc(i, j), bc(i, j + 1), u_u, thread_id);
			//		}

			//		J.AssignValue(bc(i, j), bc(i, j), u_c, thread_id);
			//	}
			//	else if (j == j_end)
			//	{
			//		// Define b_vector
			//		T H_star_n;

			//		H_star_n = max(u.grid.x_max, u.grid.y_max);
			//		
			//		T dxm = one_over_dx*(u(i, j) - u(i - 1, j)), dym = one_over_dy*(u(i, j) - u(i, j - 1));
			//		
			//		T b_vector_value = (T)1/sqrt(2)*dxm + (T)1/sqrt(2)*dym - H_star_n;

			//		b_vector[bc(i, j)] = b_vector_value;

			//		// Assign Jacobian
			//		T u_l, u_d, u_c;

			//		u_l = -(T)1/sqrt(2)*one_over_dx;
			//		u_d = -(T)1/sqrt(2)*one_over_dx;
			//		u_c = (T)2/sqrt(2)*one_over_dx;

			//		if (bc(i - 1, j) > -1)
			//		{
			//			J.AssignValue(bc(i, j), bc(i - 1, j), u_l, thread_id);
			//		}
			//		if (bc(i, j - 1) > -1)
			//		{
			//			J.AssignValue(bc(i, j), bc(i, j - 1), u_d, thread_id);
			//		}

			//		J.AssignValue(bc(i, j), bc(i, j), u_c, thread_id);
			//	}
			//	else
			//	{
			//		// Define b_vector
			//		T H_star_n;

			//		H_star_n = u.grid.x_max;
			//		
			//		T dxm = one_over_dx*(u(i, j) - u(i - 1, j));
			//		
			//		T b_vector_value = (T)1/sqrt(2)*dxm - H_star_n;

			//		b_vector[bc(i, j)] = b_vector_value;

			//		// Assign Jacobian
			//		T u_l, u_c;

			//		u_l = -one_over_dx;
			//		u_c = one_over_dx;

			//		if (bc(i - 1, j) > -1)
			//		{
			//			J.AssignValue(bc(i, j), bc(i - 1, j), u_l, thread_id);
			//		}

			//		J.AssignValue(bc(i, j), bc(i, j), u_c, thread_id);
			//	}
			//}
			//else if (j == j_start && i != i_start && i != i_end)
			//{
			//	// Define b_vector
			//	T H_star_n;

			//	H_star_n = -u.grid.y_min;
			//		
			//	T dyp = one_over_dx*(u(i, j + 1) - u(i, j));
			//		
			//	T b_vector_value = -(T)1/sqrt(2)*dyp - H_star_n;

			//	b_vector[bc(i, j)] = b_vector_value;

			//	// Assign Jacobian
			//	T u_u, u_c;

			//	u_u = -one_over_dy;
			//	u_c = one_over_dy;

			//	if (bc(i, j + 1) > -1)
			//	{
			//		J.AssignValue(bc(i, j), bc(i, j + 1), u_u, thread_id);
			//	}

			//	J.AssignValue(bc(i, j), bc(i, j), u_c, thread_id);
			//}
			//else if (j == j_end && i != i_start && i != i_end)
			//{
			//	// Define b_vector
			//	T H_star_n;

			//	H_star_n = u.grid.y_max;
			//		
			//	T dym = one_over_dx*(u(i, j) - u(i, j - 1));
			//		
			//	T b_vector_value = (T)1/sqrt(2)*dym - H_star_n;

			//	b_vector[bc(i, j)] = b_vector_value;

			//	// Assign Jacobian
			//	T u_d, u_c;

			//	u_d = -one_over_dy;
			//	u_c = one_over_dy;

			//	if (bc(i, j - 1) > -1)
			//	{
			//		J.AssignValue(bc(i, j), bc(i, j - 1), u_d, thread_id);
			//	}

			//	J.AssignValue(bc(i, j), bc(i, j), u_c, thread_id);
			//}
			/*else
			{*/
				if ((MA_array.values[1])(i, j) >= (MA_array.values[0])(i, j))
				{
					b_vector[bc(i, j)] = (MA_array.values[0])(i, j);
					
					T u_l, u_r, u_c, u_d, u_u;

					if (dxx >= delta)
					{
						if (dyy >= delta)
						{
							u_l = one_over_dx2*dyy;
							u_r = one_over_dx2*dyy;
							if (i == mid_i && j == mid_j)
							{
								u_c = -2*one_over_dx2*(dxx + dyy) - 1;
							}
							else
							{
								u_c = -2*one_over_dx2*(dxx + dyy);
							}
							u_d = one_over_dy2*dxx;
							u_u = one_over_dy2*dxx;

							if (bc(i - 1, j) > -1)
							{
								J.AssignValue(bc(i, j), bc(i - 1, j), u_l, thread_id);
							}
							if (bc(i + 1, j) > -1)
							{
								J.AssignValue(bc(i, j), bc(i + 1, j), u_r, thread_id);
							}
							if (bc(i, j - 1) > -1)
							{
								J.AssignValue(bc(i, j), bc(i, j - 1), u_d, thread_id);
							}
							if (bc(i, j + 1) > -1)
							{
								J.AssignValue(bc(i, j), bc(i, j + 1), u_u, thread_id);
							}

							J.AssignValue(bc(i, j), bc(i, j), u_c, thread_id); 
						}
						else
						{
							u_l = one_over_dx2*delta;
							u_r = one_over_dx2*delta;
							if (i == mid_i && j == mid_j)
							{
								u_c = -2*one_over_dx2*(delta - 1) - 1;
							}
							else
							{
								u_c = -2*one_over_dx2*(delta - 1);
							}
							u_d = -one_over_dy2;
							u_u = -one_over_dy2;

							if (bc(i - 1, j) > -1)
							{
								J.AssignValue(bc(i, j), bc(i - 1, j), u_l, thread_id);
							}
							if (bc(i + 1, j) > -1)
							{
								J.AssignValue(bc(i, j), bc(i + 1, j), u_r, thread_id);
							}
							if (bc(i, j - 1) > -1)
							{
								J.AssignValue(bc(i, j), bc(i, j - 1), u_d, thread_id);
							}
							if (bc(i, j + 1) > -1)
							{
								J.AssignValue(bc(i, j), bc(i, j + 1), u_u, thread_id);
							}

							J.AssignValue(bc(i, j), bc(i, j), u_c, thread_id);
						}
					}
					else
					{
						T u_l, u_r, u_c, u_d, u_u;

						if (dyy >= delta)
						{
							u_l = -one_over_dx2;
							u_r = -one_over_dx2;
							if (i == mid_i && j == mid_j)
							{
								u_c = -2*one_over_dx2*(delta - 1) - 1;
							}
							else
							{
								u_c = -2*one_over_dx2*(delta - 1);
							}
							u_d = one_over_dy2*delta;
							u_u = one_over_dy2*delta;

							if (bc(i - 1, j) > -1)
							{
								J.AssignValue(bc(i, j), bc(i - 1, j), u_l, thread_id);
							}
							if (bc(i + 1, j) > -1)
							{
								J.AssignValue(bc(i, j), bc(i + 1, j), u_r, thread_id);
							}
							if (bc(i, j - 1) > -1)
							{
								J.AssignValue(bc(i, j), bc(i, j - 1), u_d, thread_id);
							}
							if (bc(i, j + 1) > -1)
							{
								J.AssignValue(bc(i, j), bc(i, j + 1), u_u, thread_id);
							}

							J.AssignValue(bc(i, j), bc(i, j), u_c, thread_id);
						}
						else
						{
							u_l = -one_over_dx2;
							u_r = -one_over_dx2;
							if (i == mid_i && j == mid_j)
							{
								u_c = 4*one_over_dx2 - 1;
							}
							else
							{
								u_c = 4*one_over_dx2;
							}
							u_d = -one_over_dy2;
							u_u = -one_over_dy2;

							if (bc(i - 1, j) > -1)
							{
								J.AssignValue(bc(i, j), bc(i - 1, j), u_l, thread_id);
							}
							if (bc(i + 1, j) > -1)
							{
								J.AssignValue(bc(i, j), bc(i + 1, j), u_r, thread_id);
							}
							if (bc(i, j - 1) > -1)
							{
								J.AssignValue(bc(i, j), bc(i, j - 1), u_d, thread_id);
							}
							if (bc(i, j + 1) > -1)
							{
								J.AssignValue(bc(i, j), bc(i, j + 1), u_u, thread_id);
							}

							J.AssignValue(bc(i, j), bc(i, j), u_c, thread_id);
						}
					}
				}
				else
				{
					b_vector[bc(i, j)] = (MA_array.values[1])(i, j);

					T u_dl, u_ur, u_c, u_dr, u_ul;

					if (dvv >= delta)
					{
						if (dvpvp >= delta)
						{
							u_dl = (T)1/2*one_over_dx2*dvpvp;
							u_ur = (T)1/2*one_over_dx2*dvpvp;
							if (i == mid_i && j == mid_j)
							{
								u_c  = -one_over_dx2*(dvv + dvpvp) - 1;
							}
							else
							{
								u_c  = -one_over_dx2*(dvv + dvpvp);
							}
							u_dr = (T)1/2*one_over_dx2*dvv;
							u_ul = (T)1/2*one_over_dx2*dvv;

							if (bc(i - 1, j - 1) > -1)
							{
								J.AssignValue(bc(i, j), bc(i - 1, j - 1), u_dl, thread_id);
							}
							if (bc(i + 1, j + 1) > -1)
							{
								J.AssignValue(bc(i, j), bc(i + 1, j + 1), u_ur, thread_id);
							}
							if (bc(i + 1, j - 1) > -1)
							{
								J.AssignValue(bc(i, j), bc(i + 1, j - 1), u_dr, thread_id);
							}
							if (bc(i - 1, j + 1) > -1)
							{
								J.AssignValue(bc(i, j), bc(i - 1, j + 1), u_ul, thread_id);
							}

							J.AssignValue(bc(i, j), bc(i, j), u_c, thread_id);
						}
						else
						{
							u_dl = (T)1/2*one_over_dx2*delta;
							u_ur = (T)1/2*one_over_dx2*delta;
							if (i == mid_i && j == mid_j)
							{
								u_c  = -one_over_dx2*(delta - 1) - 1;
							}
							else
							{
								u_c  = -one_over_dx2*(delta - 1);
							}
							u_dr = -(T)1/2*one_over_dx2;
							u_ul = -(T)1/2*one_over_dx2;

							if (bc(i - 1, j - 1) > -1)
							{
								J.AssignValue(bc(i, j), bc(i - 1, j - 1), u_dl, thread_id);
							}
							if (bc(i + 1, j + 1) > -1)
							{
								J.AssignValue(bc(i, j), bc(i + 1, j + 1), u_ur, thread_id);
							}
							if (bc(i + 1, j - 1) > -1)
							{
								J.AssignValue(bc(i, j), bc(i + 1, j - 1), u_dr, thread_id);
							}
							if (bc(i - 1, j + 1) > -1)
							{
								J.AssignValue(bc(i, j), bc(i - 1, j + 1), u_ul, thread_id);
							}

							J.AssignValue(bc(i, j), bc(i, j), u_c, thread_id);
						}
					}
					else
					{
						T u_dl, u_ur, u_c, u_dr, u_ul;

						if (dvpvp >= delta)
						{
							u_dl = -(T)1/2*one_over_dx2;
							u_ur = -(T)1/2*one_over_dx2;
							if (i == mid_i && j == mid_j)
							{
								u_c  = -one_over_dx2*(delta - 1) - 1;
							}
							else
							{
								u_c  = -one_over_dx2*(delta - 1);
							}
							u_dr = (T)1/2*one_over_dx2*delta;
							u_ul = (T)1/2*one_over_dx2*delta;

							if (bc(i - 1, j - 1) > -1)
							{
								J.AssignValue(bc(i, j), bc(i - 1, j - 1), u_dl, thread_id);
							}
							if (bc(i + 1, j + 1) > -1)
							{
								J.AssignValue(bc(i, j), bc(i + 1, j + 1), u_ur, thread_id);
							}
							if (bc(i + 1, j - 1) > -1)
							{
								J.AssignValue(bc(i, j), bc(i + 1, j - 1), u_dr, thread_id);
							}
							if (bc(i - 1, j + 1) > -1)
							{
								J.AssignValue(bc(i, j), bc(i - 1, j + 1), u_ul, thread_id);
							}

							J.AssignValue(bc(i, j), bc(i, j), u_c, thread_id);
						}
						else
						{
							u_dl = -(T)1/2*one_over_dx2;
							u_ur = -(T)1/2*one_over_dx2;
							if (i == mid_i && j == mid_j)
							{
								u_c  = (T)2*one_over_dx2 - 1;
							}
							else
							{
								u_c  = (T)2*one_over_dx2;
							}
							u_dr = -(T)1/2*one_over_dx2;
							u_ul = -(T)1/2*one_over_dx2;

							if (bc(i - 1, j - 1) > -1)
							{
								J.AssignValue(bc(i, j), bc(i - 1, j - 1), u_dl, thread_id);
							}
							if (bc(i + 1, j + 1) > -1)
							{
								J.AssignValue(bc(i, j), bc(i + 1, j + 1), u_ur, thread_id);
							}
							if (bc(i + 1, j - 1) > -1)
							{
								J.AssignValue(bc(i, j), bc(i + 1, j - 1), u_dr, thread_id);
							}
							if (bc(i - 1, j + 1) > -1)
							{
								J.AssignValue(bc(i, j), bc(i - 1, j + 1), u_ul, thread_id);
							}

							J.AssignValue(bc(i, j), bc(i, j), u_c, thread_id);
						}
					}
				} 
			//}
		}
		END_GRID_ITERATION_2D;
	}

	int CountNonZeroElements(const FIELD_STRUCTURE_2D<int>& bc, const int& thread_id)
	{
		int nnz(0);
				
		BEGIN_GRID_ITERATION_2D(bc.partial_grids[thread_id])
		{
			if (bc(i, j) > -1)
			{
				nnz++;
				if (bc(i + 1, j) > -1)
				{
					nnz++;
				}
				if (bc(i - 1, j) > -1)
				{
					nnz++;
				}
				if (bc(i, j + 1) > -1)
				{
					nnz++;
				}
				if (bc(i, j - 1) > -1)
				{
					nnz++;
				}
			}
		}
		END_GRID_ITERATION_SUM(nnz);

		return nnz;
	}

	void NewtonMethod(CSR_MATRIX<T>& J, VECTOR_ND<T>& x, VECTOR_ND<T>& b_vector, const int& thread_id)
	{
		// For initializing Jacobian
		const int num_all_full_cells = AssignSequentialindicesToFullCells(bc, thread_id);
		const int nnz = CountNonZeroElements(bc, thread_id);
		
		HEAD_THREAD_WORK(J.Initialize(num_all_full_cells, nnz, multithreading));

		BEGIN_HEAD_THREAD_WORK
		{
			multithreading->SplitDomainIndex1D(0, J.N);
		}
		END_HEAD_THREAD_WORK;

		BEGIN_HEAD_THREAD_WORK
		{
			x.Initialize(J.N, true);	
		}
		END_HEAD_THREAD_WORK;
		
		const int N(x.num_dimension), start_ix(multithreading->start_ix_1D[thread_id]), end_ix(multithreading->end_ix_1D[thread_id]);

		BEGIN_HEAD_THREAD_WORK
		{
			num_iteration = 0;
		}
		END_HEAD_THREAD_WORK;
		
		BEGIN_HEAD_THREAD_WORK
		{
			Jb.Initialize(N, true);
		}
		END_HEAD_THREAD_WORK;
			
		T* xval(x.values);

		// For Gauss-Seidel Test
		/*CSR_MATRIX<T> A_Matrix;
		A_Matrix.Initialize(3, 9);

		A_Matrix.AssignValue(0, 0, 12);
		A_Matrix.AssignValue(0, 1, 3);
		A_Matrix.AssignValue(0, 2, -5);
		A_Matrix.AssignValue(1, 0, 1);
		A_Matrix.AssignValue(1, 1, 5);
		A_Matrix.AssignValue(1, 2, 3);
		A_Matrix.AssignValue(2, 0, 3);
		A_Matrix.AssignValue(2, 1, 7);
		A_Matrix.AssignValue(2, 2, 13);

		VECTOR_ND<T> b_vector_test;
		b_vector_test.Initialize(3, true);
		
		b_vector_test[0] = 1;
		b_vector_test[1] = 28;
		b_vector_test[2] = 76;
		
		VECTOR_ND<T> x_vector_sol;
		x_vector_sol.Initialize(3, true);

		x_vector_sol[0] = 1;
		x_vector_sol[1] = 0;
		x_vector_sol[2] = 1;

		sub_linear_solver->Solve(A_Matrix, x_vector_sol, b_vector_test, bc, thread_id);

		cout << x_vector_sol[0] << endl;
		cout << x_vector_sol[1] << endl;
		cout << x_vector_sol[2] << endl;*/

		T stopping_criterion;

		while (num_iteration < max_iteration_for_Newton_Method)
		{
			BEGIN_HEAD_THREAD_WORK
			{
				stopping_criterion = 0;
			}
			END_HEAD_THREAD_WORK;

			CalculateJacobian(J, b_vector, thread_id);
			
			sub_linear_solver->Solve(J, Jb, b_vector, bc, thread_id);

			GridToVector(u, x, bc, thread_id);

			for (int i = start_ix; i <= end_ix; i++)
			{
				xval[i] = xval[i] - alpha*Jb[i];
			}
			multithreading->Sync(thread_id);
			
			for (int i = start_ix; i <= end_ix; i++)
			{
				stopping_criterion += POW2(alpha*Jb[i]);
			}
			multithreading->SyncSum(thread_id, stopping_criterion);

			stopping_criterion = sqrt(stopping_criterion);

			VectorToGrid(x, u, bc, thread_id);

			/*ofstream fout;
			fout.open("solution");
			for (int j = u.j_start; j <= u.j_end; j++)
			{
				for (int i = u.i_start; i <= u.i_end; i++)
					fout << u(i, j) << " ";

					fout << "\n";
			}
			fout.close();*/

			if (stopping_criterion < tolerance)
			{
				cout << "--------------Newton Method--------------" << endl;
				cout << "Converge!!" << endl;
				cout << "Iteration Number : " << num_iteration << endl;
				cout << "Increment: " << stopping_criterion << endl;
				cout << "-----------------------------------------" << endl;
				
				ofstream fout;
				fout.open("solution");
				for (int j = u.j_start; j <= u.j_end; j++)
				{
					for (int i = u.i_start; i <= u.i_end; i++)
						fout << u(i, j) << " ";

						fout << "\n";
				}
				fout.close();
				break;
			}
			if (num_iteration == max_iteration_for_Newton_Method - 1)
			{
				cout << "--------------Newton Method--------------" << endl;
				cout << "Not Converge!!" << endl;
				cout << "Iteration Number : " << num_iteration << endl;
				cout << "Increment: " << stopping_criterion << endl;
				cout << "-----------------------------------------" << endl;
				ofstream fout;
				fout.open("solution");
				for (int j = u.j_start; j <= u.j_end; j++)
				{
					for (int i = u.i_start; i <= u.i_end; i++)
						fout << u(i, j) << " ";

						fout << "\n";
				}
				fout.close();
			}

			BEGIN_HEAD_THREAD_WORK
			{
				num_iteration++;
				cout << "Increment: " << stopping_criterion << endl;
			}
			END_HEAD_THREAD_WORK;
		}
	}

	int AssignSequentialindicesToFullCells(const FIELD_STRUCTURE_2D<int>& bc, const int& thread_id)
	{
		// Count number of full cells
		int full_ix(0);

		BEGIN_GRID_ITERATION_2D(bc.partial_grids[thread_id])
		{
			if (bc(i, j) > -1)
			{
				++full_ix;
			}	
		}
		END_GRID_ITERATION_2D;

		// Indexing the value of boundary condition field
		int start_full_ix, end_full_ix;
		multithreading->SyncDomainIndices1D(thread_id, full_ix, start_full_ix, end_full_ix);

		BEGIN_GRID_ITERATION_2D(bc.partial_grids[thread_id])
		{
			if (bc(i, j) > -1)
			{
				bc(i, j) = start_full_ix++;
			}		
		}
		END_GRID_ITERATION_2D;

		assert(start_full_ix - 1 == end_full_ix);

		multithreading->SyncSum(thread_id, full_ix);

		return full_ix;
	}

	void SetupBoundaryCondition(FIELD_STRUCTURE_2D<int>& bc_input, const int& thread_id)
	{
		ARRAY_2D<int>& bc_array(bc_input.array_for_this);
		GRID_STRUCTURE_2D& grid(bc_input.grid);

		int i(0), j(0);
		
		BEGIN_GRID_ITERATION_2D(bc_input.partial_grids_ghost[thread_id])
		{
			// Speed-up variable
			if (i < grid.i_start || i > grid.i_end || j < grid.j_start || j > grid.j_end)
			{
				bc_array(i, j) = BC_IMPLICIT;
			}
			else
			{
				bc_array(i, j) = BC_FULL;
			}
		}
		END_GRID_ITERATION_2D;
	}

	void SetupDensity(const int& thread_id)
	{
		if (test_number == 1)
		{
			ARRAY<T> sub(density_x.grid.i_res);
			ARRAY<T> sub_x(density_x.grid.i_res);
			ARRAY<T> sub_xx(density_x.grid.i_res);
			
			// For the integration
			ARRAY<T> sub_h(density_x.grid.i_res);
			ARRAY<T> sub_h_x(density_x.grid.i_res);
			ARRAY<T> sub_h_xx(density_x.grid.i_res);

			BEGIN_HEAD_THREAD_WORK
			{
				multithreading->SplitDomainIndex1D(0, sub.length);
			}
			END_HEAD_THREAD_WORK

			const int N(sub.length), start_ix(multithreading->start_ix_1D[thread_id]), end_ix(multithreading->end_ix_1D[thread_id]);

			for (int i = start_ix; i <= end_ix; i++)
			{
				T coor = density_x.x_min + i*density_x.dx;
				sub[i] = (-(T)1/(8*PI)*POW2(coor) + (T)1/(256*PI*PI*PI) + (T)1/(32*PI))*cos(8*PI*coor) + (T)1/(32*PI*PI)*coor*sin(8*PI*coor);
				sub_x[i] = -8*PI*(-(T)1/(8*PI)*POW2(coor) + (T)1/(256*PI*PI*PI) + (T)1/(32*PI))*sin(8*PI*coor) + (T)1/(32*PI*PI)*sin(8*PI*coor) + (T)8*PI/(32*PI*PI)*coor*cos(8*PI*coor);
				sub_xx[i] = 64*PI*PI*(-(T)1/(8*PI)*POW2(coor) + (T)1/(256*PI*PI*PI) + (T)1/(32*PI))*cos(8*PI*coor) + (T)8*PI/(32*PI*PI)*cos(8*PI*coor) + (T)8*PI/(32*PI*PI)*cos(8*PI*coor) - (T)64*PI*PI/(32*PI*PI)*coor*sin(8*PI*coor);
			}
			multithreading->Sync(thread_id);

			for (int i = start_ix; i <= end_ix; i++)
			{
				T coor = density_x.x_min + i*density_x.dx*(T)0.5;
				sub_h[i] = (-(T)1/(8*PI)*POW2(coor) + (T)1/(256*PI*PI*PI) + (T)1/(32*PI))*cos(8*PI*coor) + (T)1/(32*PI*PI)*coor*sin(8*PI*coor);
				sub_h_x[i] = -8*PI*(-(T)1/(8*PI)*POW2(coor) + (T)1/(256*PI*PI*PI) + (T)1/(32*PI))*sin(8*PI*coor) + (T)1/(32*PI*PI)*sin(8*PI*coor) + (T)8*PI/(32*PI*PI)*coor*cos(8*PI*coor);
				sub_h_xx[i] = 64*PI*PI*(-(T)1/(8*PI)*POW2(coor) + (T)1/(256*PI*PI*PI) + (T)1/(32*PI))*cos(8*PI*coor) + (T)8*PI/(32*PI*PI)*cos(8*PI*coor) + (T)8*PI/(32*PI*PI)*cos(8*PI*coor) - (T)64*PI*PI/(32*PI*PI)*coor*sin(8*PI*coor);
			}
			multithreading->Sync(thread_id);

			BEGIN_GRID_ITERATION_2D(density_x.partial_grids[thread_id])
			{
				density_x(i, j) = 1 + 4*(sub_xx[i]*sub[j] + sub[i]*sub_xx[j]) + 16*(sub[i]*sub[j]*sub_xx[i]*sub_xx[j] - POW2(sub_x[i])*POW2(sub_x[j]));
			}
			END_GRID_ITERATION_2D;

			BEGIN_GRID_ITERATION_2D(density_y.partial_grids[thread_id])
			{
				density_y(i, j) = (T)1;
			}
			END_GRID_ITERATION_2D;
			
			BEGIN_GRID_ITERATION_2D(density_x.partial_grids[thread_id])
			{
				density_h_x(i, j) = 1 + 4*(sub_h_xx[i]*sub_h[j] + sub_h[i]*sub_h_xx[j]) + 16*(sub_h[i]*sub_h[j]*sub_h_xx[i]*sub_h_xx[j] - POW2(sub_h_x[i])*POW2(sub_h_x[j]));
			}
			END_GRID_ITERATION_2D;

			BEGIN_GRID_ITERATION_2D(density_y.partial_grids[thread_id])
			{
				density_h_y(i, j) = (T)1;
			}
			END_GRID_ITERATION_2D;

			// True solution
			BEGIN_GRID_ITERATION_2D(true_solution.partial_grids[thread_id])
			{
				true_solution(i, j).x = (true_solution.grid.x_min + i*true_solution.grid.dx) + 4*sub_x[i]*sub[j];
				true_solution(i, j).y = (true_solution.grid.y_min + j*true_solution.grid.dy) + 4*sub[i]*sub_x[j];
			}
			END_GRID_ITERATION_2D;
		}
	}

	void SetupInitialForNewtonMethod(const int& thread_id)
	{
		if (test_number == 1)
		{
			BEGIN_GRID_ITERATION_2D(u.partial_grids[thread_id])
			{
				T x_coor = u.grid.x_min + i*u.grid.dx, y_coor = u.grid.y_min + j*u.grid.dy;

				u(i, j) = (T)0.5*POW2(x_coor) + (T)0.5*POW2(y_coor);
			}
			END_GRID_ITERATION_2D;
			
			//SetupBoundaryCondition(bc, thread_id);

			//BEGIN_GRID_ITERATION_2D(RHS.partial_grids[thread_id])
			//{
			//	RHS(i, j) = sqrt(density_x(i, j)/density_y(i, j));
			//}
			//END_GRID_ITERATION_2D;

			//// For integration
			//BEGIN_GRID_ITERATION_2D(RHS_h.partial_grids[thread_id])
			//{
			//	RHS_h(i, j) = sqrt(density_h_x(i, j)/density_h_y(i, j));
			//}
			//END_GRID_ITERATION_2D;

			//T integral_value(0);

			//BEGIN_GRID_ITERATION_2D(RHS_h.partial_grids[thread_id])
			//{
			//	integral_value += RHS_h(i, j)*RHS_h.grid.dx*RHS_h.grid.dy;
			//}
			//END_GRID_ITERATION_SUM(integral_value);	
			//
			//T rhs_value(0);

			//BEGIN_GRID_ITERATION_2D(RHS_h.partial_grids[thread_id])
			//{
			//	rhs_value += 2*RHS_h.grid.dx*RHS_h.grid.dy;
			//}
			//END_GRID_ITERATION_SUM(rhs_value);

			//T k = rhs_value/integral_value;

			//BEGIN_GRID_ITERATION_2D(RHS.partial_grids[thread_id])
			//{
			//	RHS(i, j) = k*RHS(i, j);
			//}
			//END_GRID_ITERATION_2D;
			//
			//poisson_solver.Solve(u, bc, RHS, thread_id);
			
				ComputeGradient();
		}
	}

	void VectorToGrid(const VECTOR_ND<T>& x, FIELD_STRUCTURE_2D<T>& solution, const FIELD_STRUCTURE_2D<int>& bc, const int& thread_id)
	{
		BEGIN_GRID_ITERATION_2D(solution.partial_grids[thread_id])
		{
			solution.array_for_this(i, j) = x[bc(i, j)];
		}
		END_GRID_ITERATION_2D;
	}

	void GridToVector(const FIELD_STRUCTURE_2D<T>& solution, VECTOR_ND<T>& x, const FIELD_STRUCTURE_2D<int>& bc, const int& thread_id)
	{
		BEGIN_GRID_ITERATION_2D(solution.partial_grids[thread_id])
		{
			x[bc(i, j)] = solution.array_for_this(i, j);
		}
		END_GRID_ITERATION_2D;
	}

	void ComputeGradient(void)
	{
		for (int i = u.grid.i_start; i <= u.grid.i_end; i++)
		{
			for (int j = u.grid.j_start; j <= u.grid.j_end; j++)
			{
				VT& grad(grad_u(i, j));

				T grad_x = (u(i + 1, j) - u(i - 1, j))*u.grid.one_over_2dx, grad_y = (u(i, j + 1) - u(i, j - 1))*u.grid.one_over_2dy;
			
				grad.x = grad_x;
				grad.y = grad_y;
			}
		}
	}

	void ComputeGradient(const FIELD_STRUCTURE_2D<T>& solution, const int& thread_id)
	{
		BEGIN_GRID_ITERATION_2D(solution.partial_grids[thread_id]);
		{
			cout << i << endl;
			cout << j << endl;

			VT& grad(grad_u(i, j));

			T grad_x = (solution(i + 1, j) - solution(i - 1, j))*u.grid.one_over_2dx, grad_y = (solution(i, j + 1) - solution(i, j - 1))*u.grid.one_over_2dy;
			
			grad.x = grad_x;
			grad.y = grad_y;
		}	
		END_GRID_ITERATION_2D;
	}
};