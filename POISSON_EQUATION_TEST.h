#pragma once

#include "COMMON_DEFINITION.h"
#include "WORLD_DISCRETIZATION_2D.h"
#include "LEVELSET_1D.h"
#include "LEVELSET_2D.h"
#include "POISSON_SOLVER.h"
#include "ADVECTION_METHOD_2D.h"

class POISSON_EQUATION_TEST
{
public: // World Discretization
	WORLD_DISCRETIZATION_2D*	world_discretization;
	
public: // For 1D problem
	// Interface Levelset
	LEVELSET_1D*				interface_levelset_1d;

	// Variables for Poisson Solver
	FIELD_STRUCTURE_1D<T>		solution_1d;
	FIELD_STRUCTURE_1D<T>		beta_1d;
	FIELD_STRUCTURE_1D<T>		rhs_1d;
	FIELD_STRUCTURE_1D<int>		boundary_condition_1d;

	// Jump Condition Field
	FIELD_STRUCTURE_1D<T>		jc_on_solution_1d;
	FIELD_STRUCTURE_1D<T>		jc_on_derivative_1d;

public: // For 2D problem
	// Interface Levelset
	LEVELSET_2D*				interface_levelset_2d;

	// Variables for Poisson Solver
	FIELD_STRUCTURE_2D<T>		solution_2d;
	FIELD_STRUCTURE_2D<T>		beta_2d;
	FIELD_STRUCTURE_2D<T>		rhs_2d;
	FIELD_STRUCTURE_2D<int>		boundary_condition_2d;

	// Jump Condition Field
	FIELD_STRUCTURE_2D<T>		jc_on_solution_2d;
	FIELD_STRUCTURE_2D<T>		jc_on_derivative_2d;

public: // Poisson Solver
	POISSON_SOLVER				poisson_solver;

public: // Control Option
	enum POISSON_SOLVER_TYPE    poisson_solver_type;

public: // Option for Boundary Condition
	bool						Dirichlet_Boundary_Condition, Neumann_Boundary_Condition;

public: // For CG/PCG
	T							tolerance;
	int							max_iteration;

public: // Convenient variables and references - For 1D
	GRID_STRUCTURE_1D			base_grid_1d;
	
public: // Convenient variables and references - For 2D
	GRID_STRUCTURE_2D			base_grid_2d;

public: // Multithreading
	MULTITHREADING*				multithreading;

public: // Options for dimension
	bool						grid_1d, grid_2d;
	
public: // Example number
	int							test_number;

public: // Constructor and Destructor
	POISSON_EQUATION_TEST(void)
		: interface_levelset_1d(0), interface_levelset_2d(0), multithreading(0), grid_1d(false), grid_2d(false)
	{}

	~POISSON_EQUATION_TEST(void)
	{}

public: // Initialization Function
	void InitializeFromBlock(const SCRIPT_BLOCK& poisson_eqn_block, MULTITHREADING* multithreading_input)
	{
		multithreading = multithreading_input;

		tolerance = poisson_eqn_block.GetFloat("tolerance", (T)1e-4);
		max_iteration = poisson_eqn_block.GetInteger("max_iteration", 30);

		Dirichlet_Boundary_Condition = poisson_eqn_block.GetBoolean("Dirichlet_Boundary_Condition", (bool)false);
		Neumann_Boundary_Condition = poisson_eqn_block.GetBoolean("Neumann_Boundary_Condition", (bool)false);

		// Poisson solver type from script
		const char* poisson_solver_type_input = poisson_eqn_block.GetString("poisson_solver_type", "Null");
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

		cout << "-------------------POISSON EQUATION TEST-------------------" << endl;
		cout << "tolerance: " << tolerance << endl;
		cout << "max_iteration: " << max_iteration << endl;

		if (Neumann_Boundary_Condition)
		{
			cout << "Neumann Boundary Condition : true" << endl;
			cout << "Dirichlet Boundary Condition : false" << endl;
		}

		if (Dirichlet_Boundary_Condition)
		{
			cout << "Neumann Boundary Condition : false" << endl;
			cout << "Dirichlet Boundary Condition : true" << endl;
		}

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

		if (grid_1d)
		{
			// Grid Initialization
			base_grid_1d.Initialize(world_discretization->world_grid_1d);
			
			// Jump Condition Field
			jc_on_solution_1d.Initialize(base_grid_1d, 1, multithreading);
			jc_on_derivative_1d.Initialize(base_grid_1d, 1, multithreading);
			jc_on_solution_1d.AssignAllValue(0);
			jc_on_derivative_1d.AssignAllValue(0); 
			
			// Initialize Fields
			solution_1d.Initialize(base_grid_1d, 1, multithreading);
			boundary_condition_1d.Initialize(base_grid_1d, 1, multithreading);
			rhs_1d.Initialize(base_grid_1d, 1, multithreading);
			beta_1d.Initialize(base_grid_1d, 1, multithreading);
			
			// Assigning initial field values
			beta_1d.array_for_this.AssignAllValues((T)1);

			// Interface Levelset
			DELETE_POINTER(interface_levelset_1d);
			interface_levelset_1d = new LEVELSET_1D();
			interface_levelset_1d->Initialize(base_grid_1d, 2, multithreading);

			// Initialize Interface Levelset
			interface_levelset_1d->AssignAllValuesLevelset(base_grid_1d.dx*(T)3);
			interface_levelset_1d->FillGhostCellsFromPointer(&(interface_levelset_1d->phi), false);
		}
		if (grid_2d)
		{
			// Grid Initialization
			base_grid_2d.Initialize(world_discretization->world_grid);

			// Jump Condition Field
			jc_on_solution_2d.Initialize(base_grid_2d, 1, multithreading);
			jc_on_derivative_2d.Initialize(base_grid_2d, 1, multithreading);
			jc_on_solution_2d.AssignAllValue(0);
			jc_on_derivative_2d.AssignAllValue(0); 
			
			// Initialize Fields
			solution_2d.Initialize(base_grid_2d, 1, multithreading);
			boundary_condition_2d.Initialize(base_grid_2d, 1, multithreading);
			rhs_2d.Initialize(base_grid_2d, 1, multithreading);
			beta_2d.Initialize(base_grid_2d, 1, multithreading);
			
			// Assigning initial field values
			beta_2d.array_for_this.AssignAllValues((T)1);
			
			// Interface Levelset
			DELETE_POINTER(interface_levelset_2d);
			interface_levelset_2d = new LEVELSET_2D();
			interface_levelset_2d->Initialize(base_grid_2d, 2, multithreading);

			// Initialize Interface Levelset
			interface_levelset_2d->AssignAllValuesLevelset(base_grid_2d.dx*(T)3);
			interface_levelset_2d->FillGhostCellsFromPointer(&(interface_levelset_2d->phi), false);
		}
		
		// Initialize Poisson solvers
		poisson_solver.Initialize(tolerance, max_iteration, 0, multithreading);
		poisson_solver.InitializeLinearSolver(poisson_solver_type);
	}

public: // Solver
	void Solve(const int& thread_id)
	{
		SetupRHS(thread_id);
		SetupBeta(thread_id);
		SetupJumpCondition(thread_id);
		DetermineSolution(thread_id);
	}

	void DetermineSolution(const int& thread_id)
	{
		if (grid_1d)
		{
			SetupBoundaryCondition(solution_1d, boundary_condition_1d, *interface_levelset_1d, thread_id);
			poisson_solver.Solve(solution_1d, beta_1d, boundary_condition_1d, rhs_1d, beta_1d, *interface_levelset_1d, jc_on_solution_1d, jc_on_derivative_1d, thread_id);
		}
		if (grid_2d)
		{
			SetupBoundaryCondition(solution_2d, boundary_condition_2d, *interface_levelset_2d, thread_id);
			interface_levelset_2d->ComputeNormals(thread_id);
			
			beta_2d.FillGhostCellsFrom(beta_2d.array_for_this, false, thread_id);
			multithreading->Sync(thread_id);
			
			poisson_solver.Solve(solution_2d,beta_2d, boundary_condition_2d, rhs_2d, *interface_levelset_2d, jc_on_solution_2d, jc_on_derivative_2d);
			//poisson_solver.Solve(solution_2d, beta_2d, boundary_condition_2d, rhs_2d, beta_2d, *interface_levelset_2d, jc_on_solution_2d, jc_on_derivative_2d, thread_id);
		}
	}

	void SetupBeta(const int& thread_id)
	{
		if (test_number == 1)
		{
			BEGIN_GRID_ITERATION_1D(beta_1d.partial_grids[thread_id])
			{
				beta_1d(i) = (T)1;
			}
			END_GRID_ITERATION_1D;
		}
		
		if (test_number == 3)
		{
			BEGIN_GRID_ITERATION_2D(beta_2d.partial_grids[thread_id])
			{
				if (interface_levelset_2d->arr(i, j) <= (T)0)
				{
					beta_2d(i, j) = (T)0.5;
				}
				else
				{
					beta_2d(i, j) = (T)1;
				}
			}
			END_GRID_ITERATION_2D;
		}

		if (test_number == 5 || test_number == 6 || test_number == 7)
		{
			BEGIN_GRID_ITERATION_2D(beta_2d.partial_grids[thread_id])
			{
				if (interface_levelset_2d->arr(i, j) <= (T)0)
				{
					beta_2d(i, j) = (T)1;
				}
				else
				{
					beta_2d(i, j) = (T)1;
				}
			}
			END_GRID_ITERATION_2D;
		}

		if (test_number == 8)
		{
			BEGIN_GRID_ITERATION_2D(beta_2d.partial_grids[thread_id])
			{
				if (interface_levelset_2d->arr(i, j) <= (T)0)
				{
					beta_2d(i, j) = (T)1;
				}
				else
				{
					beta_2d(i, j) = (T)1/10;
				}
			}
			END_GRID_ITERATION_2D;
		}
	}

	void SetupJumpCondition(const int& thread_id)
	{
		if (test_number == 1)
		{
			BEGIN_GRID_ITERATION_1D(jc_on_solution_1d.partial_grids[thread_id])
			{
				jc_on_solution_1d(i) = (T)1;
			}
			END_GRID_ITERATION_1D;

			BEGIN_GRID_ITERATION_1D(jc_on_derivative_1d.partial_grids[thread_id])
			{
				jc_on_derivative_1d(i) = (T)0;
			}
			END_GRID_ITERATION_1D;
		}
		
		if (test_number == 3)
		{
			BEGIN_GRID_ITERATION_2D(jc_on_solution_2d.partial_grids[thread_id])
			{
				T x_coor = jc_on_solution_2d.x_min + i*jc_on_solution_2d.dx, y_coor = jc_on_solution_2d.y_min + j*jc_on_solution_2d.dy;
					
				jc_on_solution_2d(i, j) = -exp(-POW2(x_coor) - POW2(y_coor));
			}
			END_GRID_ITERATION_2D;

			BEGIN_GRID_ITERATION_2D(jc_on_solution_2d.partial_grids[thread_id])
			{
				T x_coor = jc_on_derivative_2d.x_min + i*jc_on_derivative_2d.dx, y_coor = jc_on_derivative_2d.y_min + j*jc_on_derivative_2d.dy;
					
				jc_on_derivative_2d(i, j) = (T)8*(2*POW2(x_coor) + 2*POW2(y_coor) - x_coor - y_coor)*exp(-POW2(x_coor) - POW2(y_coor));
			}
			END_GRID_ITERATION_2D;
		}

		if (test_number == 5)
		{
			BEGIN_GRID_ITERATION_2D(jc_on_solution_2d.partial_grids[thread_id])
			{
				T x_coor = jc_on_solution_2d.x_min + i*jc_on_solution_2d.dx, y_coor = jc_on_solution_2d.y_min + j*jc_on_solution_2d.dy;
					
				jc_on_solution_2d(i, j) = (T)0;
			}
			END_GRID_ITERATION_2D;

			BEGIN_GRID_ITERATION_2D(jc_on_solution_2d.partial_grids[thread_id])
			{
				T x_coor = jc_on_derivative_2d.x_min + i*jc_on_derivative_2d.dx, y_coor = jc_on_derivative_2d.y_min + j*jc_on_derivative_2d.dy;
					
				jc_on_derivative_2d(i, j) = (T)2;
			}
			END_GRID_ITERATION_2D;
		}

		if (test_number == 6)
		{
			BEGIN_GRID_ITERATION_2D(jc_on_solution_2d.partial_grids[thread_id])
			{
				T x_coor = jc_on_solution_2d.x_min + i*jc_on_solution_2d.dx, y_coor = jc_on_solution_2d.y_min + j*jc_on_solution_2d.dy;
					
				jc_on_solution_2d(i, j) = -exp(x_coor)*cos(y_coor);
			}
			END_GRID_ITERATION_2D;

			BEGIN_GRID_ITERATION_2D(jc_on_solution_2d.partial_grids[thread_id])
			{
				T x_coor = jc_on_derivative_2d.x_min + i*jc_on_derivative_2d.dx, y_coor = jc_on_derivative_2d.y_min + j*jc_on_derivative_2d.dy;
					
				jc_on_derivative_2d(i, j) = (T)2*exp(x_coor)*(y_coor*sin(y_coor) - x_coor*cos(y_coor));
			}
			END_GRID_ITERATION_2D;
		}

		if (test_number == 7)
		{
			BEGIN_GRID_ITERATION_2D(jc_on_solution_2d.partial_grids[thread_id])
			{
				T x_coor = jc_on_solution_2d.x_min + i*jc_on_solution_2d.dx, y_coor = jc_on_solution_2d.y_min + j*jc_on_solution_2d.dy;
					
				jc_on_solution_2d(i, j) = POW2(y_coor) - POW2(x_coor);
			}
			END_GRID_ITERATION_2D;

			BEGIN_GRID_ITERATION_2D(jc_on_solution_2d.partial_grids[thread_id])
			{
				T x_coor = jc_on_derivative_2d.x_min + i*jc_on_derivative_2d.dx, y_coor = jc_on_derivative_2d.y_min + j*jc_on_derivative_2d.dy;
					
				jc_on_derivative_2d(i, j) = (T)4*(POW2(y_coor) - POW2(x_coor));
			}
			END_GRID_ITERATION_2D;
		}

		if (test_number == 8)
		{
			BEGIN_GRID_ITERATION_2D(jc_on_solution_2d.partial_grids[thread_id])
			{
				T x_coor = jc_on_solution_2d.x_min + i*jc_on_solution_2d.dx, y_coor = jc_on_solution_2d.y_min + j*jc_on_solution_2d.dy;
					
				jc_on_solution_2d(i, j) = (T)0.1*(POW2(x_coor) + POW2(y_coor)) - 0.01*log(2*sqrt(POW2(x_coor) + POW2(y_coor))) - (POW2(x_coor) + POW2(y_coor));
			}
			END_GRID_ITERATION_2D;

			BEGIN_GRID_ITERATION_2D(jc_on_solution_2d.partial_grids[thread_id])
			{
				T x_coor = jc_on_derivative_2d.x_min + i*jc_on_derivative_2d.dx, y_coor = jc_on_derivative_2d.y_min + j*jc_on_derivative_2d.dy;
					
				jc_on_derivative_2d(i, j) = ((T)4*(POW2(y_coor) + POW2(x_coor)) - 0.1/(POW2(y_coor) + POW2(x_coor)) - (T)2)*(x_coor*interface_levelset_2d->normal(i, j).x + y_coor*interface_levelset_2d->normal(i, j).y);
			}
			END_GRID_ITERATION_2D;
		}
	}

	void SetupRHS(const int& thread_id)
	{
		if (test_number == 1)
		{
			BEGIN_GRID_ITERATION_1D(rhs_1d.partial_grids[thread_id])
			{
				rhs_1d(i) = (T)0;
			}
			END_GRID_ITERATION_1D;
		}
		
		if (test_number == 3)
		{
			BEGIN_GRID_ITERATION_2D(rhs_2d.partial_grids[thread_id])
			{
				if (interface_levelset_2d->arr(i, j) <= (T)0)
				{
					T x_coor = rhs_2d.x_min + i*rhs_2d.dx, y_coor = rhs_2d.y_min + j*rhs_2d.dy;
					
					rhs_2d(i, j) = (T)8*(POW2(x_coor) + POW2(y_coor) - 1)*exp(-POW2(x_coor) - POW2(y_coor));
				}
				else
				{
					rhs_2d(i, j) = (T)0;
				}
			}
			END_GRID_ITERATION_2D;
		}

		if (test_number == 5 || test_number == 6 || test_number == 7)
		{
			BEGIN_GRID_ITERATION_2D(rhs_2d.partial_grids[thread_id])
			{
				rhs_2d(i, j) = (T)0;
			}
			END_GRID_ITERATION_2D;
		}

		if (test_number == 8)
		{
			BEGIN_GRID_ITERATION_2D(rhs_2d.partial_grids[thread_id])
			{
				if (interface_levelset_2d->arr(i, j) <= (T)0)
				{
					rhs_2d(i, j) = (T)4;
				}
				else
				{
					T x_coor = rhs_2d.x_min + i*rhs_2d.dx, y_coor = rhs_2d.y_min + j*rhs_2d.dy;

					rhs_2d(i, j) = (T)16*(POW2(x_coor) + POW2(y_coor));
				}
			}
			END_GRID_ITERATION_2D;
		}
	}

public: // Member Functions
	void SetupBoundaryCondition(FIELD_STRUCTURE_1D<T>& pressure_input, FIELD_STRUCTURE_1D<int>& bc_input, const LEVELSET_1D& water_levelset_input, const int& thread_id)
	{
		ARRAY_1D<int>& bc_array(bc_input.array_for_this);
		GRID_STRUCTURE_1D& grid(bc_input.grid);

		if (Dirichlet_Boundary_Condition)
		{
			BEGIN_GRID_ITERATION_1D(bc_input.partial_grids_ghost[thread_id])
			{
				// Speed-up variable
				if (i < grid.i_start || i > grid.i_end)
				{
					bc_array(i) = BC_DIR;
					pressure_input(i) = (T)0;
				}
				else
				{
					bc_array(i) = BC_FULL;
				}
			}
			END_GRID_ITERATION_1D;
		}
		
		if (Neumann_Boundary_Condition)
		{
			BEGIN_GRID_ITERATION_1D(bc_input.partial_grids_ghost[thread_id])
			{
				// Speed-up variable
				if (i < grid.i_start || i > grid.i_end)
				{
					bc_array(i) = BC_NEUM;
				}
				else
				{
					bc_array(i) = BC_FULL;
				}
			}
			END_GRID_ITERATION_1D;
		}
	}

	void SetupBoundaryCondition(FIELD_STRUCTURE_2D<T>& pressure_input, FIELD_STRUCTURE_2D<int>& bc_input, const LEVELSET_2D& water_levelset_input, const int& thread_id)
	{
		ARRAY_2D<int>& bc_array(bc_input.array_for_this);
		GRID_STRUCTURE_2D& grid(bc_input.grid);

		if (Dirichlet_Boundary_Condition)
		{
			BEGIN_GRID_ITERATION_2D(bc_input.partial_grids_ghost[thread_id])
			{
				// Speed-up variable
				if (i < grid.i_start || i > grid.i_end || j < grid.j_start || j > grid.j_end)
				{
					bc_array(i, j) = BC_DIR;
					
					if (test_number == 5)
					{
						T x_coor = bc_input.x_min + i*bc_input.dx, y_coor = bc_input.y_min + j*bc_input.dy;
						pressure_input(i, j) = 1 + log((T)2*sqrt(POW2(x_coor) + POW2(y_coor)));
					}
					else if (test_number == 8)
					{
						T x_coor = bc_input.x_min + i*bc_input.dx, y_coor = bc_input.y_min + j*bc_input.dy;
						pressure_input(i, j) = 0.1*POW2(POW2(x_coor) + POW2(y_coor)) - (T)0.01*log((T)2*sqrt(POW2(x_coor) + POW2(y_coor)));
					}
					else
					{
						pressure_input(i, j) = (T)0;
					}
				}
				else
				{
					bc_array(i, j) = BC_FULL;
				}
			}
			END_GRID_ITERATION_2D;
		}
		
		if (Neumann_Boundary_Condition)
		{
			BEGIN_GRID_ITERATION_2D(bc_input.partial_grids_ghost[thread_id])
			{
				// Speed-up variable
				if (i < grid.i_start || i > grid.i_end || j < grid.j_start || j > grid.j_end)
				{
					bc_array(i, j) = BC_NEUM;
				}
				else
				{
					bc_array(i, j) = BC_FULL;
				}
			}
			END_GRID_ITERATION_2D;
		}
	}
};