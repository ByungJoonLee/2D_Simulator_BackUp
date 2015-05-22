#pragma once

#include "LINEAR_SOLVER.h"
#include "LEVELSET_2D.h"
#include "CG_METHOD.h"
#include "PCG_METHOD.h"

class POISSON_SOLVER_2D
{
public: // Essential Data
	T									tolerance;
	T									sqr_tolerance;
	int									max_iteration;
	
	GRID_STRUCTURE_2D*					grid_ghost;
	int									ghost_width;

public: // For matrix-vector type solver
	CSR_MATRIX<T>						A;
	VECTOR_ND<T>						x;
	VECTOR_ND<T>						b;

public: // Solver
	LINEAR_SOLVER*						linear_solver;

public: // Density 
	T									density_p, density_m;

public: // Constructors and Destructor
	POISSON_SOLVER_2D(void)
		: tolerance((T)1e-4), sqr_tolerance(tolerance*tolerance), max_iteration(100), linear_solver(0)
	{}

	~POISSON_SOLVER_2D(void)
	{
		DELETE_POINTER(linear_solver);
	}

public: // Initialization Functions
	void Initialize(const T& tolerance_input, const int& max_itr_input, GRID_STRUCTURE_2D* grid_ghost_input = 0, const int ghost_width_input = 0)
	{
		ghost_width = ghost_width_input;

		max_iteration = max_itr_input; 

		InitializeLinearSolver(CG);

		SetTolerance(tolerance_input);
	}

	void InitializeLinearSolver(const POISSON_SOLVER_TYPE linear_solver_type)
	{
		DELETE_POINTER(linear_solver);
		
		switch (linear_solver_type)
		{
		case CG:
			linear_solver = new CG_METHOD();
			break;
		
		case PCG:
			linear_solver = new PCG_METHOD();
			break;
		
		case NO_SOLVER:
		
		default:
			linear_solver = 0;
			break;
		}

		linear_solver->Initialize(tolerance, max_iteration);
	}

	void InitializePressure(FIELD_STRUCTURE_2D<T>& pressure)
	{
		int i, j;
		LOOPS_2D(i, j, pressure.i_start, pressure.j_start, pressure.i_end, pressure.j_end)
		{
			pressure(i, j) = (T)0;
		}	
	}

public: // Solver
	void Solve(FIELD_STRUCTURE_2D<T>& pressure, FIELD_STRUCTURE_2D<T>& density, FIELD_STRUCTURE_2D<int>& bc, const FIELD_STRUCTURE_2D<T>& div, const LEVELSET_2D& levelset, const FIELD_STRUCTURE_2D<T>& jc_on_solution, const FIELD_STRUCTURE_2D<T>& jc_on_derivative)
	{
		assert(linear_solver != 0);

		//BuildLinearSystemNodeDirichlet(A, x, b, pressure, bc, div);	

   		//BuildLinearSystemNodeJumpCondition(A, x, b, pressure, bc, div, levelset, jc_on_solution, jc_on_derivative);
		
		FIELD_STRUCTURE_2D<T> one_over_density;
		one_over_density.Initialize(density.grid, 2);
		
		GRID_ITERATION_2D(one_over_density.grid)
		{
			one_over_density(i, j) = (T)1/density(i, j);
		}

		one_over_density.FillGhostCellsFrom(one_over_density.array_for_this, false);

		BuildLinearSystemNodeJumpConditionVaribleCoefficient(A, x, b, pressure, one_over_density, bc, div, levelset, jc_on_solution, jc_on_derivative);
		
		linear_solver->Solve(A, x, b, bc);
		
		VectorToGrid(x, pressure, bc);
	}

	void SolveForViscosity(FIELD_STRUCTURE_2D<T>& velocity, const FIELD_STRUCTURE_2D<T>& coef_1, const FIELD_STRUCTURE_2D<T>& coef_2, const FIELD_STRUCTURE_2D<T>& coef_3, const FIELD_STRUCTURE_2D<T>& coef_4, const FIELD_STRUCTURE_2D<T>& coef_5, FIELD_STRUCTURE_2D<int>& bc, const FIELD_STRUCTURE_2D<T>& explicit_term)
	{
		assert(linear_solver != 0);

		BuildLinearSystemNodeForSemiImplicitViscosity(A, x, b, coef_1, coef_2, coef_3, coef_4, coef_5, velocity, bc, explicit_term);
		
		linear_solver->Solve(A, x, b, bc);
		
		VectorToGrid(x, velocity, bc);
	}

	void SolveForAxisymmetric(FIELD_STRUCTURE_2D<T>& pressure, const FIELD_STRUCTURE_2D<T>& coef_1, const FIELD_STRUCTURE_2D<T>& coef_2, const FIELD_STRUCTURE_2D<T>& coef_3, const FIELD_STRUCTURE_2D<T>& coef_4, FIELD_STRUCTURE_2D<int>& bc, const FIELD_STRUCTURE_2D<T>& div)
	{
		assert(linear_solver != 0);

		BuildLinearSystemNodeForAxiSymmetric(A, x, b, coef_1, coef_2, coef_3, coef_4, pressure, bc, div);

		linear_solver->Solve(A, x, b, bc);
		VectorToGrid(x, pressure, bc);
	}

public: // Member Functions
	void SetTolerance(const T& tolerance_input)
	{
		tolerance = tolerance_input;
		sqr_tolerance = tolerance*tolerance;

		linear_solver->SetTolerance(tolerance);
	}

	void BuildLinearSystemNodeForAxiSymmetric(CSR_MATRIX<T>& A_matrix, VECTOR_ND<T>& x_vector, VECTOR_ND<T>& b_vector, const FIELD_STRUCTURE_2D<T>& coef_1, const FIELD_STRUCTURE_2D<T>& coef_2, const FIELD_STRUCTURE_2D<T>& coef_3, const FIELD_STRUCTURE_2D<T>& coef_4, const FIELD_STRUCTURE_2D<T>& pressure, const FIELD_STRUCTURE_2D<int>& bc, const FIELD_STRUCTURE_2D<T>& div)
	{
		const int num_all_full_cells = AssignSequentialindicesToFullCells(bc);
		const int nnz = CountNonZeroElements(bc);
		
		A_matrix.Initialize(num_all_full_cells, nnz);
		x_vector.Initialize(num_all_full_cells, true);
		b_vector.Initialize(num_all_full_cells);

		// Speed up variables 
		int i_res(pressure.grid.i_res), j_res(pressure.grid.j_res), i_start(pressure.i_start), i_end(pressure.i_end), j_start(pressure.j_start), j_end(pressure.j_end);
		T dx(pressure.dx), dy(pressure.dy), one_over_dx(pressure.one_over_dx), one_over_dy(pressure.one_over_dy), one_over_dx2(pressure.one_over_dx2), one_over_dy2(pressure.one_over_dy2);
		T x_min(pressure.grid.x_min), y_min(pressure.grid.y_min), x_max(pressure.grid.x_max), y_max(pressure.grid.y_max);

		const T dxdx = POW2(pressure.grid.dx), inv_dxdx = (T)1/dxdx, dydy = POW2(pressure.grid.dy), inv_dydy = (T)1/dydy;
		
		GRID_ITERATION_2D(pressure.grid)
		{
			if (bc(i, j) < 0)
			{
				continue;
			}

			T coef_ijk = 0;					// For optimization, inv_dxdx is multiplied at the end
				
			// If neighbor is full cell
			if (bc(i - 1, j) > -1)
			{
				coef_ijk += coef_2(i, j);
				A_matrix.AssignValue(bc(i, j), bc(i - 1, j), -coef_2(i, j)*inv_dxdx);
			}
			if (bc(i + 1, j) > -1)
			{
				coef_ijk += coef_1(i, j);
				A_matrix.AssignValue(bc(i, j), bc(i + 1, j), -coef_1(i, j)*inv_dxdx);
			}
			if (bc(i, j - 1) > -1)
			{
				coef_ijk += coef_4(i, j);
				A_matrix.AssignValue(bc(i, j), bc(i, j - 1), -coef_4(i, j)*inv_dydy);
			}
			if (bc(i, j + 1) > -1)
			{
				coef_ijk += coef_3(i, j);
				A_matrix.AssignValue(bc(i, j), bc(i, j + 1), -coef_3(i, j)*inv_dydy);
			}
			
			if (bc(i - 1, j) == BC_DIR)
			{
				coef_ijk += coef_2(i, j);
			}
			if (bc(i + 1, j) == BC_DIR)
			{
				coef_ijk += coef_1(i, j);
			}
			if (bc(i, j - 1) == BC_DIR)
			{
				coef_ijk += coef_4(i, j);
			}
			if (bc(i, j + 1) == BC_DIR)
			{
				coef_ijk += coef_3(i, j);
			}
						
			if (coef_ijk == 0)
			{
				coef_ijk = 1;
			}

			A_matrix.AssignValue(bc(i, j), bc(i, j), inv_dxdx*(T)coef_ijk);
			
			b_vector[bc(i, j)] = div(i, j);

			if (bc(i - 1, j) == BC_DIR)
			{
				b_vector[bc(i, j)] += coef_2(i, j)*inv_dxdx*pressure(i - 1, j);
			}
			if (bc(i + 1, j) == BC_DIR)
			{
				b_vector[bc(i, j)] += coef_1(i, j)*inv_dxdx*pressure(i + 1, j);
			}
			if (bc(i, j - 1) == BC_DIR)
			{
				b_vector[bc(i, j)] += coef_4(i, j)*inv_dydy*pressure(i, j - 1);
			}
			if (bc(i, j + 1) == BC_DIR)
			{
				b_vector[bc(i, j)] += coef_3(i, j)*inv_dydy*pressure(i, j + 1);
			}
								
			x_vector[bc(i, j)] = pressure(i, j);
		
		}
	}

	void BuildLinearSystemNodeForSemiImplicitViscosity(CSR_MATRIX<T>& A_matrix, VECTOR_ND<T>& x_vector, VECTOR_ND<T>& b_vector, const FIELD_STRUCTURE_2D<T>& coef_1, const FIELD_STRUCTURE_2D<T>& coef_2, const FIELD_STRUCTURE_2D<T>& coef_3, const FIELD_STRUCTURE_2D<T>& coef_4, const FIELD_STRUCTURE_2D<T>& coef_5, const FIELD_STRUCTURE_2D<T>& velocity, const FIELD_STRUCTURE_2D<int>& bc, const FIELD_STRUCTURE_2D<T>& explicit_term)
	{
		const int num_all_full_cells = AssignSequentialindicesToFullCells(bc);
		const int nnz = CountNonZeroElements(bc);
		
		A_matrix.Initialize(num_all_full_cells, nnz);
		x_vector.Initialize(num_all_full_cells, true);
		b_vector.Initialize(num_all_full_cells);

		// Speed up variables 
		int i_res(velocity.grid.i_res), j_res(velocity.grid.j_res), i_start(velocity.i_start), i_end(velocity.i_end), j_start(velocity.j_start), j_end(velocity.j_end);
		T dx(velocity.dx), dy(velocity.dy), one_over_dx(velocity.one_over_dx), one_over_dy(velocity.one_over_dy), one_over_dx2(velocity.one_over_dx2), one_over_dy2(velocity.one_over_dy2);
		T x_min(velocity.grid.x_min), y_min(velocity.grid.y_min), x_max(velocity.grid.x_max), y_max(velocity.grid.y_max);

		const T dxdx = POW2(velocity.grid.dx), inv_dxdx = (T)1/dxdx, dydy = POW2(velocity.grid.dy), inv_dydy = (T)1/dydy;
		
		GRID_ITERATION_2D(velocity.grid)
		{
			if (bc(i, j) < 0)
			{
				continue;
			}

			T coef_ijk = 0;					// For optimization, inv_dxdx is multiplied at the end
				
			// If neighbor is full cell
			if (bc(i - 1, j) > -1)
			{
				coef_ijk += coef_2(i, j);
				A_matrix.AssignValue(bc(i, j), bc(i - 1, j), -coef_2(i, j)*inv_dxdx);
			}
			if (bc(i + 1, j) > -1)
			{
				coef_ijk += coef_1(i, j);
				A_matrix.AssignValue(bc(i, j), bc(i + 1, j), -coef_1(i, j)*inv_dxdx);
			}
			if (bc(i, j - 1) > -1)
			{
				coef_ijk += coef_4(i, j);
				A_matrix.AssignValue(bc(i, j), bc(i, j - 1), -coef_4(i, j)*inv_dydy);
			}
			if (bc(i, j + 1) > -1)
			{
				coef_ijk += coef_3(i, j);
				A_matrix.AssignValue(bc(i, j), bc(i, j + 1), -coef_3(i, j)*inv_dydy);
			}
			
			if (bc(i - 1, j) == BC_DIR)
			{
				coef_ijk += coef_2(i, j);
			}
			if (bc(i + 1, j) == BC_DIR)
			{
				coef_ijk += coef_1(i, j);
			}
			if (bc(i, j - 1) == BC_DIR)
			{
				coef_ijk += coef_4(i, j);
			}
			if (bc(i, j + 1) == BC_DIR)
			{
				coef_ijk += coef_3(i, j);
			}
						
			if (coef_ijk == 0)
			{
				coef_ijk = 1;
			}

			A_matrix.AssignValue(bc(i, j), bc(i, j), (T)1 + inv_dxdx*(T)coef_ijk + coef_5(i, j));
			
			b_vector[bc(i, j)] = explicit_term(i, j);

			if (bc(i - 1, j) == BC_DIR)
			{
				b_vector[bc(i, j)] += coef_2(i, j)*inv_dxdx*velocity(i - 1, j);
			}
			if (bc(i + 1, j) == BC_DIR)
			{
				b_vector[bc(i, j)] += coef_1(i, j)*inv_dxdx*velocity(i + 1, j);
			}
			if (bc(i, j - 1) == BC_DIR)
			{
				b_vector[bc(i, j)] += coef_4(i, j)*inv_dydy*velocity(i, j - 1);
			}
			if (bc(i, j + 1) == BC_DIR)
			{
				b_vector[bc(i, j)] += coef_3(i, j)*inv_dydy*velocity(i, j + 1);
			}
								
			x_vector[bc(i, j)] = velocity(i, j);
		}
	}

	void BuildLinearSystemNodeJumpConditionVaribleCoefficient(CSR_MATRIX<T>& A_matrix, VECTOR_ND<T>& x_vector, VECTOR_ND<T>& b_vector, const FIELD_STRUCTURE_2D<T>& pressure, const FIELD_STRUCTURE_2D<T>& variable, const FIELD_STRUCTURE_2D<int>& bc, const FIELD_STRUCTURE_2D<T>& div, const LEVELSET_2D& levelset, const FIELD_STRUCTURE_2D<T>& jc_on_solution, const FIELD_STRUCTURE_2D<T>& jc_on_derivative)
	{
		const int num_all_full_cells = AssignSequentialindicesToFullCells(bc);
		const int nnz = CountNonZeroElements(bc);
		
		A_matrix.Initialize(num_all_full_cells, nnz);
		x_vector.Initialize(num_all_full_cells, true);
		b_vector.Initialize(num_all_full_cells);

		// Speed up variables 
		int i_start(pressure.i_start), i_end(pressure.i_end), j_start(pressure.j_start), j_end(pressure.j_end);
		T dx(pressure.dx), dy(pressure.dy), one_over_dx(pressure.one_over_dx), one_over_dy(pressure.one_over_dy), one_over_dx2(pressure.one_over_dx2), one_over_dy2(pressure.one_over_dy2);

		const T dxdx = POW2(pressure.grid.dx), inv_dxdx = (T)1/dxdx, dydy = POW2(pressure.grid.dy), inv_dydy = (T)1/dydy;
		
		int i, j;
		LOOPS_2D(i, j, i_start, j_start, i_end, j_end)
		{
			if (bc(i, j) < 0)
			{
				continue;
			}

			T coef_ij = 0;					// For optimization, inv_dxdx is multiplied at the end
				
			// Define betas
			T beta_l, beta_r, beta_b, beta_t;
			T levelset_c = levelset(i, j), levelset_l = levelset(i - 1, j), levelset_r = levelset(i + 1, j), levelset_b = levelset(i, j - 1), levelset_t = levelset(i, j + 1);

			beta_l = variable(i, j)*variable(i - 1, j)*(abs(levelset_l) + abs(levelset_c))/(variable(i, j)*abs(levelset_l) + variable(i - 1, j)*abs(levelset_c));
			beta_r = variable(i, j)*variable(i + 1, j)*(abs(levelset_r) + abs(levelset_c))/(variable(i + 1, j)*abs(levelset_c) + variable(i, j)*abs(levelset_r));
			beta_b = variable(i, j)*variable(i, j - 1)*(abs(levelset_b) + abs(levelset_c))/(variable(i, j)*abs(levelset_b) + variable(i, j - 1)*abs(levelset_c));
			beta_t = variable(i, j)*variable(i, j + 1)*(abs(levelset_t) + abs(levelset_c))/(variable(i, j + 1)*abs(levelset_c) + variable(i, j)*abs(levelset_t));

			// If neighbor is full cell
			if (bc(i - 1, j) > -1)
			{
				coef_ij += beta_l;
				A_matrix.AssignValue(bc(i, j), bc(i - 1, j), -beta_l*inv_dxdx);
			}
			if (bc(i + 1, j) > -1)
			{
				coef_ij += beta_r;
				A_matrix.AssignValue(bc(i, j), bc(i + 1, j), -beta_r*inv_dxdx);
			}
			if (bc(i, j - 1) > -1)
			{
				coef_ij += beta_b;
				A_matrix.AssignValue(bc(i, j), bc(i, j - 1), -beta_b*inv_dydy);
			}
			if (bc(i, j + 1) > -1)
			{
				coef_ij += beta_t;
				A_matrix.AssignValue(bc(i, j), bc(i, j + 1), -beta_t*inv_dydy);
			}
			
			// Dirichlet Boundary Condition
			if (bc(i - 1, j) == BC_DIR)
			{
				coef_ij += beta_l;
			}
			if (bc(i + 1, j) == BC_DIR)
			{
				coef_ij += beta_r;
			}
			if (bc(i, j - 1) == BC_DIR)
			{
				coef_ij += beta_b;
			}
			if (bc(i, j + 1) == BC_DIR)
			{
				coef_ij += beta_t;
			}
			
			// Neumann Boudary Condition
			if ((i == i_start) && (j == j_start))
			{
				if (bc(i - 1, j) == BC_NEUM)
				{
					coef_ij += beta_l;
				}
				if (bc(i + 1, j) == BC_NEUM)
				{
					coef_ij += beta_r;
				}
				if (bc(i, j - 1) == BC_NEUM)
				{
					coef_ij += beta_b;
				}
				if (bc(i, j + 1) == BC_NEUM)
				{
					coef_ij += beta_t;
				}
			}
			else
			{
				if (bc(i - 1, j) == BC_NEUM)
				{
					coef_ij += 0;
				}
				if (bc(i + 1, j) == BC_NEUM)
				{
					coef_ij += 0;
				}
				if (bc(i, j - 1) == BC_NEUM)
				{
					coef_ij += 0;
				}
				if (bc(i, j + 1) == BC_NEUM)
				{
					coef_ij += 0;
				}
			}
						
			if (coef_ij == 0)
			{
				coef_ij = 1;
			}

			A_matrix.AssignValue(bc(i, j), bc(i, j), inv_dxdx*(T)coef_ij);
			
			// Divide given region into different region - Boundary Capturing method for Poisson Equation on irregular domain
			T F_L, F_R, F_B ,F_T;
				
			// Left arm Stencil
			T subcell_l = abs(levelset(i-1, j))/(abs(levelset(i, j)) + abs(levelset(i-1, j)));
			T a_l = (jc_on_solution(i, j)*abs(levelset(i-1, j)) + jc_on_solution(i-1, j)*abs(levelset(i, j)))/(abs(levelset(i, j)) + abs(levelset(i-1, j)));
			T b_l = (jc_on_derivative(i, j)*levelset.normal(i, j).x*abs(levelset(i-1, j)) + jc_on_derivative(i-1, j)*levelset.normal(i-1, j).x*abs(levelset(i, j)))/(abs(levelset(i, j)) + abs(levelset(i-1, j)));
				
			if (((levelset(i, j) <= 0) && (levelset(i-1, j) <= 0)) || ((levelset(i, j) > 0) && (levelset(i-1, j) > 0)))
			{
				F_L = 0;
			}
			else if ((levelset(i, j) <= 0) && (levelset(i-1, j) > 0))
			{
				F_L = one_over_dx2*a_l*beta_l - one_over_dx*beta_l*b_l*subcell_l/variable(i - 1, j);
			}
			else if ((levelset(i, j) > 0) && (levelset(i-1, j) <= 0))
			{
				F_L = -one_over_dx2*a_l*beta_l + one_over_dx*beta_l*b_l*subcell_l/variable(i - 1, j);
			}

			// Right arm Stencil
			T subcell_r = abs(levelset(i+1, j))/(abs(levelset(i, j)) + abs(levelset(i+1, j)));
			T a_r = (jc_on_solution(i, j)*abs(levelset(i+1, j)) + jc_on_solution(i+1, j)*abs(levelset(i, j)))/(abs(levelset(i, j)) + abs(levelset(i+1, j)));
			T b_r = (jc_on_derivative(i, j)*levelset.normal(i, j).x*abs(levelset(i+1, j)) + jc_on_derivative(i+1, j)*levelset.normal(i+1, j).x*abs(levelset(i, j)))/(abs(levelset(i, j)) + abs(levelset(i+1, j)));
			
			if (((levelset(i, j) <= 0) && (levelset(i+1, j) <= 0)) || ((levelset(i, j) > 0) && (levelset(i+1, j) > 0)))
			{
				F_R = 0;
			}
			else if ((levelset(i, j) <= 0) && (levelset(i+1, j) > 0))
			{
				F_R = one_over_dx2*a_r*beta_r + one_over_dx*beta_r*b_r*subcell_r/variable(i + 1, j);
			}
			else if ((levelset(i, j) > 0) && (levelset(i+1, j) <= 0))
			{
				F_R = -one_over_dx2*a_r*beta_r - one_over_dx*beta_r*b_r*subcell_r/variable(i + 1, j);
			}

			// Bottom arm Stencil
			T subcell_b = abs(levelset(i, j-1))/(abs(levelset(i, j)) + abs(levelset(i, j-1)));
			T a_b = (jc_on_solution(i, j)*abs(levelset(i, j-1)) + jc_on_solution(i, j-1)*abs(levelset(i, j)))/(abs(levelset(i, j)) + abs(levelset(i, j-1)));
			T b_b = (jc_on_derivative(i, j)*levelset.normal(i, j).y*abs(levelset(i, j-1)) + jc_on_derivative(i, j-1)*levelset.normal(i, j-1).y*abs(levelset(i, j)))/(abs(levelset(i, j)) + abs(levelset(i, j-1)));
				
			if (((levelset(i, j) <= 0) && (levelset(i, j-1) <= 0)) || ((levelset(i, j) > 0) && (levelset(i, j-1) > 0)))
			{
				F_B = 0;
			}
			else if ((levelset(i, j) <= 0) && (levelset(i, j-1) > 0))
			{
				F_B = one_over_dy2*a_b*beta_b - one_over_dy*beta_b*b_b*subcell_b/variable(i, j - 1);
			}
			else if ((levelset(i, j) > 0) && (levelset(i, j-1) <= 0))
			{
				F_B = -one_over_dy2*a_b*beta_b + one_over_dy*beta_b*b_b*subcell_b/variable(i, j - 1);
			}

			// Top arm Stencil
			T subcell_t = abs(levelset(i, j+1))/(abs(levelset(i, j)) + abs(levelset(i, j+1)));
			T a_t = (jc_on_solution(i, j)*abs(levelset(i, j+1)) + jc_on_solution(i, j+1)*abs(levelset(i, j)))/(abs(levelset(i, j)) + abs(levelset(i, j+1)));
			T b_t = (jc_on_derivative(i, j)*levelset.normal(i, j).y*abs(levelset(i, j+1)) + jc_on_derivative(i, j+1)*levelset.normal(i, j+1).y*abs(levelset(i, j)))/(abs(levelset(i, j)) + abs(levelset(i, j+1)));
				
			if (((levelset(i, j) <= 0) && (levelset(i, j+1) <= 0)) || ((levelset(i, j) > 0) && (levelset(i, j+1) > 0)))
			{
				F_T = 0;
			}
			else if ((levelset(i, j) <= 0) && (levelset(i, j+1) > 0))
			{
				F_T = one_over_dy2*a_t*beta_t + one_over_dy*beta_t*b_t*subcell_t/variable(i, j + 1);
			}
			else if ((levelset(i, j) > 0) && (levelset(i, j+1) <= 0))
			{
				F_T = -one_over_dy2*a_t*beta_t - one_over_dy*beta_t*b_t*subcell_t/variable(i, j + 1);
			}
				
			T F_X = F_L + F_R, F_Y = F_B + F_T;

			b_vector[bc(i, j)] = -div(i, j) - F_X - F_Y;
			
			if (bc(i - 1, j) == BC_DIR)
			{
				b_vector[bc(i, j)] += beta_l*inv_dxdx*pressure(i - 1, j);
				//b_vector[bc(i, j)] += variable(i - 1, j)*inv_dxdx*pressure(i - 1, j);
			}
			if (bc(i + 1, j) == BC_DIR)
			{
				b_vector[bc(i, j)] += beta_r*inv_dxdx*pressure(i + 1, j);
				//b_vector[bc(i, j)] += variable(i + 1, j)*inv_dxdx*pressure(i + 1, j);
			}
			if (bc(i, j - 1) == BC_DIR)
			{
				b_vector[bc(i, j)] += beta_b*inv_dydy*pressure(i, j - 1);
				//b_vector[bc(i, j)] += variable(i, j - 1)*inv_dydy*pressure(i, j - 1);
			}
			if (bc(i, j + 1) == BC_DIR)
			{
				b_vector[bc(i, j)] += beta_t*inv_dydy*pressure(i, j + 1);
				//b_vector[bc(i, j)] += variable(i, j + 1)*inv_dydy*pressure(i, j + 1);
			}
					
			x_vector[bc(i, j)] = pressure(i, j);
		}	
	}

	// This is amazing technique:)
	int AssignSequentialindicesToFullCells(const FIELD_STRUCTURE_2D<int>& bc)
	{
		// Count number of full cells
		int full_ix(0);

		int i, j;
		LOOPS_2D(i, j, bc.i_start, bc.j_start, bc.i_end, bc.j_end)
		{
			if (bc(i, j) > -1)
			{
				++full_ix;
			}	
		}

		// Indexing the value of boundary condition field
		int start_full_ix(0);

		LOOPS_2D(i, j, bc.i_start, bc.j_start, bc.i_end, bc.j_end)
		{
			if (bc(i, j) > -1)
			{
				bc(i, j) = start_full_ix++;
			}		
		}

		return full_ix;
	}

	int CountNonZeroElements(const FIELD_STRUCTURE_2D<int>& bc)
	{
		int nnz(0);

		int i, j;
		LOOPS_2D(i, j, bc.i_start, bc.j_start, bc.i_end, bc.j_end)
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
		
		return nnz;
	}

	void VectorToGrid(const VECTOR_ND<T>& x, FIELD_STRUCTURE_2D<T>& pressure, const FIELD_STRUCTURE_2D<int>& bc)
	{
		int i, j;
		LOOPS_2D(i, j, pressure.i_start, pressure.j_start, pressure.i_end, pressure.j_end)
		{
			pressure.array_for_this(i, j) = x[bc(i, j)];
		}
	}
};