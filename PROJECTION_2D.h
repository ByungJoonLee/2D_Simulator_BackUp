#pragma once

#include "LEVELSET_2D.h"
#include "POISSON_SOLVER_2D.h"

enum PROJECTION_TYPE						{FREE_SURFACE_WATER, AIR, MULTIPHASE};

class PROJECTION_2D
{
public: // References of the variables
	bool							use_water_solver;
	bool							use_multiphase;

	// For MAC grid
	FIELD_STRUCTURE_2D<T>&			water_velocity_field_mac_x;
	FIELD_STRUCTURE_2D<T>&			water_velocity_field_mac_y;
	FIELD_STRUCTURE_2D<T>&			vector_field_mac_ghost_x;
	FIELD_STRUCTURE_2D<T>&			vector_field_mac_ghost_y;

	LEVELSET_2D&					water_levelset;

	// For Vortex Sheet Problem
	LEVELSET_2D&					vortex_levelset;

public: // Variables defined and used in this class
	FIELD_STRUCTURE_2D<T>			pressure_field;
	FIELD_STRUCTURE_2D<T>			projection_density_field;
	FIELD_STRUCTURE_2D<int>			boundary_condition_field;
	FIELD_STRUCTURE_2D<T>			divergence_field;

public: // Jump Condition Field
	FIELD_STRUCTURE_2D<T>			jc_on_solution;
	FIELD_STRUCTURE_2D<T>			jc_on_derivative;
	
public: // Control Options
	enum POISSON_SOLVER_TYPE        poisson_solver_type;
	enum PROJECTION_TYPE            projection_type;

public: // Subsolver
	POISSON_SOLVER_2D				poisson_solver;

public: // Convenient variables and references
	GRID_STRUCTURE_2D&				base_grid;
	const T							dx, dy;
	const T							dxdx;
	const T							one_over_dx, one_over_dy;
	const T							one_over_2dx, one_over_2dy;
	const T							half_dx;
	const T							inv_2dx;
	const T							inv_2dy;

	ARRAY_2D<T>&					pressure_array;
	ARRAY_2D<T>&					divergence_array;
	ARRAY_2D<int>&					boundary_condition_array;
	ARRAY_2D<T>&					pressure_density_array;

	T								tolerance;
	int								max_iteration;

	T								water_density, air_density, oil_density;
	T								water_viscosity, air_viscosity, oil_viscosity;
	T								surface_tension;

	// Max Velocity
	T								max_velocity_x;
	T								max_velocity_y;

	// Options For Simulation
	bool							air_water_simulation;
	bool							oil_water_simulation;
	bool							vortex_sheet_problem;
	bool							dimensionless_form;

	// Option
	bool							air_bubble_rising;
	bool							water_drop;
	
	// Option for Viscosity treatment
	bool							use_delta_function_formulation;
	bool							use_jump_condition_on_viscosity;

	// Option for continuous surface tension force
	bool							CSF_model;

	// Option for Boundary Condition
	bool							Dirichlet_Boundary_Condition, Neumann_Boundary_Condition;

public: // Dimensionless Variable
	T								We;
	T								R1;

public: // Constructors and Destructor
	PROJECTION_2D(LEVELSET_2D& water_levelset_input, LEVELSET_2D& vortex_levelset_input, FIELD_STRUCTURE_2D<T>& velocity_field_mac_x_input, FIELD_STRUCTURE_2D<T>& velocity_field_mac_y_input, FIELD_STRUCTURE_2D<T>& velocity_field_ghost_mac_x_input, FIELD_STRUCTURE_2D<T>& velocity_field_ghost_mac_y_input)
		: water_velocity_field_mac_x(velocity_field_mac_x_input),
		  water_velocity_field_mac_y(velocity_field_mac_y_input),
		  vector_field_mac_ghost_x(velocity_field_ghost_mac_x_input),
		  vector_field_mac_ghost_y(velocity_field_ghost_mac_y_input),
		  water_levelset(water_levelset_input), 
		  vortex_levelset(vortex_levelset_input), 
		  base_grid(water_levelset_input.grid),
		  dx(base_grid.dx), dy(base_grid.dy),
		  dxdx(base_grid.dx*base_grid.dx),
		  one_over_dx((T)1/dx), one_over_dy((T)1/dy),
		  one_over_2dx((T)0.5*one_over_dx), one_over_2dy((T)0.5*one_over_dy),
		  inv_2dx((T)0.5/dx), inv_2dy((T)0.5/dy),
		  half_dx((T)0.5*dx),
		  pressure_array(pressure_field.array_for_this), pressure_density_array(projection_density_field.array_for_this), divergence_array(divergence_field.array_for_this),
		  boundary_condition_array(boundary_condition_field.array_for_this), max_velocity_x(0), max_velocity_y(0), CSF_model(false), Dirichlet_Boundary_Condition(false), Neumann_Boundary_Condition(false),
		  air_water_simulation(false), oil_water_simulation(false), vortex_sheet_problem(false), dimensionless_form(false)
	{}

	~PROJECTION_2D(void)
	{}

public: // Initialization Functions
	void InitializeFromBlock(const SCRIPT_BLOCK& projection_block)
	{
		tolerance = projection_block.GetFloat("tolerance", (T)1e-4);
		max_iteration = projection_block.GetInteger("max_iteration", 30);
		
		air_bubble_rising = projection_block.GetBoolean("air_bubble_rising", false);
		water_drop = projection_block.GetBoolean("water_drop", false);

		Dirichlet_Boundary_Condition = projection_block.GetBoolean("Dirichlet_Boundary_Condition", (bool)false);
		Neumann_Boundary_Condition = projection_block.GetBoolean("Neumann_Boundary_Condition", (bool)false);

		// Poisson solver type from script
		const char* poisson_solver_type_input = projection_block.GetString("poisson_solver_type", "Null");
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
		cout << "-------------------PROJECTION-------------------" << endl;
		cout << "tolerance: " << tolerance << endl;
		cout << "max_iteration: " << max_iteration << endl;
				
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

		if (air_bubble_rising)
		{
			cout << "Air Bubble Rising Simulation is activated!" << endl;
		}
		else if (water_drop)
		{
			cout << "Water Drop Simulation is activated!" << endl;
		}

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

		// Jump Condition field
		jc_on_solution.Initialize(base_grid, 1);
		jc_on_derivative.Initialize(base_grid, 1);
		jc_on_solution.AssignAllValue(0);
		jc_on_derivative.AssignAllValue(0);

		// Initialize fields
		pressure_field.Initialize(base_grid, 1);
		//pressure_field.AssignAllValue((T)0);
		boundary_condition_field.Initialize(base_grid, 1);
		divergence_field.Initialize(base_grid, 1);
		projection_density_field.Initialize(base_grid, 1);

		// Assigning initial field values
		projection_density_field.array_for_this.AssignAllValues((T)1);

		// Initialize Poisson solvers
		poisson_solver.Initialize(tolerance, max_iteration);
		poisson_solver.InitializeLinearSolver(poisson_solver_type);
	}

public: // Solver
	void Solve(const T& dt)
	{
		DetermineProjectionDensity();
		DetermineDivergence(dt);
		DetermineJumpConditionField(dt);
		DeterminePressure();
		UpdateVelocity(dt);
	}

	void DetermineProjectionDensity(void)
	{
		int i(0), j(0);
		if (air_water_simulation)
		{
			if (air_bubble_rising)
			{
				GRID_ITERATION_2D(water_levelset.grid)
				{
					if (water_levelset(i, j) <= 0)
					{
						projection_density_field(i, j) = air_density;
					}
					else if (water_levelset(i, j) > 0)
					{
						projection_density_field(i, j) = water_density;
					}
				}	
			}
			else if (water_drop)
			{
				GRID_ITERATION_2D(water_levelset.grid)
				{
					if (water_levelset(i, j) <= 0)
					{
						projection_density_field(i, j) = water_density;
					}
					else if (water_levelset(i, j) > 0)
					{
						projection_density_field(i, j) = air_density;
					}
				}	
			} 
		}
		if (oil_water_simulation)
		{
			GRID_ITERATION_2D(water_levelset.grid)
			{
				if (dimensionless_form)
				{
					if (water_levelset(i, j) <= 0)
					{
						projection_density_field(i, j) = water_density/oil_density;
					}
					else if (water_levelset(i, j) > 0)
					{
						projection_density_field(i, j) = 1;
					} 
				}
				else
				{
					if (water_levelset(i, j) <= 0)
					{
						projection_density_field(i, j) = water_density;
					}
					else if (water_levelset(i, j) > 0)
					{
						projection_density_field(i, j) = oil_density;
					} 
				}
			}
		}
	}
			
	void DetermineDivergence(const T& dt)	
	{
		int i(0), j(0);
		if (air_water_simulation)
		{
			LOOPS_2D(i, j, divergence_field.i_start, divergence_field.j_start, divergence_field.i_end, divergence_field.j_end)
			{
				divergence_field(i, j) = one_over_dx*(water_velocity_field_mac_x(i + 1, j) - water_velocity_field_mac_x(i, j)) + one_over_dy*(water_velocity_field_mac_y(i, j + 1) - water_velocity_field_mac_y(i, j));
			}
				
			// Scaled by time step
		
			const T one_over_dt = (T)1/dt;
        
			LOOPS_2D(i, j, divergence_field.i_start, divergence_field.j_start, divergence_field.i_end, divergence_field.j_end)
			{
				divergence_field(i, j) *= one_over_dt;
			}
		}
		if (oil_water_simulation)
		{
			LOOPS_2D(i, j, divergence_field.i_start, divergence_field.j_start, divergence_field.i_end, divergence_field.j_end)
			{
				T x_min = water_velocity_field_mac_x.grid.x_min, dx = water_velocity_field_mac_x.dx;
				
				divergence_field(i, j) = one_over_dx*((x_min + (i + 1)*dx)*water_velocity_field_mac_x(i + 1, j) - (x_min + (i - 1)*dx)*water_velocity_field_mac_x(i, j)) + one_over_dy*((x_min + (i + (T)0.5)*dx)*water_velocity_field_mac_y(i, j + 1) - (x_min + (i + (T)0.5)*dx)*water_velocity_field_mac_y(i, j));
			
			}
				
			// Scaled by time step
		
			const T one_over_dt = (T)1/dt;
        
			LOOPS_2D(i, j, divergence_field.i_start, divergence_field.j_start, divergence_field.i_end, divergence_field.j_end)
			{
				divergence_field(i, j) *= one_over_dt;
			}
		}
	}

	void DeterminePressure(void)
	{
		SetupBoundaryCondition(pressure_field, boundary_condition_field, water_levelset);
		
		water_levelset.ComputeNormals();
		
		projection_density_field.FillGhostCellsFrom(projection_density_field.array_for_this, true);

		int i(0), j(0);
		LOOPS_2D(i, j, jc_on_solution.i_start, jc_on_solution.j_start, jc_on_solution.i_end, jc_on_solution.j_end)
		{
			jc_on_derivative(i, j) = (T)0;
		}
		
		if (air_water_simulation)
		{
			poisson_solver.Solve(pressure_field, projection_density_field, boundary_condition_field, divergence_field, water_levelset, jc_on_solution, jc_on_derivative);
		}
		if (oil_water_simulation)
		{
			FIELD_STRUCTURE_2D<T> coef_1, coef_2, coef_3, coef_4;
			coef_1.Initialize(pressure_field.grid, 2);
			coef_2.Initialize(pressure_field.grid, 2);
			coef_3.Initialize(pressure_field.grid, 2);
			coef_4.Initialize(pressure_field.grid, 2);

			GRID_ITERATION_2D(coef_1.grid)
			{
				T x_min = coef_1.grid.x_min, dx = coef_1.grid.dx;
				coef_1(i, j) = (x_min + (i + (T)0.5)*dx)/((T)0.5*(projection_density_field(i, j) + projection_density_field(i + 1, j)));
				coef_2(i, j) = (x_min + (i - (T)0.5)*dx)/((T)0.5*(projection_density_field(i, j) + projection_density_field(i - 1, j)));
				coef_3(i, j) = (x_min + i*dx)/((T)0.5*(projection_density_field(i, j) + projection_density_field(i, j + 1)));
				coef_4(i, j) = (x_min + i*dx)/((T)0.5*(projection_density_field(i, j) + projection_density_field(i, j - 1)));
			}

			poisson_solver.SolveForAxisymmetric(pressure_field, coef_1, coef_2, coef_3, coef_4, boundary_condition_field, divergence_field);
		}
		if (vortex_sheet_problem)
		{
			FIELD_STRUCTURE_2D<T> delta;
			delta.Initialize(divergence_field.grid, 2);

			DeltaFunction(water_levelset.signed_distance_field, (T)12*delta.grid.dx, delta);
			jc_on_solution.AssignAllValue((T)0);

			FIELD_STRUCTURE_2D<T> righthand_side;
			righthand_side.Initialize(divergence_field.grid, 2);

			GRID_ITERATION_2D(divergence_field.grid)
			{
				righthand_side(i, j) = delta(i, j)*vortex_levelset(i, j);
			}

			righthand_side.FillGhostCellsFrom(righthand_side.array_for_this, true);

			poisson_solver.Solve(pressure_field, projection_density_field, boundary_condition_field, righthand_side, water_levelset, jc_on_solution, jc_on_derivative);
		}
	}

	void UpdateVelocity(const T& dt)
	{
		UpdateVelocityByPressureGradientVariableDensity(water_velocity_field_mac_x, dt);
		UpdateVelocityByPressureGradientVariableDensity(water_velocity_field_mac_y, dt);
	}

	void DetermineJumpConditionField(const T& dt)
	{
		int i(0), j(0);
		int i_start(jc_on_solution.i_start), i_end(jc_on_solution.i_end), j_start(jc_on_solution.j_start), j_end(jc_on_solution.j_end);

		water_levelset.ComputeCurvatures();
		
		if (use_delta_function_formulation == true)
		{
			if (CSF_model)
			{
				LOOPS_2D(i, j, i_start, j_start, i_end, j_end)
				{
					jc_on_solution(i, j) = (T)0;
				}
			}
			else
			{
				LOOPS_2D(i, j, i_start, j_start, i_end, j_end)
				{
					jc_on_solution(i, j) = surface_tension*water_levelset.curvature(i, j);
				}
			}
		}
		else if (use_jump_condition_on_viscosity == true)
		{
			T viscosity_p, viscosity_m;
			if (air_bubble_rising == true)
			{
				viscosity_p = water_viscosity;
				viscosity_m = air_viscosity;
			}
			else if (water_drop == true)
			{
				viscosity_p = air_viscosity;
				viscosity_m = water_viscosity;
			}

			LOOPS_2D(i, j, i_start, j_start, i_end, j_end)
			{
				T curvature = water_levelset.curvature(i, j);
				
				T viscosity_difference = viscosity_p - viscosity_m;
				T ux = (T)0.5*(water_velocity_field_mac_x(i + 1, j) + water_velocity_field_mac_x(i + 2, j) - water_velocity_field_mac_x(i - 1, j) - water_velocity_field_mac_x(i, j))*water_velocity_field_mac_x.one_over_2dx;
				T uy = (T)0.5*(water_velocity_field_mac_x(i, j + 1) + water_velocity_field_mac_x(i + 1, j + 1) - water_velocity_field_mac_x(i, j - 1) - water_velocity_field_mac_x(i + 1, j - 1))*water_velocity_field_mac_x.one_over_2dy;
				T vx = (T)0.5*(water_velocity_field_mac_y(i + 1, j) + water_velocity_field_mac_y(i + 1, j + 1) - water_velocity_field_mac_y(i - 1, j) - water_velocity_field_mac_y(i - 1, j + 1))*water_velocity_field_mac_y.one_over_2dx;
				T vy = (T)0.5*(water_velocity_field_mac_y(i, j + 2) + water_velocity_field_mac_y(i, j + 1) - water_velocity_field_mac_y(i, j - 1) - water_velocity_field_mac_y(i, j))*water_velocity_field_mac_y.one_over_2dy;

				T add_cond = (T)2*viscosity_difference*(ux*POW2(water_levelset.normal(i, j).x) + uy*water_levelset.normal(i, j).x*water_levelset.normal(i, j).y + vx*water_levelset.normal(i, j).x*water_levelset.normal(i, j).y + vy*POW2(water_levelset.normal(i, j).y));
				jc_on_solution(i, j) = surface_tension*curvature + add_cond;
			}
		}
	}

public: // Members Functions
	void SetProjectionType(const PROJECTION_TYPE& projection_type_input)
	{
		projection_type = projection_type_input;
	}

	void SetupBoundaryCondition(FIELD_STRUCTURE_2D<T>& pressure_input, FIELD_STRUCTURE_2D<int>& bc_input, const LEVELSET_2D& water_levelset_input)
	{
		ARRAY_2D<int>& bc_array(bc_input.array_for_this);
		GRID_STRUCTURE_2D& grid(bc_input.grid);

		int i(0), j(0);
		
		if (Dirichlet_Boundary_Condition)
		{
			LOOPS_2D(i, j, bc_input.i_start_g, bc_input.j_start_g, bc_input.i_end_g, bc_input.j_end_g)
			{
				// Speed-up variable
				if (i < grid.i_start || i > grid.i_end || j < grid.j_start || j > grid.j_end)
				{
					bc_array(i, j) = BC_DIR;
					if (vortex_sheet_problem)
					{
						if (j == grid.j_start - 1)
						{
							pressure_input(i, j) = 0;
						}
						if (j == grid.j_end + 1)
						{
							pressure_input(i, j) = 0;
						}
						if (i == grid.i_start - 1)
						{
							pressure_input(i, j) = pressure_input(grid.i_end - 1, j);
						}
						if (i == grid.i_end + 1)
						{
							pressure_input(i, j) = pressure_input(grid.i_start + 1, j);
						}
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
		}
		
		if (Neumann_Boundary_Condition)
		{
			LOOPS_2D(i, j, bc_input.i_start_g, bc_input.j_start_g, bc_input.i_end_g, bc_input.j_end_g)
			{
				// Speed-up variable
				if (i < grid.i_start || i > grid.i_end || j < grid.j_start || j > grid.j_end)
				{
					bc_array(i, j) = BC_NEUM;
					//pressure_input(i, j) = (T)0;
					/*if (i < grid.i_start)
					{
						pressure_input(i, j) = pressure_input(grid.i_start, j);
					}
					if (i > grid.i_end)
					{
						pressure_input(i, j) = pressure_input(grid.i_end, j);
					}
					if (j < grid.j_start)
					{
						pressure_input(i, j) = pressure_input(i, grid.j_start);
					}
					if (j > grid.j_end)
					{
						pressure_input(i, j) = pressure_input(i, grid.j_end);
					}*/
				}
				else
				{
					bc_array(i, j) = BC_FULL;
				}
			}
		}
	}

	void UpdateVelocityByPressureGradientVariableDensity(FIELD_STRUCTURE_2D<VT>& velocity_field, const T dt)
	{
		T pressure_plus, pressure_minus;

		GRID_ITERATION_2D(velocity_field.grid)
		{
			if (boundary_condition_array(i, j) < 0)
			{
				continue;
			}

			VT& velocity_ij = velocity_field.array_for_this(i, j);

			const T pressure_ij = pressure_array(i, j);
			const T density_ij = projection_density_field(i, j);
			const T one_over_density_ij = (T)1/density_ij;

			// Note: We correct dx for better object interaction
			for (int d = 0; d < 2; d++)
			{
				VI neighbor_index(i, j);

				// Positive neighbor in d-direction
				neighbor_index.values[d]++;

				const T inv_density_plus = (T)2/(density_ij + projection_density_field(neighbor_index));

				int neighbor_type = boundary_condition_array(neighbor_index);
				
				if (neighbor_type == BC_DIR)
				{
					pressure_plus = (T)0;
				}
				else
				{
					pressure_plus = pressure_array(neighbor_index);
				}

				// Negative neighbor in d-direction
				neighbor_index.values[d] -= 2;

				const T inv_density_minus = (T)2/(density_ij + projection_density_field(neighbor_index));

				neighbor_type = boundary_condition_array(neighbor_index);

				if (neighbor_type == BC_DIR)
				{
					pressure_minus = (T)0;
				}
				else
				{
					pressure_minus = pressure_array(neighbor_index);
				}

				velocity_ij.values[d] -= dt*((pressure_plus - pressure_ij)*inv_density_plus + (pressure_ij - pressure_minus)*inv_density_minus)*(T)0.5/dx;
				//velocity_ij.values[d] -= ((pressure_plus - pressure_minus)*one_over_density_ij)/dx;

			}	
		}	
	}
	
	void UpdateVelocityByPressureGradientVariableDensity(FIELD_STRUCTURE_2D<T>& velocity_field, const T dt)
	{
		if (CSF_model)
		{
			if (oil_water_simulation)
			{
				GRID_ITERATION_2D(velocity_field.grid)
				{
					if (boundary_condition_array(i, j) < 0)
					{
						continue;
					}

					T& velocity_ij = velocity_field.array_for_this(i, j);

					const T one_over_density_half_x = (T)1/((T)0.5*(projection_density_field(i, j) + projection_density_field(i - 1, j)));
					const T one_over_density_half_y = (T)1/((T)0.5*(projection_density_field(i, j) + projection_density_field(i, j - 1)));

					if (velocity_field.is_x_component == true)
					{
						velocity_ij -= dt*(pressure_field(i, j) - pressure_field(i - 1, j))*one_over_density_half_x*one_over_dx;

						max_velocity_x = MAX(abs(velocity_ij), max_velocity_x);

					}
					if (velocity_field.is_y_component == true)
					{
						velocity_ij -= dt*(pressure_field(i, j) - pressure_field(i, j - 1))*one_over_density_half_y*one_over_dy;

						max_velocity_y = MAX(abs(velocity_ij), max_velocity_y);
					}
				}	  
			}
		}
		else
		{
			// Set up the half density field
			FIELD_STRUCTURE_2D<T> levelset_x_c, levelset_y_c;
			levelset_x_c.Initialize(velocity_field.grid, 2);
			levelset_y_c.Initialize(velocity_field.grid, 2);

			GRID_ITERATION_2D(velocity_field.grid)
			{
				levelset_x_c(i, j) = (T)0.5*(water_levelset(i, j) + water_levelset(i - 1, j));
				levelset_y_c(i, j) = (T)0.5*(water_levelset(i, j) + water_levelset(i, j - 1));
			}

			FIELD_STRUCTURE_2D<T> density_half_x, density_half_y;
			density_half_x.Initialize(velocity_field.grid, 2);
			density_half_y.Initialize(velocity_field.grid, 2);

			// Sharp density Capturing
			DetermineDensityField(levelset_x_c, (T)1.5*base_grid.dx, density_half_x);
			DetermineDensityField(levelset_y_c, (T)1.5*base_grid.dy, density_half_y);

			GRID_ITERATION_2D(velocity_field.grid)
			{
				if (boundary_condition_array(i, j) < 0)
				{
					continue;
				}

				T& velocity_ij = velocity_field.array_for_this(i, j);

				const T density_ij = projection_density_field(i, j);
				const T one_over_density_ij = (T)1/density_ij;
				const T levelset_ij = water_levelset(i, j);
				const T levelset_ij_l = water_levelset(i - 1, j);
				const T levelset_ij_b = water_levelset(i, j - 1);
				const T one_over_density_half_x = (T)1/density_half_x(i, j);
				const T one_over_density_half_y = (T)1/density_half_y(i, j);
				const T jump_condition = jc_on_solution(i, j);
				const T jump_condition_l = jc_on_solution(i - 1, j);
				const T jump_condition_b = jc_on_solution(i, j - 1);

				const T jump_condition_gamma_x = (jump_condition_l*abs(levelset_ij) + jump_condition*abs(levelset_ij_l))/(abs(levelset_ij) + abs(levelset_ij_l));
				const T jump_condition_gamma_y = (jump_condition_b*abs(levelset_ij) + jump_condition*abs(levelset_ij_b))/(abs(levelset_ij) + abs(levelset_ij_b));

				if (velocity_field.is_x_component == true)
				{
					if (levelset_ij*levelset_ij_l >= 0)
					{
						velocity_ij -= dt*(pressure_field(i, j) - pressure_field(i - 1, j))*one_over_density_ij*one_over_dx;
					}
					else
					{
						if (levelset_ij < 0)
						{
							//velocity_ij -= dt*(pressure_field(i, j) - pressure_field(i - 1, j) + jump_condition)*one_over_density_half_x*one_over_dx;
							velocity_ij -= dt*(pressure_field(i, j) - pressure_field(i - 1, j) + jump_condition_gamma_x)*one_over_density_half_x*one_over_dx;
						}
						else if (levelset_ij > 0)
						{
							//velocity_ij -= dt*(pressure_field(i, j) - pressure_field(i - 1, j) - jump_condition)*one_over_density_half_x*one_over_dx;
							velocity_ij -= dt*(pressure_field(i, j) - pressure_field(i - 1, j) - jump_condition_gamma_x)*one_over_density_half_x*one_over_dx;
						}
					}
					max_velocity_x = MAX(abs(velocity_ij), max_velocity_x);

				}
				if (velocity_field.is_y_component == true)
				{
					if (levelset_ij*levelset_ij_b >= 0)
					{
						velocity_ij -= dt*(pressure_field(i, j) - pressure_field(i, j - 1))*one_over_density_ij*one_over_dy;
					}
					else
					{
						if (levelset_ij < 0)
						{
							//velocity_ij -= dt*(pressure_field(i, j) - pressure_field(i, j - 1) + jump_condition)*one_over_density_half_y*one_over_dy;
							velocity_ij -= dt*(pressure_field(i, j) - pressure_field(i, j - 1) + jump_condition_gamma_y)*one_over_density_half_y*one_over_dy;
						}
						else if (levelset_ij > 0)
						{
							//velocity_ij -= dt*(pressure_field(i, j) - pressure_field(i, j - 1) - jump_condition)*one_over_density_half_y*one_over_dy;
							velocity_ij -= dt*(pressure_field(i, j) - pressure_field(i, j - 1) - jump_condition_gamma_y)*one_over_density_half_y*one_over_dy;
						}
					}
					max_velocity_y = MAX(abs(velocity_ij), max_velocity_y);
				}
			}	 
		}
	}

	void HeavisideFunction(FIELD_STRUCTURE_2D<T>& phi, const T& epsilon, FIELD_STRUCTURE_2D<T>& heaviside)
	{
		int i(0), j(0);
		const int i_start(phi.i_start), j_start(phi.j_start), i_end(phi.i_end), j_end(phi.j_end);

		LOOPS_2D(i, j, i_start, j_start, i_end, j_end)
		{
			if (phi(i, j) < - epsilon)
			{
				heaviside(i, j) = 0;
			}
			else if (phi(i, j) > epsilon)
			{
				heaviside(i, j) = 1;
			}
			else
			{
				T one_over_epsilon = (T)1/epsilon, one_over_pi = (T)1/PI;
                heaviside(i, j) = (T)0.5 + phi(i, j)*one_over_epsilon*(T)0.5 + (T)0.5*one_over_pi*sin(PI*phi(i, j)*one_over_epsilon);
			}
		}
	}

	void DeltaFunction(FIELD_STRUCTURE_2D<T>& phi, const T& epsilon, FIELD_STRUCTURE_2D<T>& delta)
	{
		int i(0), j(0);
		const int i_start(phi.i_start), j_start(phi.j_start), i_end(phi.i_end), j_end(phi.j_end);

		LOOPS_2D(i, j, i_start, j_start, i_end, j_end)
		{
			if (abs(phi(i, j)) < epsilon)
			{
				T one_over_epsilon = (T)1/epsilon;
				delta(i, j) = (T)0.5*one_over_epsilon*((T)1 + cos(PI*phi(i, j)*one_over_epsilon));
			}
			else
			{
				delta(i, j) = (T)0;
			}
		}
	}

	void DetermineDensityField(FIELD_STRUCTURE_2D<T>& phi, const T& epsilon, FIELD_STRUCTURE_2D<T>& density)	
	{
		int i(0), j(0);
		const int i_start(density.i_start), i_end(density.i_end), j_start(density.j_start), j_end(density.j_end);
		
		FIELD_STRUCTURE_2D<T> heaviside_phi;
		heaviside_phi.Initialize(density.grid, 2);

		HeavisideFunction(phi, epsilon, heaviside_phi);

		LOOPS_2D(i, j, i_start, j_start, i_end, j_end)
		{
			if (air_bubble_rising == true)
			{
				density(i, j) = air_density + (water_density - air_density)*heaviside_phi(i, j);
			}
			else if (water_drop == true)
			{
				density(i, j) = water_density + (air_density - water_density)*heaviside_phi(i, j);
			}
		}	
		density.FillGhostCellsFrom(density.array_for_this, true);
	}
};

		  
		 
