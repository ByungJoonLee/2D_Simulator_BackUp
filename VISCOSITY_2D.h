#pragma once

#include "POISSON_SOLVER.h"
#include "LEVELSET_2D.h"

class VISCOSITY_2D
{
public: 
	MULTITHREADING&				multithreading;

public:
	// Water Levelset
	LEVELSET_2D&				water_levelset;
	FIELD_STRUCTURE_2D<T>&		water_signed_distance_field;

	GRID_STRUCTURE_2D&			base_grid;

	// For MAC Grid
	FIELD_STRUCTURE_2D<T>&		water_velocity_field_mac_x;
	FIELD_STRUCTURE_2D<T>&		water_velocity_field_mac_y;
		
	FIELD_STRUCTURE_2D<T>&		scalar_field_ghost;
	FIELD_STRUCTURE_2D<T>&		velocity_field_mac_ghost_x;
	FIELD_STRUCTURE_2D<T>&		velocity_field_mac_ghost_y;
	
	// Viscosity and Density
	T							water_viscosity, air_viscosity, oil_viscosity;
	T							water_density, air_density, oil_density;

	// Options for Bubble
	bool						air_bubble_rising, water_drop;
	
	// Options for Viscosity
	bool						use_delta_function_formulation;
	
    // Options For Simulation
    bool						air_water_simulation;
    bool						oil_water_simulation;
	bool						dimensionless_form;

public: // Sub Solver - For the semi-implicit method
	bool						semi_implicit_approach;

	enum POISSON_SOLVER_TYPE	poisson_solver_type;
	POISSON_SOLVER				poisson_solver;

	FIELD_STRUCTURE_2D<int>		boundary_condition_field_x;
	FIELD_STRUCTURE_2D<int>		boundary_condition_field_y;
	
	T							tolerance;
	int							max_iteration;

	
public: // Sub Variables for Delta Function approach
	// x-component 
	// Density
	FIELD_STRUCTURE_2D<T>		density_half_x;
	
	// Viscosity
	FIELD_STRUCTURE_2D<T>		viscosity_half_x;

	// Levelset
	FIELD_STRUCTURE_2D<T>		levelset_x_half;
	FIELD_STRUCTURE_2D<T>		levelset_cu_x;
	FIELD_STRUCTURE_2D<T>		levelset_cd_x;
	FIELD_STRUCTURE_2D<T>		levelset_r_x;
	FIELD_STRUCTURE_2D<T>		levelset_c_x;
	
	// Viscosity
	FIELD_STRUCTURE_2D<T>		viscosity_cu_x;
	FIELD_STRUCTURE_2D<T>		viscosity_cd_x;
	FIELD_STRUCTURE_2D<T>		viscosity_r_x;
	FIELD_STRUCTURE_2D<T>		viscosity_c_x;
	
	// y-component
	// Density
	FIELD_STRUCTURE_2D<T>		density_half_y;
	
	// Viscosity
	FIELD_STRUCTURE_2D<T>		viscosity_half_y;

	// Levelset
	FIELD_STRUCTURE_2D<T>		levelset_y_half;
	FIELD_STRUCTURE_2D<T>		levelset_cu_y;
	FIELD_STRUCTURE_2D<T>		levelset_cl_y;
	FIELD_STRUCTURE_2D<T>		levelset_u_y;
	FIELD_STRUCTURE_2D<T>		levelset_c_y;

	// Viscosity
	FIELD_STRUCTURE_2D<T>		viscosity_cu_y;
	FIELD_STRUCTURE_2D<T>		viscosity_cl_y;
	FIELD_STRUCTURE_2D<T>		viscosity_u_y;
	FIELD_STRUCTURE_2D<T>		viscosity_c_y;

	// For Semi-Implicit Method
	FIELD_STRUCTURE_2D<T>		explicit_term_x, temp_x, coef_1_x, coef_2_x, coef_3_x, coef_4_x, coef_5_x;
	FIELD_STRUCTURE_2D<T>		explicit_term_y, temp_y, coef_1_y, coef_2_y, coef_3_y, coef_4_y, coef_5_y;
	
	// Speedup Variable
	T							one_over_dx, one_over_dy;

public: // Dimensionless Variable
	T							Re;

public: // Constructor and Destructor
	VISCOSITY_2D(MULTITHREADING& multithreading_input, LEVELSET_2D& water_levelset_input, FIELD_STRUCTURE_2D<T>& water_velocity_field_mac_x_input, FIELD_STRUCTURE_2D<T>& water_velocity_field_mac_y_input,
				 FIELD_STRUCTURE_2D<T>& scalar_field_ghost_input, FIELD_STRUCTURE_2D<T>& velocity_field_mac_ghost_x_input, FIELD_STRUCTURE_2D<T>& velocity_field_mac_ghost_y_input)
				 : multithreading(multithreading_input), water_levelset(water_levelset_input), water_signed_distance_field(water_levelset_input.signed_distance_field), base_grid(water_levelset.grid), 
				 water_velocity_field_mac_x(water_velocity_field_mac_x_input), water_velocity_field_mac_y(water_velocity_field_mac_y_input),
				 scalar_field_ghost(scalar_field_ghost_input), velocity_field_mac_ghost_x(velocity_field_mac_ghost_x_input), velocity_field_mac_ghost_y(velocity_field_mac_ghost_y_input), 
                 air_density((T)0), water_density((T)0), oil_density((T)0), air_viscosity((T)0), water_viscosity((T)0), oil_viscosity((T)0), air_bubble_rising(false), water_drop(false), air_water_simulation(false), oil_water_simulation(false), semi_implicit_approach(false), dimensionless_form(false)
	{}

	~VISCOSITY_2D(void)
	{}

public: // Initialization Function
	void InitializeFromScriptBlock(const SCRIPT_BLOCK& script_block)
	{
		// Viscosity 
		if (air_water_simulation)
		{
            water_viscosity = script_block.GetFloat("water_viscosity", (T)1);
		    air_viscosity = script_block.GetFloat("air_viscosity", (T)1);
		}
		if (oil_water_simulation)
		{
            water_viscosity = script_block.GetFloat("water_viscosity", (T)1);
		    oil_viscosity = script_block.GetFloat("oil_viscosity", (T)1);
		}
		
		// Density
		if (air_water_simulation)
		{
            water_density = script_block.GetFloat("water_density", (T)1);
		    air_density = script_block.GetFloat("air_density", (T)1);
		}
		if (oil_water_simulation)
		{
            water_density = script_block.GetFloat("water_density", (T)1);
		    oil_density = script_block.GetFloat("oil_density", (T)1);
		}
		
		// Options for Bubble
		if (air_water_simulation)
		{
            air_bubble_rising = script_block.FindBlock("PROJECTION_WATER").GetBoolean("air_bubble_rising", (bool)false);
		    water_drop = script_block.FindBlock("PROJECTION_WATER").GetBoolean("water_drop", (bool)false);
		}
		
		// Options for Viscosity
		use_delta_function_formulation = script_block.GetBoolean("use_delta_function_formulation", (bool)false);
		
		// Semi-Implicit Approach
		semi_implicit_approach = script_block.GetBoolean("semi_implicit_approach", (bool)false);

		// Dimensionless Variable
		Re = script_block.GetFloat("Re", (T)1);

		// x-component 
		// Density
		density_half_x.Initialize(water_velocity_field_mac_x.grid, 2, &multithreading);

		// Viscosity
		viscosity_half_x.Initialize(water_velocity_field_mac_x.grid, 2, &multithreading);

		// Levelset
		levelset_x_half.Initialize(water_velocity_field_mac_x.grid, 2, &multithreading);
		levelset_cu_x.Initialize(water_velocity_field_mac_x.grid, 2, &multithreading);
		levelset_cd_x.Initialize(water_velocity_field_mac_x.grid, 2, &multithreading);
		levelset_r_x.Initialize(water_velocity_field_mac_x.grid, 2, &multithreading);
		levelset_c_x.Initialize(water_velocity_field_mac_x.grid, 2, &multithreading);

		// Viscosity
		viscosity_cu_x.Initialize(water_velocity_field_mac_x.grid, 2, &multithreading);
		viscosity_cd_x.Initialize(water_velocity_field_mac_x.grid, 2, &multithreading);
		viscosity_r_x.Initialize(water_velocity_field_mac_x.grid, 2, &multithreading);
		viscosity_c_x.Initialize(water_velocity_field_mac_x.grid, 2, &multithreading);

		// y-component
		// Density
		density_half_y.Initialize(water_velocity_field_mac_y.grid, 2, &multithreading);

		// Viscosity
		viscosity_half_y.Initialize(water_velocity_field_mac_y.grid, 2, &multithreading);

		// Levelset
		levelset_y_half.Initialize(water_velocity_field_mac_y.grid, 2, &multithreading);
		levelset_cu_y.Initialize(water_velocity_field_mac_y.grid, 2, &multithreading);
		levelset_cl_y.Initialize(water_velocity_field_mac_y.grid, 2, &multithreading);
		levelset_u_y.Initialize(water_velocity_field_mac_y.grid, 2, &multithreading);
		levelset_c_y.Initialize(water_velocity_field_mac_y.grid, 2, &multithreading);

		// Viscosity
		viscosity_cu_y.Initialize(water_velocity_field_mac_y.grid, 2, &multithreading);
		viscosity_cl_y.Initialize(water_velocity_field_mac_y.grid, 2, &multithreading);
		viscosity_u_y.Initialize(water_velocity_field_mac_y.grid, 2, &multithreading);
		viscosity_c_y.Initialize(water_velocity_field_mac_y.grid, 2, &multithreading);

		// Speedup variables
		one_over_dx = base_grid.one_over_dx;
		one_over_dy = base_grid.one_over_dy;

		// For the semi-implicit solver
		tolerance = script_block.FindBlock("PROJECTION").GetFloat("tolerance", (T)1e-4);
		max_iteration = script_block.FindBlock("PROJECTION").GetInteger("max_iteration", 30);

		const char* poisson_solver_type_input = script_block.FindBlock("PROJECTION").GetString("poisson_solver_type", "NULL");
		if (!strcmp(poisson_solver_type_input, "CG"))
		{
			poisson_solver_type = CG;
		}
		else if (!strcmp(poisson_solver_type_input, "PCG"))
		{
			poisson_solver_type = PCG;
		}
		else
		{
			poisson_solver_type = CG;
		}

		// Initialize the Poisson Solver
		poisson_solver.Initialize(tolerance, max_iteration, 0, &multithreading);
		poisson_solver.InitializeLinearSolver(poisson_solver_type);
		
		boundary_condition_field_x.Initialize(water_velocity_field_mac_x.grid, 1, &multithreading);
		boundary_condition_field_y.Initialize(water_velocity_field_mac_y.grid, 1, &multithreading);

		// Initialize the Semi-Implicit solver
		// x-component
		explicit_term_x.Initialize(water_velocity_field_mac_x.grid, 2, &multithreading);
		temp_x.Initialize(water_velocity_field_mac_x.grid, 2, &multithreading);
		coef_1_x.Initialize(water_velocity_field_mac_x.grid, 2, &multithreading);
		coef_2_x.Initialize(water_velocity_field_mac_x.grid, 2, &multithreading);
		coef_3_x.Initialize(water_velocity_field_mac_x.grid, 2, &multithreading);
		coef_4_x.Initialize(water_velocity_field_mac_x.grid, 2, &multithreading);
		coef_5_x.Initialize(water_velocity_field_mac_x.grid, 2, &multithreading);

		// y-component
		explicit_term_y.Initialize(water_velocity_field_mac_y.grid, 2, &multithreading);
		temp_y.Initialize(water_velocity_field_mac_y.grid, 2, &multithreading);
		coef_1_y.Initialize(water_velocity_field_mac_y.grid, 2, &multithreading);
		coef_2_y.Initialize(water_velocity_field_mac_y.grid, 2, &multithreading);
		coef_3_y.Initialize(water_velocity_field_mac_y.grid, 2, &multithreading);
		coef_4_y.Initialize(water_velocity_field_mac_y.grid, 2, &multithreading);
		coef_5_y.Initialize(water_velocity_field_mac_y.grid, 2, &multithreading);

		if (semi_implicit_approach)
		{
			cout << "--------------SEMI-IMPLICIT--------------" << endl;
			cout << "tolerance: " << tolerance << endl;
			cout << "max iteration: " << max_iteration << endl;
		}
		else
		{
			cout << "Semi-implicit for Viscosity: false" << endl;
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
	}

public: // Main Function
	void ApplyViscosity(const T& epsilon_for_mollification, const T& dt)
	{
		if (use_delta_function_formulation == true)
		{
			if (air_water_simulation)
			{
				// Speed Up Variable
				int i(0), j(0);
				const T one_over_dx = base_grid.one_over_dx, one_over_dy = base_grid.one_over_dy;

				const ARRAY_2D<T>& water_array(water_levelset.arr);
				const ARRAY_2D<T>& velocity_x(water_velocity_field_mac_x.array_for_this);
				const ARRAY_2D<T>& velocity_y(water_velocity_field_mac_y.array_for_this);
				
				FIELD_STRUCTURE_2D<T> temp_x, temp_y;
				temp_x.Initialize(water_velocity_field_mac_x.grid, 2);
				temp_y.Initialize(water_velocity_field_mac_y.grid, 2);

				GRID_ITERATION_2D(water_velocity_field_mac_x.grid)
				{
					temp_x(i, j) = velocity_x(i, j);
				}

				GRID_ITERATION_2D(water_velocity_field_mac_y.grid)
				{
					temp_y(i, j) = velocity_y(i, j);
				}

				GRID_ITERATION_2D(water_velocity_field_mac_x.grid)
				{
					levelset_x_half(i, j) = (T)0.5*(water_array(i - 1, j) + water_array(i, j));
					levelset_cu_x(i, j) = (T)0.25*(water_array(i - 1, j) + water_array(i, j) + water_array(i - 1, j + 1) + water_array(i, j + 1));
					levelset_cd_x(i, j) = (T)0.25*(water_array(i - 1, j - 1) + water_array(i, j - 1) + water_array(i - 1, j) + water_array(i, j));
					levelset_r_x(i, j) = water_array(i, j);
					levelset_c_x(i, j) = water_array(i - 1, j);
				}

				DetermineDensityField(levelset_x_half, epsilon_for_mollification, density_half_x);
				DetermineViscosityField(levelset_cu_x, epsilon_for_mollification, viscosity_cu_x);
				DetermineViscosityField(levelset_cd_x, epsilon_for_mollification, viscosity_cd_x);
				DetermineViscosityField(levelset_r_x, epsilon_for_mollification, viscosity_r_x);
				DetermineViscosityField(levelset_c_x, epsilon_for_mollification, viscosity_c_x);

				GRID_ITERATION_2D(water_velocity_field_mac_y.grid)
				{
					levelset_y_half(i, j) = (T)0.5*(water_array(i, j - 1) + water_array(i, j));
					levelset_cu_y(i, j) = (T)0.25*(water_array(i, j - 1) + water_array(i + 1, j - 1) + water_array(i, j) + water_array(i + 1, j));
					levelset_cl_y(i, j) = (T)0.25*(water_array(i - 1, j - 1) + water_array(i, j - 1) + water_array(i - 1, j) + water_array(i, j));
					levelset_u_y(i, j) = water_array(i, j);
					levelset_c_y(i, j) = water_array(i, j - 1);
				}

				DetermineDensityField(levelset_y_half, epsilon_for_mollification, density_half_y);
				DetermineViscosityField(levelset_cu_y, epsilon_for_mollification, viscosity_cu_y);
				DetermineViscosityField(levelset_cl_y, epsilon_for_mollification, viscosity_cl_y);
				DetermineViscosityField(levelset_u_y, epsilon_for_mollification, viscosity_u_y);
				DetermineViscosityField(levelset_c_y, epsilon_for_mollification, viscosity_c_y);

				GRID_ITERATION_2D(water_velocity_field_mac_x.grid)
				{
					T one_over_density_half_x = (T)1/density_half_x(i, j);
					T coef = dt*one_over_density_half_x;
					T first_update = (T)2*coef*one_over_dx*(viscosity_r_x(i, j)*one_over_dx*(temp_x(i + 1, j) - temp_x(i, j)) - viscosity_c_x(i, j)*one_over_dx*(temp_x(i, j) - temp_x(i - 1, j)));
					T second_update = coef*one_over_dy*(viscosity_cu_x(i, j)*(one_over_dy*(temp_x(i, j + 1) - temp_x(i, j)) + one_over_dx*(temp_y(i, j + 1) - temp_y(i - 1, j + 1))) - viscosity_cd_x(i, j)*(one_over_dy*(temp_x(i, j) - temp_x(i, j - 1)) + one_over_dx*(temp_y(i, j) - temp_y(i - 1, j))));
					velocity_x(i, j) += first_update + second_update;
				}

				GRID_ITERATION_2D(water_velocity_field_mac_y.grid)
				{
					T one_over_density_half_y = (T)1/density_half_y(i, j);
					T coef = dt*one_over_density_half_y;
					T first_update = (T)2*coef*one_over_dy*(viscosity_u_y(i, j)*one_over_dy*(temp_y(i, j + 1) - temp_y(i, j)) - viscosity_c_y(i, j)*one_over_dy*(temp_y(i, j) - temp_y(i, j - 1)));
					T second_update = coef*one_over_dx*(viscosity_cu_y(i, j)*(one_over_dy*(temp_x(i + 1, j) - temp_x(i + 1, j - 1)) + one_over_dx*(temp_y(i + 1, j) - temp_y(i, j))) - viscosity_cl_y(i, j)*(one_over_dy*(temp_x(i, j) - temp_x(i, j - 1)) + one_over_dx*(temp_y(i, j) - temp_y(i - 1, j))));
					velocity_y(i, j) += first_update + second_update;
				} 
			}
		}
		if (oil_water_simulation)
		{
			if (semi_implicit_approach)
			{
				// Speed Up Variable
				const ARRAY_2D<T>& water_array(water_levelset.arr);
			
				const ARRAY_2D<T>& velocity_x(water_velocity_field_mac_x.array_for_this);
				const ARRAY_2D<T>& velocity_y(water_velocity_field_mac_y.array_for_this);
					
				if (dimensionless_form)
				{
					// Define the coefficient function -- coefficients of x-velocity
					GRID_ITERATION_2D(water_velocity_field_mac_x.grid)
					{
						levelset_x_half(i, j) = (T)0.5*(water_array(i - 1, j) + water_array(i, j));
						levelset_cu_x(i, j) = (T)0.25*(water_array(i - 1, j) + water_array(i, j) + water_array(i - 1, j + 1) + water_array(i, j + 1));
						levelset_cd_x(i, j) = (T)0.25*(water_array(i - 1, j) + water_array(i, j) + water_array(i - 1, j - 1) + water_array(i, j - 1));
						levelset_r_x(i, j) = water_array(i, j);
						levelset_c_x(i, j) = water_array(i - 1, j);
					}
				
					DetermineDensityField(levelset_x_half, epsilon_for_mollification, density_half_x);
					DetermineViscosityField(levelset_x_half, epsilon_for_mollification, viscosity_half_x);
					DetermineViscosityField(levelset_cu_x, epsilon_for_mollification, viscosity_cu_x);
					DetermineViscosityField(levelset_cd_x, epsilon_for_mollification, viscosity_cd_x);
					DetermineViscosityField(levelset_r_x, epsilon_for_mollification, viscosity_r_x);
					DetermineViscosityField(levelset_c_x, epsilon_for_mollification, viscosity_c_x);
					
					// x-component
					GRID_ITERATION_2D(temp_x.grid)
					{
						temp_x(i, j) = water_velocity_field_mac_x(i, j);
					}

					GRID_ITERATION_2D(coef_1_x.grid)
					{
						T one_over_r = 1/(coef_1_x.grid.x_min + (i + (T)0.5)*coef_1_x.grid.dx);
						T one_over_r2 = POW2(one_over_r);
						coef_1_x(i, j) = one_over_r*(coef_1_x.grid.x_min + (i + 1)*coef_1_x.grid.dx)*(T)2*dt/(Re*density_half_x(i, j))*viscosity_r_x(i, j);	
						coef_2_x(i, j) = one_over_r*(coef_1_x.grid.x_min + i*coef_1_x.grid.dx)*(T)2*dt/(Re*density_half_x(i, j))*viscosity_c_x(i, j);
						coef_3_x(i, j) = dt/(Re*density_half_x(i, j))*viscosity_cu_x(i, j);
						coef_4_x(i, j) = dt/(Re*density_half_x(i, j))*viscosity_cd_x(i, j);
						coef_5_x(i, j) = (T)2*dt/(Re*density_half_x(i, j))*viscosity_half_x(i, j)*one_over_r2;
					}

					GRID_ITERATION_2D(explicit_term_x.grid)
					{
						explicit_term_x(i, j) = temp_x(i, j) + dt/(Re*density_half_x(i, j))*(one_over_dy*(viscosity_cu_x(i, j)*one_over_dx*(velocity_y(i, j + 1) - velocity_y(i - 1, j + 1)) - viscosity_cd_x(i, j)*one_over_dx*(velocity_y(i, j) - velocity_y(i - 1, j))));
					}

					SetupBoundaryConditionsForVelocity(water_velocity_field_mac_x, boundary_condition_field_x);

					poisson_solver.SolveForViscosity(water_velocity_field_mac_x, coef_1_x, coef_2_x, coef_3_x, coef_4_x, coef_5_x, boundary_condition_field_x, explicit_term_x);
					
					// Define for coefficient function -- coefficients of y-velocity
					GRID_ITERATION_2D(water_velocity_field_mac_y.grid)
					{
						levelset_y_half(i, j) = (T)0.5*(water_array(i, j - 1) + water_array(i, j));
						levelset_cu_y(i, j) = (T)0.25*(water_array(i, j - 1) + water_array(i + 1, j - 1) + water_array(i, j) + water_array(i + 1, j));
						levelset_cl_y(i, j) = (T)0.25*(water_array(i - 1, j - 1) + water_array(i, j - 1) + water_array(i - 1, j) + water_array(i, j));
						levelset_u_y(i, j) = water_array(i, j);
						levelset_c_y(i, j) = water_array(i, j - 1);
					}
				
					DetermineDensityField(levelset_y_half, epsilon_for_mollification, density_half_y);
					DetermineViscosityField(levelset_y_half, epsilon_for_mollification, viscosity_half_y);
					DetermineViscosityField(levelset_cu_y, epsilon_for_mollification, viscosity_cu_y);
					DetermineViscosityField(levelset_cl_y, epsilon_for_mollification, viscosity_cl_y);
					DetermineViscosityField(levelset_u_y, epsilon_for_mollification, viscosity_u_y);
					DetermineViscosityField(levelset_c_y, epsilon_for_mollification, viscosity_c_y);
										
					// y-component
					GRID_ITERATION_2D(temp_y.grid)
					{
						temp_y(i, j) = water_velocity_field_mac_y(i, j);
					}
	
					GRID_ITERATION_2D(coef_1_y.grid)
					{
						T r_coor = coef_1_y.x_min + i*coef_1_y.dx;
						T one_over_r = 1/r_coor;
						coef_1_y(i, j) = one_over_r*dt/(Re*density_half_y(i, j))*viscosity_cu_y(i, j)*(r_coor + (T)0.5*coef_1_y.dx);
						coef_2_y(i, j) = one_over_r*dt/(Re*density_half_y(i, j))*viscosity_cl_y(i, j)*(r_coor - (T)0.5*coef_1_y.dx);
						coef_3_y(i, j) = (T)2*dt/(Re*density_half_y(i, j))*viscosity_u_y(i, j);
						coef_4_y(i, j) = (T)2*dt/(Re*density_half_y(i, j))*viscosity_c_y(i, j);
						coef_5_y(i, j) = (T)0;
					}
						
					GRID_ITERATION_2D(explicit_term_y.grid)
					{
						T r_coor = coef_1_y.x_min + i*coef_1_y.dx;
						T one_over_r = 1/r_coor;
						explicit_term_y(i, j) = temp_y(i, j) + dt/(Re*density_half_y(i, j))*((one_over_dx*(viscosity_cu_y(i, j)*one_over_dy*(temp_x(i + 1, j) - temp_x(i + 1, j - 1)) - viscosity_cl_y(i, j)*one_over_dy*(temp_x(i, j) - temp_x(i, j - 1)))) + one_over_r*viscosity_half_y(i, j)*one_over_dy*((T)0.5*(water_velocity_field_mac_x(i, j) + water_velocity_field_mac_x(i + 1, j)) - (T)0.5*(water_velocity_field_mac_x(i, j - 1) + water_velocity_field_mac_x(i + 1, j - 1))));
					}
	
					SetupBoundaryConditionsForVelocity(water_velocity_field_mac_y, boundary_condition_field_y);
						
					poisson_solver.SolveForViscosity(water_velocity_field_mac_y, coef_1_y, coef_2_y, coef_3_y, coef_4_y, coef_5_y, boundary_condition_field_y, explicit_term_y);
				}
			}
		}
	}

	void ApplyViscosity(const T& epsilon_for_mollification, const T& dt, const int& thread_id)
	{
		if (use_delta_function_formulation)
		{
			if (air_water_simulation)
			{
				// Speed Up Variable
				const T one_over_dx = base_grid.one_over_dx, one_over_dy = base_grid.one_over_dy;

				const ARRAY_2D<T>& water_array(water_levelset.arr);
				const ARRAY_2D<T>& velocity_x(water_velocity_field_mac_x.array_for_this);
				const ARRAY_2D<T>& velocity_y(water_velocity_field_mac_y.array_for_this);
				
				FIELD_STRUCTURE_2D<T> temp_x, temp_y;
				temp_x.Initialize(water_velocity_field_mac_x.grid, 2, &multithreading);
				temp_y.Initialize(water_velocity_field_mac_y.grid, 2, &multithreading);

				GRID_ITERATION_2D(water_velocity_field_mac_x.partial_grids[thread_id])
				{
					temp_x(i, j) = velocity_x(i, j);
				}
				multithreading.Sync(thread_id);

				GRID_ITERATION_2D(water_velocity_field_mac_y.partial_grids[thread_id])
				{
					temp_y(i, j) = velocity_y(i, j);
				}
				multithreading.Sync(thread_id);

				GRID_ITERATION_2D(water_velocity_field_mac_x.partial_grids[thread_id])
				{
					levelset_x_half(i, j) = (T)0.5*(water_array(i - 1, j) + water_array(i, j));
					levelset_cu_x(i, j) = (T)0.25*(water_array(i - 1, j) + water_array(i, j) + water_array(i - 1, j + 1) + water_array(i, j + 1));
					levelset_cd_x(i, j) = (T)0.25*(water_array(i - 1, j - 1) + water_array(i, j - 1) + water_array(i - 1, j) + water_array(i, j));
					levelset_r_x(i, j) = water_array(i, j);
					levelset_c_x(i, j) = water_array(i - 1, j);
				}
				multithreading.Sync(thread_id);

				DetermineDensityField(levelset_x_half, epsilon_for_mollification, density_half_x, thread_id);
				DetermineViscosityField(levelset_cu_x, epsilon_for_mollification, viscosity_cu_x, thread_id);
				DetermineViscosityField(levelset_cd_x, epsilon_for_mollification, viscosity_cd_x, thread_id);
				DetermineViscosityField(levelset_r_x, epsilon_for_mollification, viscosity_r_x, thread_id);
				DetermineViscosityField(levelset_c_x, epsilon_for_mollification, viscosity_c_x, thread_id);

				GRID_ITERATION_2D(water_velocity_field_mac_y.partial_grids[thread_id])
				{
					levelset_y_half(i, j) = (T)0.5*(water_array(i, j - 1) + water_array(i, j));
					levelset_cu_y(i, j) = (T)0.25*(water_array(i, j - 1) + water_array(i + 1, j - 1) + water_array(i, j) + water_array(i + 1, j));
					levelset_cl_y(i, j) = (T)0.25*(water_array(i - 1, j - 1) + water_array(i, j - 1) + water_array(i - 1, j) + water_array(i, j));
					levelset_u_y(i, j) = water_array(i, j);
					levelset_c_y(i, j) = water_array(i, j - 1);
				}
				multithreading.Sync(thread_id);

				DetermineDensityField(levelset_y_half, epsilon_for_mollification, density_half_y, thread_id);
				DetermineViscosityField(levelset_cu_y, epsilon_for_mollification, viscosity_cu_y, thread_id);
				DetermineViscosityField(levelset_cl_y, epsilon_for_mollification, viscosity_cl_y, thread_id);
				DetermineViscosityField(levelset_u_y, epsilon_for_mollification, viscosity_u_y, thread_id);
				DetermineViscosityField(levelset_c_y, epsilon_for_mollification, viscosity_c_y, thread_id);

				GRID_ITERATION_2D(water_velocity_field_mac_x.partial_grids[thread_id])
				{
					T one_over_density_half_x = (T)1/density_half_x(i, j);
					T coef = dt*one_over_density_half_x;
					T first_update = (T)2*coef*one_over_dx*(viscosity_r_x(i, j)*one_over_dx*(temp_x(i + 1, j) - temp_x(i, j)) - viscosity_c_x(i, j)*one_over_dx*(temp_x(i, j) - temp_x(i - 1, j)));
					T second_update = coef*one_over_dy*(viscosity_cu_x(i, j)*(one_over_dy*(temp_x(i, j + 1) - temp_x(i, j)) + one_over_dx*(temp_y(i, j + 1) - temp_y(i - 1, j + 1))) - viscosity_cd_x(i, j)*(one_over_dy*(temp_x(i, j) - temp_x(i, j - 1)) + one_over_dx*(temp_y(i, j) - temp_y(i - 1, j))));
					velocity_x(i, j) += first_update + second_update;
				}
				multithreading.Sync(thread_id);

				GRID_ITERATION_2D(water_velocity_field_mac_y.partial_grids[thread_id])
				{
					T one_over_density_half_y = (T)1/density_half_y(i, j);
					T coef = dt*one_over_density_half_y;
					T first_update = (T)2*coef*one_over_dy*(viscosity_u_y(i, j)*one_over_dy*(temp_y(i, j + 1) - temp_y(i, j)) - viscosity_c_y(i, j)*one_over_dy*(temp_y(i, j) - temp_y(i, j - 1)));
					T second_update = coef*one_over_dx*(viscosity_cu_y(i, j)*(one_over_dy*(temp_x(i + 1, j) - temp_x(i + 1, j - 1)) + one_over_dx*(temp_y(i + 1, j) - temp_y(i, j))) - viscosity_cl_y(i, j)*(one_over_dy*(temp_x(i, j) - temp_x(i, j - 1)) + one_over_dx*(temp_y(i, j) - temp_y(i - 1, j))));
					velocity_y(i, j) += first_update + second_update;
				} 
				multithreading.Sync(thread_id);
			}
		}
		if (oil_water_simulation)
		{
			if (semi_implicit_approach)
			{
				// Speed Up Variable
				const ARRAY_2D<T>& water_array(water_levelset.arr);
			
				const ARRAY_2D<T>& velocity_x(water_velocity_field_mac_x.array_for_this);
				const ARRAY_2D<T>& velocity_y(water_velocity_field_mac_y.array_for_this);
					
				if (dimensionless_form)
				{
					// Define the coefficient function -- coefficients of x-velocity
					GRID_ITERATION_2D(water_velocity_field_mac_x.grid)
					{
						levelset_x_half(i, j) = (T)0.5*(water_array(i - 1, j) + water_array(i, j));
						levelset_cu_x(i, j) = (T)0.25*(water_array(i - 1, j) + water_array(i, j) + water_array(i - 1, j + 1) + water_array(i, j + 1));
						levelset_cd_x(i, j) = (T)0.25*(water_array(i - 1, j) + water_array(i, j) + water_array(i - 1, j - 1) + water_array(i, j - 1));
						levelset_r_x(i, j) = water_array(i, j);
						levelset_c_x(i, j) = water_array(i - 1, j);
					}
				
					DetermineDensityField(levelset_x_half, epsilon_for_mollification, density_half_x);
					DetermineViscosityField(levelset_x_half, epsilon_for_mollification, viscosity_half_x);
					DetermineViscosityField(levelset_cu_x, epsilon_for_mollification, viscosity_cu_x);
					DetermineViscosityField(levelset_cd_x, epsilon_for_mollification, viscosity_cd_x);
					DetermineViscosityField(levelset_r_x, epsilon_for_mollification, viscosity_r_x);
					DetermineViscosityField(levelset_c_x, epsilon_for_mollification, viscosity_c_x);
					
					// x-component
					GRID_ITERATION_2D(temp_x.grid)
					{
						temp_x(i, j) = water_velocity_field_mac_x(i, j);
					}

					GRID_ITERATION_2D(coef_1_x.grid)
					{
						T one_over_r = 1/(coef_1_x.grid.x_min + (i + (T)0.5)*coef_1_x.grid.dx);
						T one_over_r2 = POW2(one_over_r);
						coef_1_x(i, j) = one_over_r*(coef_1_x.grid.x_min + (i + 1)*coef_1_x.grid.dx)*(T)2*dt/(Re*density_half_x(i, j))*viscosity_r_x(i, j);	
						coef_2_x(i, j) = one_over_r*(coef_1_x.grid.x_min + i*coef_1_x.grid.dx)*(T)2*dt/(Re*density_half_x(i, j))*viscosity_c_x(i, j);
						coef_3_x(i, j) = dt/(Re*density_half_x(i, j))*viscosity_cu_x(i, j);
						coef_4_x(i, j) = dt/(Re*density_half_x(i, j))*viscosity_cd_x(i, j);
						coef_5_x(i, j) = (T)2*dt/(Re*density_half_x(i, j))*viscosity_half_x(i, j)*one_over_r2;
					}

					GRID_ITERATION_2D(explicit_term_x.grid)
					{
						explicit_term_x(i, j) = temp_x(i, j) + dt/(Re*density_half_x(i, j))*(one_over_dy*(viscosity_cu_x(i, j)*one_over_dx*(velocity_y(i, j + 1) - velocity_y(i - 1, j + 1)) - viscosity_cd_x(i, j)*one_over_dx*(velocity_y(i, j) - velocity_y(i - 1, j))));
					}

					SetupBoundaryConditionsForVelocity(water_velocity_field_mac_x, boundary_condition_field_x);

					poisson_solver.SolveForViscosity(water_velocity_field_mac_x, coef_1_x, coef_2_x, coef_3_x, coef_4_x, coef_5_x, boundary_condition_field_x, explicit_term_x);
					
					// Define for coefficient function -- coefficients of y-velocity
					GRID_ITERATION_2D(water_velocity_field_mac_y.grid)
					{
						levelset_y_half(i, j) = (T)0.5*(water_array(i, j - 1) + water_array(i, j));
						levelset_cu_y(i, j) = (T)0.25*(water_array(i, j - 1) + water_array(i + 1, j - 1) + water_array(i, j) + water_array(i + 1, j));
						levelset_cl_y(i, j) = (T)0.25*(water_array(i - 1, j - 1) + water_array(i, j - 1) + water_array(i - 1, j) + water_array(i, j));
						levelset_u_y(i, j) = water_array(i, j);
						levelset_c_y(i, j) = water_array(i, j - 1);
					}
				
					DetermineDensityField(levelset_y_half, epsilon_for_mollification, density_half_y);
					DetermineViscosityField(levelset_y_half, epsilon_for_mollification, viscosity_half_y);
					DetermineViscosityField(levelset_cu_y, epsilon_for_mollification, viscosity_cu_y);
					DetermineViscosityField(levelset_cl_y, epsilon_for_mollification, viscosity_cl_y);
					DetermineViscosityField(levelset_u_y, epsilon_for_mollification, viscosity_u_y);
					DetermineViscosityField(levelset_c_y, epsilon_for_mollification, viscosity_c_y);
										
					// y-component
					GRID_ITERATION_2D(temp_y.grid)
					{
						temp_y(i, j) = water_velocity_field_mac_y(i, j);
					}
	
					GRID_ITERATION_2D(coef_1_y.grid)
					{
						T r_coor = coef_1_y.x_min + i*coef_1_y.dx;
						T one_over_r = 1/r_coor;
						coef_1_y(i, j) = one_over_r*dt/(Re*density_half_y(i, j))*viscosity_cu_y(i, j)*(r_coor + (T)0.5*coef_1_y.dx);
						coef_2_y(i, j) = one_over_r*dt/(Re*density_half_y(i, j))*viscosity_cl_y(i, j)*(r_coor - (T)0.5*coef_1_y.dx);
						coef_3_y(i, j) = (T)2*dt/(Re*density_half_y(i, j))*viscosity_u_y(i, j);
						coef_4_y(i, j) = (T)2*dt/(Re*density_half_y(i, j))*viscosity_c_y(i, j);
						coef_5_y(i, j) = (T)0;
					}
						
					GRID_ITERATION_2D(explicit_term_y.grid)
					{
						T r_coor = coef_1_y.x_min + i*coef_1_y.dx;
						T one_over_r = 1/r_coor;
						explicit_term_y(i, j) = temp_y(i, j) + dt/(Re*density_half_y(i, j))*((one_over_dx*(viscosity_cu_y(i, j)*one_over_dy*(temp_x(i + 1, j) - temp_x(i + 1, j - 1)) - viscosity_cl_y(i, j)*one_over_dy*(temp_x(i, j) - temp_x(i, j - 1)))) + one_over_r*viscosity_half_y(i, j)*one_over_dy*((T)0.5*(water_velocity_field_mac_x(i, j) + water_velocity_field_mac_x(i + 1, j)) - (T)0.5*(water_velocity_field_mac_x(i, j - 1) + water_velocity_field_mac_x(i + 1, j - 1))));
					}
	
					SetupBoundaryConditionsForVelocity(water_velocity_field_mac_y, boundary_condition_field_y);
						
					poisson_solver.SolveForViscosity(water_velocity_field_mac_y, coef_1_y, coef_2_y, coef_3_y, coef_4_y, coef_5_y, boundary_condition_field_y, explicit_term_y);
				}
			}
		}
	}

public: // Member Function
	void SetupBoundaryConditionsForVelocity(FIELD_STRUCTURE_2D<T>& velocity_input, FIELD_STRUCTURE_2D<int>& bc_input)
	{
		ARRAY_2D<int>& bc_array(bc_input.array_for_this);
		GRID_STRUCTURE_2D& grid(bc_input.grid);

		GRID_ITERATION_2D(bc_input.grid)
		{
			if (i < grid.i_start || i > grid.i_end || j < grid.j_start || j > grid.j_end)
			{
				bc_array(i, j) = BC_DIR;
				velocity_input(i, j) = 0;
			}
			else
			{
				bc_array(i, j) = BC_FULL;
			}
		}
	}

public: // Helpful Functions
	void HeavisideFunction(FIELD_STRUCTURE_2D<T>& phi, const T& epsilon, FIELD_STRUCTURE_2D<T>& heaviside)
	{
		T one_over_epsilon = (T)1/epsilon, one_over_pi = (T)1/PI;

		GRID_ITERATION_2D(phi.grid)
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
				heaviside(i, j) = (T)0.5 + phi(i, j)*one_over_epsilon*(T)0.5 + (T)0.5*one_over_pi*sin(PI*phi(i, j)*one_over_epsilon);
			}
		}
	}
	
	void HeavisideFunction(FIELD_STRUCTURE_2D<T>& phi, const T& epsilon, FIELD_STRUCTURE_2D<T>& heaviside, const int& thread_id)
	{
		T one_over_epsilon = (T)1/epsilon, one_over_pi = (T)1/PI;

		GRID_ITERATION_2D(phi.partial_grids[thread_id])
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
				heaviside(i, j) = (T)0.5 + phi(i, j)*one_over_epsilon*(T)0.5 + (T)0.5*one_over_pi*sin(PI*phi(i, j)*one_over_epsilon);
			}
		}
		multithreading.Sync(thread_id);
	}

	void DetermineViscosityField(FIELD_STRUCTURE_2D<T>& phi, const T& epsilon, FIELD_STRUCTURE_2D<T>& viscosity)
	{
		FIELD_STRUCTURE_2D<T> heaviside_phi;
		heaviside_phi.Initialize(viscosity.grid, 2);

		HeavisideFunction(phi, epsilon, heaviside_phi);

		if (air_water_simulation)
		{
            GRID_ITERATION_2D(phi.grid)
		    {
		    	if (air_bubble_rising)
		    	{
		    		viscosity(i, j) = air_viscosity + (water_viscosity - air_viscosity)*heaviside_phi(i, j);
		    	}
		    	if (water_drop)
		    	{
		    		viscosity(i, j) = water_viscosity + (air_viscosity - water_viscosity)*heaviside_phi(i, j);
		    	}
		    }
		}
		
		if (oil_water_simulation)
		{
            GRID_ITERATION_2D(phi.grid)
		    {
				if (dimensionless_form)
				{
					viscosity(i, j) = oil_viscosity/water_viscosity + ((T)1 - oil_viscosity/water_viscosity)*heaviside_phi(i, j);
				}
				else
				{
					viscosity(i, j) = oil_viscosity + (water_viscosity - oil_viscosity)*heaviside_phi(i, j);
				}
			}
		}

		viscosity.FillGhostCellsFrom(viscosity.array_for_this, true);
	}

	void DetermineViscosityField(FIELD_STRUCTURE_2D<T>& phi, const T& epsilon, FIELD_STRUCTURE_2D<T>& viscosity, const int& thread_id)
	{
		FIELD_STRUCTURE_2D<T> heaviside_phi;
		heaviside_phi.Initialize(viscosity.grid, 2, &multithreading);

		HeavisideFunction(phi, epsilon, heaviside_phi, thread_id);
		multithreading.Sync(thread_id);

		if (air_water_simulation)
		{
            GRID_ITERATION_2D(phi.partial_grids[thread_id])
		    {
		    	if (air_bubble_rising)
		    	{
		    		viscosity(i, j) = air_viscosity + (water_viscosity - air_viscosity)*heaviside_phi(i, j);
		    	}
		    	if (water_drop)
		    	{
		    		viscosity(i, j) = water_viscosity + (air_viscosity - water_viscosity)*heaviside_phi(i, j);
		    	}
		    }
			multithreading.Sync(thread_id);
		}
		
		if (oil_water_simulation)
		{
            GRID_ITERATION_2D(phi.partial_grids[thread_id])
		    {
				if (dimensionless_form)
				{
					viscosity(i, j) = oil_viscosity/water_viscosity + ((T)1 - oil_viscosity/water_viscosity)*heaviside_phi(i, j);
				}
				else
				{
					viscosity(i, j) = oil_viscosity + (water_viscosity - oil_viscosity)*heaviside_phi(i, j);
				}
			}
			multithreading.Sync(thread_id);
		}

		viscosity.FillGhostCellsFrom(viscosity.array_for_this, true);
	}

	void DetermineDensityField(FIELD_STRUCTURE_2D<T>& phi, const T& epsilon, FIELD_STRUCTURE_2D<T>& density)
	{
		FIELD_STRUCTURE_2D<T> heaviside_phi;
		heaviside_phi.Initialize(density.grid, 2);

		HeavisideFunction(phi, epsilon, heaviside_phi);
		
		if (air_water_simulation)
		{
            GRID_ITERATION_2D(phi.grid)
		    {
		    	if (air_bubble_rising)
		    	{
		    		density(i, j) = air_density + (water_density - air_density)*heaviside_phi(i, j);
		    	}
		    	if (water_drop)
		    	{
		    		density(i, j) = water_density + (air_density - water_density)*heaviside_phi(i, j);
		    	}
		    }
		}
		
		if (oil_water_simulation)
		{
            GRID_ITERATION_2D(phi.grid)
		    {
				if (dimensionless_form)
				{
					density(i, j) = oil_density/water_density + ((T)1 - oil_density/water_density)*heaviside_phi(i, j);
				}
				else
				{
					density(i, j) = oil_density + (water_density - oil_density)*heaviside_phi(i, j);
				}
				
			}
		}

		density.FillGhostCellsFrom(density.array_for_this, true);
	}
	
	void DetermineDensityField(FIELD_STRUCTURE_2D<T>& phi, const T& epsilon, FIELD_STRUCTURE_2D<T>& density, const int& thread_id)
	{
		FIELD_STRUCTURE_2D<T> heaviside_phi;
		heaviside_phi.Initialize(density.grid, 2, &multithreading);

		HeavisideFunction(phi, epsilon, heaviside_phi, thread_id);
		
		multithreading.Sync(thread_id);

		if (air_water_simulation)
		{
            GRID_ITERATION_2D(phi.partial_grids[thread_id])
		    {
		    	if (air_bubble_rising)
		    	{
		    		density(i, j) = air_density + (water_density - air_density)*heaviside_phi(i, j);
		    	}
		    	if (water_drop)
		    	{
		    		density(i, j) = water_density + (air_density - water_density)*heaviside_phi(i, j);
		    	}
		    }
			multithreading.Sync(thread_id);
		}
		
		if (oil_water_simulation)
		{
            GRID_ITERATION_2D(phi.partial_grids[thread_id])
		    {
				if (dimensionless_form)
				{
					density(i, j) = oil_density/water_density + ((T)1 - oil_density/water_density)*heaviside_phi(i, j);
				}
				else
				{
					density(i, j) = oil_density + (water_density - oil_density)*heaviside_phi(i, j);
				}
			}
			multithreading.Sync(thread_id);
		}

		density.FillGhostCellsFrom(density.array_for_this, true, thread_id);
	}
};