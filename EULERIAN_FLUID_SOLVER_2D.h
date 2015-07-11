#pragma once

#include "COMMON_DEFINITION.h"
#include "LEVELSET_2D.h"
#include "WORLD_DISCRETIZATION_2D.h"

// SubSolvers
#include "ADVECTION_2D.h"
#include "PROJECTION_2D.h"
#include "VISCOSITY_2D.h"

class EULERIAN_FLUID_SOLVER_2D
{
public: // Grids
	GRID_STRUCTURE_2D				base_grid;
	GRID_STRUCTURE_2D				ghost_grid;

	ARRAY<GRID_STRUCTURE_2D>		partial_base_grids;
	ARRAY<GRID_STRUCTURE_2D>		partial_ghost_grids;

public: // Dynamic Data
	// Velocity Field for General Grid
	FIELD_STRUCTURE_2D<VT>*			water_velocity_field;

	// Velocity Field for MAC Grid
	FIELD_STRUCTURE_2D<T>*			water_velocity_field_mac_x;
	FIELD_STRUCTURE_2D<T>*			water_velocity_field_mac_y;

	// Velocity Field for RK2
	FIELD_STRUCTURE_2D<T>			water_velocity_field_x_rk_3rd;
	FIELD_STRUCTURE_2D<T>			water_velocity_field_y_rk_3rd;

	// Water Levelset
	LEVELSET_2D*					water_levelset;

	// Vortex Levelset
	LEVELSET_2D*					vortex_levelset;
	FIELD_STRUCTURE_2D<T>           magnitude_of_gradient_old;
	FIELD_STRUCTURE_2D<T>			magnitude_of_gradient_new;
	FIELD_STRUCTURE_2D<T>           temp_for_vortex;
	
	WORLD_DISCRETIZATION_2D*		world_discretization;

	// Levelset For Runge-Kutta Method
	FIELD_STRUCTURE_2D<T>			temp_for_levelset;

	// Sub-Function for Reinitialization
	FIELD_STRUCTURE_2D<T>           sign_function;
	FIELD_STRUCTURE_2D<T>           phi_0;
public: // Ghost Fields
	FIELD_STRUCTURE_2D<T>			scalar_field_ghost;
	FIELD_STRUCTURE_2D<VT>			vector_field_ghost;
	FIELD_STRUCTURE_2D<T>			vector_field_mac_ghost_x;
	FIELD_STRUCTURE_2D<T>			vector_field_mac_ghost_y;

public: // Simulation Options
	bool							air_water_simulation;
	bool							oil_water_simulation;
	bool							vortex_sheet_problem;
	bool							dimensionless_form;
	bool							is_vertical, is_parallel;

public: // Viscosity and Surface Tension, Density Field - Set 1 as inside and 2 as outside
	T								water_viscosity, air_viscosity, oil_viscosity;
	T								epsilon_for_mollification;
	T								water_density, air_density, oil_density;
	T								beta_p, beta_m;
	VT								gravity;
	T								surface_tension;
	FIELD_STRUCTURE_2D<T>			viscosity_field;
	FIELD_STRUCTURE_2D<T>			density_field;

public: // Dimensionless Number
	T								R1, R2, m, K, g, f, a, A, V0, We;
	T								one_over_We;

public: // CFL number
	T								cfl_number;

public: // For CFL TimeStep Debugging
	T								c_f, g_f, s_f, v_f;

public: // Scaling number for Reinitialization
	int								scaling_number_for_reinitialization;
	int								iteration_number_for_reinitialization;

public: // SubSolvers
	ADVECTION_2D*					advecting_field_variables;
	PROJECTION_2D*					water_projection;
	VISCOSITY_2D*					viscosity_solver;

public: // For debugging
	bool							is_velocity_advection_active;
	bool							is_sourcing_active;
	bool							is_projection_active;
	bool							is_levelset_advection_active;
	bool							is_viscosity_active;

public: // Reinitialization
	bool							fastsweeping_reinitialization;
	bool							sussmanstyle_reinitialization;

public: // Option for viscosity term -- Need to be fixed as for another header
	bool							use_delta_function_formulation;
	bool							use_jump_condition_on_viscosity;

public: // CSF Model
	bool							CSF_model;

public: // MAC grid
	bool							use_mac_grid;

public: // Option for time marching
	bool							use_rk_3rd, use_rk_3rd_for_reinitialization;
	
public: // Speed up variable
	int								i_res_x, j_res_x, i_res_y, j_res_y;
	int								i_start_x, j_start_x, i_end_x, j_end_x, i_start_y, j_start_y, i_end_y, j_end_y;
	
public: // Multithreading
	MULTITHREADING*					multithreading;

public: // Option for Time Marching
	int								order_for_time_advancing;

public: // Constructor and Destructor
	EULERIAN_FLUID_SOLVER_2D(void)
		: advecting_field_variables(0), water_projection(0), water_levelset(0), viscosity_solver(0), vortex_levelset(0), water_velocity_field(0), water_velocity_field_mac_x(0), water_velocity_field_mac_y(0),
		air_water_simulation(false), oil_water_simulation(false), dimensionless_form(false), vortex_sheet_problem(false),
		is_velocity_advection_active(false), is_sourcing_active(false), is_projection_active(false), is_levelset_advection_active(false), is_viscosity_active(false),
		fastsweeping_reinitialization(false), sussmanstyle_reinitialization(false), 
		use_delta_function_formulation(false), use_jump_condition_on_viscosity(false), 
		CSF_model(false), is_vertical(false), is_parallel(false), 
		use_rk_3rd(false), use_rk_3rd_for_reinitialization(false), use_mac_grid(false),
		multithreading(0), order_for_time_advancing((int)1)
	{}

	~EULERIAN_FLUID_SOLVER_2D(void)
	{
		DELETE_POINTER(advecting_field_variables);
		DELETE_POINTER(water_projection);
		DELETE_POINTER(water_velocity_field);
		DELETE_POINTER(water_velocity_field_mac_x);
		DELETE_POINTER(water_velocity_field_mac_y);
		DELETE_POINTER(water_levelset);
		DELETE_POINTER(vortex_levelset);
		DELETE_POINTER(viscosity_solver);
	}

public: // Initialization Functions
	void InitializeFromScriptBlock(const SCRIPT_BLOCK& script_block, MULTITHREADING* multithreading_input)
	{
		// Initialize grids from outside
		SCRIPT_BLOCK grid_sb = script_block.FindBlock("GRID_STRUCTURE_2D");

		base_grid.InitializeFromBlock(grid_sb);

		InitializeFromScriptBlock(base_grid, script_block, multithreading_input);
	}

	void InitializeFromScriptBlock(const GRID_STRUCTURE_2D& world_grid, const SCRIPT_BLOCK& fluid_solver_block, MULTITHREADING* multithreading_input)
	{
		// Simulation properties from script
		int ghost_width = fluid_solver_block.GetInteger("ghost_width", 3);
		
		// Reinitialization
		fastsweeping_reinitialization = fluid_solver_block.GetBoolean("fastsweeping_reinitialization", (bool)false);
		sussmanstyle_reinitialization = fluid_solver_block.GetBoolean("sussmanstyle_reinitialization", (bool)false);
		scaling_number_for_reinitialization = fluid_solver_block.GetInteger("scaling_number_for_reinitialization", (int)5);
		iteration_number_for_reinitialization = fluid_solver_block.GetInteger("iteration_number_for_reinitialization", (int)10);

		// Viscosity
		if (air_water_simulation)
		{
			water_viscosity = fluid_solver_block.GetFloat("water_viscosity", (T)1);
			air_viscosity = fluid_solver_block.GetFloat("air_viscosity", (T)0);
		}
		if (oil_water_simulation)
		{
			water_viscosity = fluid_solver_block.GetFloat("water_viscosity", (T)1);
			oil_viscosity = fluid_solver_block.GetFloat("oil_viscosity", (T)1);
		}

		// Gravity
		gravity = fluid_solver_block.GetVector2("gravity", VT((T)0, (T)-9.8));
		
		// Density 
		if (air_water_simulation)
		{
			water_density = fluid_solver_block.GetFloat("water_density", 1);
			air_density = fluid_solver_block.GetFloat("air_density", 1);
		}
		if (oil_water_simulation)
		{
			water_density = fluid_solver_block.GetFloat("water_density", 1);
			oil_density = fluid_solver_block.GetFloat("oil_density", 1);
		}
		if (vortex_sheet_problem)
		{
			water_density = 1;
			oil_density = 1;
		}
		surface_tension = fluid_solver_block.GetFloat("surface_tension", (T)1);

		// Dimensionless Form
		if (oil_water_simulation)
		{
			dimensionless_form = fluid_solver_block.GetBoolean("dimensionless_form", (bool)false);
		}

		// Time marching order
		order_for_time_advancing = fluid_solver_block.GetInteger("order_for_time_advancing", (int)1);

		// For the pipe location
		if (oil_water_simulation)
		{
			is_vertical = world_discretization->is_vertical;
			is_parallel = world_discretization->is_parallel;
		}
		// CFL number
		cfl_number = fluid_solver_block.GetFloat("cfl_number", (T)0.5);

		// MAC grid
		use_mac_grid = fluid_solver_block.GetBoolean("use_mac_grid", (T)false);

		// Debugging
		is_velocity_advection_active = fluid_solver_block.GetBoolean("is_velocity_advection_active", (bool)false);
		is_sourcing_active = fluid_solver_block.GetBoolean("is_sourcing_active", (bool)false);
		is_projection_active = fluid_solver_block.GetBoolean("is_projection_active", (bool)false);
		is_levelset_advection_active = fluid_solver_block.GetBoolean("is_levelset_advection_active", (bool)false);
		is_viscosity_active = fluid_solver_block.GetBoolean("is_viscosity_active", (bool)false);

		// Option for viscosity term
		use_delta_function_formulation = fluid_solver_block.GetBoolean("use_delta_function_formulation", (bool)false);
		use_jump_condition_on_viscosity = fluid_solver_block.GetBoolean("use_jump_condition_on_viscosity", (bool)false);

		// Option for surface tension term
		CSF_model = fluid_solver_block.GetBoolean("CSF_model", (bool)false);

		// Time marching
		use_rk_3rd = fluid_solver_block.GetBoolean("use_rk_3rd", (bool)false);
		use_rk_3rd_for_reinitialization = fluid_solver_block.GetBoolean("use_rk_3rd_for_reinitialization", (bool)false);

		cout << "---------------FLUID SOLVER VARIABLES---------------" << endl;
		cout << "Ghost width : " << ghost_width << endl;
		cout << "Gravity: (" << gravity.i << " ," << gravity.j << ") " << endl;
		
		if (air_water_simulation)
		{
			cout << "Water density: " << water_density << endl;
			cout << "Air density: " << air_density << endl;
			cout << "Water viscosity: " << water_viscosity << endl;
			cout << "Air viscosity: " << air_viscosity << endl;
			cout << "Surface tension: " << surface_tension << endl;
		}
		if (oil_water_simulation)
		{
			cout << "Water density: " << water_density << endl;
			cout << "Oil density: " << oil_density << endl;
			cout << "Water viscosity: " << water_viscosity << endl;
			cout << "Oil viscosity: " << air_viscosity << endl;
			cout << "Surface tension: " << surface_tension << endl;
			
			if (dimensionless_form)
			{
				cout << "Dimensionless Form: true" << endl;
			}
			else
			{
				cout << "Dimensionless Form: false" << endl;
			}
		}
				
		if (order_for_time_advancing == 1)
		{
			cout << "Time advancing: Forward Euler" << endl;
		}
		if (order_for_time_advancing == 2)
		{
			cout << "Time advancing: Runge-Kutta 2nd" << endl;
		}
		if (order_for_time_advancing == 3)
		{
			cout << "Time advancing: Runge-Kutta 3rd" << endl;
		}

		if (use_rk_3rd)
		{
			cout << "Use RK 3rd : true" << endl;
		}
		else
		{
			cout << "Use RK 3rd : false" << endl;
		}

		if (use_rk_3rd_for_reinitialization)
		{
			cout << "Use RK 3rd for Reinitialization: true" << endl;
		}
		else
		{
			cout << "Use RK 3rd for Reinitialization: false" << endl;
		}

		if (fastsweeping_reinitialization)
		{
			cout << "Reinitialized by Fast Sweeping Method: true" << endl;
		}
		else
		{
			cout << "Reinitialized by Fast Sweeping Method: false" << endl;
		}

		if (sussmanstyle_reinitialization)
		{
			cout << "Reinitialized by SussmanStyle: true" << endl;
		}
		else
		{
			cout << "Reinitialized by SussmanStyle: false" << endl;
		}

		if (air_water_simulation)
		{
			if (world_discretization->large_bubble)
			{
				cout << "Large Bubble Simulation is activated" << endl;
			}
			else if (world_discretization->small_bubble)
			{
				cout << "Small Bubble Simulation is activated" << endl;
			}
		}
		
		// Grid
		base_grid.Initialize(world_grid);
		ghost_grid.Initialize(base_grid.Enlarged(ghost_width));

		// Multithreading
		multithreading = multithreading_input;
		base_grid.SplitInYDirection(multithreading->num_threads, partial_base_grids);
		ghost_grid.SplitInYDirection(multithreading->num_threads, partial_ghost_grids);
		
		// Epsilon for Mollification
		epsilon_for_mollification = (T)1.5*base_grid.dx;
		
		// Initialize Data Fields
		DELETE_POINTER(water_velocity_field);
		water_velocity_field = new FIELD_STRUCTURE_2D<VT>();
		water_velocity_field->Initialize(base_grid.i_res, base_grid.j_res, base_grid.i_start, base_grid.j_start, base_grid.x_min, base_grid.y_min, base_grid.x_max, base_grid.y_max, 2, false, true, multithreading);
		
		// Water Velocity Field for MAC grid
		DELETE_POINTER(water_velocity_field_mac_x);
		water_velocity_field_mac_x = new FIELD_STRUCTURE_2D<T>();
		water_velocity_field_mac_x->Initialize(base_grid.i_res + 1, base_grid.j_res, base_grid.i_start, base_grid.j_start, base_grid.x_min - (T)0.5*base_grid.dx, base_grid.y_min, base_grid.x_max + (T)0.5*base_grid.dx, base_grid.y_max, 2, false, true, multithreading);
		water_velocity_field_mac_x->is_x_component = true;
		water_velocity_field_mac_x->is_y_component = false;

		DELETE_POINTER(water_velocity_field_mac_y);
		water_velocity_field_mac_y = new FIELD_STRUCTURE_2D<T>();
		water_velocity_field_mac_y->Initialize(base_grid.i_res, base_grid.j_res + 1, base_grid.i_start, base_grid.j_start, base_grid.x_min, base_grid.y_min - (T)0.5*base_grid.dy, base_grid.x_max, base_grid.y_max + (T)0.5*base_grid.dy, 2, false, true, multithreading);
		water_velocity_field_mac_y->is_y_component = true;
		water_velocity_field_mac_y->is_x_component = false;
		
		// Velocity Field For RK 3rd
		water_velocity_field_x_rk_3rd.Initialize(base_grid.i_res + 1, base_grid.j_res, base_grid.i_start, base_grid.j_start, base_grid.x_min - (T)0.5*base_grid.dx, base_grid.y_min, base_grid.x_max + (T)0.5*base_grid.dx, base_grid.y_max, 2, false, true, multithreading);
		water_velocity_field_y_rk_3rd.Initialize(base_grid.i_res, base_grid.j_res + 1, base_grid.i_start, base_grid.j_start, base_grid.x_min, base_grid.y_min - (T)0.5*base_grid.dy, base_grid.x_max, base_grid.y_max + (T)0.5*base_grid.dy, 2, false, true, multithreading);

		// Ghost Field
		scalar_field_ghost.Initialize(base_grid, 3, true, false, multithreading);
		vector_field_ghost.Initialize(water_velocity_field->grid, 3, false, true, multithreading);
		vector_field_mac_ghost_x.Initialize(water_velocity_field_mac_x->grid, 3, false, true, multithreading);
		vector_field_mac_ghost_y.Initialize(water_velocity_field_mac_y->grid, 3, false, true, multithreading);
		
		// Viscosity Field and Density Field
		viscosity_field.Initialize(base_grid, 2, true, false, multithreading);
		density_field.Initialize(base_grid, 2, true, false, multithreading);

		// Water Levelset
		DELETE_POINTER(water_levelset);
		water_levelset = new LEVELSET_2D();
		water_levelset->Initialize(base_grid, 2, multithreading);

		// Vortex Levelset
		DELETE_POINTER(vortex_levelset);
		vortex_levelset = new LEVELSET_2D();
		vortex_levelset->Initialize(base_grid, 2, multithreading);

		// Initialize the values for MAC grid
		if (air_water_simulation)
		{
			water_velocity_field_mac_x->array_for_this.AssignAllValues(T());
			water_velocity_field_mac_y->array_for_this.AssignAllValues(T());
		}
		if (oil_water_simulation)
		{
			water_velocity_field_mac_x->array_for_this.AssignAllValues(T());
			water_velocity_field_mac_y->array_for_this.AssignAllValues(T());

			a = fluid_solver_block.GetFloat("a", (T)1);
			R2 = fluid_solver_block.GetFloat("R2", (T)1);
			R1 = R2/a;
			m = fluid_solver_block.GetFloat("m", (T)1);
			K = fluid_solver_block.GetFloat("K", (T)1);
			g = fluid_solver_block.GetFloat("g", (T)1);
			f = 1/(K-1)*g*(oil_density - K*water_density);
			A = m*K + POW2(a) - (T)1 + (T)2*(K - 1)*log(a);
			V0 = (f + water_density*g)*POW2(R1)/((T)4*water_viscosity)*A;
			We = oil_density*R1*POW2(V0)/surface_tension;
			one_over_We = 1/We;

			cout << "--------------PROPETIES--------------" << endl;
			cout << "a (R2/R1) = " << a << endl;
			cout << "R2 = " << R2 << endl;
			cout << "R1 = " << R1 << endl;
			cout << "m (mu2/mu1) = " << m << endl;
			cout << "K (f + rho_1*g)/(f + rho_2*g) = " << K << endl;
			cout << "A = " << A << endl;
			cout << "g = " << g << endl;
			cout << "V0 = " << V0 << endl;

			int ghost_width_y;
			ghost_width_y = water_velocity_field_mac_y->ghost_width;

			// Speed up variables
			i_res_x = water_velocity_field_mac_x->grid.i_res, j_res_x = water_velocity_field_mac_x->grid.j_res, i_res_y = water_velocity_field_mac_y->grid.i_res, j_res_y = water_velocity_field_mac_y->grid.j_res;
			i_start_x = water_velocity_field_mac_x->i_start, i_end_x = water_velocity_field_mac_x->i_end , j_start_x = water_velocity_field_mac_x->j_start, j_end_x = water_velocity_field_mac_x->j_end;
			i_start_y = water_velocity_field_mac_y->i_start, i_end_y = water_velocity_field_mac_y->i_end , j_start_y = water_velocity_field_mac_y->j_start, j_end_y = water_velocity_field_mac_y->j_end;

			if (dimensionless_form)
			{
				if (is_vertical)
				{
					for (int j = j_start_y; j <= j_end_y; j++)
					{
						for (int i = i_start_y; i <= i_end_y; i++)
						{
							T r_coor = water_velocity_field_mac_y->grid.x_min + i*water_velocity_field_mac_y->grid.dx;
							
							if (r_coor >= 1 && r_coor <= a)
							{
								water_velocity_field_mac_y->array_for_this(i, j) = 1/A*(POW2(a) - POW2(r_coor) - (T)2*(K - 1)*log(r_coor/a));
							}
							else if (r_coor < 1)
							{
								water_velocity_field_mac_y->array_for_this(i, j) = (T)1 - 1/A*(m*POW2(r_coor)*K);
							}
						}
					}
					for (int j = j_start_y; j <= j_end_y; j++)
					{
						water_velocity_field_mac_y->array_for_this(i_start_y - 1, j) = V0;
					}
				}
			}
		}
		
		// Initialize Advection Solver
		DELETE_POINTER(advecting_field_variables);
		advecting_field_variables = new ADVECTION_2D(*water_levelset, scalar_field_ghost, *vortex_levelset, *water_velocity_field, vector_field_ghost, *water_velocity_field_mac_x, *water_velocity_field_mac_y, vector_field_mac_ghost_x, vector_field_mac_ghost_y, *multithreading); 
		advecting_field_variables->use_mac_grid = use_mac_grid;
		advecting_field_variables->InitializeFromBlock(fluid_solver_block);
		
		if (oil_water_simulation)
		{
			advecting_field_variables->is_vertical = is_vertical;
			advecting_field_variables->is_parallel = is_parallel;
		}
		
		// Initialize Projection Solver
		DELETE_POINTER(water_projection);
		water_projection = new PROJECTION_2D(*water_levelset, *vortex_levelset, *water_velocity_field_mac_x, *water_velocity_field_mac_y, vector_field_mac_ghost_x, vector_field_mac_ghost_y, multithreading);
		
		if (air_water_simulation)
		{
			water_projection->air_water_simulation = true;
			water_projection->oil_water_simulation = false;
			water_projection->InitializeFromBlock(fluid_solver_block.FindBlock("PROJECTION"));
			water_projection->water_density = water_density;
			water_projection->air_density = air_density;
			water_projection->surface_tension = fluid_solver_block.GetFloat("surface_tension", (T)1);
			water_projection->use_delta_function_formulation = use_delta_function_formulation;
			water_projection->use_jump_condition_on_viscosity = use_jump_condition_on_viscosity;
			water_projection->water_viscosity = water_viscosity;
			water_projection->air_viscosity = air_viscosity;
		}
		if (oil_water_simulation)
		{
			water_projection->air_water_simulation = false;
			water_projection->oil_water_simulation = true;
			water_projection->InitializeFromBlock(fluid_solver_block.FindBlock("PROJECTION"));
			if (dimensionless_form)
			{
				water_projection->dimensionless_form = true;
				water_projection->We = We;
				water_projection->R1 = R1;
			}
			else
			{
				water_projection->dimensionless_form = false;
				water_projection->We = (T)0;
			}
			
			water_projection->water_density = water_density;
			water_projection->oil_density = oil_density;
			water_projection->surface_tension = fluid_solver_block.GetFloat("surface_tension", (T)1);
			water_projection->use_delta_function_formulation = use_delta_function_formulation;
			water_projection->use_jump_condition_on_viscosity = use_jump_condition_on_viscosity;
			water_projection->water_viscosity = water_viscosity;
			water_projection->oil_viscosity = oil_viscosity;
		}
		
		if (vortex_sheet_problem)
		{
			water_projection->vortex_sheet_problem = true;
			water_projection->InitializeFromBlock(fluid_solver_block.FindBlock("PROJECTION"));
			water_projection->water_density = water_density;
			water_projection->oil_density = water_density;
		}

		if (air_water_simulation)
		{
			if (water_projection->air_bubble_rising)
			{
				water_projection->poisson_solver.density_p = water_density;
				water_projection->poisson_solver.density_m = air_density;
			}
			if (water_projection->water_drop)
			{
				water_projection->poisson_solver.density_p = air_density;
				water_projection->poisson_solver.density_m = water_density;
			}
		}
		
		// Initialize viscosity solver
		DELETE_POINTER(viscosity_solver);
		viscosity_solver = new VISCOSITY_2D(*multithreading, *water_levelset, *water_velocity_field_mac_x, *water_velocity_field_mac_y, scalar_field_ghost, vector_field_mac_ghost_x, vector_field_mac_ghost_y);

		if (air_water_simulation)
		{
			viscosity_solver->air_water_simulation = true;
			viscosity_solver->oil_water_simulation = false;
		}
		if (oil_water_simulation)
		{
			viscosity_solver->air_water_simulation = false;
			viscosity_solver->oil_water_simulation = true;
			
			if (dimensionless_form)
			{
				viscosity_solver->dimensionless_form = true;
			}
			else
			{
				viscosity_solver->dimensionless_form = false;
			}

			viscosity_solver->InitializeFromScriptBlock(fluid_solver_block);
		}

		if (oil_water_simulation)
		{
			water_levelset->is_axisymmetric = true;
		}

		// Initialize water levelset
		water_levelset->AssignAllValuesLevelset(base_grid.dx*(T)3);
		water_levelset->FillGhostCellsFromPointer(&(water_levelset->phi), false);
		water_levelset->curvature_by_normal_vector = fluid_solver_block.GetBoolean("curvature_by_normal_vector", false);
		water_levelset->curvature_by_levelset = fluid_solver_block.GetBoolean("curvature_by_levelset", false);

		// Initialize vortex levelset
		vortex_levelset->AssignAllValuesLevelset((T)1);
		vortex_levelset->FillGhostCellsFromPointer(&(vortex_levelset->phi), false);
		
		sign_function.Initialize(water_levelset->signed_distance_field.grid, 3, multithreading);
	}

public: // Advancing
	void AdvanceOneTimeStep(const T& dt)
	{
		if (order_for_time_advancing == 1)
		{
			cout << "Time step : " << dt << endl;

			if (is_levelset_advection_active)
			{
				if (vortex_sheet_problem)
				{
					magnitude_of_gradient_old.Initialize(base_grid, 2);

					water_levelset->ComputeGradient();

					GRID_ITERATION_2D(base_grid)
					{
						magnitude_of_gradient_old(i, j) = sqrt(POW2(water_levelset->gradient(i, j).x) + POW2(water_levelset->gradient(i, j).y));
					}

					Levelset_Advection(dt);
						
					water_levelset->FillGhostCellsContinuousDerivativesFromPointer(&(water_levelset->phi), true);
					
					Vortex_Advection(dt);
						
					vortex_levelset->FillGhostCellsContinuousDerivativesFromPointer(&(vortex_levelset->phi), true);
				
					Reinitialization(dt);
					
					water_levelset->FillGhostCellsContinuousDerivativesFromPointer(&(water_levelset->phi), true);
										
					magnitude_of_gradient_new.Initialize(base_grid, 2);

					water_levelset->ComputeGradient();

					GRID_ITERATION_2D(base_grid)
					{
						magnitude_of_gradient_new(i, j) = sqrt(POW2(water_levelset->gradient(i, j).x) + POW2(water_levelset->gradient(i, j).y));
					}
				
					GRID_ITERATION_2D(base_grid)
					{
						T ratio_of_gradient = magnitude_of_gradient_new(i, j)/magnitude_of_gradient_old(i, j);
						vortex_levelset->arr(i, j) = ratio_of_gradient*vortex_levelset->arr(i, j);
					}
				}
				else
				{
					Levelset_Advection(dt);
					
					water_levelset->FillGhostCellsContinuousDerivativesFromPointer(&(water_levelset->phi), true);

					Reinitialization(dt);
				
					water_levelset->FillGhostCellsContinuousDerivativesFromPointer(&(water_levelset->phi), true);
				}
			}
			
			int i(0), j(0);
			
			SolveForRK(dt);
		}
		
		if (order_for_time_advancing == 2)
		{
			cout << "Time step : " << dt << endl;

			if (is_levelset_advection_active)
			{
				if (vortex_sheet_problem)
				{
					magnitude_of_gradient_old.Initialize(base_grid, 2);

					water_levelset->ComputeGradient();

					GRID_ITERATION_2D(base_grid)
					{
						magnitude_of_gradient_old(i, j) = sqrt(POW2(water_levelset->gradient(i, j).x) + POW2(water_levelset->gradient(i, j).y));
					}
					
					temp_for_levelset.Initialize(base_grid, 2);
			
					GRID_ITERATION_2D(base_grid)
					{
						T temp = water_levelset->arr(i, j);
						temp_for_levelset(i, j) = temp;
					}

					Levelset_Advection(dt);
					Levelset_Advection(dt);
	
					GRID_ITERATION_2D(base_grid)
					{
						water_levelset->arr(i, j) = (T)0.5*temp_for_levelset(i, j) + (T)0.5*water_levelset->arr(i, j);
					}
						
					water_levelset->FillGhostCellsContinuousDerivativesFromPointer(&(water_levelset->phi), true);
					
					temp_for_vortex.Initialize(base_grid, 2);
			
					GRID_ITERATION_2D(base_grid)
					{
						T temp = vortex_levelset->arr(i, j);
						temp_for_vortex(i, j) = temp;
					}

					Vortex_Advection(dt);
					Vortex_Advection(dt);
	
					GRID_ITERATION_2D(base_grid)
					{
						vortex_levelset->arr(i, j) = (T)0.5*temp_for_vortex(i, j) + (T)0.5*vortex_levelset->arr(i, j);
					}
					
					vortex_levelset->FillGhostCellsContinuousDerivativesFromPointer(&(vortex_levelset->phi), true);
				
					temp_for_levelset.Initialize(base_grid, 2);
			
					GRID_ITERATION_2D(base_grid)
					{
						T temp = water_levelset->arr(i, j);
						temp_for_levelset(i, j) = temp;
					}

					Reinitialization(dt);
					Reinitialization(dt);
	
					GRID_ITERATION_2D(base_grid)
					{
						water_levelset->arr(i, j) = (T)0.5*temp_for_levelset(i, j) + (T)0.5*water_levelset->arr(i, j);
					}
					
					water_levelset->FillGhostCellsContinuousDerivativesFromPointer(&(water_levelset->phi), true);
										
					magnitude_of_gradient_new.Initialize(base_grid, 2);

					water_levelset->ComputeGradient();

					GRID_ITERATION_2D(base_grid)
					{
						magnitude_of_gradient_new(i, j) = sqrt(POW2(water_levelset->gradient(i, j).x) + POW2(water_levelset->gradient(i, j).y));
					}
				
					GRID_ITERATION_2D(base_grid)
					{
						T ratio_of_gradient = magnitude_of_gradient_new(i, j)/magnitude_of_gradient_old(i, j);
						vortex_levelset->arr(i, j) = ratio_of_gradient*vortex_levelset->arr(i, j);
					}
				}
				else
				{
					temp_for_levelset.Initialize(base_grid, 2);

					GRID_ITERATION_2D(base_grid)
					{
						T temp = water_levelset->arr(i, j);
						temp_for_levelset(i, j) = temp;
					}

					Levelset_Advection(dt);
					Levelset_Advection(dt);

					GRID_ITERATION_2D(base_grid)
					{
						water_levelset->arr(i, j) = (T)0.5*temp_for_levelset(i, j) + (T)0.5*water_levelset->arr(i, j);
					}

					water_levelset->FillGhostCellsContinuousDerivativesFromPointer(&(water_levelset->phi), true);

					temp_for_levelset.Initialize(base_grid, 2);

					GRID_ITERATION_2D(base_grid)
					{
						T temp = water_levelset->arr(i, j);
						temp_for_levelset(i, j) = temp;
					}

					Reinitialization(dt);
					Reinitialization(dt);

					GRID_ITERATION_2D(base_grid)
					{
						water_levelset->arr(i, j) = (T)0.5*temp_for_levelset(i, j) + (T)0.5*water_levelset->arr(i, j);
					}

					water_levelset->FillGhostCellsContinuousDerivativesFromPointer(&(water_levelset->phi), true);
				}
			}
			
			int i(0), j(0);
			
			GRID_ITERATION_2D(water_velocity_field_mac_x->grid)
			{
				T Temp_x = water_velocity_field_mac_x->array_for_this(i, j);
				water_velocity_field_x_rk_3rd(i, j) = Temp_x;
			}
			   GRID_ITERATION_2D(water_velocity_field_mac_y->grid)
			{
				T Temp_y = water_velocity_field_mac_y->array_for_this(i, j);
				water_velocity_field_y_rk_3rd(i, j) = Temp_y;
			}

			SolveForRK(dt);
			SolveForRK(dt);

			GRID_ITERATION_2D(water_velocity_field_mac_x->grid)
			{
				water_velocity_field_mac_x->array_for_this(i ,j) = (T)0.5*water_velocity_field_x_rk_3rd(i, j) + (T)0.5*water_velocity_field_mac_x->array_for_this(i, j);
			}
			GRID_ITERATION_2D(water_velocity_field_mac_y->grid)
			{
				water_velocity_field_mac_y->array_for_this(i, j) = (T)0.5*water_velocity_field_y_rk_3rd(i, j) + (T)0.5*water_velocity_field_mac_y->array_for_this(i, j);
			}
		}

		if (order_for_time_advancing == 3)
		{
			cout << "Time step : " << dt << endl;
			
			if (is_levelset_advection_active)
			{
				if (vortex_sheet_problem)
				{
					magnitude_of_gradient_old.Initialize(base_grid, 2);

					water_levelset->ComputeGradient();

					GRID_ITERATION_2D(base_grid)
					{
						magnitude_of_gradient_old(i, j) = sqrt(POW2(water_levelset->gradient(i, j).x) + POW2(water_levelset->gradient(i, j).y));
					}
					
					temp_for_levelset.Initialize(base_grid, 2);
			
					GRID_ITERATION_2D(base_grid)
					{
						T temp = water_levelset->arr(i, j);
						temp_for_levelset(i, j) = temp;
					}

					Levelset_Advection(dt);
					Levelset_Advection(dt);
	
					GRID_ITERATION_2D(base_grid)
					{
						water_levelset->arr(i, j) = (T)0.75*temp_for_levelset(i, j) + (T)0.25*water_levelset->arr(i, j);
					}
						
					Levelset_Advection(dt);
	
					GRID_ITERATION_2D(base_grid)
					{
						water_levelset->arr(i, j) = 1/(T)3*temp_for_levelset(i, j) + 2/(T)3*water_levelset->arr(i, j);
					}
					water_levelset->FillGhostCellsContinuousDerivativesFromPointer(&(water_levelset->phi), true);
					
					temp_for_vortex.Initialize(base_grid, 2);
			
					GRID_ITERATION_2D(base_grid)
					{
						T temp = vortex_levelset->arr(i, j);
						temp_for_vortex(i, j) = temp;
					}

					Vortex_Advection(dt);
					Vortex_Advection(dt);
	
					GRID_ITERATION_2D(base_grid)
					{
						vortex_levelset->arr(i, j) = (T)0.75*temp_for_vortex(i, j) + (T)0.25*vortex_levelset->arr(i, j);
					}
					
					Vortex_Advection(dt);
	
					GRID_ITERATION_2D(base_grid)
					{
						vortex_levelset->arr(i, j) = 1/(T)3*temp_for_vortex(i, j) + 2/(T)3*vortex_levelset->arr(i, j);
					}
					vortex_levelset->FillGhostCellsContinuousDerivativesFromPointer(&(vortex_levelset->phi), true);
				
					temp_for_levelset.Initialize(base_grid, 2);
			
					GRID_ITERATION_2D(base_grid)
					{
						T temp = water_levelset->arr(i, j);
						temp_for_levelset(i, j) = temp;
					}

					Reinitialization(dt);
					Reinitialization(dt);
	
					GRID_ITERATION_2D(base_grid)
					{
						water_levelset->arr(i, j) = (T)0.75*temp_for_levelset(i, j) + (T)0.25*water_levelset->arr(i, j);
					}
					
					Reinitialization(dt);
	
					GRID_ITERATION_2D(base_grid)
					{
						water_levelset->arr(i, j) = 1/(T)3*temp_for_levelset(i, j) + 2/(T)3*water_levelset->arr(i, j);
					}
					water_levelset->FillGhostCellsContinuousDerivativesFromPointer(&(water_levelset->phi), true);
										
					FIELD_STRUCTURE_2D<T> magnitude_of_gradient_new;
					magnitude_of_gradient_new.Initialize(base_grid, 2);

					water_levelset->ComputeGradient();

					GRID_ITERATION_2D(base_grid)
					{
						magnitude_of_gradient_new(i, j) = sqrt(POW2(water_levelset->gradient(i, j).x) + POW2(water_levelset->gradient(i, j).y));
					}
				
					GRID_ITERATION_2D(base_grid)
					{
						T ratio_of_gradient = magnitude_of_gradient_new(i, j)/magnitude_of_gradient_old(i, j);
						vortex_levelset->arr(i, j) = ratio_of_gradient*vortex_levelset->arr(i, j);
					}
				}
				else
				{
					temp_for_levelset.Initialize(base_grid, 2);

					GRID_ITERATION_2D(base_grid)
					{
						T temp = water_levelset->arr(i, j);
						temp_for_levelset(i, j) = temp;
					}

					Levelset_Advection(dt);
					Levelset_Advection(dt);

					GRID_ITERATION_2D(base_grid)
					{
						water_levelset->arr(i, j) = (T)0.75*temp_for_levelset(i, j) + (T)0.25*water_levelset->arr(i, j);
					}

					Levelset_Advection(dt);

					GRID_ITERATION_2D(base_grid)
					{
						water_levelset->arr(i, j) = 1/(T)3*temp_for_levelset(i, j) + 2/(T)3*water_levelset->arr(i, j);
					}
					water_levelset->FillGhostCellsContinuousDerivativesFromPointer(&(water_levelset->phi), true);

					temp_for_levelset.Initialize(base_grid, 2);

					GRID_ITERATION_2D(base_grid)
					{
						T temp = water_levelset->arr(i, j);
						temp_for_levelset(i, j) = temp;
					}

					Reinitialization(dt);
					Reinitialization(dt);

					GRID_ITERATION_2D(base_grid)
					{
						water_levelset->arr(i, j) = (T)0.75*temp_for_levelset(i, j) + (T)0.25*water_levelset->arr(i, j);
					}

					Reinitialization(dt);

					GRID_ITERATION_2D(base_grid)
					{
						water_levelset->arr(i, j) = 1/(T)3*temp_for_levelset(i, j) + 2/(T)3*water_levelset->arr(i, j);
					}
					water_levelset->FillGhostCellsContinuousDerivativesFromPointer(&(water_levelset->phi), true);
				}
			}
			
			int i(0), j(0);
			
			GRID_ITERATION_2D(water_velocity_field_mac_x->grid)
			{
				T Temp_x = water_velocity_field_mac_x->array_for_this(i, j);
				water_velocity_field_x_rk_3rd(i, j) = Temp_x;
			}
			   GRID_ITERATION_2D(water_velocity_field_mac_y->grid)
			{
				T Temp_y = water_velocity_field_mac_y->array_for_this(i, j);
				water_velocity_field_y_rk_3rd(i, j) = Temp_y;
			}

			SolveForRK(dt);
			SolveForRK(dt);

			GRID_ITERATION_2D(water_velocity_field_mac_x->grid)
			{
				water_velocity_field_mac_x->array_for_this(i ,j) = (T)0.75*water_velocity_field_x_rk_3rd(i, j) + (T)0.25*water_velocity_field_mac_x->array_for_this(i, j);
			}
			GRID_ITERATION_2D(water_velocity_field_mac_y->grid)
			{
				water_velocity_field_mac_y->array_for_this(i, j) = (T)0.75*water_velocity_field_y_rk_3rd(i, j) + (T)0.25*water_velocity_field_mac_y->array_for_this(i, j);
			}

			SolveForRK(dt);

			GRID_ITERATION_2D(water_velocity_field_mac_x->grid)
			{
				water_velocity_field_mac_x->array_for_this(i, j) = 1/(T)3*water_velocity_field_x_rk_3rd(i, j) + 2/(T)3*water_velocity_field_mac_x->array_for_this(i, j);
			}
			GRID_ITERATION_2D(water_velocity_field_mac_y->grid)
			{
			    water_velocity_field_mac_y->array_for_this(i, j) = 1/(T)3*water_velocity_field_y_rk_3rd(i, j) + 2/(T)3*water_velocity_field_mac_y->array_for_this(i, j);
			}
		}
	}

	void AdvanceOneTimeStepThread(const T& dt, const int& thread_id = 0)
	{
		if (order_for_time_advancing == 1)
		{
			BEGIN_HEAD_THREAD_WORK
			{
				cout << "Time step : " << dt << endl;
			}
			END_HEAD_THREAD_WORK;

			if (is_levelset_advection_active)
			{
				if (vortex_sheet_problem)
				{
					magnitude_of_gradient_old.Initialize(base_grid, 2, multithreading);

					water_levelset->ComputeGradient(thread_id);

					BEGIN_GRID_ITERATION_2D(partial_base_grids[thread_id])
					{
						magnitude_of_gradient_old(i, j) = sqrt(POW2(water_levelset->gradient(i, j).x) + POW2(water_levelset->gradient(i, j).y));
					}
					END_GRID_ITERATION_2D;

					Levelset_Advection(dt, thread_id);
						
					Vortex_Advection(dt, thread_id);
						
					Reinitialization(dt, thread_id);
					
					water_levelset->FillGhostCellsContinuousDerivativesFromPointer(&(water_levelset->phi), true);
										
					magnitude_of_gradient_new.Initialize(base_grid, 2, multithreading);

					water_levelset->ComputeGradient(thread_id);

					BEGIN_GRID_ITERATION_2D(partial_base_grids[thread_id])
					{
						magnitude_of_gradient_new(i, j) = sqrt(POW2(water_levelset->gradient(i, j).x) + POW2(water_levelset->gradient(i, j).y));
					}
					END_GRID_ITERATION_2D;
				
					BEGIN_GRID_ITERATION_2D(partial_base_grids[thread_id])
					{
						T ratio_of_gradient = magnitude_of_gradient_new(i, j)/magnitude_of_gradient_old(i, j);
						vortex_levelset->arr(i, j) = ratio_of_gradient*vortex_levelset->arr(i, j);
					}
					END_GRID_ITERATION_2D;
				}
				else
				{
					Levelset_Advection(dt, thread_id);
					
					Reinitialization(dt, thread_id);
				}	
			}
			
			SolveForRK(dt, thread_id);
		}
		
		if (order_for_time_advancing == 2)
		{
			BEGIN_HEAD_THREAD_WORK
			{
				cout << "Time step : " << dt << endl;
			}
			END_HEAD_THREAD_WORK;

			if (is_levelset_advection_active)
			{
				if (vortex_sheet_problem)
				{
					magnitude_of_gradient_old.Initialize(base_grid, 2, multithreading);

					water_levelset->ComputeGradient(thread_id);

					BEGIN_GRID_ITERATION_2D(partial_base_grids[thread_id])
					{
						magnitude_of_gradient_old(i, j) = sqrt(POW2(water_levelset->gradient(i, j).x) + POW2(water_levelset->gradient(i, j).y));
					}
					END_GRID_ITERATION_2D;

					temp_for_levelset.Initialize(base_grid, 2, multithreading);
			
					BEGIN_GRID_ITERATION_2D(partial_base_grids[thread_id])
					{
						T temp = water_levelset->arr(i, j);
						temp_for_levelset(i, j) = temp;
					}
					END_GRID_ITERATION_2D;

					Levelset_Advection(dt, thread_id);
					Levelset_Advection(dt, thread_id);
	
					BEGIN_GRID_ITERATION_2D(partial_base_grids[thread_id])
					{
						water_levelset->arr(i, j) = (T)0.5*temp_for_levelset(i, j) + (T)0.5*water_levelset->arr(i, j);
					}
					END_GRID_ITERATION_2D;

					temp_for_vortex.Initialize(base_grid, 2);
			
					BEGIN_GRID_ITERATION_2D(partial_base_grids[thread_id])
					{
						T temp = vortex_levelset->arr(i, j);
						temp_for_vortex(i, j) = temp;
					}
					END_GRID_ITERATION_2D;
					
					Vortex_Advection(dt, thread_id);
					Vortex_Advection(dt, thread_id);
	
					BEGIN_GRID_ITERATION_2D(partial_base_grids[thread_id])
					{
						vortex_levelset->arr(i, j) = (T)0.5*temp_for_vortex(i, j) + (T)0.5*vortex_levelset->arr(i, j);
					}
					END_GRID_ITERATION_2D;

					temp_for_levelset.Initialize(base_grid, 2);
			
					BEGIN_GRID_ITERATION_2D(partial_base_grids[thread_id])
					{
						T temp = water_levelset->arr(i, j);
						temp_for_levelset(i, j) = temp;
					}
					END_GRID_ITERATION_2D;

					Reinitialization(dt, thread_id);
					Reinitialization(dt, thread_id);
	
					BEGIN_GRID_ITERATION_2D(partial_base_grids[thread_id])
					{
						water_levelset->arr(i, j) = (T)0.5*temp_for_levelset(i, j) + (T)0.5*water_levelset->arr(i, j);
					}
					END_GRID_ITERATION_2D;

					magnitude_of_gradient_new.Initialize(base_grid, 2, multithreading);

					water_levelset->ComputeGradient(thread_id);

					BEGIN_GRID_ITERATION_2D(partial_base_grids[thread_id])
					{
						magnitude_of_gradient_new(i, j) = sqrt(POW2(water_levelset->gradient(i, j).x) + POW2(water_levelset->gradient(i, j).y));
					}
					END_GRID_ITERATION_2D;

					BEGIN_GRID_ITERATION_2D(partial_base_grids[thread_id])
					{
						T ratio_of_gradient = magnitude_of_gradient_new(i, j)/magnitude_of_gradient_old(i, j);
						vortex_levelset->arr(i, j) = ratio_of_gradient*vortex_levelset->arr(i, j);
					}
					END_GRID_ITERATION_2D;
				}
				else
				{
					BEGIN_HEAD_THREAD_WORK
					{
						temp_for_levelset.Initialize(base_grid, 2, multithreading);
					}
					END_HEAD_THREAD_WORK;

					BEGIN_GRID_ITERATION_2D(partial_base_grids[thread_id])
					{
						T temp = water_levelset->arr(i, j);
						temp_for_levelset(i, j) = temp;
					}
					END_GRID_ITERATION_2D;

					Levelset_Advection(dt, thread_id);
					Levelset_Advection(dt, thread_id);

					BEGIN_GRID_ITERATION_2D(partial_base_grids[thread_id])
					{
						water_levelset->arr(i, j) = (T)0.5*temp_for_levelset(i, j) + (T)0.5*water_levelset->arr(i, j);
					}
					END_GRID_ITERATION_2D;

					BEGIN_GRID_ITERATION_2D(partial_base_grids[thread_id])
					{
						T temp = water_levelset->arr(i, j);
						temp_for_levelset(i, j) = temp;
					}
					END_GRID_ITERATION_2D;

					Reinitialization(dt, thread_id);
					Reinitialization(dt, thread_id);

					BEGIN_GRID_ITERATION_2D(partial_base_grids[thread_id])
					{
						water_levelset->arr(i, j) = (T)0.5*temp_for_levelset(i, j) + (T)0.5*water_levelset->arr(i, j);
					}
					END_GRID_ITERATION_2D;
				}
			}
			
			BEGIN_GRID_ITERATION_2D(water_velocity_field_mac_x->partial_grids[thread_id])
			{
				T Temp_x = water_velocity_field_mac_x->array_for_this(i, j);
				water_velocity_field_x_rk_3rd(i, j) = Temp_x;
			}
			END_GRID_ITERATION_2D;

			BEGIN_GRID_ITERATION_2D(water_velocity_field_mac_y->partial_grids[thread_id])
			{
				T Temp_y = water_velocity_field_mac_y->array_for_this(i, j);
				water_velocity_field_y_rk_3rd(i, j) = Temp_y;
			}
			END_GRID_ITERATION_2D;

			SolveForRK(dt, thread_id);
			SolveForRK(dt, thread_id);

			BEGIN_GRID_ITERATION_2D(water_velocity_field_mac_x->partial_grids[thread_id])
			{
				water_velocity_field_mac_x->array_for_this(i ,j) = (T)0.5*water_velocity_field_x_rk_3rd(i, j) + (T)0.5*water_velocity_field_mac_x->array_for_this(i, j);
			}
			END_GRID_ITERATION_2D;

			BEGIN_GRID_ITERATION_2D(water_velocity_field_mac_y->partial_grids[thread_id])
			{
				water_velocity_field_mac_y->array_for_this(i, j) = (T)0.5*water_velocity_field_y_rk_3rd(i, j) + (T)0.5*water_velocity_field_mac_y->array_for_this(i, j);
			}
			END_GRID_ITERATION_2D;
		}

		if (order_for_time_advancing == 3)
		{
			BEGIN_HEAD_THREAD_WORK
			{
				cout << "Time step : " << dt << endl;
			}
			END_HEAD_THREAD_WORK;
			
			if (is_levelset_advection_active)
			{
				if (vortex_sheet_problem)
				{
					magnitude_of_gradient_old.Initialize(base_grid, 2, multithreading);

					water_levelset->ComputeGradient(thread_id);

					BEGIN_GRID_ITERATION_2D(partial_base_grids[thread_id])
					{
						magnitude_of_gradient_old(i, j) = sqrt(POW2(water_levelset->gradient(i, j).x) + POW2(water_levelset->gradient(i, j).y));
					}
					END_GRID_ITERATION_2D;

					temp_for_levelset.Initialize(base_grid, 2, multithreading);
			
					BEGIN_GRID_ITERATION_2D(partial_base_grids[thread_id])
					{
						T temp = water_levelset->arr(i, j);
						temp_for_levelset(i, j) = temp;
					}
					END_GRID_ITERATION_2D;

					Levelset_Advection(dt, thread_id);
					Levelset_Advection(dt, thread_id);
	
					BEGIN_GRID_ITERATION_2D(partial_base_grids[thread_id])
					{
						water_levelset->arr(i, j) = (T)0.75*temp_for_levelset(i, j) + (T)0.25*water_levelset->arr(i, j);
					}
					END_GRID_ITERATION_2D;

					Levelset_Advection(dt, thread_id);
	
					BEGIN_GRID_ITERATION_2D(partial_base_grids[thread_id])
					{
						water_levelset->arr(i, j) = 1/(T)3*temp_for_levelset(i, j) + 2/(T)3*water_levelset->arr(i, j);
					}
					END_GRID_ITERATION_2D;

					temp_for_vortex.Initialize(base_grid, 2, multithreading);
			
					BEGIN_GRID_ITERATION_2D(partial_base_grids[thread_id])
					{
						T temp = vortex_levelset->arr(i, j);
						temp_for_vortex(i, j) = temp;
					}
					END_GRID_ITERATION_2D;

					Vortex_Advection(dt, thread_id);
					Vortex_Advection(dt, thread_id);
	
					BEGIN_GRID_ITERATION_2D(partial_base_grids[thread_id])
					{
						vortex_levelset->arr(i, j) = (T)0.75*temp_for_vortex(i, j) + (T)0.25*vortex_levelset->arr(i, j);
					}
					END_GRID_ITERATION_2D;

					Vortex_Advection(dt, thread_id);
	
					BEGIN_GRID_ITERATION_2D(partial_base_grids[thread_id])
					{
						vortex_levelset->arr(i, j) = 1/(T)3*temp_for_vortex(i, j) + 2/(T)3*vortex_levelset->arr(i, j);
					}
					END_GRID_ITERATION_2D;
				
					temp_for_levelset.Initialize(base_grid, 2, multithreading);
			
					BEGIN_GRID_ITERATION_2D(partial_base_grids[thread_id])
					{
						T temp = water_levelset->arr(i, j);
						temp_for_levelset(i, j) = temp;
					}
					END_GRID_ITERATION_2D;

					Reinitialization(dt, thread_id);
					Reinitialization(dt, thread_id);
	
					BEGIN_GRID_ITERATION_2D(partial_base_grids[thread_id])
					{
						water_levelset->arr(i, j) = (T)0.75*temp_for_levelset(i, j) + (T)0.25*water_levelset->arr(i, j);
					}
					END_GRID_ITERATION_2D;
					
					Reinitialization(dt, thread_id);
	
					BEGIN_GRID_ITERATION_2D(partial_base_grids[thread_id])
					{
						water_levelset->arr(i, j) = 1/(T)3*temp_for_levelset(i, j) + 2/(T)3*water_levelset->arr(i, j);
					}
					END_GRID_ITERATION_2D;
										
					magnitude_of_gradient_new.Initialize(base_grid, 2, multithreading);

					water_levelset->ComputeGradient(thread_id);

					BEGIN_GRID_ITERATION_2D(partial_base_grids[thread_id])
					{
						magnitude_of_gradient_new(i, j) = sqrt(POW2(water_levelset->gradient(i, j).x) + POW2(water_levelset->gradient(i, j).y));
					}
					END_GRID_ITERATION_2D;

					BEGIN_GRID_ITERATION_2D(partial_base_grids[thread_id])
					{
						T ratio_of_gradient = magnitude_of_gradient_new(i, j)/magnitude_of_gradient_old(i, j);
						vortex_levelset->arr(i, j) = ratio_of_gradient*vortex_levelset->arr(i, j);
					}
					END_GRID_ITERATION_2D;
				}
				else
				{
					BEGIN_HEAD_THREAD_WORK
					{
						temp_for_levelset.Initialize(base_grid, 2, multithreading);
					}
					END_HEAD_THREAD_WORK;

					BEGIN_GRID_ITERATION_2D(partial_base_grids[thread_id])
					{
						T temp = water_levelset->arr(i, j);
						temp_for_levelset(i, j) = temp;
					}
					END_GRID_ITERATION_2D;

					Levelset_Advection(dt, thread_id);
					Levelset_Advection(dt, thread_id);

					BEGIN_GRID_ITERATION_2D(partial_base_grids[thread_id])
					{
						water_levelset->arr(i, j) = (T)0.75*temp_for_levelset(i, j) + (T)0.25*water_levelset->arr(i, j);
					}
					END_GRID_ITERATION_2D;

					Levelset_Advection(dt, thread_id);

					BEGIN_GRID_ITERATION_2D(partial_base_grids[thread_id])
					{
						water_levelset->arr(i, j) = 1/(T)3*temp_for_levelset(i, j) + 2/(T)3*water_levelset->arr(i, j);
					}
					END_GRID_ITERATION_2D;

					BEGIN_GRID_ITERATION_2D(partial_base_grids[thread_id])
					{
						T temp = water_levelset->arr(i, j);
						temp_for_levelset(i, j) = temp;
					}
					END_GRID_ITERATION_2D;

					Reinitialization(dt, thread_id);
					Reinitialization(dt, thread_id);

					BEGIN_GRID_ITERATION_2D(partial_base_grids[thread_id])
					{
						water_levelset->arr(i, j) = (T)0.75*temp_for_levelset(i, j) + (T)0.25*water_levelset->arr(i, j);
					}
					END_GRID_ITERATION_2D;

					Reinitialization(dt, thread_id);

					BEGIN_GRID_ITERATION_2D(partial_base_grids[thread_id])
					{
						water_levelset->arr(i, j) = 1/(T)3*temp_for_levelset(i, j) + 2/(T)3*water_levelset->arr(i, j);
					}
					END_GRID_ITERATION_2D;
				}
			}
			
			int i(0), j(0);
			
			BEGIN_GRID_ITERATION_2D(water_velocity_field_mac_x->partial_grids[thread_id])
			{
				T Temp_x = water_velocity_field_mac_x->array_for_this(i, j);
				water_velocity_field_x_rk_3rd(i, j) = Temp_x;
			}
			END_GRID_ITERATION_2D;

			BEGIN_GRID_ITERATION_2D(water_velocity_field_mac_y->partial_grids[thread_id])
			{
				T Temp_y = water_velocity_field_mac_y->array_for_this(i, j);
				water_velocity_field_y_rk_3rd(i, j) = Temp_y;
			}
			END_GRID_ITERATION_2D;

			SolveForRK(dt, thread_id);
			SolveForRK(dt, thread_id);

			BEGIN_GRID_ITERATION_2D(water_velocity_field_mac_x->partial_grids[thread_id])
			{
				water_velocity_field_mac_x->array_for_this(i ,j) = (T)0.75*water_velocity_field_x_rk_3rd(i, j) + (T)0.25*water_velocity_field_mac_x->array_for_this(i, j);
			}
			END_GRID_ITERATION_2D;

			BEGIN_GRID_ITERATION_2D(water_velocity_field_mac_y->partial_grids[thread_id])
			{
				water_velocity_field_mac_y->array_for_this(i, j) = (T)0.75*water_velocity_field_y_rk_3rd(i, j) + (T)0.25*water_velocity_field_mac_y->array_for_this(i, j);
			}
			END_GRID_ITERATION_2D;

			SolveForRK(dt, thread_id);

			BEGIN_GRID_ITERATION_2D(water_velocity_field_mac_x->partial_grids[thread_id])
			{
				water_velocity_field_mac_x->array_for_this(i, j) = 1/(T)3*water_velocity_field_x_rk_3rd(i, j) + 2/(T)3*water_velocity_field_mac_x->array_for_this(i, j);
			}
			END_GRID_ITERATION_2D;
			
			BEGIN_GRID_ITERATION_2D(water_velocity_field_mac_y->partial_grids[thread_id])
			{
			    water_velocity_field_mac_y->array_for_this(i, j) = 1/(T)3*water_velocity_field_y_rk_3rd(i, j) + 2/(T)3*water_velocity_field_mac_y->array_for_this(i, j);
			}
			END_GRID_ITERATION_2D;
		}
	}
	
	void Solve(const T& dt)
	{
		cout << "Time step : " << dt << endl;

		if (is_levelset_advection_active)
		{
			Levelset_Advection(dt);
			
			if (vortex_sheet_problem)
			{
				Vortex_Advection(dt);
			}
		}
				
		if (is_velocity_advection_active)
		{
			Velocity_Advection(dt);
		}

		Reinitialization(dt);

		ApplyGraivity(dt);

		if (is_sourcing_active)
		{
			Sourcing(dt);
		}
		
		SetupBoundaryConditionForVelocity();
		
		if (is_projection_active)
		{
			Projection(dt);
		}
	}

	void SolveForRK(const T& dt)
	{
		if (is_velocity_advection_active)
		{
			Velocity_Advection(dt);
		}
		
		ApplyGraivity(dt);
		
		if (is_sourcing_active)
		{
			Sourcing(dt);
		}
		
		SetupBoundaryConditionForVelocity();
		
		if (is_projection_active)
		{
			Projection(dt);
		}
	}

	void SolveForRK(const T& dt, const int& thread_id)
	{
		if (is_velocity_advection_active)
		{
			Velocity_Advection(dt, thread_id);
		}
		
		ApplyGraivity(dt, thread_id);
		
		if (is_sourcing_active)
		{
			Sourcing(dt, thread_id);
		}
		
		SetupBoundaryConditionForVelocity();
		
		if (is_projection_active)
		{
			Projection(dt, thread_id);
		}
		
	}

public: // Simulation Steps
	void Velocity_Advection(const T& dt)
	{
		advecting_field_variables->Solve_Velocity(dt);
	}

	void Velocity_Advection(const T& dt, const int& thread_id)
	{
		advecting_field_variables->Solve_Velocity(dt, thread_id);
	}

	void Levelset_Advection(const T& dt)
	{
		water_levelset->ComputeNormals();
		advecting_field_variables->Solve_Levelset(dt);
		water_levelset->FillGhostCellsContinuousDerivativesFrom(water_levelset->phi, true);
	}

	void Levelset_Advection(const T& dt, const int& thread_id)
	{
		water_levelset->ComputeNormals();
		advecting_field_variables->Solve_Levelset(dt, thread_id);
		water_levelset->FillGhostCellsContinuousDerivativesFrom(water_levelset->phi, true);
	}

	void Vortex_Advection(const T& dt)
	{
		advecting_field_variables->Solve_Vortex(dt);
		vortex_levelset->FillGhostCellsFrom(vortex_levelset->phi, true);
	}

	void Vortex_Advection(const T& dt, const int& thread_id)
	{
		advecting_field_variables->Solve_Vortex(dt, thread_id);
		vortex_levelset->FillGhostCellsFrom(vortex_levelset->phi, true);
	}

	void Sourcing(const T& dt)
	{
		DetermineDensity();
				
		if (is_viscosity_active)
		{
			viscosity_solver->ApplyViscosity(epsilon_for_mollification, dt);
			//ApplyViscosity(dt);
		}
						
		if (CSF_model)
		{
			ApplySurfaceTension(dt);
			water_projection->CSF_model = true;
		}
	}
	
	void Sourcing(const T& dt, const int& thread_id)
	{
		DetermineDensity(thread_id);
				
		if (is_viscosity_active)
		{
			viscosity_solver->ApplyViscosity(epsilon_for_mollification, dt, thread_id);
			//ApplyViscosity(dt);
		}
						
		if (CSF_model)
		{
			ApplySurfaceTension(dt);
			water_projection->CSF_model = true;
		}
	}

	void DetermineDensity(void)
	{
		if (air_water_simulation)
		{
			if (water_projection->air_bubble_rising)
			{
				GRID_ITERATION_2D(water_levelset->grid)
				{
					if (water_levelset->arr(i, j) <= 0)
					{
						density_field(i, j) = air_density;
					}
					else if (water_levelset->arr(i, j) > 0)
					{
						density_field(i, j) = water_density;
					}
				}
			}
			else if (water_projection->water_drop)
			{
				GRID_ITERATION_2D(water_levelset->grid)
				{
					if (water_levelset->arr(i, j) <= 0)
					{
						density_field(i, j) = water_density;
					}
					else if (water_levelset->arr(i, j) > 0)
					{
						density_field(i, j) = air_density;
					}
				}
			} 
		}
		if (oil_water_simulation)
		{
			GRID_ITERATION_2D(density_field.grid)
			{
				if (dimensionless_form)
				{
					if (water_levelset->arr(i, j) <= 0)
					{
						density_field(i, j) = water_density/oil_density;
					}
					else
					{
						density_field(i, j) = 1;
					}
				}
				else
				{
					if (water_levelset->arr(i, j) <= 0)
					{
						density_field(i, j) = oil_density;
					}
					else
					{
						density_field(i, j) = water_density;
					}
				}
			}
		}
	}

	void DetermineDensity(const int& thread_id)
	{
		if (air_water_simulation)
		{
			if (water_projection->air_bubble_rising)
			{
				BEGIN_GRID_ITERATION_2D(water_levelset->partial_grids[thread_id])
				{
					if (water_levelset->arr(i, j) <= 0)
					{
						density_field(i, j) = air_density;
					}
					else if (water_levelset->arr(i, j) > 0)
					{
						density_field(i, j) = water_density;
					}
				}
				END_GRID_ITERATION_2D;
			}
			else if (water_projection->water_drop)
			{
				BEGIN_GRID_ITERATION_2D(water_levelset->partial_grids[thread_id])
				{
					if (water_levelset->arr(i, j) <= 0)
					{
						density_field(i, j) = water_density;
					}
					else if (water_levelset->arr(i, j) > 0)
					{
						density_field(i, j) = air_density;
					}
				}
				END_GRID_ITERATION_2D;
			} 
		}
		if (oil_water_simulation)
		{
			GRID_ITERATION_2D(density_field.grid)
			{
				if (dimensionless_form)
				{
					if (water_levelset->arr(i, j) <= 0)
					{
						density_field(i, j) = water_density/oil_density;
					}
					else
					{
						density_field(i, j) = 1;
					}
				}
				else
				{
					if (water_levelset->arr(i, j) <= 0)
					{
						density_field(i, j) = oil_density;
					}
					else
					{
						density_field(i, j) = water_density;
					}
				}
			}
		}
	}

	void SetupBoundaryConditionForVelocity(void)
	{
		int i_start_x(water_velocity_field_mac_x->i_start), i_end_x(water_velocity_field_mac_x->i_end), j_start_x(water_velocity_field_mac_x->j_start), j_end_x(water_velocity_field_mac_x->j_end);
		int i_start_y(water_velocity_field_mac_y->i_start), i_end_y(water_velocity_field_mac_y->i_end), j_start_y(water_velocity_field_mac_y->j_start), j_end_y(water_velocity_field_mac_y->j_end);

		if (water_projection->air_bubble_rising)
		{
			// Boundary condition for x velocity - Wall Boundary Condition : This must equal to zero always
			/*for (int j = j_start_x; j <= j_end_x; j++)
			{
				water_velocity_field_mac_x->array_for_this(i_start_x , j) = -water_velocity_field_mac_x->array_for_this(i_start_x + 1, j);
				water_velocity_field_mac_x->array_for_this(i_end_x, j) = -water_velocity_field_mac_x->array_for_this(i_end_x - 1, j);
			}*/
			for (int j = j_start_x; j <= j_end_x; j++)
			{
				water_velocity_field_mac_x->array_for_this(i_start_x , j) = (T)0;
				water_velocity_field_mac_x->array_for_this(i_end_x, j) = (T)0;
			}
			/*for (int j = j_start_x + 1; j <= j_end_x - 1; j++)
			{
				water_velocity_field_mac_x->array_for_this(i_start_x + 1, j) = -water_velocity_field_mac_x->array_for_this(i_start_x, j);
				water_velocity_field_mac_x->array_for_this(i_end_x - 1, j) = -water_velocity_field_mac_x->array_for_this(i_end_x, j);
			}*/
		
			// Boundary condition for y velocity - Note that the flux through the boundary must be zero
			/*for (int i = i_start_y; i <= i_end_y; i++)
			{
				water_velocity_field_mac_y->array_for_this(i, j_start_y) = -water_velocity_field_mac_y->array_for_this(i, j_start_y + 1);
				water_velocity_field_mac_y->array_for_this(i, j_end_y) = -water_velocity_field_mac_y->array_for_this(i, j_end_y - 1);
			}*/
			for (int i = i_start_y; i <= i_end_y; i++)
			{
				water_velocity_field_mac_y->array_for_this(i, j_start_y) = (T)0;
				water_velocity_field_mac_y->array_for_this(i, j_end_y) = (T)0;
			}
			/*for (int i = i_start_y + 1; i <= i_end_y - 1; i++)
			{
				water_velocity_field_mac_y->array_for_this(i, j_start_y + 1) = -water_velocity_field_mac_y->array_for_this(i, j_start_y);
				water_velocity_field_mac_y->array_for_this(i, j_end_y - 1) = -water_velocity_field_mac_y->array_for_this(i, j_end_y);
			}*/

			//// Inflow and outflow boundary condition
			//for (int i = i_start_y; i <= i_end_y; i++)
			//{
			//	water_velocity_field_mac_y->array_for_this(i, j_start_y) = -water_velocity_field_mac_y->array_for_this(i, j_start_y + 1);
			//	water_velocity_field_mac_y->array_for_this(i, j_end_y) = 2*water_velocity_field_mac_y->array_for_this(i, j_end_y - 1) - water_velocity_field_mac_y->array_for_this(i, j_end_y - 2);
			//	//water_velocity_field_mac_y->array_for_this(i, j_end_y) = water_velocity_field_mac_y->array_for_this(i, j_end_y - 1);
			//}

			//for (int j = j_start_y; j <= j_end_y; j++)
			//{
			//	water_velocity_field_mac_y->array_for_this(i_start_y - 1, j) = -water_velocity_field_mac_y->array_for_this(i_start_y, j);
			//	water_velocity_field_mac_y->array_for_this(i_end_y + 1, j) = -water_velocity_field_mac_y->array_for_this(i_end_y, j);
			//}
			
			// Calculate the flux and modify the outflow boundary condition
			/*T integration_1(0), integration_2(0);
			for (int i = i_start_y; i < i_end_y; i++)
			{
				integration_1 += water_velocity_field_mac_y->dx*0.5*(water_velocity_field_mac_y->array_for_this(i, j_start_y) + water_velocity_field_mac_y->array_for_this(i + 1, j_start_y));	
				integration_2 += water_velocity_field_mac_y->dx*0.5*(water_velocity_field_mac_y->array_for_this(i, j_end_y) + water_velocity_field_mac_y->array_for_this(i + 1, j_end_y));	
			}
			
			T integration_constant = -integration_1/integration_2;
			for (int i = i_start_y; i <= i_end_y; i++)
			{
				water_velocity_field_mac_y->array_for_this(i, j_end_y) *= integration_constant;
			}*/
		}
		if (water_projection->water_drop)
		{
			// Boundary condition for x velocity - Wall Boundary Condition : This must equal to zero always
			for (int j = j_start_x; j <= j_end_x; j++)
			{
				water_velocity_field_mac_x->array_for_this(i_start_x , j) = 0;
				water_velocity_field_mac_x->array_for_this(i_end_x, j) = 0;
			}
			for (int j = j_start_x + 1; j <= j_end_x - 1; j++)
			{
				water_velocity_field_mac_x->array_for_this(i_start_x + 1, j) = -water_velocity_field_mac_x->array_for_this(i_start_x, j);
				water_velocity_field_mac_x->array_for_this(i_end_x - 1, j) = -water_velocity_field_mac_x->array_for_this(i_end_x, j);
			}
		
			// Boundary condition for y velocity - Note that the flux through the boundary must be zero
			for (int i = i_start_y; i <= i_end_y; i++)
			{
				water_velocity_field_mac_y->array_for_this(i, j_start_y) = 0;
				water_velocity_field_mac_y->array_for_this(i, j_end_y) = 0;
			}
			for (int i = i_start_y + 1; i <= i_end_y - 1; i++)
			{
				water_velocity_field_mac_y->array_for_this(i, j_start_y + 1) = -water_velocity_field_mac_y->array_for_this(i, j_start_y);
				water_velocity_field_mac_y->array_for_this(i, j_end_y - 1) = -water_velocity_field_mac_y->array_for_this(i, j_end_y);
			}
			
			// Calculate the flux and modify the outflow boundary condition
			/*T integration_1(0), integration_2(0);
			for (int i = i_start_y; i < i_end_y; i++)
			{
				integration_1 += water_velocity_field_mac_y->dx*0.5*(water_velocity_field_mac_y->array_for_this(i, j_start_y) + water_velocity_field_mac_y->array_for_this(i + 1, j_start_y));	
				integration_2 += water_velocity_field_mac_y->dx*0.5*(water_velocity_field_mac_y->array_for_this(i, j_end_y) + water_velocity_field_mac_y->array_for_this(i + 1, j_end_y));	
			}
			
			T integration_constant = -integration_2/integration_1;
			for (int i = i_start_y; i <= i_end_y; i++)
			{
				water_velocity_field_mac_y->array_for_this(i, j_end_y) *= integration_constant;
			}*/
		}
		if (oil_water_simulation)
		{
			for (int j = j_start_x; j <= j_end_x; j++)
			{
				water_velocity_field_mac_x->array_for_this(i_end_x, j) = (T)0;
			}
			for (int i = i_start_y; i <= i_end_y; i++)
			{
				water_velocity_field_mac_y->array_for_this(i, j_end_y) = (T)0;
			}
			for (int j = j_start_y; j <= j_end_y; j++)
			{
				water_velocity_field_mac_y->array_for_this(i_start_y, j) = V0;
			}
		}
	}

	void Projection(const T& dt)
	{
		water_projection->surface_tension = surface_tension;
		water_projection->Solve(dt);
	}

	void Projection(const T& dt, const int& thread_id)
	{
		water_projection->surface_tension = surface_tension;
		water_projection->Solve(dt, thread_id);
	}

	void Reinitialization(const T& dt)
	{
		if (fastsweeping_reinitialization)
		{
			water_levelset->FastSweepingMethodOriginal(4);
			//water_levelset->FastSweepingMethod();
		}
		if (sussmanstyle_reinitialization)
		{
			T fictious_dt = base_grid.dx/(T)scaling_number_for_reinitialization;
			
			int iter(0);
			
			sign_function.Initialize(water_levelset->signed_distance_field.grid, 3);

			GRID_ITERATION_2D(sign_function.grid)
			{
				sign_function(i, j) = water_levelset->arr(i, j)/sqrt(POW2(water_levelset->arr(i, j)) + POW2(water_levelset->grid.dx));
			}
			
			phi_0.Initialize(water_levelset->grid, 2);

			GRID_ITERATION_2D(phi_0.grid)
			{
				phi_0(i, j) = water_levelset->signed_distance_field(i, j);
			}
			phi_0.FillGhostCellsContinuousDerivativesFrom(phi_0.array_for_this, true);
			
			while (iter < iteration_number_for_reinitialization)
			{
				advecting_field_variables->ReinitializationBySussman(fictious_dt, sign_function);
				iter = iter + 1;
			}
			
		}
	}

	void Reinitialization(const T& dt, const int& thread_id)
	{
		if (fastsweeping_reinitialization)
		{
			water_levelset->FastSweepingMethodOriginal(4);
			//water_levelset->FastSweepingMethod();
		}
		if (sussmanstyle_reinitialization)
		{
			T fictious_dt = base_grid.dx/(T)scaling_number_for_reinitialization;
			
			int iter(0);
			
			BEGIN_GRID_ITERATION_2D(sign_function.partial_grids[thread_id])
			{
				sign_function(i, j) = water_levelset->arr(i, j)/sqrt(POW2(water_levelset->arr(i, j)) + POW2(water_levelset->grid.dx));
			}
			END_GRID_ITERATION_2D;

			while (iter < iteration_number_for_reinitialization)
			{
				advecting_field_variables->ReinitializationBySussman(fictious_dt, sign_function, thread_id);
				iter = iter + 1;
			}
			
		}
	}

public:	// Sourcing Functions
	void HeavisideFunction(FIELD_STRUCTURE_2D<T>& phi, const T& epsilon, FIELD_STRUCTURE_2D<T>& heaviside)
	{
		int i(0), j(0);
		const int i_start(phi.i_start), j_start(phi.j_start), i_end(phi.i_end), j_end(phi.j_end);

		T one_over_epsilon = 1/(T)epsilon, one_over_pi = 1/(T)PI;
		
		LOOPS_2D(i, j, i_start, j_start, i_end, j_end)
		{
			if (phi(i, j) < -epsilon)
			{
				heaviside(i, j) = 0;
			}
			else if (phi(i, j) > epsilon)
			{
				heaviside(i, j) = 1;
			}
			else
			{
				heaviside(i, j) = (T)0.5 + (T)0.5*phi(i, j)*one_over_epsilon + (T)0.5*one_over_pi*sin((T)PI*phi(i, j)*one_over_epsilon);
			}
		}
	}

	void DeltaFunction(FIELD_STRUCTURE_2D<T>& phi, const T& epsilon, FIELD_STRUCTURE_2D<T>& delta)
	{
		int i(0), j(0);
		const int i_start(phi.i_start), j_start(phi.j_start), i_end(phi.i_end), j_end(phi.j_end);

		T one_over_epsilon = 1/(T)epsilon;

		LOOPS_2D(i, j, i_start, j_start, i_end, j_end)
		{
			if (phi(i, j) < -epsilon)
			{
				delta(i, j) = 0;
			}
			else if (phi(i, j) >= -epsilon && phi(i, j) <= epsilon)
			{
				delta(i, j) = (T)0.5*one_over_epsilon + (T)0.5*one_over_epsilon*cos(PI*phi(i, j)*one_over_epsilon);
			}
			else
			{
				delta(i, j) = 0;
			}
		}
	}
	
	void DetermineViscosityField(FIELD_STRUCTURE_2D<T>& phi, const T& epsilon, FIELD_STRUCTURE_2D<T>& viscosity)	
	{
		int i(0), j(0);
		const int i_start(viscosity.i_start), i_end(viscosity.i_end), j_start(viscosity.j_start), j_end(viscosity.j_end);
		
		FIELD_STRUCTURE_2D<T> heaviside_phi;
		heaviside_phi.Initialize(viscosity.grid, 2);

		HeavisideFunction(phi, epsilon, heaviside_phi);

		LOOPS_2D(i, j, i_start, j_start, i_end, j_end)
		{
			if (water_projection->air_bubble_rising == true)
			{
				viscosity(i, j) = air_viscosity + (water_viscosity - air_viscosity)*heaviside_phi(i, j);
			}
			if (water_projection->water_drop == true)
			{
				viscosity(i, j) = water_viscosity + (air_viscosity - water_viscosity)*heaviside_phi(i, j);
			}
		}
		viscosity.FillGhostCellsFrom(viscosity.array_for_this, true);
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
			if (water_projection->air_bubble_rising == true)
			{
				density(i, j) = air_density + (water_density - air_density)*heaviside_phi(i, j);
			}
			if (water_projection->water_drop == true)
			{
				density(i, j) = water_density + (air_density - water_density)*heaviside_phi(i, j);
			}
		}	
		density.FillGhostCellsFrom(density.array_for_this, true);
	}

	void ComputeJumpconitionMatrixForViscosity(FIELD_STRUCTURE_2D<T>& J11, FIELD_STRUCTURE_2D<T>& J12, FIELD_STRUCTURE_2D<T>& J21, FIELD_STRUCTURE_2D<T>& J22)
	{
		water_levelset->ComputeTangential();
		
		T viscosity_p, viscosity_m;

		if (water_projection->air_bubble_rising == true)
		{
			viscosity_p = water_viscosity;
			viscosity_m = air_viscosity;
		}
		else if (water_projection->water_drop == true)
		{
			viscosity_p = air_viscosity;
			viscosity_m = water_viscosity;
		}
		
		T viscosity_difference = viscosity_p - viscosity_m;
		
		FIELD_STRUCTURE_2D<T> smooth_velocity_x, smooth_velocity_y;
		smooth_velocity_x.Initialize(base_grid, 2);
		smooth_velocity_y.Initialize(base_grid, 2);

		GRID_ITERATION_2D(base_grid)
		{
			smooth_velocity_x(i, j) = (T)0.5*(water_velocity_field_mac_x->array_for_this(i, j) + water_velocity_field_mac_x->array_for_this(i + 1, j));
			smooth_velocity_y(i, j) = (T)0.5*(water_velocity_field_mac_y->array_for_this(i, j) + water_velocity_field_mac_y->array_for_this(i, j + 1));
		}

		smooth_velocity_x.FillGhostCellsFrom(smooth_velocity_x.array_for_this, true);
		smooth_velocity_y.FillGhostCellsFrom(smooth_velocity_y.array_for_this, true);

		FIELD_STRUCTURE_2D<T> smooth_velocity_x_d_x, smooth_velocity_x_d_y, smooth_velocity_y_d_x, smooth_velocity_y_d_y;
		smooth_velocity_x_d_x.Initialize(base_grid, 2);
		smooth_velocity_x_d_y.Initialize(base_grid, 2);
		smooth_velocity_y_d_x.Initialize(base_grid, 2);
		smooth_velocity_y_d_y.Initialize(base_grid, 2);

		GRID_ITERATION_2D(base_grid)
		{
			smooth_velocity_x_d_x(i, j) = base_grid.one_over_2dx*(smooth_velocity_x(i + 1, j) - smooth_velocity_x(i - 1, j));
			smooth_velocity_x_d_y(i, j) = base_grid.one_over_2dy*(smooth_velocity_x(i, j + 1) - smooth_velocity_x(i, j - 1));
			smooth_velocity_y_d_x(i, j) = base_grid.one_over_2dx*(smooth_velocity_y(i + 1, j) - smooth_velocity_y(i - 1, j));
			smooth_velocity_y_d_y(i, j) = base_grid.one_over_2dy*(smooth_velocity_y(i, j + 1) - smooth_velocity_y(i, j - 1));
		}
		
		smooth_velocity_x_d_x.FillGhostCellsFrom(smooth_velocity_x_d_x.array_for_this, true);
		smooth_velocity_x_d_y.FillGhostCellsFrom(smooth_velocity_x_d_y.array_for_this, true);
		smooth_velocity_y_d_x.FillGhostCellsFrom(smooth_velocity_y_d_x.array_for_this, true);
		smooth_velocity_y_d_y.FillGhostCellsFrom(smooth_velocity_y_d_y.array_for_this, true);

		// Main Calculation of entries in Matrix
		GRID_ITERATION_2D(base_grid)
		{
			T sqr_t1(POW2(water_levelset->tangential(i, j).x)), sqr_t2(POW2(water_levelset->tangential(i, j).y)), t1t2(water_levelset->tangential(i, j).x*water_levelset->tangential(i, j).y);
			T sqr_n1(POW2(water_levelset->normal(i, j).x)), sqr_n2(POW2(water_levelset->normal(i, j).y)), n1n2(water_levelset->normal(i, j).x*water_levelset->normal(i, j).y);
			T xdx(smooth_velocity_x_d_x(i, j)), xdy(smooth_velocity_x_d_y(i, j)), ydx(smooth_velocity_y_d_x(i, j)), ydy(smooth_velocity_y_d_y(i, j));

			J11(i, j) = sqr_t1*xdx + t1t2*xdy;
			J11(i, j) += sqr_n1*(sqr_n1*xdx + n1n2*ydx) + n1n2*(sqr_n1*xdy + n1n2*ydy);
			J11(i, j) -= sqr_n1*(sqr_t1*xdx + t1t2*xdy) + n1n2*(sqr_t1*ydx + t1t2*ydy);
			J11(i, j) *= viscosity_difference;

			J12(i, j) = t1t2*xdx + sqr_t2*xdy;
			J12(i, j) += n1n2*(sqr_n1*xdx + n1n2*ydx) + sqr_n2*(sqr_n1*xdy + n1n2*ydy);
			J12(i, j) -= n1n2*(sqr_t1*xdx + t1t2*xdy) + sqr_n2*(sqr_t1*ydx + t1t2*ydy);
			J12(i, j) *= viscosity_difference;

			J21(i, j) = sqr_t1*ydx + t1t2*ydy;
			J21(i, j) += sqr_n1*(n1n2*xdx + sqr_n2*ydx) + n1n2*(n1n2*xdy + sqr_n2*ydy);
			J21(i, j) -= sqr_n1*(t1t2*xdx + sqr_t2*xdy) + n1n2*(t1t2*ydx + sqr_t2*ydy);
			J21(i, j) *= viscosity_difference;

			J22(i, j) = t1t2*ydx + sqr_t2*ydy;
			J22(i, j) += n1n2*(n1n2*xdx + sqr_n2*ydx) + sqr_n2*(n1n2*xdy + sqr_n2*ydy);
			J22(i, j) -= n1n2*(t1t2*xdx + sqr_t2*xdy) + sqr_n2*(t1t2*ydx + sqr_t2*ydy);
			J22(i, j) *= viscosity_difference;
		}		
	}

	void ApplyViscosity(const T& dt)
	{
		if (use_delta_function_formulation == true)
		{
			// Speed Up Variable
			int i(0), j(0);
			const T one_over_dx = base_grid.one_over_dx, one_over_dy = base_grid.one_over_dy;

			const ARRAY_2D<T>& water_array(water_levelset->arr);
			const ARRAY_2D<T>& velocity_x(water_velocity_field_mac_x->array_for_this);
			const ARRAY_2D<T>& velocity_y(water_velocity_field_mac_y->array_for_this);
			FIELD_STRUCTURE_2D<T> temp_x, temp_y;
			temp_x.Initialize(water_velocity_field_mac_x->grid, 2);
			temp_y.Initialize(water_velocity_field_mac_y->grid, 2);

			LOOPS_2D(i, j, i_start_x, j_start_x, i_end_x, j_end_x)
			{
				temp_x(i, j) = velocity_x(i, j);
			}
			
			LOOPS_2D(i, j, i_start_y, j_start_y, i_end_y, j_end_y)
			{
				temp_y(i, j) = velocity_y(i, j);
			}

			// Define for coefficient function -- coefficients of x-velocity
			FIELD_STRUCTURE_2D<T> density_half_x, viscosity_cu_x, viscosity_cd_x, viscosity_r_x, viscosity_c_x, levelset_x_half, levelset_cu_x, levelset_cd_x, levelset_r_x, levelset_c_x;
			density_half_x.Initialize(water_velocity_field_mac_x->grid, 2);
			viscosity_cu_x.Initialize(water_velocity_field_mac_x->grid, 2);
			viscosity_cd_x.Initialize(water_velocity_field_mac_x->grid, 2);
			viscosity_r_x.Initialize(water_velocity_field_mac_x->grid, 2);
			viscosity_c_x.Initialize(water_velocity_field_mac_x->grid, 2);
			levelset_x_half.Initialize(water_velocity_field_mac_x->grid, 2);
			levelset_cu_x.Initialize(water_velocity_field_mac_x->grid, 2);
			levelset_cd_x.Initialize(water_velocity_field_mac_x->grid, 2);
			levelset_r_x.Initialize(water_velocity_field_mac_x->grid, 2);
			levelset_c_x.Initialize(water_velocity_field_mac_x->grid, 2);
	
			LOOPS_2D(i, j, i_start_x, j_start_x, i_end_x, j_end_x)
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
				
			// Define for coefficient function -- coefficients of y-velocity
			FIELD_STRUCTURE_2D<T> density_half_y, viscosity_cu_y, viscosity_cl_y, viscosity_u_y, viscosity_c_y, levelset_y_half, levelset_cu_y, levelset_cl_y, levelset_u_y, levelset_c_y;
			density_half_y.Initialize(water_velocity_field_mac_y->grid, 2);
			viscosity_cu_y.Initialize(water_velocity_field_mac_y->grid, 2);
			viscosity_cl_y.Initialize(water_velocity_field_mac_y->grid, 2);
			viscosity_u_y.Initialize(water_velocity_field_mac_y->grid, 2);
			viscosity_c_y.Initialize(water_velocity_field_mac_y->grid, 2);
			levelset_y_half.Initialize(water_velocity_field_mac_y->grid, 2);
			levelset_cu_y.Initialize(water_velocity_field_mac_y->grid, 2);
			levelset_cl_y.Initialize(water_velocity_field_mac_y->grid, 2);
			levelset_u_y.Initialize(water_velocity_field_mac_y->grid, 2);
			levelset_c_y.Initialize(water_velocity_field_mac_y->grid, 2);
	
			LOOPS_2D(i, j, i_start_y, j_start_y, i_end_y, j_end_y)
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
	
			LOOPS_2D(i, j, i_start_x, j_start_x, i_end_x, j_end_x)
			{
				T one_over_density_half_x = (T)1/density_half_x(i, j);
				T coef = dt*one_over_density_half_x;
				T first_update = (T)2*coef*one_over_dx*(viscosity_r_x(i, j)*one_over_dx*(temp_x(i + 1, j) - temp_x(i, j)) - viscosity_c_x(i, j)*one_over_dx*(temp_x(i, j) - temp_x(i - 1, j)));
				T second_update = coef*one_over_dy*(viscosity_cu_x(i, j)*(one_over_dy*(temp_x(i, j + 1) - temp_x(i, j)) + one_over_dx*(temp_y(i, j + 1) - temp_y(i - 1, j + 1))) - viscosity_cd_x(i, j)*(one_over_dy*(temp_x(i, j) - temp_x(i, j - 1)) + one_over_dx*(temp_y(i, j) - temp_y(i - 1, j))));
				velocity_x(i, j) += first_update + second_update;
			}

			LOOPS_2D(i, j, i_start_y, j_start_y, i_end_y, j_end_y)
			{
				T one_over_density_half_y = (T)1/density_half_y(i, j);
				T coef = dt*one_over_density_half_y;
				T first_update = (T)2*coef*one_over_dy*(viscosity_u_y(i, j)*one_over_dy*(temp_y(i, j + 1) - temp_y(i, j)) - viscosity_c_y(i, j)*one_over_dy*(temp_y(i, j) - temp_y(i, j - 1)));
				T second_update = coef*one_over_dx*(viscosity_cu_y(i, j)*(one_over_dy*(temp_x(i + 1, j) - temp_x(i + 1, j - 1)) + one_over_dx*(temp_y(i + 1, j) - temp_y(i, j))) - viscosity_cl_y(i, j)*(one_over_dy*(temp_x(i, j) - temp_x(i, j - 1)) + one_over_dx*(temp_y(i, j) - temp_y(i - 1, j))));
				velocity_y(i, j) += first_update + second_update;
			}
		}
		else if (use_jump_condition_on_viscosity == true)
		{
			// Set up the viscosity and density
			T viscosity_p, viscosity_m;
			T density_p, density_m;

			if (water_projection->air_bubble_rising == true)
			{
				viscosity_p = water_viscosity;
				viscosity_m = air_viscosity;
				density_p = water_density;
				density_m = air_density;
			}
			else if (water_projection->water_drop == true)
			{
				viscosity_p = air_viscosity;
				viscosity_m = water_viscosity;
				density_p = air_density;
				density_m = water_density;
			}

			FIELD_STRUCTURE_2D<T> j11, j12, j21, j22;
			j11.Initialize(base_grid, 2);
			j12.Initialize(base_grid, 2);
			j21.Initialize(base_grid, 2);
			j22.Initialize(base_grid, 2);

			ComputeJumpconitionMatrixForViscosity(j11, j12, j21, j22);
			
			// Offsets for Elements of Jump Condition Matrix
			
			// Offsets for Levelset
			FIELD_STRUCTURE_2D<T> levelset_m_x, levelset_r_x, levelset_l_x, levelset_u_x, levelset_d_x, levelset_m_y, levelset_r_y, levelset_l_y, levelset_u_y, levelset_d_y;
			// x-components
			levelset_m_x.Initialize(water_velocity_field_mac_x->grid, 2);
			levelset_r_x.Initialize(water_velocity_field_mac_x->grid, 2);
			levelset_l_x.Initialize(water_velocity_field_mac_x->grid, 2);
			levelset_u_x.Initialize(water_velocity_field_mac_x->grid, 2);
			levelset_d_x.Initialize(water_velocity_field_mac_x->grid, 2);
			// y-components
			levelset_m_y.Initialize(water_velocity_field_mac_y->grid, 2);
			levelset_r_y.Initialize(water_velocity_field_mac_y->grid, 2);
			levelset_l_y.Initialize(water_velocity_field_mac_y->grid, 2);
			levelset_u_y.Initialize(water_velocity_field_mac_y->grid, 2);
			levelset_d_y.Initialize(water_velocity_field_mac_y->grid, 2);

			GRID_ITERATION_2D(water_velocity_field_mac_x->grid)
			{
				levelset_m_x(i, j) = (T)0.5*(water_levelset->arr(i - 1, j) + water_levelset->arr(i, j));
				levelset_r_x(i, j) = (T)0.5*(water_levelset->arr(i, j) + water_levelset->arr(i + 1, j));
				levelset_l_x(i, j) = (T)0.5*(water_levelset->arr(i - 2, j) + water_levelset->arr(i - 1, j));
				levelset_u_x(i, j) = (T)0.5*(water_levelset->arr(i - 1, j + 1) + water_levelset->arr(i, j + 1));
				levelset_d_x(i, j) = (T)0.5*(water_levelset->arr(i - 1, j - 1) + water_levelset->arr(i, j - 1));
			}

			GRID_ITERATION_2D(water_velocity_field_mac_y->grid)
			{
				levelset_m_y(i, j) = (T)0.5*(water_levelset->arr(i, j - 1) + water_levelset->arr(i, j));
				levelset_r_y(i, j) = (T)0.5*(water_levelset->arr(i + 1, j - 1) + water_levelset->arr(i + 1, j));
				levelset_l_y(i, j) = (T)0.5*(water_levelset->arr(i - 1, j - 1) + water_levelset->arr(i - 1, j));
				levelset_u_y(i, j) = (T)0.5*(water_levelset->arr(i, j) + water_levelset->arr(i, j + 1));
				levelset_d_y(i, j) = (T)0.5*(water_levelset->arr(i, j - 2) + water_levelset->arr(i, j - 1));
			}
			
			
			// Offsets for Velocity
			FIELD_STRUCTURE_2D<T> u_m, u_r, u_l, u_u, u_d, v_m, v_r, v_l, v_u, v_d;
			u_m.Initialize(water_velocity_field_mac_x->grid, 2);
			u_r.Initialize(water_velocity_field_mac_x->grid, 2);
			u_l.Initialize(water_velocity_field_mac_x->grid, 2);
			u_u.Initialize(water_velocity_field_mac_x->grid, 2);
			u_d.Initialize(water_velocity_field_mac_x->grid, 2);
			v_m.Initialize(water_velocity_field_mac_y->grid, 2);
			v_r.Initialize(water_velocity_field_mac_y->grid, 2);
			v_l.Initialize(water_velocity_field_mac_y->grid, 2);
			v_u.Initialize(water_velocity_field_mac_y->grid, 2);
			v_d.Initialize(water_velocity_field_mac_y->grid, 2);

			GRID_ITERATION_2D(water_velocity_field_mac_x->grid)
			{
				u_m(i, j) = water_velocity_field_mac_x->array_for_this(i, j);
				u_r(i, j) = water_velocity_field_mac_x->array_for_this(i + 1, j);
				u_l(i, j) = water_velocity_field_mac_x->array_for_this(i - 1, j);
				u_u(i, j) = water_velocity_field_mac_x->array_for_this(i, j + 1);
				u_d(i, j) = water_velocity_field_mac_x->array_for_this(i, j - 1);
			}

			GRID_ITERATION_2D(water_velocity_field_mac_y->grid)
			{
				v_m(i, j) = water_velocity_field_mac_y->array_for_this(i, j);
				v_r(i, j) = water_velocity_field_mac_y->array_for_this(i + 1, j);
				v_l(i, j) = water_velocity_field_mac_y->array_for_this(i - 1, j);
				v_u(i, j) = water_velocity_field_mac_y->array_for_this(i, j + 1);
				v_d(i, j) = water_velocity_field_mac_y->array_for_this(i, j - 1);
			}
			
			GRID_ITERATION_2D(water_velocity_field_mac_x->grid)
			{
				T uxl, uxr, uxx;
				T j11m, j11r, j11l, j11i;
				T effective_viscosity_x;

				if (levelset_m_x(i, j) > 0)
				{
					uxl = viscosity_p*(water_velocity_field_mac_x->one_over_dx*(u_m(i, j) - u_l(i, j)));
					uxr = viscosity_p*(water_velocity_field_mac_x->one_over_dx*(u_r(i, j) - u_m(i, j)));
					if (levelset_l_x(i, j) <= 0)
					{
						T subcell_l = abs(levelset_l_x(i, j))/(abs(levelset_l_x(i, j)) + abs(levelset_m_x(i, j)));
						j11m = (T)0.5*(j11(i - 1, j) + j11(i, j)), j11r = (T)0.5*(j11(i, j) + j11(i + 1, j)), j11l = (T)0.5*(j11(i - 2, j) + j11(i - 1, j));
						j11i = subcell_l*j11m + (1 - subcell_l)*j11l;
						
						effective_viscosity_x = viscosity_p*viscosity_m/(viscosity_p*subcell_l + viscosity_m*(1 - subcell_l));
						
						uxl = effective_viscosity_x*(water_velocity_field_mac_x->one_over_dx*(u_m(i, j) - u_l(i, j))) + effective_viscosity_x*j11i*subcell_l/viscosity_m; 
					}
					if (levelset_r_x(i, j) <= 0)
					{
						T subcell_r = abs(levelset_r_x(i, j))/(abs(levelset_r_x(i, j)) + abs(levelset_m_x(i, j)));
						j11m = (T)0.5*(j11(i - 1, j) + j11(i, j)), j11r = (T)0.5*(j11(i, j) + j11(i + 1, j)), j11l = (T)0.5*(j11(i - 2, j) + j11(i - 1, j));
						j11i = subcell_r*j11m + (1 - subcell_r)*j11r;
						
						effective_viscosity_x = viscosity_p*viscosity_m/(viscosity_m*subcell_r + viscosity_p*(1 - subcell_r));
						
						uxr = effective_viscosity_x*(water_velocity_field_mac_x->one_over_dx*(u_r(i, j) - u_m(i, j))) + effective_viscosity_x*j11i*subcell_r/viscosity_m; 
					}
					uxx = water_velocity_field_mac_x->one_over_dx*(uxr - uxl);
				}
				else
				{
					uxl = viscosity_m*(water_velocity_field_mac_x->one_over_dx*(u_m(i, j) - u_l(i, j)));
					uxr = viscosity_m*(water_velocity_field_mac_x->one_over_dx*(u_r(i, j) - u_m(i, j)));
						
					if (levelset_l_x(i, j) > 0)
					{
						T subcell_l = abs(levelset_l_x(i, j))/(abs(levelset_l_x(i, j)) + abs(levelset_m_x(i, j)));
						j11m = (T)0.5*(j11(i - 1, j) + j11(i, j)), j11r = (T)0.5*(j11(i, j) + j11(i + 1, j)), j11l = (T)0.5*(j11(i - 2, j) + j11(i - 1, j));
						j11i = subcell_l*j11m + (1 - subcell_l)*j11l;
						
						effective_viscosity_x = viscosity_p*viscosity_m/(viscosity_m*subcell_l + viscosity_p*(1 - subcell_l));
						
						uxl = effective_viscosity_x*(water_velocity_field_mac_x->one_over_dx*(u_m(i, j) - u_l(i, j))) - effective_viscosity_x*j11i*subcell_l/viscosity_p; 
					}
					if (levelset_r_x(i, j) > 0)
					{
						T subcell_r = abs(levelset_r_x(i, j))/(abs(levelset_r_x(i, j)) + abs(levelset_m_x(i, j)));
						j11m = (T)0.5*(j11(i - 1, j) + j11(i, j)), j11r = (T)0.5*(j11(i, j) + j11(i + 1, j)), j11l = (T)0.5*(j11(i - 2, j) + j11(i - 1, j));
						j11i = subcell_r*j11m + (1 - subcell_r)*j11r;
						
						effective_viscosity_x = viscosity_p*viscosity_m/(viscosity_m*subcell_r + viscosity_p*(1 - subcell_r));
						
						uxr = effective_viscosity_x*(water_velocity_field_mac_x->one_over_dx*(u_r(i, j) - u_m(i, j))) - effective_viscosity_x*j11i*subcell_r/viscosity_p; 
					}
					uxx = water_velocity_field_mac_x->one_over_dx*(uxr - uxl);
				}

				T uyu, uyd, uyy;
				T j12m, j12u, j12d, j12i;
				T effective_viscosity_y;

				if (levelset_m_x(i, j) > 0)
				{
					uyu = viscosity_p*(water_velocity_field_mac_x->one_over_dy*(u_u(i, j) - u_m(i, j)));
					uyd = viscosity_p*(water_velocity_field_mac_x->one_over_dy*(u_m(i, j) - u_l(i, j)));
						
					if (levelset_d_x(i, j) <= 0)
					{
						T subcell_d = abs(levelset_d_x(i, j))/(abs(levelset_d_x(i, j)) + abs(levelset_m_x(i, j)));
						j12m = (T)0.5*(j12(i - 1, j) + j12(i, j)), j12u = (T)0.5*(j12(i - 1, j + 1) + j12(i, j + 1)), j12d = (T)0.5*(j11(i - 1, j - 1) + j11(i, j - 1));
						j12i = subcell_d*j12m + (1 - subcell_d)*j12d;
						
						effective_viscosity_y = viscosity_p*viscosity_m/(viscosity_p*subcell_d + viscosity_m*(1 - subcell_d));
						
						uyd = effective_viscosity_y*(water_velocity_field_mac_x->one_over_dy*(u_m(i, j) - u_d(i, j))) + effective_viscosity_y*j12i*subcell_d/viscosity_m; 
					}
					if (levelset_u_x(i, j) <= 0)
					{
						T subcell_u = abs(levelset_u_x(i, j))/(abs(levelset_u_x(i, j)) + abs(levelset_m_x(i, j)));
						j12m = (T)0.5*(j12(i - 1, j) + j12(i, j)), j12u = (T)0.5*(j12(i - 1, j + 1) + j12(i, j + 1)), j12d = (T)0.5*(j11(i - 1, j - 1) + j11(i, j - 1));
						j12i = subcell_u*j12m + (1 - subcell_u)*j12u;
						
						effective_viscosity_y = viscosity_p*viscosity_m/(viscosity_m*subcell_u + viscosity_p*(1 - subcell_u));
						
						uyu = effective_viscosity_y*(water_velocity_field_mac_x->one_over_dy*(u_u(i, j) - u_m(i, j))) + effective_viscosity_y*j11i*subcell_u/viscosity_m; 
					}
					uyy = water_velocity_field_mac_x->one_over_dy*(uyu - uyd);
				}
				else
				{
					uyu = viscosity_m*(water_velocity_field_mac_x->one_over_dy*(u_u(i, j) - u_m(i, j)));
					uyd = viscosity_m*(water_velocity_field_mac_x->one_over_dy*(u_m(i, j) - u_d(i, j)));
						
					if (levelset_d_x(i, j) > 0)
					{
						T subcell_d = abs(levelset_d_x(i, j))/(abs(levelset_d_x(i, j)) + abs(levelset_m_x(i, j)));
						j12m = (T)0.5*(j12(i - 1, j) + j12(i, j)), j12u = (T)0.5*(j12(i - 1, j + 1) + j12(i, j + 1)), j12d = (T)0.5*(j11(i - 1, j - 1) + j11(i, j - 1));
						j12i = subcell_d*j12m + (1 - subcell_d)*j12d;
						
						effective_viscosity_y = viscosity_p*viscosity_m/(viscosity_m*subcell_d + viscosity_p*(1 - subcell_d));
						
						uyd = effective_viscosity_y*(water_velocity_field_mac_x->one_over_dy*(u_m(i, j) - u_d(i, j))) - effective_viscosity_y*j12i*subcell_d/viscosity_p; 
					}
					if (levelset_u_x(i, j) > 0)
					{
						T subcell_u = abs(levelset_u_x(i, j))/(abs(levelset_u_x(i, j)) + abs(levelset_m_x(i, j)));
						j12m = (T)0.5*(j12(i - 1, j) + j12(i, j)), j12u = (T)0.5*(j12(i - 1, j + 1) + j12(i, j + 1)), j12d = (T)0.5*(j11(i - 1, j - 1) + j11(i, j - 1));
						j12i = subcell_u*j12m + (1 - subcell_u)*j12u;
						
						effective_viscosity_y = viscosity_p*viscosity_m/(viscosity_m*subcell_u + viscosity_p*(1 - subcell_u));
						
						uyu = effective_viscosity_y*(water_velocity_field_mac_x->one_over_dy*(u_u(i, j) - u_m(i, j))) - effective_viscosity_y*j12i*subcell_u/viscosity_p; 
					}
					uyy = water_velocity_field_mac_x->one_over_dy*(uyu - uyd);
				}

				if (levelset_m_x(i, j) <= 0)
				{
					T one_over_density = (T)1/density_m;
					water_velocity_field_mac_x->array_for_this(i, j) += dt*one_over_density*(uxx + uyy);
				}
				else if (levelset_m_x(i, j) > 0)
				{
					T one_over_density = (T)1/density_p;
					water_velocity_field_mac_x->array_for_this(i, j) += dt*one_over_density*(uxx + uyy);
				}
			}
			
			GRID_ITERATION_2D(water_velocity_field_mac_y->grid)
			{
				T vxl, vxr, vxx;
				T j21m, j21r, j21l, j21i;
				T effective_viscosity_x;

				if (levelset_m_y(i, j) > 0)
				{
					vxl = viscosity_p*(water_velocity_field_mac_y->one_over_dx*(u_m(i, j) - u_l(i, j)));
					vxr = viscosity_p*(water_velocity_field_mac_y->one_over_dx*(u_r(i, j) - u_m(i, j)));
						
					if (levelset_l_y(i, j) <= 0)
					{
						T subcell_l = abs(levelset_l_y(i, j))/(abs(levelset_l_y(i, j)) + abs(levelset_m_y(i, j)));
						j21m = (T)0.5*(j21(i, j - 1) + j21(i, j)), j21r = (T)0.5*(j21(i + 1, j - 1) + j21(i + 1, j)), j21l = (T)0.5*(j11(i - 1, j - 1) + j11(i - 1, j));
						j21i = subcell_l*j21m + (1 - subcell_l)*j21l;
						
						effective_viscosity_x = viscosity_p*viscosity_m/(viscosity_p*subcell_l + viscosity_m*(1 - subcell_l));
						
						vxl = effective_viscosity_x*(water_velocity_field_mac_y->one_over_dx*(u_m(i, j) - u_l(i, j))) + effective_viscosity_x*j21i*subcell_l/viscosity_m; 
					}
					if (levelset_r_y(i, j) <= 0)
					{
						T subcell_r = abs(levelset_r_y(i, j))/(abs(levelset_r_y(i, j)) + abs(levelset_m_y(i, j)));
						j21m = (T)0.5*(j21(i, j - 1) + j21(i, j)), j21r = (T)0.5*(j21(i + 1, j - 1) + j21(i + 1, j)), j21l = (T)0.5*(j11(i - 1, j - 1) + j11(i - 1, j));
						j21i = subcell_r*j21m + (1 - subcell_r)*j21r;
						
						effective_viscosity_x = viscosity_p*viscosity_m/(viscosity_m*subcell_r + viscosity_p*(1 - subcell_r));
						
						vxr = effective_viscosity_x*(water_velocity_field_mac_y->one_over_dx*(u_r(i, j) - u_m(i, j))) + effective_viscosity_x*j21i*subcell_r/viscosity_m; 
					}
					vxx = water_velocity_field_mac_y->one_over_dx*(vxr - vxl);
				}
				else
				{
					vxl = viscosity_m*(water_velocity_field_mac_y->one_over_dx*(u_m(i, j) - u_l(i, j)));
					vxr = viscosity_m*(water_velocity_field_mac_y->one_over_dx*(u_r(i, j) - u_m(i, j)));
						
					if (levelset_l_y(i, j) > 0)
					{
						T subcell_l = abs(levelset_l_y(i, j))/(abs(levelset_l_y(i, j)) + abs(levelset_m_y(i, j)));
						j21m = (T)0.5*(j21(i, j - 1) + j21(i, j)), j21r = (T)0.5*(j21(i + 1, j - 1) + j21(i + 1, j)), j21l = (T)0.5*(j11(i - 1, j - 1) + j11(i - 1, j));
						j21i = subcell_l*j21m + (1 - subcell_l)*j21l;
						
						effective_viscosity_x = viscosity_p*viscosity_m/(viscosity_m*subcell_l + viscosity_p*(1 - subcell_l));
						
						vxl = effective_viscosity_x*(water_velocity_field_mac_y->one_over_dx*(u_m(i, j) - u_l(i, j))) - effective_viscosity_x*j21i*subcell_l/viscosity_p; 
					}
					if (levelset_r_y(i, j) > 0)
					{
						T subcell_r = abs(levelset_r_y(i, j))/(abs(levelset_r_y(i, j)) + abs(levelset_m_y(i, j)));
						j21m = (T)0.5*(j21(i, j - 1) + j21(i, j)), j21r = (T)0.5*(j21(i + 1, j - 1) + j21(i + 1, j)), j21l = (T)0.5*(j11(i - 1, j - 1) + j11(i - 1, j));
						j21i = subcell_r*j21m + (1 - subcell_r)*j21r;
						
						effective_viscosity_x = viscosity_p*viscosity_m/(viscosity_m*subcell_r + viscosity_p*(1 - subcell_r));
						
						vxr = effective_viscosity_x*(water_velocity_field_mac_y->one_over_dx*(u_r(i, j) - u_m(i, j))) - effective_viscosity_x*j21i*subcell_r/viscosity_p; 
					}
					vxx = water_velocity_field_mac_y->one_over_dx*(vxr - vxl);
				}

				T vyu, vyd, vyy;
				T j22m, j22u, j22d, j22i;
				T effective_viscosity_y;

				if (levelset_m_y(i, j) > 0)
				{
					vyu = viscosity_p*(water_velocity_field_mac_y->one_over_dy*(u_u(i, j) - u_m(i, j)));
					vyd = viscosity_p*(water_velocity_field_mac_y->one_over_dy*(u_m(i, j) - u_l(i, j)));
						
					if (levelset_d_y(i, j) <= 0)
					{
						T subcell_d = abs(levelset_d_y(i, j))/(abs(levelset_d_y(i, j)) + abs(levelset_m_y(i, j)));
						j22m = (T)0.5*(j22(i, j - 1) + j12(i, j)), j22u = (T)0.5*(j22(i, j) + j22(i, j + 1)), j22d = (T)0.5*(j22(i, j - 2) + j22(i, j - 1));
						j22i = subcell_d*j22m + (1 - subcell_d)*j22d;
						
						effective_viscosity_y = viscosity_p*viscosity_m/(viscosity_p*subcell_d + viscosity_m*(1 - subcell_d));
						
						vyd = effective_viscosity_y*(water_velocity_field_mac_y->one_over_dy*(u_m(i, j) - u_d(i, j))) + effective_viscosity_y*j22i*subcell_d/viscosity_m; 
					}
					if (levelset_u_y(i, j) <= 0)
					{
						T subcell_u = abs(levelset_u_y(i, j))/(abs(levelset_u_y(i, j)) + abs(levelset_m_y(i, j)));
						j22m = (T)0.5*(j22(i, j - 1) + j12(i, j)), j22u = (T)0.5*(j22(i, j) + j22(i, j + 1)), j22d = (T)0.5*(j22(i, j - 2) + j22(i, j - 1));
						j22i = subcell_u*j22m + (1 - subcell_u)*j22u;
						
						effective_viscosity_y = viscosity_p*viscosity_m/(viscosity_m*subcell_u + viscosity_p*(1 - subcell_u));
						
						vyu = effective_viscosity_y*(water_velocity_field_mac_y->one_over_dy*(u_u(i, j) - u_m(i, j))) + effective_viscosity_y*j22i*subcell_u/viscosity_m; 
					}
					vyy = water_velocity_field_mac_y->one_over_dy*(vyu - vyd);
				}
				else
				{
					vyu = viscosity_m*(water_velocity_field_mac_y->one_over_dy*(u_u(i, j) - u_m(i, j)));
					vyd = viscosity_m*(water_velocity_field_mac_y->one_over_dy*(u_m(i, j) - u_d(i, j)));
					
					if (levelset_d_y(i, j) > 0)
					{
						T subcell_d = abs(levelset_d_y(i, j))/(abs(levelset_d_y(i, j)) + abs(levelset_m_y(i, j)));
						j22m = (T)0.5*(j22(i, j - 1) + j12(i, j)), j22u = (T)0.5*(j22(i, j) + j22(i, j + 1)), j22d = (T)0.5*(j22(i, j - 2) + j22(i, j - 1));
						j22i = subcell_d*j22m + (1 - subcell_d)*j22d;
						
						effective_viscosity_y = viscosity_p*viscosity_m/(viscosity_m*subcell_d + viscosity_p*(1 - subcell_d));
						
						vyd = effective_viscosity_y*(water_velocity_field_mac_y->one_over_dy*(u_m(i, j) - u_d(i, j))) - effective_viscosity_y*j22i*subcell_d/viscosity_p; 
					}
					else if (levelset_u_y(i, j) > 0)
					{
						T subcell_u = abs(levelset_u_y(i, j))/(abs(levelset_u_y(i, j)) + abs(levelset_m_y(i, j)));
						j22m = (T)0.5*(j22(i, j - 1) + j12(i, j)), j22u = (T)0.5*(j22(i, j) + j22(i, j + 1)), j22d = (T)0.5*(j22(i, j - 2) + j22(i, j - 1));
						j22i = subcell_u*j22m + (1 - subcell_u)*j22u;
						
						effective_viscosity_y = viscosity_p*viscosity_m/(viscosity_m*subcell_u + viscosity_p*(1 - subcell_u));
						
						vyu = effective_viscosity_y*(water_velocity_field_mac_y->one_over_dy*(u_u(i, j) - u_m(i, j))) - effective_viscosity_y*j22i*subcell_u/viscosity_p; 
					}
					vyy = water_velocity_field_mac_y->one_over_dy*(vyu - vyd);
				}

				if (levelset_m_y(i, j) <= 0)
				{
					T one_over_density = (T)1/density_m;
					water_velocity_field_mac_y->array_for_this(i, j) += dt*one_over_density*(vxx + vyy);
				}
				else if (levelset_m_y(i, j) > 0)
				{
					T one_over_density = (T)1/density_p;
					water_velocity_field_mac_y->array_for_this(i, j) += dt*one_over_density*(vxx + vyy);
				}
			}
		}
	}

	void ApplyGraivity(const T& dt)
	{
		const VT gravity_dt(dt*gravity);
		
		if (air_water_simulation)
		{
			ARRAY_2D<T>& water_velocity_x(water_velocity_field_mac_x->array_for_this);
			ARRAY_2D<T>& water_velocity_y(water_velocity_field_mac_y->array_for_this);

			int i, j;
			LOOPS_2D(i, j, water_velocity_field_mac_x->i_start, water_velocity_field_mac_x->j_start, water_velocity_field_mac_x->i_end, water_velocity_field_mac_x->j_end)
			{
				water_velocity_x(i, j) += gravity_dt.x;
			}	
			LOOPS_2D(i, j, water_velocity_field_mac_y->i_start, water_velocity_field_mac_y->j_start, water_velocity_field_mac_y->i_end, water_velocity_field_mac_y->j_end)
			{
				water_velocity_y(i, j) += gravity_dt.y;
			} 
		}
		if (oil_water_simulation)
		{
			if (dimensionless_form)
			{
				if (is_vertical)
				{
					ARRAY_2D<T>& water_velocity_x(water_velocity_field_mac_x->array_for_this);
					ARRAY_2D<T>& water_velocity_y(water_velocity_field_mac_y->array_for_this);

					int i, j;
					LOOPS_2D(i, j, water_velocity_field_mac_x->i_start, water_velocity_field_mac_x->j_start, water_velocity_field_mac_x->i_end, water_velocity_field_mac_x->j_end)
					{
						water_velocity_x(i, j) += R1/POW2(V0)*gravity_dt.x;
					}	
					LOOPS_2D(i, j, water_velocity_field_mac_y->i_start, water_velocity_field_mac_y->j_start, water_velocity_field_mac_y->i_end, water_velocity_field_mac_y->j_end)
					{
						water_velocity_y(i, j) += R1/POW2(V0)*gravity_dt.y;
						// Driving Pressure
						water_velocity_y(i, j) += dt*R1/(POW2(V0)*oil_density)*f;
					}   
				}
			}
		}

		SetupBoundaryConditionForVelocity();
	}

	void ApplyGraivity(const T& dt, const int& thread_id)
	{
		const VT gravity_dt(dt*gravity);
		
		if (air_water_simulation)
		{
			ARRAY_2D<T>& water_velocity_x(water_velocity_field_mac_x->array_for_this);
			ARRAY_2D<T>& water_velocity_y(water_velocity_field_mac_y->array_for_this);

			BEGIN_GRID_ITERATION_2D(water_velocity_field_mac_x->partial_grids[thread_id])
			{
				water_velocity_x(i, j) += gravity_dt.x;
			}	
			END_GRID_ITERATION_2D;
			
			BEGIN_GRID_ITERATION_2D(water_velocity_field_mac_y->partial_grids[thread_id])
			{
				water_velocity_y(i, j) += gravity_dt.y;
			} 
			END_GRID_ITERATION_2D;
		}
		if (oil_water_simulation)
		{
			if (dimensionless_form)
			{
				if (is_vertical)
				{
					ARRAY_2D<T>& water_velocity_x(water_velocity_field_mac_x->array_for_this);
					ARRAY_2D<T>& water_velocity_y(water_velocity_field_mac_y->array_for_this);

					BEGIN_GRID_ITERATION_2D(water_velocity_field_mac_x->partial_grids[thread_id])
					{
						water_velocity_x(i, j) += R1/POW2(V0)*gravity_dt.x;
					}
					END_GRID_ITERATION_2D;

					BEGIN_GRID_ITERATION_2D(water_velocity_field_mac_y->partial_grids[thread_id])
					{
						water_velocity_y(i, j) += R1/POW2(V0)*gravity_dt.y;
						// Driving Pressure
						water_velocity_y(i, j) += dt*R1/(POW2(V0)*oil_density)*f;
					}   
					END_GRID_ITERATION_2D;
				}
			}
		}
	}

	void ApplySurfaceTension(const T& dt)
	{
		if (air_water_simulation)
		{
			
			FIELD_STRUCTURE_2D<T> density_x, density_y, normal_x, normal_y;
			density_x.Initialize(water_velocity_field_mac_x->grid, 2);
			density_y.Initialize(water_velocity_field_mac_y->grid, 2);
			normal_x.Initialize(water_velocity_field_mac_x->grid, 2);
			normal_y.Initialize(water_velocity_field_mac_y->grid, 2);
			
			GRID_ITERATION_2D(water_velocity_field_mac_x->grid)
			{
				normal_x(i, j) = (water_levelset->phi(i, j) - water_levelset->phi(i - 1, j))*water_levelset->grid.one_over_dx;	
			}
	
			GRID_ITERATION_2D(water_velocity_field_mac_y->grid)
			{
				normal_y(i, j) = (water_levelset->phi(i, j) - water_levelset->phi(i, j - 1))*water_levelset->grid.one_over_dy;
			}
			
			FIELD_STRUCTURE_2D<T> delta_c, delta_l, delta_d;
			delta_c.Initialize(water_levelset->grid, 2);
			delta_l.Initialize(water_levelset->grid, 2);
			delta_d.Initialize(water_levelset->grid, 2);

			FIELD_STRUCTURE_2D<T> phi_c, phi_l, phi_d, phi_x, phi_y;
			phi_c.Initialize(water_levelset->grid, 2);
			phi_l.Initialize(water_levelset->grid, 2);
			phi_d.Initialize(water_levelset->grid, 2);
			phi_x.Initialize(water_velocity_field_mac_x->grid, 2);
			phi_y.Initialize(water_velocity_field_mac_y->grid, 2);
			
			GRID_ITERATION_2D(water_velocity_field_mac_x->grid)
			{
				phi_x(i, j) = (T)0.5*(water_levelset->arr(i, j) + water_levelset->arr(i - 1, j));
			}

			GRID_ITERATION_2D(water_velocity_field_mac_y->grid)
			{
				phi_y(i, j) = (T)0.5*(water_levelset->arr(i, j) + water_levelset->arr(i, j - 1));
			}

			DetermineDensityField(phi_x, epsilon_for_mollification, density_x);
			DetermineDensityField(phi_y, epsilon_for_mollification, density_y);

			water_levelset->FillGhostCellsContinuousDerivativesFrom(water_levelset->arr, true);

			GRID_ITERATION_2D(water_levelset->grid)
			{
				phi_c(i, j) = water_levelset->arr(i, j);
				phi_l(i, j) = water_levelset->arr(i - 1, j);
				phi_d(i, j) = water_levelset->arr(i, j - 1);
			}
			
			DeltaFunction(phi_c, epsilon_for_mollification, delta_c);
			DeltaFunction(phi_l, epsilon_for_mollification, delta_l);
			DeltaFunction(phi_d, epsilon_for_mollification, delta_d);
	
			water_levelset->ComputeCurvatures();
	
			water_levelset->curvature.FillGhostCellsContinuousDerivativesFrom(water_levelset->curvature.array_for_this, false);
	
			FIELD_STRUCTURE_2D<T> curvature_x, curvature_y;
			curvature_x.Initialize(water_velocity_field_mac_x->grid, 2);
			curvature_y.Initialize(water_velocity_field_mac_y->grid, 2);
	
			GRID_ITERATION_2D(curvature_x.grid)
			{
				curvature_x(i, j) = (T)0.5*(water_levelset->curvature(i, j) + water_levelset->curvature(i - 1, j));
			}
	
			GRID_ITERATION_2D(curvature_y.grid)
			{
				curvature_y(i, j) = (T)0.5*(water_levelset->curvature(i, j) + water_levelset->curvature(i, j - 1));
			}
	
			GRID_ITERATION_2D(water_velocity_field_mac_x->grid)
			{
				T delta = (T)0.5*(delta_c(i, j) + delta_l(i, j));
				T surf_x = dt*surface_tension*curvature_x(i, j)*delta*normal_x(i, j)/density_x(i, j);

				water_velocity_field_mac_x->array_for_this(i, j) = water_velocity_field_mac_x->array_for_this(i, j) + surf_x;
			}
			
			GRID_ITERATION_2D(water_velocity_field_mac_y->grid)
			{
				T delta = (T)0.5*(delta_c(i, j) + delta_d(i, j));

				water_velocity_field_mac_y->array_for_this(i, j) = water_velocity_field_mac_y->array_for_this(i, j) + dt*surface_tension*curvature_y(i, j)*delta*normal_y(i, j)/density_y(i, j); 
			}
		}
		if (oil_water_simulation)
		{
			T one_over_weber = 1/(T)We;

			FIELD_STRUCTURE_2D<T> curv_r, delta_c, delta_l, density_c, density_l, phi_c, phi_l;
			curv_r.Initialize(water_velocity_field_mac_x->grid, 2);
			phi_c.Initialize(water_levelset->grid, 2);
			phi_l.Initialize(water_levelset->grid, 2);
			delta_c.Initialize(water_levelset->grid, 2);
			delta_l.Initialize(water_levelset->grid, 2);
			density_c.Initialize(water_levelset->grid, 2);
			density_l.Initialize(water_levelset->grid, 2);

			GRID_ITERATION_2D(water_velocity_field_mac_x->grid)
			{
				curv_r(i, j) = (T)0.5*(water_levelset->curvature(i, j) + water_levelset->curvature(i - 1, j));
			}

			GRID_ITERATION_2D(water_levelset->grid)
			{
				phi_c(i, j) = water_levelset->arr(i, j);
				phi_l(i, j) = water_levelset->arr(i - 1, j);
				density_c(i, j) = density_field(i, j);
				density_l(i, j) = density_field(i - 1, j);
			}

			DeltaFunction(phi_c, epsilon_for_mollification, delta_c);
			DeltaFunction(phi_l, epsilon_for_mollification, delta_l);

			GRID_ITERATION_2D(water_velocity_field_mac_x->grid)
			{
				T delta = (T)0.5*(delta_c(i, j) + delta_l(i, j));
				T density = (T)0.5*(density_c(i, j) + density_l(i, j));
				
				T surf = one_over_weber*curv_r(i, j)*delta*base_grid.one_over_dx*(phi_c(i, j) - phi_l(i, j))/density;

				water_velocity_field_mac_x->array_for_this(i, j) += dt*surf;
			}

			FIELD_STRUCTURE_2D<T> curv_x, phi_d, density_d, delta_d;
			curv_x.Initialize(water_velocity_field_mac_y->grid, 2);
			phi_d.Initialize(water_levelset->grid, 2);
			density_d.Initialize(water_levelset->grid, 2);
			delta_d.Initialize(water_levelset->grid, 2);

			GRID_ITERATION_2D(water_velocity_field_mac_y->grid)
			{
				curv_x(i, j) = (T)0.5*(water_levelset->curvature(i, j) + water_levelset->curvature(i, j - 1));
			}

			GRID_ITERATION_2D(water_levelset->grid)
			{
				phi_d(i, j) = water_levelset->arr(i, j - 1);
				density_d(i, j) = density_field(i, j - 1);
			}

			DeltaFunction(phi_d, epsilon_for_mollification, delta_d);

			GRID_ITERATION_2D(water_velocity_field_mac_y->grid)
			{
				T delta = (T)0.5*(delta_c(i, j) + delta_d(i, j));
				T density = (T)0.5*(density_c(i, j) + density_d(i, j));

				T surf= one_over_weber*curv_x(i, j)*delta*base_grid.one_over_dx*(phi_c(i, j) - phi_d(i, j))/density;

				water_velocity_field_mac_y->array_for_this(i, j) += dt*surf;
			}
		}	
	}
	
public: // Functions related to the time
	T CFLOneTimeStep(void)
	{
		T u_max(0), v_max(0);
		
		int i(0), j(0);
		GRID_ITERATION_2D(water_velocity_field_mac_x->grid)
		{
			if (u_max <= abs(water_velocity_field_mac_x->array_for_this(i, j)))
			{
				u_max = abs(water_velocity_field_mac_x->array_for_this(i, j));
			}
		}
		GRID_ITERATION_2D(water_velocity_field_mac_y->grid)
		{
			if (v_max <= abs(water_velocity_field_mac_y->array_for_this(i, j)))
			{
				v_max = abs(water_velocity_field_mac_y->array_for_this(i, j));
			}
		}
		
		// CFL condition For Advection Term
		c_f = (u_max*water_velocity_field_mac_x->grid.one_over_dx + v_max*water_velocity_field_mac_y->grid.one_over_dy);
		
		g_f = 0;
		
		// CFL condition For Gravity Term
		if (air_water_simulation)
		{
			g_f = sqrt((T)9.8*base_grid.one_over_dy);
		}
		if (oil_water_simulation)
		{
			if (dimensionless_form)
			{
				g_f = sqrt(abs(R1/POW2(V0)*gravity.y)*water_velocity_field_mac_y->one_over_dy);
			}
			else
			{
				g_f = sqrt(abs(gravity.y)*water_velocity_field_mac_y->one_over_dy);
			}
		}

		// CFL condition For Viscosity Term
		v_f = 0;
		if (air_water_simulation)
		{
			v_f = max(water_viscosity/water_density, air_viscosity/air_density)*((T)2*base_grid.one_over_dx2 + (T)2*base_grid.one_over_dy2);
		}
		if (oil_water_simulation)
		{
			if (viscosity_solver->semi_implicit_approach)
			{
				v_f = 0;
			}
			else
			{
				v_f = max(water_viscosity/water_density, oil_viscosity/oil_density)*((T)2*base_grid.one_over_dx2 + (T)2*base_grid.one_over_dy2);
			}
		}

		// CFL condition for Curvature Term
		T max_abs_cur(0);
		GRID_ITERATION_2D(water_levelset->curvature.grid)
		{
			if (abs(water_levelset->curvature.array_for_this(i, j)) >= max_abs_cur )
			{
				max_abs_cur = abs(water_levelset->curvature.array_for_this(i, j));
			}
		}

		s_f = 0;
		if (air_water_simulation)
		{
			s_f = sqrt(surface_tension*max_abs_cur/(min(water_density, air_density)*(POW2(min(base_grid.dx, base_grid.dy)))));
		}
		if (oil_water_simulation)
		{
			if (dimensionless_form)
			{
				s_f = sqrt((R1/We*max_abs_cur)/(MIN(water_density, oil_density)*POW2(MIN(water_velocity_field_mac_x->grid.dx, water_velocity_field_mac_y->grid.dy))));
			}
			else
			{
				s_f = sqrt((surface_tension*max_abs_cur)/(MIN(water_density, oil_density)*POW2(MIN(water_velocity_field_mac_x->grid.dx, water_velocity_field_mac_y->grid.dy))));
			}	
		}
		
		T CFL = (T)0.5*(c_f + v_f + sqrt(POW2(c_f + v_f) + (T)4*POW2(g_f) + (T)4*POW2(s_f)));

		T CFL_dt = (T)cfl_number/CFL;

		return CFL_dt;
	}	
};

	
		
		


	

