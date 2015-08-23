#pragma once

#include "EULERIAN_FLUID_SOLVER_2D.h"
#include "NUMERICAL_INTEGRATION.h"
#include "MAXIMUM_ENTROPY_ANALYSIS.h"
#include "WORLD_DISCRETIZATION_2D.h"
#include "POISSON_EQUATION_TEST.h"
#include "MONGE_AMPERE_SOLVER.h"
#include <ctype.h>
#include <boost/chrono.hpp>

class SIMULATION_WORLD
{
public: // Primitive Solvers which are updated
	EULERIAN_FLUID_SOLVER_2D		eulerian_solver;
	NUMERICAL_INTEGRATION			numerical_integration;
    POISSON_EQUATION_TEST			poisson_equation_test;
	MONGE_AMPERE_SOLVER				monge_ampere_solver;
	// Temporary
    MAXIMUM_ENTROPY_ANALYSIS        maximum_entropy_analysis;

public: 
	WORLD_DISCRETIZATION_2D			world_discretization;

public: // Essential for script update
	std::string						script_abs_dir;
	std::string						app_abs_dir;
	std::string						script_file_path;
	std::string						script_base_name;

public: // Simulation Properties
	T								frame_rate;
	T								dt, accu_dt, max_dt;
	T								CFL;
	int								num_current_frame;
	
public: // Options for Simulation
	bool							fluid_solver;
	bool							numerical_test_solver;

	// Fluid Solver Options
	bool							air_water_simulation;
	bool							oil_water_simulation;
	bool							vortex_sheet_problem;
	
	// Numerical Test Options
	bool							numerical_integration_test;
    bool                            BJ_model_test;
    bool                            signal_processing_test;
    bool							poisson_equation_with_jump_condition;
	bool							monge_ampere_solver_test;

	// Option for Air-Water Simulation
	bool							large_bubble, small_bubble;

	// Option
	int								last_frame;
	bool							auto_run;
	bool							is_polygonize_levelset;

	// Option for Poisson Solver Test
	int								test_number;

	// Grid stuff
	T								grid_scale;

public: // Options for Oil Water Simulation
	bool							is_vertical, is_parallel;

public: // Wave length and Wave number for Oil Water Simulation
	T								A_0, alpha;

public: // Multithreading
	MULTITHREADING*					multithreading;

public: // Constructor and Destructor
	SIMULATION_WORLD(void)
		: multithreading(0), fluid_solver(false), numerical_test_solver(false), air_water_simulation(false), large_bubble(false), small_bubble(false), oil_water_simulation(false), numerical_integration_test(false), poisson_equation_with_jump_condition(false), is_vertical(false), is_parallel(false), num_current_frame(0), last_frame((int)1000000), dt((T)0), max_dt((T)0), test_number(0)
	{}

	~SIMULATION_WORLD(void)
	{}

public: // Initialization Functions
	void InitializeFromScript(const char* script)
	{
		DELETE_POINTER(multithreading);
		multithreading = new MULTITHREADING;

		SetCurrentDirectory(script_abs_dir.c_str());

		SCRIPT_READER script_reader(script);

		num_current_frame = 0;

		SCRIPT_BLOCK script_block_for_this = script_reader.FindBlock("SIMULATION_WORLD");

		dt = script_block_for_this.GetFloat("dt", (T)0.01);
		max_dt = script_block_for_this.GetFloat("max_dt", (T)100);
		CFL = script_block_for_this.GetFloat("CFL", (T)0);
		frame_rate = script_block_for_this.GetFloat("frame_rate", (T)24);
		accu_dt = (T)0;

		// Simulation Options
		fluid_solver = script_block_for_this.GetBoolean("fluid_solver", false);
		numerical_test_solver = script_block_for_this.GetBoolean("numerical_test_solver", false);
		
		if (fluid_solver)
		{
			air_water_simulation = script_block_for_this.FindBlock("FLUID_SOLVER_OPTIONS").GetBoolean("air_water_simulation", false);
			oil_water_simulation = script_block_for_this.FindBlock("FLUID_SOLVER_OPTIONS").GetBoolean("oil_water_simulation", false);
			vortex_sheet_problem = script_block_for_this.FindBlock("FLUID_SOLVER_OPTIONS").GetBoolean("vortex_sheet_problem", false);
		}
		if (numerical_test_solver)
		{
			numerical_integration_test = script_block_for_this.FindBlock("NUMERICAL_TEST_OPTIONS").GetBoolean("numerical_integration_test", false);
			poisson_equation_with_jump_condition = script_block_for_this.FindBlock("NUMERICAL_TEST_OPTIONS").GetBoolean("poisson_equation_with_jump_condition", false);
			signal_processing_test = script_block_for_this.FindBlock("NUMERICAL_TEST_OPTIONS").GetBoolean("signal_processing_test", false);
			monge_ampere_solver_test = script_block_for_this.FindBlock("NUMERICAL_TEST_OPTIONS").GetBoolean("monge_ampere_solver_test", false);

			if (numerical_integration_test)
			{
                BJ_model_test = script_block_for_this.FindBlock("NUMERICAL_TEST_OPTIONS").FindBlock("NUMERICAL_INTEGRATION_OPTION").GetBoolean("BJ_model_test", false);    
			}
			if (poisson_equation_with_jump_condition)
			{
				test_number = script_block_for_this.FindBlock("NUMERICAL_TEST_OPTIONS").FindBlock("POISSON_EQUATION_WITH_JUMP_CONDITION_TEST_NUMBER").GetInteger("test_number", (int)1);
			} // Test number 1~2 gives 1d example, and the others will be 2d
			if (monge_ampere_solver_test)
			{
				test_number = script_block_for_this.FindBlock("NUMERICAL_TEST_OPTIONS").FindBlock("MONGE_AMPERE_EQUATION_TEST_NUMBER").GetInteger("test_number", (int)1);
			}
		}
		
		// Multithreading
		multithreading->Initialize(script_block_for_this.GetInteger("number_of_threads"));

		last_frame = script_block_for_this.GetInteger("last_frame");
		auto_run = script_block_for_this.GetBoolean("auto_run");

		// Display
		cout << "---------------SIMULATION WORLD VARIABLES---------------" << endl; 
		if (fluid_solver)
		{
			cout << "dt = " << dt << endl;
			cout << "CFL: " << CFL << endl;
			cout << "frame rate: " << frame_rate << endl;

			cout << "               <Fluid Simulation>          " << endl;
			if (air_water_simulation)
			{
				cout << "Air-Water Simulation is activated!" << endl;
			}
			if (oil_water_simulation)
			{
				cout << "Oil-Water Simulation is activated!" << endl;
			}
			if (vortex_sheet_problem)
			{
				cout << "Vortex Sheet Problem is activated!" << endl;
			} 
		}
		if (numerical_test_solver)
		{
			cout << "             <Numerical Test>          " << endl;
			if (numerical_integration_test)
			{
				cout << "Numerical Integration Test is activated!" << endl;
				
				if (BJ_model_test)
				{
                    cout << "BJ model Test is activated!" << endl;
				}
			}
			if (poisson_equation_with_jump_condition)
			{
				cout << "Poisson Equation with Jump Condition Test is activated!" << endl;
				cout << "Test Number : #" << test_number << endl;
			}
			if (signal_processing_test)
			{
				cout << "Signal Processing Test is activated!" << endl;
			}
			if (monge_ampere_solver_test)
			{
				cout << "Monge-Ampere Solver Test is activated!" << endl;
				cout << "Test Number : #" << test_number << endl;
			}
		}

		if (fluid_solver)
		{
			if (oil_water_simulation)
			{
				A_0 = script_block_for_this.GetFloat("A_0", (T)1);
				alpha = script_block_for_this.GetFloat("alpha", (T)1);
			} 
		}

		cout << "Number of threads: " << multithreading->num_threads << endl;

		// Initialize world discretization
		if (fluid_solver)
		{
			// Falsify the numerical test options
			world_discretization.poisson_equation_with_jump_condition_test = false;
			
			if (air_water_simulation)
			{
				world_discretization.air_water_simulation = true;
				world_discretization.oil_water_simulation = false;
				world_discretization.Initialize(script_block_for_this.FindBlock("WORLD_DISCRETIZATION").FindBlock("AIR_WATER_SIMULATION"));
			}
			if (oil_water_simulation)
			{
				world_discretization.air_water_simulation = false;
				world_discretization.oil_water_simulation = true;
				world_discretization.Initialize(script_block_for_this.FindBlock("WORLD_DISCRETIZATION").FindBlock("OIL_WATER_SIMULATION"));
				is_vertical = world_discretization.is_vertical;
				is_parallel = world_discretization.is_parallel;
			}
			if (vortex_sheet_problem)
			{
				world_discretization.Initialize(script_block_for_this.FindBlock("WORLD_DISCRETIZATION").FindBlock("VORTEX_SHEET_PROBLEM"));
			} 
		}
		if (numerical_test_solver)
		{
			// Falsify the fluid solver options
			world_discretization.air_water_simulation = false;
			world_discretization.oil_water_simulation = false;

			if (numerical_integration_test)
			{
				world_discretization.poisson_equation_with_jump_condition_test = false;

				world_discretization.Initialize(script_block_for_this.FindBlock("WORLD_DISCRETIZATION").FindBlock("NUMERICAL_INTEGRATION_TEST"));
			} 
			if (poisson_equation_with_jump_condition)
			{
				world_discretization.poisson_equation_with_jump_condition_test = true;
				
				if (test_number == 1 || test_number == 2)
				{
					world_discretization.grid_1d = true;
					world_discretization.grid_2d = false;
				}
				else
				{
					world_discretization.grid_1d = false;
					world_discretization.grid_2d = true;	
				}

				world_discretization.Initialize(script_block_for_this.FindBlock("WORLD_DISCRETIZATION").FindBlock("POISSON_EQUATION_WITH_JUMP_CONDITION_TEST"));
			}
			if (monge_ampere_solver_test)
			{
				world_discretization.poisson_equation_with_jump_condition_test = false;

				world_discretization.Initialize(script_block_for_this.FindBlock("WORLD_DISCRETIZATION").FindBlock("MONGE_AMPERE_EQUATION_TEST"));
			}
			if (signal_processing_test)
			{
				world_discretization.poisson_equation_with_jump_condition_test = false;

				maximum_entropy_analysis.Nsamps = script_block_for_this.FindBlock("WORLD_DISCRETIZATION").FindBlock("MAXIMUM_ENTROPY_METHOD_TEST").GetInteger("number_of_samples", (int)1000);
				maximum_entropy_analysis.min_f = script_block_for_this.FindBlock("WORLD_DISCRETIZATION").FindBlock("MAXIMUM_ENTROPY_METHOD_TEST").GetFloat("minimum_frequency", (float)2.0)/maximum_entropy_analysis.Nsamps;
				maximum_entropy_analysis.max_f = script_block_for_this.FindBlock("WORLD_DISCRETIZATION").FindBlock("MAXIMUM_ENTROPY_METHOD_TEST").GetFloat("maximum_frequency", (float)0.5);
				maximum_entropy_analysis.ihr = script_block_for_this.FindBlock("WORLD_DISCRETIZATION").FindBlock("MAXIMUM_ENTROPY_METHOD_TEST").GetInteger("ihr", (int)1);
				maximum_entropy_analysis.nihr = maximum_entropy_analysis.Nsamps*maximum_entropy_analysis.ihr;

				world_discretization.Initialize(maximum_entropy_analysis.nihr, 0, maximum_entropy_analysis.min_f, maximum_entropy_analysis.max_f);
			}
		}

		// Option for Air-Water
		if (fluid_solver)
		{
			if (air_water_simulation)
			{
				large_bubble = world_discretization.large_bubble;
				small_bubble = world_discretization.small_bubble;
			} 
		}

		InitializeNumericalSolversFromScript(script_reader);
		InitializeObjectListFromScript(script_reader);

		grid_scale = 1;

		SetCurrentDirectory(app_abs_dir.c_str());
	}

	void InitializeNumericalSolversFromScript(SCRIPT_READER& script_reader)
	{
		if (fluid_solver)
		{
			eulerian_solver.world_discretization = &world_discretization;
		}
		if (numerical_test_solver)
		{
			if (numerical_integration_test)
			{
				numerical_integration.world_discretization = &world_discretization;
			}
			if (poisson_equation_with_jump_condition)
			{
				poisson_equation_test.world_discretization = &world_discretization;
			}
			if (monge_ampere_solver_test)
			{
				monge_ampere_solver.world_discretization = &world_discretization;
			}
			if (signal_processing_test)
			{
				maximum_entropy_analysis.world_discretization = &world_discretization;
			}
		}
		
		if (fluid_solver)
		{
			if (air_water_simulation)
			{
				eulerian_solver.air_water_simulation = air_water_simulation;
				eulerian_solver.oil_water_simulation = oil_water_simulation;
				eulerian_solver.InitializeFromScriptBlock(world_discretization.world_grid, script_reader.FindBlock("FLUID_SOLVER_UNIFORM").FindBlock("AIR_WATER_SIMULATION"), multithreading);
			}
			if (oil_water_simulation)
			{
				eulerian_solver.air_water_simulation = air_water_simulation;
				eulerian_solver.oil_water_simulation = oil_water_simulation;
				eulerian_solver.InitializeFromScriptBlock(world_discretization.world_grid, script_reader.FindBlock("FLUID_SOLVER_UNIFORM").FindBlock("OIL_WATER_SIMULATION"), multithreading);
			}
			if (vortex_sheet_problem)
			{
				eulerian_solver.vortex_sheet_problem = vortex_sheet_problem;
				eulerian_solver.InitializeFromScriptBlock(world_discretization.world_grid, script_reader.FindBlock("FLUID_SOLVER_UNIFORM").FindBlock("VORTEX_SHEET_PROBLEM"), multithreading);
			} 
		}
		if (numerical_test_solver)
		{
			if (numerical_integration_test)
			{
				numerical_integration.InitializeFromScriptBlock(world_discretization.world_grid, script_reader.FindBlock("NUMERICAL_INTEGRATION_TEST"), multithreading);
			}
			if (signal_processing_test)
			{
                maximum_entropy_analysis.InitializeFromScriptBlock(script_reader.FindBlock("SIGNAL_PROCESSING_TEST"));
			}
			if (signal_processing_test)
			{
				cout << "bins_per_time_unit : " << maximum_entropy_analysis.bptu << endl;
				cout << "minimum_filter_length : " << maximum_entropy_analysis.mfl << endl;
				cout << "number_of_sample : " << maximum_entropy_analysis.Nsamps << endl;
				cout << "number_of_coefficients : " << maximum_entropy_analysis.Ncoef << endl;
			}
			if (poisson_equation_with_jump_condition)
			{
				poisson_equation_test.test_number = test_number;
				poisson_equation_test.grid_1d = world_discretization.grid_1d;
				poisson_equation_test.grid_2d = world_discretization.grid_2d;
				poisson_equation_test.InitializeFromBlock(script_reader.FindBlock("POISSON_EQUATION_TEST"), multithreading);
			}
			if (monge_ampere_solver_test)
			{
				monge_ampere_solver.test_number = test_number;
				monge_ampere_solver.InitializeFromBlock(script_reader.FindBlock("MONGE_AMPERE_EQUATION_TEST"), multithreading);
			}
		}
	}

	void InitializeObjectListFromScript(SCRIPT_READER& script_reader)
	{
		GRID_STRUCTURE_2D& base_grid = world_discretization.world_grid;

		//SCRIPT_BLOCK script_block = script_reader.FindBlock("SIMULATION_WORLD");

		if (fluid_solver)
		{
			if (air_water_simulation)
			{
				LEVELSET_2D& water_levelset = *eulerian_solver.water_levelset;
				GRID_STRUCTURE_2D& water_grid = water_levelset.grid;

				int i(0), j(0);

				if (large_bubble)
				{
					LOOPS_2D(i, j, water_grid.i_start, water_grid.j_start, water_grid.i_end, water_grid.j_end)
					{
						T x_coor = water_grid.x_min + i*water_grid.dx, y_coor = water_grid.y_min + j*water_grid.dy;
						water_levelset(i, j) = sqrt(POW2(x_coor) + POW2(y_coor + (T)0.5)) - (T)1/(T)3;
					}
					water_levelset.FillGhostCellsFromThreaded(&(water_levelset.phi), false);
				}
				if (small_bubble)
				{
					LOOPS_2D(i, j, water_grid.i_start, water_grid.j_start, water_grid.i_end, water_grid.j_end)
					{
						T x_coor = water_grid.x_min + i*water_grid.dx, y_coor = water_grid.y_min + j*water_grid.dy;
						water_levelset(i, j) = sqrt(POW2(x_coor) + POW2(y_coor + (T)0.005)) - (T)1/(T)300;
					}
					water_levelset.FillGhostCellsFromThreaded(&(water_levelset.phi), false);
				}
			} 
		}
		else if (numerical_test_solver)
		{
			if (numerical_integration_test)
			{
				LEVELSET_2D& object_levelset = *numerical_integration.object_levelset;
				GRID_STRUCTURE_2D& object_grid = object_levelset.grid;

				SCRIPT_BLOCK script_block = script_reader.FindBlock("NUMERICAL_INTEGRATION_TEST");
				
				bool test_1 = script_block.GetBoolean("test_1", false);
				bool test_2 = script_block.GetBoolean("test_2", false);
				bool test_3 = script_block.GetBoolean("test_3", false);
				bool test_4 = script_block.GetBoolean("test_4", false);

				int i(0), j(0);

                if (test_1)
				{
					/*T dev_x = (double)rand()/(double)RAND_MAX - (double)0.5;
                    T dev_y = (double)rand()/(double)RAND_MAX - (double)0.5;*/

                    T dev_x = (T)0;
                    T dev_y = (T)0;

                    LOOPS_2D(i, j, object_grid.i_start, object_grid.j_start, object_grid.i_end, object_grid.j_end)
					{
						T x_coor = object_grid.x_min + i*object_grid.dx, y_coor = object_grid.y_min + j*object_grid.dy;
						//object_levelset(i, j) = POW2((x_coor - dev_x)/(T)1.5) + POW2((y_coor - dev_y)/(T)0.75) - (T)1;
                        object_levelset(i, j) = sqrt(POW2(x_coor - dev_x) + POW2(y_coor - dev_y)) - (T)1;
                    }
					object_levelset.FillGhostCellsFromThreaded(&(object_levelset.phi), false); 

                    //numerical_integration.true_integral_value = (T)7.266336165;
                    numerical_integration.true_integral_value = (T)6.283185307179586;
				}
			}

			if (poisson_equation_with_jump_condition)
			{
				if (test_number == 1)
				{
					LEVELSET_1D& interface_levelset = *poisson_equation_test.interface_levelset_1d;
					GRID_STRUCTURE_1D& interface_grid = interface_levelset.grid;

					int i(0);

					LOOPS_1D(i, interface_grid.i_start, interface_grid.i_end)
					{
						T x_coor = interface_grid.x_min + i*interface_grid.dx;

						interface_levelset(i) = x_coor - (T)0.5;
					}

					interface_levelset.FillGhostCellsFromThreaded(&(interface_levelset.phi), false);
				}
				
				if (test_number == 3)
				{
					LEVELSET_2D& interface_levelset = *poisson_equation_test.interface_levelset_2d;
					GRID_STRUCTURE_2D& interface_grid = interface_levelset.grid;

					int i(0), j(0);

					LOOPS_2D(i, j, interface_grid.i_start, interface_grid.j_start, interface_grid.i_end, interface_grid.j_end)
					{
						T x_coor = interface_grid.x_min + i*interface_grid.dx, y_coor = interface_grid.y_min + j*interface_grid.dy;
						
						interface_levelset(i, j) = sqrt(POW2(x_coor - (T)0.5) + POW2(y_coor - (T)0.5)) - (T)0.25;
					}
					interface_levelset.FillGhostCellsFromThreaded(&(interface_levelset.phi), false);
				}

				if (test_number == 5 || test_number == 6 || test_number == 7)
				{
					LEVELSET_2D& interface_levelset = *poisson_equation_test.interface_levelset_2d;
					GRID_STRUCTURE_2D& interface_grid = interface_levelset.grid;

					int i(0), j(0);

					LOOPS_2D(i, j, interface_grid.i_start, interface_grid.j_start, interface_grid.i_end, interface_grid.j_end)
					{
						T x_coor = interface_grid.x_min + i*interface_grid.dx, y_coor = interface_grid.y_min + j*interface_grid.dy;
						
						interface_levelset(i, j) = sqrt(POW2(x_coor) + POW2(y_coor)) - (T)0.5;
					}
					interface_levelset.FillGhostCellsFromThreaded(&(interface_levelset.phi), false);
				}

				if (test_number == 8)
				{
					LEVELSET_2D& interface_levelset = *poisson_equation_test.interface_levelset_2d;
					GRID_STRUCTURE_2D& interface_grid = interface_levelset.grid;

					int N_exact(8000);

					ARRAY<VT> exact_interface(N_exact);

					for (int p = 0; p < N_exact; p++)
					{
						T theta = (T)2*PI*p/N_exact;
						exact_interface[p].x = (T)0.02*sqrt(5) + (0.5 + 0.2*sin(5*theta))*cos(theta);
						exact_interface[p].y = (T)0.02*sqrt(5) + (0.5 + 0.2*sin(5*theta))*sin(theta);
					}

					int i(0), j(0);
					LOOPS_2D(i, j, interface_grid.i_start, interface_grid.j_start, interface_grid.i_end, interface_grid.j_end)
					{
						T x_coor = interface_grid.x_min + i*interface_grid.dx, y_coor = interface_grid.y_min + j*interface_grid.dy;
						
						T min_val = sqrt(POW2(x_coor - exact_interface[0].x) + POW2(y_coor - exact_interface[0].y));
						
						for (int k = 1; k < N_exact; k++)
						{
							T temp = sqrt(POW2(x_coor - exact_interface[k].x) + POW2(y_coor - exact_interface[k].y));
							if (temp <= min_val)
							{
								min_val = temp;
							}
						}

						int p_index = 0;
						bool inside = false;

						/*while (p_index < (int)10000)
						{
							T theta = (T)2*PI*p_index/N_exact;
							if (POW2(x_coor - 0.02*sqrt(5)) + POW2(y_coor - 0.02*sqrt(5)) - POW2(0.5 + 0.2*sin(5*theta)) <= 0)
							{
								inside = true;
								break;
							}
							if (POW2(x_coor - 0.02*sqrt(5)) + POW2(y_coor - 0.02*sqrt(5)) - POW2(0.5 + 0.2*sin(5*theta)) > 0)
							{
								p_index += 1;
								if (p_index == N_exact)
								{
									interface_levelset(i, j) = min_val;
									inside = false;
									break;
								}
								continue;
							}
						}*/

						while (p_index < (int)4000)
						{
							T theta = (T)2*PI*p_index/N_exact, theta_p = (T)2*PI*(p_index + 1)/N_exact;
							T co_theta(0);
							T co_ratio = y_coor/x_coor;

							// first quater
							if (x_coor >= (T)0 && y_coor >= (T)0)
							{
								co_theta = atan(co_ratio);
							}
							else if (x_coor <= (T)0 && y_coor >= (T)0)
							{
								co_theta = PI + atan(co_ratio);
							}
							else if (x_coor <= (T)0 && y_coor <= (T)0)
							{
								co_theta = PI + atan(co_ratio);
							}
							else
							{
								co_theta = 2*PI + atan(co_ratio);
							}
							
							if (POW2(x_coor - 0.02*sqrt(5)) + POW2(y_coor - 0.02*sqrt(5)) - POW2(0.5 + 0.2*sin(5*co_theta)) > 0)
							{
							interface_levelset(i, j) = min_val;
							inside = false;
							break;
							}
							else if (POW2(x_coor - 0.02*sqrt(5)) + POW2(y_coor - 0.02*sqrt(5)) - POW2(0.5 + 0.2*sin(5*co_theta)) < 0)
							{
								interface_levelset(i, j) = -min_val;
								inside = true;
								break;
							}
							else
							{
								interface_levelset(i, j) = 0;
							}
														
							p_index += 1;
						}

						//while (p_index < (int)10000)
						//{
						//	T theta = (T)2*PI*p_index/N_exact;
						//	if (POW2(x_coor - 0.02*sqrt(5)) + POW2(y_coor - 0.02*sqrt(5)) - POW2(0.5 + 0.2*sin(5*theta)) > 0)
						//	{
						//		interface_levelset(i, j) = -min_val;
						//		inside = false;
						//		break;
						//	}
						//	if (POW2(x_coor - 0.02*sqrt(5)) + POW2(y_coor - 0.02*sqrt(5)) - POW2(0.5 + 0.2*sin(5*theta)) < 0)
						//	{
						//		p_index += 1;
						//		if (p_index == N_exact/40)
						//		{
						//			interface_levelset(i, j) = -min_val;
						//			inside = true;
						//			break;
						//		}
						//		/*if (p_index == N_exact)
						//		{
						//			interface_levelset(i, j) = -min_val;
						//			inside = true;
						//			break;
						//		}*/
						//		continue;
						//	}
						//}
						/*if (inside == true)
						{
							interface_levelset(i, j) = 0;
						}*/
					}
					interface_levelset.FillGhostCellsFromThreaded(&(interface_levelset.phi), false);
				}
			}
		}
	}

	// For the Fluid Solver
	void AdvanceOneTimeStep(const T& dt_input)
	{
		eulerian_solver.AdvanceOneTimeStep(dt_input);
	}
	
	void AdvanceOneTimeStepThread(const int& thread_id, const T& dt_input)
	{
		eulerian_solver.AdvanceOneTimeStepThread(dt_input, thread_id);
	}

	void AdvanceOneFrameThread(const int& thread_id)
	{
		BEGIN_HEAD_THREAD_WORK
		{
			dt = DetermineTimeStep()*CFL;
			accu_dt += dt;
			cout << "Max x-velocity   : " << eulerian_solver.water_projection->max_velocity_x << endl;
			cout << "Max y-velocity   : " << eulerian_solver.water_projection->max_velocity_y << endl;
			cout << "Time Step        : " << dt << endl;
			cout << "Accumulated Time : " << accu_dt << endl;
			cout << "-------CFL Time Step-------" << endl;
			cout << "c_f              : " << eulerian_solver.c_f << endl;
			cout << "g_f              : " << eulerian_solver.g_f << endl;
			cout << "s_f              : " << eulerian_solver.s_f << endl;
			cout << "v_f              : " << eulerian_solver.v_f << endl;
			
			if (accu_dt > max_dt)
			{
				cout << "Time step reaches at max!:)" << endl;
				exit(0);
			}
		}
		END_HEAD_THREAD_WORK;

		AdvanceOneTimeStepThread(thread_id, dt);
		
		multithreading->Sync(thread_id);
	}

	void AdvanceOneFrame()
	{
		eulerian_solver.water_levelset->ComputeCurvaturesThreaded();
		multithreading->RunThreads(&SIMULATION_WORLD::AdvanceOneFrameThread, this);
		num_current_frame++;
	}

	// For Numerical Integration
	void CalculateIntegralValue()
	{
		numerical_integration.object_levelset->ComputeNormals();
		multithreading->RunThreads(&SIMULATION_WORLD::CalculateIntegralValueThreaded, this);
	    numerical_integration.CalculateRelativeError();
	}

	// For Poisson Equation Test
	void DeterminePoissonSolution()
	{
		multithreading->RunThreads(&SIMULATION_WORLD::DeterminePoissonSolutionThreaded, this);
	}
	
	// For Monge Ampere Equation
	void SolveMongeAmpereEquation()
	{
		multithreading->RunThreads(&SIMULATION_WORLD::SolveMongeAmpereEquationThreaded, this);
	}

	void DeterminePoissonSolutionThreaded(const int& thread_id)
	{
		poisson_equation_test.Solve(thread_id);
	}

    void EstimateSignal()
	{
	    maximum_entropy_analysis.Estimate();    
	}

	void CalculateIntegralValueThreaded(const int& thread_id)
	{
		if (BJ_model_test)
		{
            numerical_integration.SurfaceIntegralThreadedBJ(thread_id);
		}
		else
		{
            numerical_integration.SurfaceIntegralThreaded(thread_id);
		}
	}

	void SolveMongeAmpereEquationThreaded(const int& thread_id)
	{
		monge_ampere_solver.SolveThreaded(thread_id);
	}
	
	T DetermineTimeStep()
	{
		//T dt_for_this = (T)1/frame_rate;
		T dt_for_this(0);

		dt_for_this = eulerian_solver.CFLOneTimeStep();
		
		return dt_for_this;
	}
};
		

