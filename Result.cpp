#pragma once

#include "CSR_MATRIX.h"
#include "CG_METHOD.h"
#include "PROJECTION_2D.h"
#include "SIMULATION_MANAGER.h"
#include "PROJECT_INFO.h"

SIMULATION_MANAGER simulation;

int main(int argc, char** argv)
{
	cout << "Script File: " << argv[1] << endl;
	string script_path = argv[1];

	if (false == simulation.Initialize(script_path))
	{
		return 2;
	}
	
	int i(0), j(0);
	
	if (simulation.simulation->air_water_simulation)
	{
		bool large_bubble(simulation.simulation->eulerian_solver.world_discretization->large_bubble), small_bubble(simulation.simulation->eulerian_solver.world_discretization->small_bubble);
		bool air_bubble_rising(simulation.simulation->eulerian_solver.water_projection->air_bubble_rising), water_drop(simulation.simulation->eulerian_solver.water_projection->water_drop);

		if (large_bubble)
		{
			GRID_ITERATION_2D(simulation.simulation->eulerian_solver.world_discretization->world_grid)
			{
				T x_min = simulation.simulation->eulerian_solver.world_discretization->world_grid.x_min, y_min = simulation.simulation->eulerian_solver.world_discretization->world_grid.y_min;
				T dx = simulation.simulation->eulerian_solver.world_discretization->world_grid.dx, dy = simulation.simulation->eulerian_solver.world_discretization->world_grid.dy;

				simulation.simulation->eulerian_solver.water_levelset->arr(i, j) = sqrt(POW2(x_min + i*dx) + POW2(y_min + j*dy)) - (T)1/(T)3;
			}
			simulation.simulation->eulerian_solver.water_levelset->FillGhostCellsContinuousDerivativesFromPointer(&(simulation.simulation->eulerian_solver.water_levelset->phi), true);
		}
		if (small_bubble)
		{
			GRID_ITERATION_2D(simulation.simulation->eulerian_solver.world_discretization->world_grid)
			{
				T x_min = simulation.simulation->eulerian_solver.world_discretization->world_grid.x_min, y_min = simulation.simulation->eulerian_solver.world_discretization->world_grid.y_min;
				T dx = simulation.simulation->eulerian_solver.world_discretization->world_grid.dx, dy = simulation.simulation->eulerian_solver.world_discretization->world_grid.dy;

				simulation.simulation->eulerian_solver.water_levelset->arr(i, j) = sqrt(POW2(x_min + i*dx) + POW2(y_min + j*dy)) - (T)1/(T)300;
				/*if (y_min + j*dy >= (T)1/(T)200)
				{
				simulation.simulation->eulerian_solver.water_levelset->arr(i, j) = abs(y_min + j*dy - (T)1/(T)200);
				}
				else if ((y_min + j*dy >= 0) && (y_min + j*dy < (T)1/(T)200))
				{
				simulation.simulation->eulerian_solver.water_levelset->arr(i, j) = y_min + j*dy - (T)1/(T)200;
				}
				else if ((y_min + j*dy < 0) && (y_min + j*dy > -(T)1/(T)200))
				{
				simulation.simulation->eulerian_solver.water_levelset->arr(i, j) = -(y_min + j*dy + (T)1/(T)200);
				}
				else if (y_min + j*dy <= -(T)1/(T)200)
				{
				simulation.simulation->eulerian_solver.water_levelset->arr(i, j) = abs(y_min + j*dy + (T)1/(T)200);
				}*/
			}
			simulation.simulation->eulerian_solver.water_levelset->FillGhostCellsContinuousDerivativesFromPointer(&(simulation.simulation->eulerian_solver.water_levelset->phi), true);
		}	

		int i_res = simulation.simulation->eulerian_solver.base_grid.i_res - 1, j_res = simulation.simulation->eulerian_solver.base_grid.j_res - 1;

		ofstream fout;

		//simulation.simulation->eulerian_solver.water_levelset->ComputeCurvatures();

		if (large_bubble)
		{
			if (air_bubble_rising)
			{
				char filename[50];
				sprintf(filename, "Test/Large_Air_Bubble_%dX%d/%dX%d_initial", i_res, j_res, i_res, j_res);
				fout.open(filename);
				for (int j = simulation.simulation->eulerian_solver.world_discretization->world_grid.j_start; j <= simulation.simulation->eulerian_solver.world_discretization->world_grid.j_end; j++)
				{
					for (int i = simulation.simulation->eulerian_solver.world_discretization->world_grid.i_start; i <= simulation.simulation->eulerian_solver.world_discretization->world_grid.i_end; i++)
						fout << simulation.simulation->eulerian_solver.water_levelset->arr(i, j) << " ";

					fout << "\n";
				}
				fout.close();
			}
			if (water_drop)
			{
				char filename[50];
				sprintf(filename, "Test/Large_Water_Drop_%dX%d/%dX%d_initial", i_res, j_res, i_res, j_res);
				fout.open(filename);
				for (int j = simulation.simulation->eulerian_solver.world_discretization->world_grid.j_start; j <= simulation.simulation->eulerian_solver.world_discretization->world_grid.j_end; j++)
				{
					for (int i = simulation.simulation->eulerian_solver.world_discretization->world_grid.i_start; i <= simulation.simulation->eulerian_solver.world_discretization->world_grid.i_end; i++)
						fout << simulation.simulation->eulerian_solver.water_levelset->arr(i, j) << " ";

					fout << "\n";
				}
				fout.close();
			}
		}
		if (small_bubble)
		{
			if (air_bubble_rising)
			{
				char filename[50];
				sprintf(filename, "Test/Small_Air_Bubble_%dX%d/%dX%d_initial", i_res, j_res, i_res, j_res);
				fout.open(filename);
				for (int j = simulation.simulation->eulerian_solver.world_discretization->world_grid.j_start; j <= simulation.simulation->eulerian_solver.world_discretization->world_grid.j_end; j++)
				{
					for (int i = simulation.simulation->eulerian_solver.world_discretization->world_grid.i_start; i <= simulation.simulation->eulerian_solver.world_discretization->world_grid.i_end; i++)
						fout << simulation.simulation->eulerian_solver.water_levelset->arr(i, j) << " ";

					fout << "\n";
				}
				fout.close();

				sprintf(filename, "Test/Small_Air_Bubble_%dX%d/Velocity/x_%dX%d_initial", i_res, j_res, i_res, j_res);
				fout.open(filename);
				for (int j = simulation.simulation->eulerian_solver.water_velocity_field_mac_x->j_start; j <= simulation.simulation->eulerian_solver.water_velocity_field_mac_x->j_end; j++)
				{
					for (int i = simulation.simulation->eulerian_solver.water_velocity_field_mac_x->i_start; i <= simulation.simulation->eulerian_solver.water_velocity_field_mac_x->i_end; i++)
						fout << simulation.simulation->eulerian_solver.water_velocity_field_mac_x->array_for_this(i, j) << " ";

					fout << "\n";
				}
				fout.close();

				sprintf(filename, "Test/Small_Air_Bubble_%dX%d/Velocity/y_%dX%d_initial", i_res, j_res, i_res, j_res);
				fout.open(filename);
				for (int j = simulation.simulation->eulerian_solver.water_velocity_field_mac_y->j_start; j <= simulation.simulation->eulerian_solver.water_velocity_field_mac_y->j_end; j++)
				{
					for (int i = simulation.simulation->eulerian_solver.water_velocity_field_mac_y->i_start; i <= simulation.simulation->eulerian_solver.water_velocity_field_mac_y->i_end; i++)
						fout << simulation.simulation->eulerian_solver.water_velocity_field_mac_y->array_for_this(i, j) << " ";

					fout << "\n";
				}
				fout.close();
			}
			if (water_drop)
			{
				char filename[50];
				sprintf(filename, "Test/Small_Water_Drop_%dX%d/%dX%d_initial", i_res, j_res, i_res, j_res);
				fout.open(filename);
				for (int j = simulation.simulation->eulerian_solver.world_discretization->world_grid.j_start; j <= simulation.simulation->eulerian_solver.world_discretization->world_grid.j_end; j++)
				{
					for (int i = simulation.simulation->eulerian_solver.world_discretization->world_grid.i_start; i <= simulation.simulation->eulerian_solver.world_discretization->world_grid.i_end; i++)
						fout << simulation.simulation->eulerian_solver.water_levelset->arr(i, j) << " ";

					fout << "\n";
				}
				fout.close();
			}
		}

		int iter(0);

		// Calculate the time
		clock_t start_time;
		clock_t end_time;
		clock_t middle_time;

		start_time = clock();

		T dt = simulation.simulation->dt;
		//T accu_dt(dt);

		int index(0);

		if (large_bubble)
		{
			if (air_bubble_rising)
			{
				while(simulation.simulation->accu_dt < 2.0)
				{
					cout << "---------------------" << iter << "---------------------" << endl;

					// Time Calculation
					middle_time = clock();

					// Simulating one time step
					simulation.simulation->AdvanceOneFrame();
					//simulation.simulation->eulerian_solver.AdvanceOneTimeStep(dt);
					//dt = simulation.simulation->eulerian_solver.CFLOneTimeStep();

					//simulation.simulation->accu_dt+= dt;
					
					iter = iter + 1;

					end_time = clock();

					double elapsed_sec = (double)(end_time - middle_time) / CLOCKS_PER_SEC;
					double total_sec = (double)(end_time - start_time) / CLOCKS_PER_SEC;

					cout << "Time in terms of given simulation: " << simulation.simulation->accu_dt<< endl;
					cout << "////////////////// Simulating Time //////////////////" << endl;
					cout << "Simulation time per each iteration: " << elapsed_sec << endl;
					cout << "Accumulated time: " << total_sec << endl;
					cout << "/////////////////////////////////////////////////////" << endl;
					cout << "Max x-Velocity: " << simulation.simulation->eulerian_solver.water_projection->max_velocity_x << endl;								
					cout << "Max y-Velocity: " << simulation.simulation->eulerian_solver.water_projection->max_velocity_y << endl;								

					if (iter % (i_res/20) == 0)
					{
						index += 1;
						char filename[50];
						sprintf(filename, "Test/Large_Air_Bubble_%dX%d/%dX%d_%d", i_res, j_res, i_res, j_res, index);
						fout.open(filename);
						for (int j = simulation.simulation->eulerian_solver.world_discretization->world_grid.j_start; j <= simulation.simulation->eulerian_solver.world_discretization->world_grid.j_end; j++)
						{
							for (int i = simulation.simulation->eulerian_solver.world_discretization->world_grid.i_start; i <= simulation.simulation->eulerian_solver.world_discretization->world_grid.i_end; i++)
								fout << simulation.simulation->eulerian_solver.water_levelset->arr(i, j) << " ";

							fout << "\n";
						}
						fout.close();
					}
				}
			}
			if (water_drop)
			{
				while(simulation.simulation->accu_dt< 0.5)
				{
					cout << "---------------------" << iter << "---------------------" << endl;

					// Time Calculation
					middle_time = clock();

					// Simulating one time step
					simulation.simulation->eulerian_solver.AdvanceOneTimeStep(dt);
					dt = simulation.simulation->eulerian_solver.CFLOneTimeStep();

					simulation.simulation->accu_dt+= dt;
					iter = iter + 1;

					end_time = clock();

					double elapsed_sec = (double)(end_time - middle_time) / CLOCKS_PER_SEC;
					double total_sec = (double)(end_time - start_time) / CLOCKS_PER_SEC;

					cout << "Time in terms of given simulation: " << simulation.simulation->accu_dt<< endl;
					cout << "////////////////// Simulating Time //////////////////" << endl;
					cout << "Simulation time per each iteration: " << elapsed_sec << endl;
					cout << "Accumulated time: " << total_sec << endl;
					cout << "/////////////////////////////////////////////////////" << endl;

					if (iter % (i_res/20) == 0)
					{
						index += 1;
						char filename[50];
						sprintf(filename, "Test/Large_Water_Drop_%dX%d/%dX%d_%d", i_res, j_res, i_res, j_res, index);
						fout.open(filename);
						for (int j = simulation.simulation->eulerian_solver.world_discretization->world_grid.j_start; j <= simulation.simulation->eulerian_solver.world_discretization->world_grid.j_end; j++)
						{
							for (int i = simulation.simulation->eulerian_solver.world_discretization->world_grid.i_start; i <= simulation.simulation->eulerian_solver.world_discretization->world_grid.i_end; i++)
								fout << simulation.simulation->eulerian_solver.water_levelset->arr(i, j) << " ";

							fout << "\n";
						}
						fout.close();
					}

				}
			}
		}
		if (small_bubble)
		{
			if (air_bubble_rising)
			{
				while(simulation.simulation->accu_dt < 2.0)
				{
					cout << "---------------------" << iter << "---------------------" << endl;

					// Time Calculation
					middle_time = clock();

					// Simulating one time step
					simulation.simulation->AdvanceOneFrame();
					//simulation.simulation->eulerian_solver.AdvanceOneTimeStep(dt);
					//dt = simulation.simulation->eulerian_solver.CFLOneTimeStep();

					//simulation.simulation->accu_dt+= dt;
					iter = iter + 1;

					end_time = clock();

					double elapsed_sec = (double)(end_time - middle_time) / CLOCKS_PER_SEC;
					double total_sec = (double)(end_time - start_time) / CLOCKS_PER_SEC;

					cout << "Time in terms of given simulation: " << simulation.simulation->accu_dt<< endl;
					cout << "////////////////// Simulating Time //////////////////" << endl;
					cout << "Simulation time per each iteration: " << elapsed_sec << endl;
					cout << "Accumulated time: " << total_sec << endl;
					cout << "/////////////////////////////////////////////////////" << endl;

					if (iter % i_res == 0)
					{
						index += 1;
						char filename[50];
						sprintf(filename, "Test/Small_Air_Bubble_%dX%d/%dX%d_%d", i_res, j_res, i_res, j_res, index);
						fout.open(filename);
						for (int j = simulation.simulation->eulerian_solver.world_discretization->world_grid.j_start; j <= simulation.simulation->eulerian_solver.world_discretization->world_grid.j_end; j++)
						{
							for (int i = simulation.simulation->eulerian_solver.world_discretization->world_grid.i_start; i <= simulation.simulation->eulerian_solver.world_discretization->world_grid.i_end; i++)
								fout << simulation.simulation->eulerian_solver.water_levelset->arr(i, j) << " ";

							fout << "\n";
						}
						fout.close();

						sprintf(filename, "Test/Small_Air_Bubble_%dX%d/Velocity/x_%dX%d_%d", i_res, j_res, i_res, j_res, index);
						fout.open(filename);
						for (int j = simulation.simulation->eulerian_solver.water_velocity_field_mac_x->j_start; j <= simulation.simulation->eulerian_solver.water_velocity_field_mac_x->j_end; j++)
						{
							for (int i = simulation.simulation->eulerian_solver.water_velocity_field_mac_x->i_start; i <= simulation.simulation->eulerian_solver.water_velocity_field_mac_x->i_end; i++)
								fout << simulation.simulation->eulerian_solver.water_velocity_field_mac_x->array_for_this(i, j) << " ";

							fout << "\n";
						}
						fout.close();

						sprintf(filename, "Test/Small_Air_Bubble_%dX%d/Velocity/y_%dX%d_%d", i_res, j_res, i_res, j_res, index);
						fout.open(filename);
						for (int j = simulation.simulation->eulerian_solver.water_velocity_field_mac_y->j_start; j <= simulation.simulation->eulerian_solver.water_velocity_field_mac_y->j_end; j++)
						{
							for (int i = simulation.simulation->eulerian_solver.water_velocity_field_mac_y->i_start; i <= simulation.simulation->eulerian_solver.water_velocity_field_mac_y->i_end; i++)
								fout << simulation.simulation->eulerian_solver.water_velocity_field_mac_y->array_for_this(i, j) << " ";

							fout << "\n";
						}
						fout.close();

						sprintf(filename, "Test/Small_Air_Bubble_%dX%d/Pressure/%d", i_res, j_res, index);
						fout.open(filename);
						for (int j = simulation.simulation->eulerian_solver.water_projection->pressure_field.j_start; j <= simulation.simulation->eulerian_solver.water_projection->pressure_field.j_end; j++)
						{
							for (int i = simulation.simulation->eulerian_solver.water_projection->pressure_field.i_start; i <= simulation.simulation->eulerian_solver.water_projection->pressure_field.i_end; i++)
								fout << simulation.simulation->eulerian_solver.water_projection->pressure_field(i, j) << " ";

							fout << "\n";
						}
						fout.close();
					}
				}
			}
			if (water_drop)
			{
				while(simulation.simulation->accu_dt< 0.05)
				{
					cout << "---------------------" << iter << "---------------------" << endl;

					// Time Calculation
					middle_time = clock();

					// Simulating one time step
					simulation.simulation->eulerian_solver.AdvanceOneTimeStep(dt);
					dt = simulation.simulation->eulerian_solver.CFLOneTimeStep();

					simulation.simulation->accu_dt+= dt;
					iter = iter + 1;

					end_time = clock();

					double elapsed_sec = (double)(end_time - middle_time) / CLOCKS_PER_SEC;
					double total_sec = (double)(end_time - start_time) / CLOCKS_PER_SEC;

					cout << "Time in terms of given simulation: " << simulation.simulation->accu_dt<< endl;
					cout << "////////////////// Simulating Time //////////////////" << endl;
					cout << "Simulation time per each iteration: " << elapsed_sec << endl;
					cout << "Accumulated time: " << total_sec << endl;
					cout << "/////////////////////////////////////////////////////" << endl;

					if (iter % i_res == 0)
					{
						index += 1;
						char filename[50];
						sprintf(filename, "Test/Small_Water_Drop_%dX%d/%dX%d_%d", i_res, j_res, i_res, j_res, index);
						fout.open(filename);
						for (int j = simulation.simulation->eulerian_solver.world_discretization->world_grid.j_start; j <= simulation.simulation->eulerian_solver.world_discretization->world_grid.j_end; j++)
						{
							for (int i = simulation.simulation->eulerian_solver.world_discretization->world_grid.i_start; i <= simulation.simulation->eulerian_solver.world_discretization->world_grid.i_end; i++)
								fout << simulation.simulation->eulerian_solver.water_levelset->arr(i, j) << " ";

							fout << "\n";
						}
						fout.close();
					}
				}
			}
		} 
	}
	
	if (simulation.simulation->oil_water_simulation)
	{
		// Interface
		GRID_ITERATION_2D(simulation.simulation->eulerian_solver.world_discretization->world_grid)
		{
			T x_min = simulation.simulation->eulerian_solver.world_discretization->world_grid.x_min, y_min = simulation.simulation->eulerian_solver.world_discretization->world_grid.y_min;
			T dx = simulation.simulation->eulerian_solver.world_discretization->world_grid.dx, dy = simulation.simulation->eulerian_solver.world_discretization->world_grid.dy;
		
			if (x_min + i*dx < simulation.simulation->A_0*cos(simulation.simulation->alpha*(y_min + j*dy)) + simulation.simulation->eulerian_solver.R1)
			{
				simulation.simulation->eulerian_solver.water_levelset->arr(i,j) = -100;
			}
			else if (x_min + i*dx > simulation.simulation->A_0*cos(simulation.simulation->alpha*(y_min + j*dy)) + simulation.simulation->eulerian_solver.R1)
			{
				simulation.simulation->eulerian_solver.water_levelset->arr(i,j) = 100;
			}
			else
			{
				simulation.simulation->eulerian_solver.water_levelset->arr(i,j) = 0;
			}
		}

		int number_for_interface(3001);

		ARRAY<VT> point_for_interface;
		point_for_interface.Initialize(number_for_interface);

		int ind(0);
		
		T step = 3.14/(number_for_interface - 1);

		for (ind = 0; ind < number_for_interface; ind++)
		{
			point_for_interface[ind] = VT(simulation.simulation->A_0*cos(simulation.simulation->alpha*(step*ind)) + simulation.simulation->eulerian_solver.R1, step*ind);
		}

		GRID_ITERATION_2D(simulation.simulation->eulerian_solver.water_levelset->exact_solution.grid)
		{
			T x_min = simulation.simulation->eulerian_solver.world_discretization->world_grid.x_min, y_min = simulation.simulation->eulerian_solver.world_discretization->world_grid.y_min;
			T dx = simulation.simulation->eulerian_solver.world_discretization->world_grid.dx, dy = simulation.simulation->eulerian_solver.world_discretization->world_grid.dy;
		
			T min_value = sqrt(POW2(x_min + i*dx - point_for_interface[0].x) + POW2(y_min + j*dy - point_for_interface[0].y));
			for (int k = 0; k < point_for_interface.length; k++)
			{
				min_value = min(min_value, sqrt(POW2(x_min + i*dx - point_for_interface[k].x) + POW2(y_min + j*dy - point_for_interface[k].y)));
			}
		
			if (simulation.simulation->eulerian_solver.water_levelset->arr(i, j) > (T)0)
			{
				simulation.simulation->eulerian_solver.water_levelset->arr(i, j) = min_value;
			}
			else if (simulation.simulation->eulerian_solver.water_levelset->arr(i, j) < (T)0)
			{
				simulation.simulation->eulerian_solver.water_levelset->arr(i, j) = -min_value;
			}
			else
			{
				simulation.simulation->eulerian_solver.water_levelset->arr(i, j) = (T)0;
			}
		}

		/*GRID_ITERATION_2D(simulation.simulation->eulerian_solver.world_discretization->world_grid)
		{
			T x_min = simulation.simulation->eulerian_solver.world_discretization->world_grid.x_min, y_min = simulation.simulation->eulerian_solver.world_discretization->world_grid.y_min;
			T dx = simulation.simulation->eulerian_solver.world_discretization->world_grid.dx, dy = simulation.simulation->eulerian_solver.world_discretization->world_grid.dy;

			simulation.simulation->eulerian_solver.water_levelset->arr(i, j) = (x_min + i*dx) - (simulation.simulation->A_0*cos(simulation.simulation->alpha*(y_min + j*dy)) + simulation.simulation->eulerian_solver.R1);
		}*/

		simulation.simulation->eulerian_solver.water_levelset->FillGhostCellsContinuousDerivativesFromPointer(&(simulation.simulation->eulerian_solver.water_levelset->phi), true);
		
		int i_res = simulation.simulation->eulerian_solver.base_grid.i_res - 1, j_res = simulation.simulation->eulerian_solver.base_grid.j_res - 1;

		ofstream fout;

		fout.open("Test/Oil_Water/Example_1/initial");
		for (int j = simulation.simulation->eulerian_solver.world_discretization->world_grid.j_start; j <= simulation.simulation->eulerian_solver.world_discretization->world_grid.j_end; j++)
		{
			for (int i = simulation.simulation->eulerian_solver.world_discretization->world_grid.i_start; i <= simulation.simulation->eulerian_solver.world_discretization->world_grid.i_end; i++)
				fout << simulation.simulation->eulerian_solver.water_levelset->arr(i, j) << " ";
			
				fout << "\n";
		}
		fout.close();
		
		fout.open("Test/Oil_Water/Example_1/Velocity/velocity_x_initial");
		for (int j = simulation.simulation->eulerian_solver.water_velocity_field_mac_x->j_start; j <= simulation.simulation->eulerian_solver.water_velocity_field_mac_x->j_end; j++)
		{
			for (int i = simulation.simulation->eulerian_solver.water_velocity_field_mac_x->i_start; i <= simulation.simulation->eulerian_solver.water_velocity_field_mac_x->i_end; i++)
				fout << simulation.simulation->eulerian_solver.water_velocity_field_mac_x->array_for_this(i, j) << " ";
			
				fout << "\n";
		}
		fout.close();

		fout.open("Test/Oil_Water/Example_1/Velocity/velocity_y_initial");
		for (int j = simulation.simulation->eulerian_solver.water_velocity_field_mac_y->j_start; j <= simulation.simulation->eulerian_solver.water_velocity_field_mac_y->j_end; j++)
		{
			for (int i = simulation.simulation->eulerian_solver.water_velocity_field_mac_y->i_start; i <= simulation.simulation->eulerian_solver.water_velocity_field_mac_y->i_end; i++)
				fout << simulation.simulation->eulerian_solver.water_velocity_field_mac_y->array_for_this(i, j) << " ";
			
				fout << "\n";
		}
		fout.close();

		int iter(0);

		// Calculate the time
		clock_t start_time;
		clock_t end_time;
		clock_t middle_time;

		start_time = clock();

		T dt = simulation.simulation->dt;
		T accu_dt(dt);

		int index(0);
		
		while(simulation.simulation->accu_dt< 0.5)
		{
			cout << "---------------------" << iter << "---------------------" << endl;

			// Time Calculation
			middle_time = clock();

			// Simulating one time step
			simulation.simulation->eulerian_solver.AdvanceOneTimeStep(dt);
			dt = simulation.simulation->eulerian_solver.CFLOneTimeStep();

			simulation.simulation->accu_dt+= dt;
			iter = iter + 1;

			end_time = clock();

			double elapsed_sec = (double)(end_time - middle_time) / CLOCKS_PER_SEC;
			double total_sec = (double)(end_time - start_time) / CLOCKS_PER_SEC;

			cout << "Time in terms of given simulation: " << simulation.simulation->accu_dt<< endl;
			cout << "////////////////// Simulating Time //////////////////" << endl;
			cout << "Simulation time per each iteration: " << elapsed_sec << endl;
			cout << "Accumulated time: " << total_sec << endl;
			cout << "/////////////////////////////////////////////////////" << endl;
			cout << "Max x-Velocity: " << simulation.simulation->eulerian_solver.water_projection->max_velocity_x << endl;								
			cout << "Max y-Velocity: " << simulation.simulation->eulerian_solver.water_projection->max_velocity_y << endl;								

			if (iter % 1 == 0)
			{
				index += 1;
				char filename[50];
				sprintf(filename, "Test/Oil_Water/Example_1/%d", index);
				fout.open(filename);
				for (int j = simulation.simulation->eulerian_solver.world_discretization->world_grid.j_start; j <= simulation.simulation->eulerian_solver.world_discretization->world_grid.j_end; j++)
				{
					for (int i = simulation.simulation->eulerian_solver.world_discretization->world_grid.i_start; i <= simulation.simulation->eulerian_solver.world_discretization->world_grid.i_end; i++)
						fout << simulation.simulation->eulerian_solver.water_levelset->arr(i, j) << " ";

						fout << "\n";
				}
				fout.close();

				char filename1[50];
				sprintf(filename1, "Test/Oil_Water/Example_1/Velocity/x_%d", index);
				fout.open(filename1);
				for (int j = simulation.simulation->eulerian_solver.water_velocity_field_mac_x->j_start; j <= simulation.simulation->eulerian_solver.water_velocity_field_mac_x->j_end; j++)
				{
					for (int i = simulation.simulation->eulerian_solver.water_velocity_field_mac_x->i_start; i <= simulation.simulation->eulerian_solver.water_velocity_field_mac_x->i_end; i++)
						fout << simulation.simulation->eulerian_solver.water_velocity_field_mac_x->array_for_this(i, j) << " ";
			
						fout << "\n";
				}
				fout.close();

				char filename2[50];
				sprintf(filename2, "Test/Oil_Water/Example_1/Velocity/y_%d", index);
				fout.open(filename2);
				for (int j = simulation.simulation->eulerian_solver.water_velocity_field_mac_y->j_start; j <= simulation.simulation->eulerian_solver.water_velocity_field_mac_y->j_end; j++)
				{
					for (int i = simulation.simulation->eulerian_solver.water_velocity_field_mac_y->i_start; i <= simulation.simulation->eulerian_solver.water_velocity_field_mac_y->i_end; i++)
						fout << simulation.simulation->eulerian_solver.water_velocity_field_mac_y->array_for_this(i, j) << " ";
			
						fout << "\n";
				}
				fout.close();
			}
		}
	}

	return 0;
}
