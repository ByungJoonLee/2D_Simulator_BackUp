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
	
	T theta = asin((T)2.5/15);
		
	GRID_ITERATION_2D(simulation.simulation->eulerian_solver.world_discretization->world_grid)
	{
		T x_min = simulation.simulation->eulerian_solver.world_discretization->world_grid.x_min, y_min = simulation.simulation->eulerian_solver.world_discretization->world_grid.y_min;
		T dx = simulation.simulation->eulerian_solver.world_discretization->world_grid.dx, dy = simulation.simulation->eulerian_solver.world_discretization->world_grid.dy;
		
		if ((POW2(x_min + i*dx - 50.0) + POW2(y_min + j*dy - 75.0) < 15*15) && ((x_min + i*dx < 47.5) || (x_min + i*dx > 52.5) || (y_min + j*dy > 75 + 15*sin(theta - PI/2) + 25)))
		{
			simulation.simulation->eulerian_solver.water_levelset->arr(i,j) = -100;
		}
		else if ((POW2(x_min + i*dx - 50.0) + POW2(y_min + j*dy - 75.0) == 15*15) && ((x_min + i*dx == 47.5) || (x_min + i*dx == 52.5)) && (y_min + j*dy == 75 + 15*sin(theta - PI/2) + 25))
		{
			simulation.simulation->eulerian_solver.water_levelset->arr(i,j) = 0;
		}
		else
		{
			simulation.simulation->eulerian_solver.water_levelset->arr(i,j) = 100;
		}
	}
	
	// Numerically Exact Solution
	int number_for_exact_point(3000);

	ARRAY<VT> exact_point;
	exact_point.Initialize(number_for_exact_point);

	int index(0);
	
	T step = ((T)2*PI - (T)2*theta)/2000;
	
	for (index = 0; index <= 2000; index++)
	{
		exact_point[index] = VT(50 + 15*cos(theta - PI/2 + index*step), 75 + 15*sin(theta - PI/2 + index*step));
	}
	for (index = 2001; index <= 2400; index++)
	{
		exact_point[index] = VT(47.5, 75 + 15*sin(theta - PI/2 + 2000*step) + (index - 2000)*(T)25/400);
	}
	for (index = 2401; index <= 2600; index++)
	{
		exact_point[index] = VT(47.5 + (index - 2400)*(T)5/200, 75 + 15*sin(theta - PI/2 + 2000*step) + 25);
	}
	for (index = 2601; index < 3000; index++)
	{
		exact_point[index] = VT(52.5, 75 + 15*sin(theta - PI/2 + 2000*step) + 25 - (index - 2600)*(T)25/400);
	}

	GRID_ITERATION_2D(simulation.simulation->eulerian_solver.water_levelset->exact_solution.grid)
	{
		T x_min = simulation.simulation->eulerian_solver.world_discretization->world_grid.x_min, y_min = simulation.simulation->eulerian_solver.world_discretization->world_grid.y_min;
		T dx = simulation.simulation->eulerian_solver.world_discretization->world_grid.dx, dy = simulation.simulation->eulerian_solver.world_discretization->world_grid.dy;
		
		T min_value = sqrt(POW2(x_min + i*dx - exact_point[0].x) + POW2(y_min + j*dy - exact_point[0].y));
		for (int k = 0; k < exact_point.length; k++)
		{
			min_value = min(min_value, sqrt(POW2(x_min + i*dx - exact_point[k].x) + POW2(y_min + j*dy - exact_point[k].y)));
		}
		
		if (simulation.simulation->eulerian_solver.water_levelset->arr(i, j) > (T)0)
		{
			simulation.simulation->eulerian_solver.water_levelset->exact_solution(i, j) = min_value;
		}
		else if (simulation.simulation->eulerian_solver.water_levelset->arr(i, j) < (T)0)
		{
			simulation.simulation->eulerian_solver.water_levelset->exact_solution(i, j) = -min_value;
		}
		else
		{
			simulation.simulation->eulerian_solver.water_levelset->exact_solution(i, j) = (T)0;
		}
	}

	GRID_ITERATION_2D(simulation.simulation->eulerian_solver.water_levelset->grid)
	{
		simulation.simulation->eulerian_solver.water_levelset->arr(i, j) = simulation.simulation->eulerian_solver.water_levelset->exact_solution(i, j);
	}

	simulation.simulation->eulerian_solver.water_levelset->FillGhostCellsContinuousDerivativesFromPointer(&(simulation.simulation->eulerian_solver.water_levelset->phi), true);
	
	simulation.simulation->eulerian_solver.water_levelset->ComputeGradient();

	T dt = (T)1;

	//simulation.simulation->eulerian_solver.Reinitialization(dt);

	ofstream fout;
	fout.open("A Gradient Augmented Levelset/Example 1/magnitude_of_gradient"); 
	for (int j = simulation.simulation->eulerian_solver.world_discretization->world_grid.j_start; j <= simulation.simulation->eulerian_solver.world_discretization->world_grid.j_end; j++)
	{
		for (int i = simulation.simulation->eulerian_solver.world_discretization->world_grid.i_start; i <= simulation.simulation->eulerian_solver.world_discretization->world_grid.i_end; i++)
			fout << sqrt(POW2(simulation.simulation->eulerian_solver.water_levelset->gradient(i, j).x) + POW2(simulation.simulation->eulerian_solver.water_levelset->gradient(i, j).y)) << " ";
	
		fout << "\n";
	}
	fout.close();

	fout.open("A Gradient Augmented Levelset/Example 1/gradient_x"); 
	for (int j = simulation.simulation->eulerian_solver.world_discretization->world_grid.j_start; j <= simulation.simulation->eulerian_solver.world_discretization->world_grid.j_end; j++)
	{
		for (int i = simulation.simulation->eulerian_solver.world_discretization->world_grid.i_start; i <= simulation.simulation->eulerian_solver.world_discretization->world_grid.i_end; i++)
			fout << simulation.simulation->eulerian_solver.water_levelset->gradient(i, j).x << " ";
	
		fout << "\n";
	}
	fout.close();

	fout.open("A Gradient Augmented Levelset/Example 1/gradient_y"); 
	for (int j = simulation.simulation->eulerian_solver.world_discretization->world_grid.j_start; j <= simulation.simulation->eulerian_solver.world_discretization->world_grid.j_end; j++)
	{
		for (int i = simulation.simulation->eulerian_solver.world_discretization->world_grid.i_start; i <= simulation.simulation->eulerian_solver.world_discretization->world_grid.i_end; i++)
			fout << simulation.simulation->eulerian_solver.water_levelset->gradient(i, j).y << " ";
	
		fout << "\n";
	}
	fout.close();

	fout.open("A Gradient Augmented Levelset/Example 1/Example_1_0"); 
	for (int j = simulation.simulation->eulerian_solver.world_discretization->world_grid.j_start; j <= simulation.simulation->eulerian_solver.world_discretization->world_grid.j_end; j++)
	{
		for (int i = simulation.simulation->eulerian_solver.world_discretization->world_grid.i_start; i <= simulation.simulation->eulerian_solver.world_discretization->world_grid.i_end; i++)
			fout << simulation.simulation->eulerian_solver.water_levelset->arr(i, j) << " ";
	
		fout << "\n";
	}
	fout.close();
	
	fout.open("A Gradient Augmented Levelset/Example 1/Exact"); 
	for (int j = simulation.simulation->eulerian_solver.world_discretization->world_grid.j_start; j <= simulation.simulation->eulerian_solver.world_discretization->world_grid.j_end; j++)
	{
		for (int i = simulation.simulation->eulerian_solver.world_discretization->world_grid.i_start; i <= simulation.simulation->eulerian_solver.world_discretization->world_grid.i_end; i++)
			fout << simulation.simulation->eulerian_solver.water_levelset->exact_solution(i, j) << " ";
	
		fout << "\n";
	}
	fout.close();
	
	GRID_ITERATION_2D(simulation.simulation->eulerian_solver.water_velocity_field->grid)
	{
		T x_min = simulation.simulation->eulerian_solver.water_velocity_field->x_min, y_min = simulation.simulation->eulerian_solver.water_velocity_field->y_min;
		T dx = simulation.simulation->eulerian_solver.water_velocity_field->dx, dy = simulation.simulation->eulerian_solver.water_velocity_field->dy;
		
		simulation.simulation->eulerian_solver.water_velocity_field->array_for_this(i, j).x = PI/314*(50 - (y_min + j*dy));
		simulation.simulation->eulerian_solver.water_velocity_field->array_for_this(i, j).y = PI/314*((x_min + i*dx) - 50);
	}

	simulation.simulation->eulerian_solver.vector_field_ghost.FillGhostCellsContinuousDerivativesFrom(simulation.simulation->eulerian_solver.water_velocity_field->array_for_this, true);
	
	GRID_ITERATION_2D(simulation.simulation->eulerian_solver.water_velocity_field_mac_x->grid)
	{
		T x_min = simulation.simulation->eulerian_solver.water_velocity_field_mac_x->x_min, y_min = simulation.simulation->eulerian_solver.water_velocity_field_mac_x->y_min;
		T dx = simulation.simulation->eulerian_solver.water_velocity_field_mac_x->dx, dy = simulation.simulation->eulerian_solver.water_velocity_field_mac_x->dy;
		
		simulation.simulation->eulerian_solver.water_velocity_field_mac_x->array_for_this(i, j) = PI/314*(50 - (y_min + j*dy));
	}

	simulation.simulation->eulerian_solver.vector_field_mac_ghost_x.FillGhostCellsFrom(simulation.simulation->eulerian_solver.water_velocity_field_mac_x->array_for_this, true);

	GRID_ITERATION_2D(simulation.simulation->eulerian_solver.water_velocity_field_mac_y->grid)
	{
		T x_min = simulation.simulation->eulerian_solver.water_velocity_field_mac_y->x_min, y_min = simulation.simulation->eulerian_solver.water_velocity_field_mac_y->y_min;
		T dx = simulation.simulation->eulerian_solver.water_velocity_field_mac_y->dx, dy = simulation.simulation->eulerian_solver.water_velocity_field_mac_y->dy;
		
		simulation.simulation->eulerian_solver.water_velocity_field_mac_y->array_for_this(i, j) = PI/314*((x_min + i*dx) - 50);
	}

	simulation.simulation->eulerian_solver.vector_field_mac_ghost_y.FillGhostCellsFrom(simulation.simulation->eulerian_solver.water_velocity_field_mac_y->array_for_this, true);

	int iter(0);
		
	while(iter < 628)
	{
		cout << "---------------------" << iter << "---------------------" << endl;
		
		// Simulating one time step
		simulation.simulation->eulerian_solver.AdvanceOneTimeStep(dt);
		
		iter = iter + 1;
		
		char filename[50];
		sprintf(filename, "A Gradient Augmented Levelset/Example 1/Example_1_%d", iter);
		fout.open(filename);
		for (int j = simulation.simulation->eulerian_solver.world_discretization->world_grid.j_start; j <= simulation.simulation->eulerian_solver.world_discretization->world_grid.j_end; j++)
		{
			for (int i = simulation.simulation->eulerian_solver.world_discretization->world_grid.i_start; i <= simulation.simulation->eulerian_solver.world_discretization->world_grid.i_end; i++)
				fout << simulation.simulation->eulerian_solver.water_levelset->arr(i, j) << " ";
			
			fout << "\n";
		}
		fout.close();

	}
		
	return 0;
}
