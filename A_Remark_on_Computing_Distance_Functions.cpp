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
	
	// Reinitialization Test
	GRID_ITERATION_2D(simulation.simulation->eulerian_solver.world_discretization->world_grid)
	{
		T x_min = simulation.simulation->eulerian_solver.world_discretization->world_grid.x_min, y_min = simulation.simulation->eulerian_solver.world_discretization->world_grid.y_min;
		T dx = simulation.simulation->eulerian_solver.world_discretization->world_grid.dx, dy = simulation.simulation->eulerian_solver.world_discretization->world_grid.dy;
		
		/*if (simulation.simulation->eulerian_solver.water_levelset->arr(i, j) = sqrt(POW2(x_min + i*dx) + POW2(y_min + j*dy)) < (T)2)
		{
			simulation.simulation->eulerian_solver.water_levelset->arr(i, j) = -100;
		}
		else if (simulation.simulation->eulerian_solver.water_levelset->arr(i, j) = sqrt(POW2(x_min + i*dx) + POW2(y_min + j*dy)) > (T)2)
		{
			simulation.simulation->eulerian_solver.water_levelset->arr(i, j) = 100;
		}
		else
		{
			simulation.simulation->eulerian_solver.water_levelset->arr(i, j) = 0;
		}*/
		/*if (i == 24 && j == 24)
		{
			simulation.simulation->eulerian_solver.water_levelset->arr(i, j) = 0;
		}
		else
		{
			simulation.simulation->eulerian_solver.water_levelset->arr(i, j) = 100;
		}*/
		simulation.simulation->eulerian_solver.water_levelset->arr(i, j) = ((T)0.1 + POW2(x_min + i*dx - (T)3.5) + POW2(y_min + j*dy - (T)2))*(sqrt(POW2(x_min + i*dx)/16 + POW2(y_min + j*dy)/4) - (T)1);
	}
	
	simulation.simulation->eulerian_solver.water_levelset->FillGhostCellsContinuousDerivativesFromPointer(&(simulation.simulation->eulerian_solver.water_levelset->phi), true);

	// Numerically Exact Solution
	int number_for_exact_point(2000);

	ARRAY<VT> exact_point;
	exact_point.Initialize(number_for_exact_point);

	int index(0);
	
	for (index = 0; index < exact_point.length; index++)
	{
		exact_point[index] = VT((T)4*cos((T)2*PI*index/exact_point.length), (T)2*sin((T)2*PI*index/exact_point.length));
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
		
	T dt = (T)1;

	// Reinitialization
	simulation.simulation->eulerian_solver.Reinitialization(dt);
	
	ofstream fout;
	char filename[100];
	sprintf(filename, "A Remark on Computing Distance Functions/2D Example/Result_2D_%d_%d_w_test", simulation.simulation->eulerian_solver.base_grid.i_res - 1, simulation.simulation->eulerian_solver.base_grid.j_res - 1);
	fout.open(filename); 
	for (int j = simulation.simulation->eulerian_solver.world_discretization->world_grid.j_start; j <= simulation.simulation->eulerian_solver.world_discretization->world_grid.j_end; j++)
	{
		for (int i = simulation.simulation->eulerian_solver.world_discretization->world_grid.i_start; i <= simulation.simulation->eulerian_solver.world_discretization->world_grid.i_end; i++)
			fout << simulation.simulation->eulerian_solver.water_levelset->arr(i, j) << " ";
	
		fout << "\n";
	}
	fout.close();
	
	char filename1[100];
	sprintf(filename1, "A Remark on Computing Distance Functions/2D Example/Exact_%d_%d", simulation.simulation->eulerian_solver.base_grid.i_res - 1, simulation.simulation->eulerian_solver.base_grid.j_res - 1);
	fout.open(filename1); 
	for (int j = simulation.simulation->eulerian_solver.world_discretization->world_grid.j_start; j <= simulation.simulation->eulerian_solver.world_discretization->world_grid.j_end; j++)
	{
		for (int i = simulation.simulation->eulerian_solver.world_discretization->world_grid.i_start; i <= simulation.simulation->eulerian_solver.world_discretization->world_grid.i_end; i++)
			fout << simulation.simulation->eulerian_solver.water_levelset->exact_solution(i, j) << " ";
	
		fout << "\n";
	}
	fout.close();
	
	return 0;
}
