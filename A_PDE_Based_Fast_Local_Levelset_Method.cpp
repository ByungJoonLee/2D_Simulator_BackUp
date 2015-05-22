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
	
	GRID_ITERATION_2D(simulation.simulation->eulerian_solver.world_discretization->world_grid)
	{
		T x_min = simulation.simulation->eulerian_solver.world_discretization->world_grid.x_min, y_min = simulation.simulation->eulerian_solver.world_discretization->world_grid.y_min;
		T dx = simulation.simulation->eulerian_solver.world_discretization->world_grid.dx, dy = simulation.simulation->eulerian_solver.world_discretization->world_grid.dy;
		
		simulation.simulation->eulerian_solver.water_levelset->arr(i, j) = y_min + j*dy + (T)0.05*sin(PI*(x_min + i*dx));
	}
	
	simulation.simulation->eulerian_solver.water_levelset->FillGhostCellsContinuousDerivativesFromPointer(&(simulation.simulation->eulerian_solver.water_levelset->phi), true);

	T dt = (T)0.1*simulation.simulation->eulerian_solver.base_grid.dx;

	ofstream fout;
	fout.open("A PDE Based Fast Local Level Set Method/Example3/initial"); 
	for (int j = simulation.simulation->eulerian_solver.world_discretization->world_grid.j_start; j <= simulation.simulation->eulerian_solver.world_discretization->world_grid.j_end; j++)
	{
		for (int i = simulation.simulation->eulerian_solver.world_discretization->world_grid.i_start; i <= simulation.simulation->eulerian_solver.world_discretization->world_grid.i_end; i++)
			fout << simulation.simulation->eulerian_solver.water_levelset->arr(i, j) << " ";
	
		fout << "\n";
	}
	fout.close();

	int iter(0);
	
	EULERIAN_FLUID_SOLVER_2D& euler(simulation.simulation->eulerian_solver);

	int index(0);

	while(iter < 2880)
	{
		cout << "---------------------" << iter << "---------------------" << endl;
		
		// Simulating one time step
		euler.water_projection->DeterminePressure();		
		
		GRID_ITERATION_2D(euler.water_velocity_field->grid)
		{
			euler.water_velocity_field->array_for_this(i, j).x = euler.water_velocity_field->one_over_2dy*(euler.water_projection->pressure_field(i, j + 1) - euler.water_projection->pressure_field(i, j - 1));
			euler.water_velocity_field->array_for_this(i, j).y = -euler.water_velocity_field->one_over_2dx*(euler.water_projection->pressure_field(i + 1, j) - euler.water_projection->pressure_field(i - 1, j));
		}

		euler.vector_field_ghost.FillGhostCellsFrom(euler.water_velocity_field->array_for_this, true);

		/*GRID_ITERATION_2D(euler.water_velocity_field_mac_x->grid)
		{
			euler.water_velocity_field_mac_x->array_for_this(i, j) = euler.water_velocity_field->one_over_dx*(euler.water_projection->pressure_field(i, j) - euler.water_projection->pressure_field(i - 1, j));
		}
		euler.vector_field_mac_ghost_x.FillGhostCellsContinuousDerivativesFrom(euler.water_velocity_field_mac_x->array_for_this, true);

		GRID_ITERATION_2D(euler.water_velocity_field_mac_y->grid)
		{
			euler.water_velocity_field_mac_y->array_for_this(i, j) = -euler.water_velocity_field->one_over_dy*(euler.water_projection->pressure_field(i, j) - euler.water_projection->pressure_field(i, j - 1));
		}
		euler.vector_field_mac_ghost_y.FillGhostCellsContinuousDerivativesFrom(euler.water_velocity_field_mac_y->array_for_this, true);*/

		euler.AdvanceOneTimeStep(dt);

		iter = iter + 1;
		
		if ((int)(iter%((simulation.simulation->eulerian_solver.base_grid.i_res - 1)/2)) == 0)
		{
			index = index + 1;

			char filename[50];
			sprintf(filename, "A PDE Based Fast Local Level Set Method/Example3/%d", index);
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

	return 0;
}
