#pragma once

#include "CSR_MATRIX.h"
#include "CG_METHOD.h"
#include "PROJECTION_2D.h"
#include "SIMULATION_MANAGER.h"
#include "PROJECT_INFO.h"

int main()
{
	LEVELSET_2D levelset;
	levelset.Initialize(0, 0, 160, 160, 0, 0, 1, 1, 2);
	
	T x_min = levelset.grid.x_min, y_min = levelset.grid.y_min;
	T dx = levelset.grid.dx, dy = levelset.grid.dy;

	GRID_ITERATION_2D(levelset.grid)
	{
		levelset(i, j) = sqrt(POW2(x_min + i*dx - 0.5) + POW2(y_min + j*dy - 0.5)) - 0.25;
	}
	
	FIELD_STRUCTURE_2D<T> f;
	f.Initialize(levelset.grid, 2);
	
	GRID_ITERATION_2D(f.grid)
	{
		if (levelset(i, j) <= 0)
		{
			f(i, j) = 8*(POW2(x_min + i*dx) + POW2(y_min + j*dy) - 1)*exp(-POW2(x_min + i*dx) - POW2(y_min + j*dy));
		}
		else
		{
			f(i, j) = 0;
		}
	}

	FIELD_STRUCTURE_2D<T> a, b;
	a.Initialize(levelset.grid, 2);
	b.Initialize(levelset.grid, 2);

	GRID_ITERATION_2D(a.grid)
	{
		a(i, j) = -exp(-POW2(x_min + i*dx) - POW2(y_min + j*dy));
		b(i, j) = 8*(2*POW2(x_min + i*dx) + 2*POW2(y_min + j*dy) - (x_min + i*dx + y_min + j*dy))*exp(-POW2(x_min + i*dx) - POW2(y_min + j*dy));
	}

	FIELD_STRUCTURE_2D<T> beta;
	beta.Initialize(levelset.grid, 2);

	GRID_ITERATION_2D(beta.grid)
	{
		if (levelset(i, j) <= 0)
		{
			beta(i, j) = 0.5;
		}
		else
		{
			beta(i, j) = 1;
		}
	}

	POISSON_SOLVER_2D poisson_solver_test;
	poisson_solver_test.Initialize(10e-6, 10000);
	poisson_solver_test.InitializeLinearSolver(PCG);

	FIELD_STRUCTURE_2D<T> pressure;
	pressure.Initialize(levelset.grid, 2);
	pressure.AssignAllValue((T)0);

	FIELD_STRUCTURE_2D<int> bc;
	bc.Initialize(levelset.grid, 2);
	bc.AssignAllValue((int)0);

	levelset.ComputeNormals();
			
	int i, j;
	LOOPS_2D(i, j, bc.i_start_g, bc.j_start_g, bc.i_end_g, bc.j_end_g)
	{
		if (i < bc.i_start || i > bc.i_end || j < bc.j_start || j > bc.j_end)
		{
			bc(i, j) = BC_DIR;
			pressure(i, j) = 0;
		}
	}

	poisson_solver_test.Solve(pressure, beta, bc, f, levelset, a, b);
	
	FIELD_STRUCTURE_2D<T> true_solution;
	true_solution.Initialize(levelset.grid, 2);

	GRID_ITERATION_2D(true_solution.grid)
	{
		if (levelset(i, j) <= 0)
		{
			true_solution(i, j) = exp(-POW2(x_min + i*dx) - POW2(y_min + j*dy));
		}
		else
		{
			true_solution(i, j) = 0;
		}
	}

	ofstream fout;
	fout.open("pressure");
	for (int j = levelset.grid.j_start; j <= levelset.grid.j_end; j++)
	{
		for (int i = levelset.grid.i_start; i <= levelset.grid.i_end; i++)
			fout << pressure(i,j) << " ";

			fout << "\n";
	}
	fout.close();

	fout.open("true");
	for (int j = levelset.grid.j_start; j <= levelset.grid.j_end; j++)
	{
		for (int i = levelset.grid.i_start; i <= levelset.grid.i_end; i++)
			fout << true_solution(i,j) << " ";

			fout << "\n";
	}
	fout.close();

	return 0;
}
