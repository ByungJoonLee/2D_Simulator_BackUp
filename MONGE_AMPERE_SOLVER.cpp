#include "MONGE_AMPERE_SOLVER.h"

void MONGE_AMPERE_SOLVER::InitializeLinearSolver(const POISSON_SOLVER_TYPE linear_solver_type)
{
	DELETE_POINTER(sub_linear_solver);

	switch (linear_solver_type)
	{
	case GS:
		sub_linear_solver = new GAUSS_SEIDEL_METHOD();
		break;
	case BICG:
		sub_linear_solver = new BICGSTAB_METHOD();
		break;

	case NO_SOLVER:

	default:
		sub_linear_solver = 0;
		break;
	}

	sub_linear_solver->Initialize(tolerance, max_iteration_for_sub_linear_solver, multithreading);
}

