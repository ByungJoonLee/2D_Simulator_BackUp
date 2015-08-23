#include "LINEAR_SOLVER.h"

// Initialization Funciton
void LINEAR_SOLVER::Initialize(const T& tolerance_input, const int& max_iteration_input, MULTITHREADING* multithreading_input)
{
	multithreading = multithreading_input;
	max_iteration = max_iteration_input;
	SetTolerance(tolerance_input);
}

// Member Functions
void LINEAR_SOLVER::Solve(const CSR_MATRIX<T>& A, VECTOR_ND<T>& x, const VECTOR_ND<T>& b, const FIELD_STRUCTURE_1D<int>& bc, const int& thread_id)
{
	cout << "virtual LINEAR_SOLVER::Solve" << endl;
	exit(1);
}

void LINEAR_SOLVER::Solve(const CSR_MATRIX<T>& A, VECTOR_ND<T>& x, const VECTOR_ND<T>& b, const FIELD_STRUCTURE_1D<int>& bc)
{
	cout << "virtual LINEAR_SOLVER::Solve" << endl;
	exit(1);
}

void LINEAR_SOLVER::Solve(const CSR_MATRIX<T>& A, VECTOR_ND<T>& x, const VECTOR_ND<T>& b, const FIELD_STRUCTURE_2D<int>& bc, const int& thread_id)
{
	cout << "virtual LINEAR_SOLVER::Solve" << endl;
	exit(1);
}

void LINEAR_SOLVER::Solve(const CSR_MATRIX<T>& A, VECTOR_ND<T>& x, const VECTOR_ND<T>& b, const FIELD_STRUCTURE_2D<int>& bc)
{
	cout << "virtual LINEAR_SOLVER::Solve" << endl;
	exit(1);
}

void LINEAR_SOLVER::SetTolerance(const T& tolerance_input)
{
	tolerance = tolerance_input;
	sqr_tolerance = POW2(tolerance);
}