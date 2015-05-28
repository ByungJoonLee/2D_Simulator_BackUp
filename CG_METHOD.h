#pragma once

#include "LINEAR_SOLVER.h"

class CG_METHOD : public LINEAR_SOLVER
{
public: // Typedef
	typedef LINEAR_SOLVER BASE;

public: // Using keyword
	using BASE::tolerance;
	using BASE::sqr_tolerance;
	using BASE::residual;
	using BASE::max_iteration;
	using BASE::num_iteration;
	using BASE::multithreading;

public: // Essential Data
	VECTOR_ND<T> res, p, Ap;

public: // Constructor and Destructor
	CG_METHOD(void)
	{}

	~CG_METHOD(void)
	{}

public: // Solver
	void Solve(const CSR_MATRIX<T>& A, VECTOR_ND<T>& x, const VECTOR_ND<T>& b, const FIELD_STRUCTURE_2D<int>& bc)
	{
		CGMethod(A, x, b);
	}

	void Solve(const CSR_MATRIX<T>& A, VECTOR_ND<T>& x, const VECTOR_ND<T>& b, const FIELD_STRUCTURE_2D<int>& bc, const int& thread_id)
	{
		CGMethod(A, x, b, thread_id);
	}

	void CGMethod(const CSR_MATRIX<T>& A, VECTOR_ND<T>& x, const VECTOR_ND<T>& b);
	void CGMethod(const CSR_MATRIX<T>& A, VECTOR_ND<T>& x, const VECTOR_ND<T>& b, const int& thread_id);
};


