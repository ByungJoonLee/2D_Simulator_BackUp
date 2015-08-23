#pragma once

#include "NONLINEAR_SOLVER.h"

class NEWTONS_METHOD : public NONLINEAR_SOLVER
{
public: // Typedef
	typedef NONLINEAR_SOLVER BASE;

public: // Using Keyword
	using BASE::tolerance;
	using BASE::sqr_tolerance;
	using BASE::residual;
	using BASE::max_iteration;
	using BASE::num_iteration;
	using BASE::multithreading;

public: // Constructor and Destructor
	NEWTONS_METHOD(void)
	{}

	~NEWTONS_METHOD(void)
	{}

public: // Solver
	void OneStepUpdate(const CSR_MATRIX<T>& A, VECTOR_ND<T>& x, const VECTOR_ND<T>& b)
	{
		
	}
};