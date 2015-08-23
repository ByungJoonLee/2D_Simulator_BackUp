#pragma once

#include "LINEAR_SOLVER.h"
#include "CSR_MATRIX.h"
#include "VECTOR_ND.h"
#include "LINEAR_SOLVER.h"
#include "FIELD_STRUCTURE_2D.h"

class GAUSS_SEIDEL_METHOD : public LINEAR_SOLVER
{
public: // Typedef
	typedef LINEAR_SOLVER BASE;

public: // Using Keyword
	using BASE::tolerance;
	using BASE::sqr_tolerance;
	using BASE::residual;
	using BASE::max_iteration;
	using BASE::num_iteration;
	using BASE::multithreading;

public: // Constructor and Destructor
	GAUSS_SEIDEL_METHOD(void)
	{}

	~GAUSS_SEIDEL_METHOD(void)
	{}

public: 
	void Solve(const CSR_MATRIX<T>& A, VECTOR_ND<T>& x, const VECTOR_ND<T>& b, const FIELD_STRUCTURE_2D<int>& bc, const int& thread_id);

	void GaussSeidelMethod(const CSR_MATRIX<T>& A, VECTOR_ND<T>& x, const VECTOR_ND<T>& b, const int& thread_id);
	T    GaussSeidelStep(const CSR_MATRIX<T>& A, VECTOR_ND<T>& x, const VECTOR_ND<T>& b, const int& thread_id);
};