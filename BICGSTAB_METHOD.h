#pragma once

#include "LINEAR_SOLVER.h"
#include "VECTOR_ND.h"

class BICGSTAB_METHOD : public LINEAR_SOLVER
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

public: // Essential Data
	VECTOR_ND<T> res, rtilde, p, Ap, s, As;

public: // Constructor and Destructor
	BICGSTAB_METHOD(void)
	{}

	~BICGSTAB_METHOD(void)
	{}

public: // Solver
	void Solve(const CSR_MATRIX<T>& A, VECTOR_ND<T>& x, const VECTOR_ND<T>& b, const FIELD_STRUCTURE_2D<int>& bc);
	void Solve(const CSR_MATRIX<T>& A, VECTOR_ND<T>& x, const VECTOR_ND<T>& b, const FIELD_STRUCTURE_2D<int>& bc, const int& thread_id);
	
	void BICGSTABMethod(const CSR_MATRIX<T>& A, VECTOR_ND<T>& x, const VECTOR_ND<T>& b);
	void BICGSTABMethod(const CSR_MATRIX<T>& A, VECTOR_ND<T>& x, const VECTOR_ND<T>& b, const int& thread_id);
};