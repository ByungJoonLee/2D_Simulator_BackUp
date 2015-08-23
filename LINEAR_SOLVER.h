#pragma once

#include "FIELD_STRUCTURE_1D.h"
#include "FIELD_STRUCTURE_2D.h"
#include "CSR_MATRIX.h"

class LINEAR_SOLVER
{
public: // Essential Data
	T					tolerance, sqr_tolerance;
	T					residual;							// Residual of previous solving
	int					max_iteration;
	int					num_iteration;						// Number of iterations of pervious solving

	MULTITHREADING*	    multithreading;

public: // Constructors and Destructor
	LINEAR_SOLVER(void)
		: tolerance((T)1), sqr_tolerance(tolerance*tolerance), residual((T)1e8), max_iteration(10), num_iteration(0), multithreading(0)
	{}

	~LINEAR_SOLVER(void)
	{}

public: // Initialization Function
	void Initialize(const T& tolerance_input, const int& max_iteration_input, MULTITHREADING* multithreading_input);

public: // Member Functions
	virtual void Solve(const CSR_MATRIX<T>& A, VECTOR_ND<T>& x, const VECTOR_ND<T>& b, const FIELD_STRUCTURE_1D<int>& bc, const int& thread_id);
	virtual void Solve(const CSR_MATRIX<T>& A, VECTOR_ND<T>& x, const VECTOR_ND<T>& b, const FIELD_STRUCTURE_1D<int>& bc);
	virtual void Solve(const CSR_MATRIX<T>& A, VECTOR_ND<T>& x, const VECTOR_ND<T>& b, const FIELD_STRUCTURE_2D<int>& bc, const int& thread_id);
	virtual void Solve(const CSR_MATRIX<T>& A, VECTOR_ND<T>& x, const VECTOR_ND<T>& b, const FIELD_STRUCTURE_2D<int>& bc);
	//virtual void Solve(const CSR_MATRIX<T>& A, VECTOR_ND<T>& x, VECTOR_ND<T>& b, const FIELD_STRUCTURE_2D<int>& bc, const int& thread_id);
	void SetTolerance(const T& tolerance_input);
};

