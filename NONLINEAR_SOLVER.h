#pragma once

#include "FIELD_STRUCTURE_1D.h"
#include "FIELD_STRUCTURE_2D.h"
#include "CSR_MATRIX.h"

class NONLINEAR_SOLVER
{
public: // Essential Data
	T					tolerance, sqr_tolerance;
	T					residual;
	int					max_iteration;
	int					num_iteration;

	MULTITHREADING*		multithreading;

public: // Constructor and Destructor
	NONLINEAR_SOLVER(void)
		: tolerance((T)1), sqr_tolerance(tolerance*tolerance), residual((T)1e8), max_iteration(10), num_iteration(0), multithreading(0)
	{}

	NONLINEAR_SOLVER(void)
	{}

public: // Initialization Function
	void Initialize(const T& tolerance_input, const int& max_iteration_input, MULTITHREADING* multithreading_input)
	{
		multithreading = multithreading_input;
		max_iteration = max_iteration_input;
		SetTolerance(tolerance_input);
	}

public: // Member Function
	void SetTolerance(const T& tolerance_input)
	{
		tolerance = tolerance_input;
		sqr_tolerance = POW2(tolerance);
	}
};