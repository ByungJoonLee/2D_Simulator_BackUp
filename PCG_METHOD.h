#pragma once

#include "LINEAR_SOLVER.h"
#include "CG_METHOD.h"

class PCG_METHOD : public LINEAR_SOLVER
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

public: // Essential data
	VECTOR_ND<T>  res, p, Ap, s;

public: // Subdata for Multithreading
	CSR_MATRIX<T> M;
	VECTOR_ND<T>  y;
	ARRAY<int>	  start_ix_1D_x, end_ix_1D_x, start_ix_1D_y, end_ix_1D_y;

public: // Constructor and Destructor
	PCG_METHOD(void)
	{}

	~PCG_METHOD(void)
	{}

public:	// Initialization Function
	void Intialize(const T& tolerance_input, const int& max_iteration_input, MULTITHREADING* multithreading_input);

public:	// Solver
	void Solve(const CSR_MATRIX<T>& A, VECTOR_ND<T>& x, const VECTOR_ND<T>& b, const FIELD_STRUCTURE_2D<int>& bc);
	void Solve(const CSR_MATRIX<T>& A, VECTOR_ND<T>& x, const VECTOR_ND<T>& b, const FIELD_STRUCTURE_2D<int>& bc, const int& thread_id);
	
	// Preconditioner as a diagonal matrix - easiest one but mediocre 
	void MultiplicationByMinverseAsDiagonal(const CSR_MATRIX<T>& A, VECTOR_ND<T>& x, const VECTOR_ND<T>& b);
	void MultiplicationByMinverseAsDiagonal(const CSR_MATRIX<T>& A, VECTOR_ND<T>& x, const VECTOR_ND<T>& b, const int& thread_id);
	void MultiplicationByMinverse(const int& i_res_input, const int& j_res_input, const CSR_MATRIX<T>& M, VECTOR_ND<T>& x, const VECTOR_ND<T>& b);
	void MultiplicationByMinverse(const int& i_res_input, const int& j_res_input, const CSR_MATRIX<T>& M, VECTOR_ND<T>& x, const VECTOR_ND<T>& b, const int& thread_id);
	
	void PCGMethod(const int& i_res_input, const int& j_res_input, const CSR_MATRIX<T>& A, VECTOR_ND<T>& x, const VECTOR_ND<T>& b);
	void PCGMethod(const int& i_res_input, const int& j_res_input, const CSR_MATRIX<T>& A, VECTOR_ND<T>& x, const VECTOR_ND<T>& b, const int& thread_id);
	
	void SplitDomainIndex1DinXDirection(const int& i_start, const int& i_res);
	void SplitDomainIndex1DinYDirection(const int& j_start, const int& j_res);
	
	int GlobalCoordinateForRedBlockLeftLowerCorner(const int& thread_id, const int& order, const int& direction);
	int GlobalCoordinateForBlackBlockLeftLowerCorner(const int& thread_id, const int& order, const int& direction);
	int GlobalCoordinateForRedBlockLeftUpperCorner(const int& thread_id, const int& order, const int& direction);
	int GlobalCoordinateForRedBlockRightLowerCorner(const int& thread_id, const int& order, const int& direction);
	int GlobalCoordinateForRedBlockRightUpperCorner(const int& thread_id, const int& order, const int& direction);
	int GlobalCoordinateForBlackBlockRightUpperCorner(const int& thread_id, const int& order, const int& direction);
	void ForwardSubstitution(const int& thread_id, const int& i_res_input, const CSR_MATRIX<T>& M, VECTOR_ND<T>& x, const VECTOR_ND<T>& b);
};