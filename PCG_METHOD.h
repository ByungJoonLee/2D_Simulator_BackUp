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

public: // Essential data
	VECTOR_ND<T>  res, p, Ap, s;

public: // Constructor and Destructor
	PCG_METHOD(void)
	{}

	~PCG_METHOD(void)
	{}

public:	// Initialization Function
	void Intialize(const T& tolerance_input, const int& max_iteration_input)
	{
		BASE::Initialize(tolerance_input, max_iteration_input);
	}

public:	// Solver
	void Solve(const CSR_MATRIX<T>& A, VECTOR_ND<T>& x, const VECTOR_ND<T>& b, const FIELD_STRUCTURE_2D<int>& bc)
	{
		PCGMethod(bc.grid.i_res, bc.grid.j_res, A, x, b);
	}

	// Preconditioner as a diagonal matrix - easiest one but mediocre 
	void MultiplicationByMinverseAsDiagonal(const CSR_MATRIX<T>& A, VECTOR_ND<T>& x, const VECTOR_ND<T>& b)
	{
		const int N = A.N;
				
		CSR_MATRIX<T> M;
		M.Initialize(N, N);
		
		M = DiagonalPreconditioner(A);
		
		for (int i = 0; i < x.num_dimension; i++)
		{
			x[i] = b[i]/M(i, i);	
		}	
	}

	void MultiplicationByMinverse(const int& i_res_input, const int& j_res_input, const CSR_MATRIX<T>& M, VECTOR_ND<T>& x, const VECTOR_ND<T>& b)
	{
		VECTOR_ND<T> y;
		y.Initialize(x.num_dimension, true);

		T one_over_M_start = 1/M(0, 0);
		
		y[0] = b[0]*one_over_M_start;
		
		int number(0), num_2(0);
		for (int i = 1; i < y.num_dimension; i++)
		{
			T sum(0);
			
			for (int k = M.row_ptr[i]; k < (M.row_ptr[i + 1] - 1); k++)
			{
				sum += M.values[k]*y[M.column_index[k]];
			}
			
			T one_over_Mii = 1/M(i, i);
			y[i] = (b[i] - sum)*one_over_Mii;
		}

		T* summation = new T[x.num_dimension];

		for (int i = 0; i < x.num_dimension; i++)
		{
			summation[i] = 0;
		}

		// Matrix-Transpose-Vector Multiplication
		// Parallel Sparse Matrix-Vector and Matrix-Trnaspose-Vector Multiplication Using Compressed Sparse Blocks
		T one_over_M_end = 1/M(x.num_dimension - 1, x.num_dimension - 1);
		x[x.num_dimension - 1] = y[x.num_dimension - 1]*one_over_M_end;
		for (int i = x.num_dimension - 1; i >= 1; i--)
		{
			for (int k = M.row_ptr[i + 1] - 1; k >= M.row_ptr[i]; k--)
			{
				if (k != M.row_ptr[i + 1] - 1)
				{
					summation[M.column_index[k]] += M.values[k]*x[i];
				}
			}
			
			T one_over_M = 1/M(i - 1, i - 1);
			x[i - 1] = (y[i - 1] - summation[i - 1])*one_over_M;
		}
		
		delete summation;
	}
	
	void PCGMethod(const int& i_res_input, const int& j_res_input, const CSR_MATRIX<T>& A, VECTOR_ND<T>& x, const VECTOR_ND<T>& b)
	{
		const int N(x.num_dimension);

		res.Initialize(N);

		s.Initialize(N);
		
		p.Initialize(N);

		Ap.Initialize(N);
		
		num_iteration = 0;
		
		T *rval(res.values), *sval(s.values), *pval(p.values), *Apval(Ap.values), *xval(x.values);

		T alpha, beta, res_old, res_new, dot_result;

		A.ComputeResidual(x, b, res);
		
		const int nz = A.nz;
		
		const int num_of_nz = (nz - N)/2 + N;
		
		CSR_MATRIX<T> M;
		M.Initialize(N, num_of_nz);
		
		IncompleteCholeskyDecomposition(i_res_input, j_res_input, A, M);

		//MultiplicationByMinverseAsDiagonal(A, p, res);
		MultiplicationByMinverse(i_res_input, j_res_input, M, p, res);
		
        DotProduct(res, p, res_new);

		while (true)
		{
			A.Multiply(p, Ap);

            DotProduct(p, Ap, dot_result);

            alpha = res_new/dot_result;

            for (int i = 0; i < N; i++)
			{
				xval[i] += alpha*pval[i];
				rval[i] -= alpha*Apval[i];
			}
			
			/*if (num_iteration % (int)sqrt(A.N) == 0)
			{
				A.ComputeResidual(x, b, res);
			}
			else
			{
				for (int i = 0; i < N; i++)
				{
					rval[i] -= alpha*Apval[i];
				}
			}*/
						
			//MultiplicationByMinverseAsDiagonal(A, s, res);
			MultiplicationByMinverse(i_res_input, j_res_input, M, s, res);
			
			res_old = res_new;

			DotProduct(res, s, res_new);
            
			beta = res_new/res_old;

            for (int i = 0; i < N; i++)
			{
				pval[i] *= beta;
				pval[i] += sval[i];
			}
			
			num_iteration++;
			
			if(num_iteration > max_iteration) break;
			if (res_new < sqr_tolerance)
			{
				cout << "Converge!!" << endl;
				cout << "Iteration Number : " << num_iteration << endl;
				break;
			}
		}
        
		residual = sqrt(res_new);
		cout << "Residual: " << residual << endl;
	}
};