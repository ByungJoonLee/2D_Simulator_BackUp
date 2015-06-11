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
	void Intialize(const T& tolerance_input, const int& max_iteration_input, MULTITHREADING* multithreading_input)
	{
		BASE::Initialize(tolerance_input, max_iteration_input, multithreading_input);
	}

public:	// Solver
	void Solve(const CSR_MATRIX<T>& A, VECTOR_ND<T>& x, const VECTOR_ND<T>& b, const FIELD_STRUCTURE_2D<int>& bc)
	{
		PCGMethod(bc.grid.i_res, bc.grid.j_res, A, x, b);
	}

	void Solve(const CSR_MATRIX<T>& A, VECTOR_ND<T>& x, const VECTOR_ND<T>& b, const FIELD_STRUCTURE_2D<int>& bc, const int& thread_id)
	{
		PCGMethod(bc.grid.i_res, bc.grid.j_res, A, x, b, thread_id);
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

	void MultiplicationByMinverseAsDiagonal(const CSR_MATRIX<T>& A, VECTOR_ND<T>& x, const VECTOR_ND<T>& b, const int& thread_id)
	{
		const int N = A.N;

		BEGIN_HEAD_THREAD_WORK
		{
			M.Initialize(N, N, multithreading);
			multithreading->SplitDomainIndex1D(0, N);
		
			M.start_ix[0] = 0;
			M.end_ix[0] = multithreading->sync_value_int[0] - 1;
			M.prev_row_array[0] = -1;
			M.values_ix_array[0] = 0;
			for (int id = 1; id < multithreading->num_threads; id++)
			{
				M.start_ix[id] = M.end_ix[id - 1] + 1;
				M.end_ix[id] = M.end_ix[id - 1] + multithreading->sync_value_int[id];
				M.prev_row_array[id] = -1;
				M.values_ix_array[id] = M.start_ix[id];
			}
		}
		END_HEAD_THREAD_WORK

		DiagonalPreconditioner(A, M, multithreading, thread_id);
		
		const int start_ix(multithreading->start_ix_1D[thread_id]), end_ix(multithreading->end_ix_1D[thread_id]);

		for (int i = start_ix; i <= end_ix; i++)
		{
			x[i] = b[i]/M(i, i);
		}
		multithreading->Sync(thread_id);
	}

	void MultiplicationByMinverse(const int& i_res_input, const int& j_res_input, const CSR_MATRIX<T>& M, VECTOR_ND<T>& x, const VECTOR_ND<T>& b)
	{
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
	
	void MultiplicationByMinverse(const int& i_res_input, const int& j_res_input, const CSR_MATRIX<T>& M, VECTOR_ND<T>& x, const VECTOR_ND<T>& b, const int& thread_id)
	{
		BEGIN_HEAD_THREAD_WORK
		{
			y.Initialize(x.num_dimension, true);
		}
		END_HEAD_THREAD_WORK;

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

		//ForwardSubstitution(thread_id, i_res_input, M, y, b);

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

	void PCGMethod(const int& i_res_input, const int& j_res_input, const CSR_MATRIX<T>& A, VECTOR_ND<T>& x, const VECTOR_ND<T>& b, const int& thread_id)
	{
		BEGIN_HEAD_THREAD_WORK
		{
			multithreading->SplitDomainIndex1D(0, A.N);
			SplitDomainIndex1DinYDirection(0, j_res_input);
			SplitDomainIndex1DinXDirection(0, i_res_input);
		}
		END_HEAD_THREAD_WORK
	
		const int N(x.num_dimension), start_ix(multithreading->start_ix_1D[thread_id]), end_ix(multithreading->end_ix_1D[thread_id]), num_of_red_block(2*multithreading->num_threads), num_of_black_block(2*multithreading->num_threads);

		BEGIN_HEAD_THREAD_WORK
		{
			res.Initialize(N);
		}
		END_HEAD_THREAD_WORK

		BEGIN_HEAD_THREAD_WORK
		{
			s.Initialize(N);
		}
		END_HEAD_THREAD_WORK

		BEGIN_HEAD_THREAD_WORK
		{
			p.Initialize(N);
		}
		END_HEAD_THREAD_WORK
	
		BEGIN_HEAD_THREAD_WORK
		{
			Ap.Initialize(N);
		}
		END_HEAD_THREAD_WORK

		BEGIN_HEAD_THREAD_WORK
		{
			num_iteration = 0;
		}
		END_HEAD_THREAD_WORK

		T *rval(res.values), *pval(p.values), *Apval(Ap.values), *xval(x.values), *sval(s.values);

		T alpha, beta, res_old, res_new, dot_result;

		A.ComputeResidual(x, b, res, thread_id);

		/*BEGIN_HEAD_THREAD_WORK
		{
			IncompleteCholeskyDecomposition(i_res_input, j_res_input, A, M);
		}
		END_HEAD_THREAD_WORK;*/

		//IncompleteCholeskyDecomposition(multithreading, i_res_input, j_res_input, k_res_input, A, M);

		//MultiplicationByMinverse(i_res_input, j_res_input, M, p, res, thread_id);

		MultiplicationByMinverseAsDiagonal(A, p, res, thread_id);

		DotProduct(res, p, res_new, multithreading, thread_id);
		
		while (num_iteration < max_iteration)
		{
			A.Multiply(p, Ap, thread_id);

			DotProduct(p, Ap, dot_result, multithreading, thread_id);
						
			alpha = res_new/dot_result;

			for (int i = start_ix; i <= end_ix; i++)
			{
				xval[i] += alpha*pval[i];
			}
			multithreading->Sync(thread_id);
		
			for (int i = start_ix; i <= end_ix; i++)
			{
				rval[i] -= alpha*Apval[i];
			}
			multithreading->Sync(thread_id);

			//MultiplicationByMinverse(i_res_input, j_res_input, M, p, res, thread_id);
			MultiplicationByMinverseAsDiagonal(A, s, res, thread_id);

			res_old = res_new;

			DotProduct(res, s, res_new, multithreading, thread_id);

			if (res_new < sqr_tolerance)
			{
				break;
			}

			beta = res_new/res_old;

			for (int i = start_ix; i <= end_ix; i++)
			{
				pval[i] *= beta;
				pval[i] += sval[i];
			}
			multithreading->Sync(thread_id);

			BEGIN_HEAD_THREAD_WORK
			{
				num_iteration++;
			}
			END_HEAD_THREAD_WORK;
		}
		multithreading->Sync(thread_id);
        
		BEGIN_HEAD_THREAD_WORK
		{
			residual = sqrt(res_new);
			cout << "Residual: " << residual << endl;
		}
		END_HEAD_THREAD_WORK;
	}

	void SplitDomainIndex1DinXDirection(const int& i_start, const int& i_res)
	{
		const int num_threads = multithreading->num_threads;
		const int i_end = i_start + i_res - 1;
		const int quotient = i_res / num_threads;
		const int remainder = i_res % num_threads;

		start_ix_1D_x.Initialize(num_threads);
		end_ix_1D_x.Initialize(num_threads);

		int i_start_p = i_start;

		for(int i = 0; i < num_threads; i++)
		{
			int i_depth = i < remainder ? (quotient + 1) : quotient;
			start_ix_1D_x[i] = i_start_p;
			end_ix_1D_x[i] = i_start_p + i_depth - 1;

			i_start_p += i_depth;
		}
	}

	void SplitDomainIndex1DinYDirection(const int& j_start, const int& j_res)
	{
		const int num_threads = multithreading->num_threads;
		const int j_end = j_start + j_res - 1;
		const int quotient = j_res / num_threads;
		const int remainder = j_res % num_threads;

		start_ix_1D_y.Initialize(num_threads);
		end_ix_1D_y.Initialize(num_threads);

		int j_start_p = j_start;

		for(int i = 0; i < num_threads; i++)
		{
			int j_depth = i < remainder ? (quotient + 1) : quotient;
			start_ix_1D_y[i] = j_start_p;
			end_ix_1D_y[i] = j_start_p + j_depth - 1;

			j_start_p += j_depth;
		}
	}

	int GlobalCoordinateForRedBlockLeftLowerCorner(const int& thread_id, const int& order, const int& direction)
	{
		// x-direction
		if (direction == (int)1)
		{
			if (thread_id % 2 == 0)
			{
				return start_ix_1D_x[2*order];
			}
			else
			{
				return start_ix_1D_x[2*order + 1];
			}
		}
		else if (direction == (int)2)
		{
			if (thread_id % 2 == 0)
			{
				return start_ix_1D_y[2*order];
			}
			else
			{
				return start_ix_1D_y[2*order + 1];
			}
		}
	}

	int GlobalCoordinateForBlackBlockLeftLowerCorner(const int& thread_id, const int& order, const int& direction)
	{
		// x-direction
		if (direction == (int)1)
		{
			if (thread_id % 2 == 0)
			{
				return start_ix_1D_x[2*order + 1];
			}
			else
			{
				return start_ix_1D_x[2*order];
			}
		}
		else if (direction == (int)2)
		{
			if (thread_id % 2 == 0)
			{
				return start_ix_1D_y[2*order + 1];
			}
			else
			{
				return start_ix_1D_y[2*order];
			}
		}
	}

	int GlobalCoordinateForRedBlockLeftUpperCorner(const int& thread_id, const int& order, const int& direction)
	{
		// x-direction
		if (direction == (int)1)
		{
			if (thread_id % 2 == 0)
			{
				return start_ix_1D_x[2*order];
			}
			else
			{
				return start_ix_1D_x[2*order + 1];
			}
		}
		else if (direction == (int)2)
		{
			if (thread_id % 2 == 0)
			{
				return end_ix_1D_y[2*order];
			}
			else
			{
				return end_ix_1D_y[2*order + 1];
			}
		}
	}

	int GlobalCoordinateForRedBlockRightLowerCorner(const int& thread_id, const int& order, const int& direction)
	{
		// x-direction
		if (direction == (int)1)
		{
			if (thread_id % 2 == 0)
			{
				return end_ix_1D_x[2*order];
			}
			else
			{
				return end_ix_1D_x[2*order + 1];
			}
		}
		else if (direction == (int)2)
		{
			if (thread_id % 2 == 0)
			{
				return start_ix_1D_y[2*order];
			}
			else
			{
				return start_ix_1D_y[2*order + 1];
			}
		}
	}

	int GlobalCoordinateForRedBlockRightUpperCorner(const int& thread_id, const int& order, const int& direction)
	{
		// x-direction
		if (direction == (int)1)
		{
			if (thread_id % 2 == 0)
			{
				return end_ix_1D_x[2*order];
			}
			else
			{
				return end_ix_1D_x[2*order + 1];
			}
		}
		else if (direction == (int)2)
		{
			if (thread_id % 2 == 0)
			{
				return end_ix_1D_y[2*order];
			}
			else
			{
				return end_ix_1D_y[2*order + 1];
			}
		}
	}

	int GlobalCoordinateForBlackBlockRightUpperCorner(const int& thread_id, const int& order, const int& direction)
	{
		// x-direction
		if (direction == (int)1)
		{
			if (thread_id % 2 == 0)
			{
				return end_ix_1D_x[2*order + 1];
			}
			else
			{
				return end_ix_1D_x[2*order];
			}
		}
		else if (direction == (int)2)
		{
			if (thread_id % 2 == 0)
			{
				return end_ix_1D_y[2*order + 1];
			}
			else
			{
				return end_ix_1D_y[2*order];
			}
		}
	}

	void ForwardSubstitution(const int& thread_id, const int& i_res_input, const CSR_MATRIX<T>& M, VECTOR_ND<T>& x, const VECTOR_ND<T>& b)
	{
		// Note that we assume number of threads must be power of 2
		const int num_of_red_block(multithreading->num_threads/2), num_of_black_block(multithreading->num_threads/2);
		int llx, lly, rux, ruy;

		// Red Block
		for (int i = 0; i < num_of_red_block; i++)
		{
			llx = GlobalCoordinateForRedBlockLeftLowerCorner(thread_id, i, 1);
			lly = GlobalCoordinateForRedBlockLeftLowerCorner(thread_id, i, 2);
			rux = GlobalCoordinateForRedBlockRightUpperCorner(thread_id, i, 1);
			ruy = GlobalCoordinateForRedBlockRightUpperCorner(thread_id, i, 2);

			// Domain 1
			int i_1D;
			i_1D = lly*i_res_input + llx;
			
			T one_over_M_start = 1/M(i_1D, i_1D);
		
			y[i_1D] = b[i_1D]*one_over_M_start;

			// Domain 2
			int number(0), num_2(0);
			for (int j = i_1D + 1; j <= lly*i_res_input + rux; j++)
			{
				T sum(0);
			
				for (int k = M.row_ptr[j]; k < (M.row_ptr[j + 1] - 1); k++)
				{
					sum += M.values[k]*y[M.column_index[k]];
				}
			    
				T one_over_Mii = 1/M(j, j);
				y[j] = (b[j] - sum)*one_over_Mii;
			}

			// Domain 3
			number = 0, num_2 = 0;
			for (int j = (lly+1)*i_res_input + llx; j <= lly*i_res_input + rux; j += i_res_input)
			{
				T sum(0);
			
				for (int k = M.row_ptr[j]; k < (M.row_ptr[j + 1] - 1); k++)
				{
					sum += M.values[k]*y[M.column_index[k]];
				}
			
				T one_over_Mii = 1/M(j, j);
				y[j] = (b[j] - sum)*one_over_Mii;
			}

			// Domain 4
			number = 0, num_2 = 0;
			for (int j = lly + 1; j <= ruy; j++)
			{
				for (int k = llx + 1; k <= rux; k++)
				{
					T sum(0);

					int t = j*i_res_input + k;
					
					for (int ks = M.row_ptr[t]; ks < (M.row_ptr[t + 1] - 1); ks++)
					{
						sum += M.values[ks]*y[M.column_index[ks]];
					}

					T one_over_Mii = 1/M(t, t);
					y[t] = (b[t] - sum)*one_over_Mii; 
				}
			}
		}

		// Black Block
		for (int i = 0; i < num_of_black_block; i++)
		{
			llx = GlobalCoordinateForBlackBlockLeftLowerCorner(thread_id, i, 1);
			lly = GlobalCoordinateForBlackBlockLeftLowerCorner(thread_id, i, 2);
			rux = GlobalCoordinateForBlackBlockRightUpperCorner(thread_id, i, 1);
			ruy = GlobalCoordinateForBlackBlockRightUpperCorner(thread_id, i, 2);

			// Domain 5
			int number = 0, num_2 = 0;
			for (int j = lly; j <= ruy - 1; j++)
			{
				for (int k = llx; k <= rux - 1; k++)
				{
					T sum(0);

					int t = j*i_res_input + k;
					
					for (int ks = M.row_ptr[t]; ks < (M.row_ptr[t + 1] - 1); ks++)
					{
						sum += M.values[ks]*y[M.column_index[ks]];
					}

					T one_over_Mii = 1/M(t, t);
					y[t] = (b[t] - sum)*one_over_Mii; 
				}
			}

			// Domain 6
			number = 0, num_2 = 0;
			for (int j = ruy*i_res_input + llx; j <= ruy*i_res_input + rux - 1; j++)
			{
				T sum(0);
			
				for (int k = M.row_ptr[j]; k < (M.row_ptr[j + 1] - 1); k++)
				{
					sum += M.values[k]*y[M.column_index[k]];
				}
			
				T one_over_Mii = 1/M(j, j);
				y[j] = (b[j] - sum)*one_over_Mii;
			}

			// Domain 7
			number = 0, num_2 = 0;
			for (int j = lly*i_res_input + rux; j <= ruy*i_res_input + rux - i_res_input; j += i_res_input)
			{
				T sum(0);
			
				for (int k = M.row_ptr[j]; k < (M.row_ptr[j + 1] - 1); k++)
				{
					sum += M.values[k]*y[M.column_index[k]];
				}
			    
				T one_over_Mii = 1/M(j, j);
				y[j] = (b[j] - sum)*one_over_Mii;
			}

			// Domain 8
			int ii = ruy*i_res_input + rux;
			T one_over_M_start = 1/M(ii, ii);
		
			T sum(0);
			for (int k = M.row_ptr[ii]; k < (M.row_ptr[ii + 1] - 1); k++)
			{
				sum += M.values[k]*y[M.column_index[k]];
			}
			    
			T one_over_Mii = 1/M(ii, ii);
			y[ii] = (b[ii] - sum)*one_over_Mii;
		}
	}
};