#include "GAUSS_SEIDEL_METHOD.h"

void GAUSS_SEIDEL_METHOD::Solve(const CSR_MATRIX<T>& A, VECTOR_ND<T>& x, const VECTOR_ND<T>& b, const FIELD_STRUCTURE_2D<int>& bc, const int& thread_id)
{
	GaussSeidelMethod(A, x, b, thread_id);
}

void GAUSS_SEIDEL_METHOD::GaussSeidelMethod(const CSR_MATRIX<T>& A, VECTOR_ND<T>& x, const VECTOR_ND<T>& b, const int& thread_id)
{
	BEGIN_HEAD_THREAD_WORK
	{
		multithreading->SplitDomainIndex1D(0, A.N);
	}
	END_HEAD_THREAD_WORK

	BEGIN_HEAD_THREAD_WORK
	{
		num_iteration = 0;
	}
	END_HEAD_THREAD_WORK

	residual = (T)0;

	for(int i = 0; i < max_iteration; i ++)	//NOTE: do not use member variable as an iterator when multithreaded.
	{
		residual = GaussSeidelStep(A, x, b, thread_id);
		
		BEGIN_HEAD_THREAD_WORK
		{
			num_iteration ++;
		}
		END_HEAD_THREAD_WORK

		if(residual < tolerance) 
		{
			/*cout << "--------------Gauss-Seidel Iteration--------------" << endl;
			cout << "Converge!!" << endl;
			cout << "Iteration Number : " << num_iteration << endl;
			cout << "Residual: " << residual << endl;
			cout << "--------------------------------------------------" << endl;*/
			break;
		}
		if (i == max_iteration - 1)
		{
			cout << "--------------Gauss-Seidel Iteration--------------" << endl;
			cout << "Not Converge T^T" << endl;
			cout << "Residual: " << residual << endl;
			cout << "--------------------------------------------------" << endl;
		}
	}

	multithreading->Sync(thread_id);
}

T GAUSS_SEIDEL_METHOD::GaussSeidelStep(const CSR_MATRIX<T>& A, VECTOR_ND<T>& x, const VECTOR_ND<T>& div, const int& thread_id) // this_matrix * x -> b
{
	assert(A.N == x.num_dimension);
	assert(x.num_dimension == div.num_dimension);
	
	const int k_start(multithreading->start_ix_1D[thread_id]), k_end(multithreading->end_ix_1D[thread_id]);

	T *divval(div.values), *xval(x.values);
	
	int v_start, v_end, vix;
	T v, A_ii, residual_temp, residual_sum(0);
	
	// For check diagonal dominancy
	T dd_check(0);

	for(int row = k_start; row <= k_end; row++) // multiply 'row'th row of this matrix and vector x to compute b[row]
	{
		v_start = A.row_ptr[row];
		v_end = A.row_ptr[row+1] - 1;

		v = (T)0;
		//A_ii = (A.values[A.row_ptr[row+1] - 1]);	//NOTE: the last value of each row is the diagonal term in our POISSON_SOLVER_UNIFORM::BuilidLinearSystem function.
		A_ii = A(row, row);
		
		/*dd_check = 0;

		for (vix = v_start; vix <= v_end; vix++)
		{
			if (A.column_index[vix] != row)
			{
				dd_check += abs(A.values[vix]);
			}
		}
						
		if (abs(A_ii) < dd_check)
		{
			cout << "A_ii = " << A_ii << endl;
			cout << "Sum of others = " << dd_check << endl;
			cout << "Error: Strict Diagonally Dominancy fails" << endl;
			exit(-1);
		}*/

		for(vix = v_start; vix <= v_end; vix ++) // iterate all components of 'row'th row of this matrix
		{
			v += A.values[vix]*xval[A.column_index[vix]];
		}

		residual_temp = divval[row] - v;

		if(A_ii != 0) xval[row] += residual_temp/A_ii;
		// this cell is surrounded by Neumann cells if A_ii = 0.

		residual_sum += POW2(residual_temp);
	}

	multithreading->SyncSum(thread_id, residual_sum);

	residual_sum = sqrt(residual_sum);

	return residual_sum;
}
