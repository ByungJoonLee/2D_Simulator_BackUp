#pragma once

#include "COMMON_DEFINITION.h"
#include "VECTOR_ND.h"

template<class TT>
class CSR_MATRIX
{
public: // Essential Data
	int						N;					// the number of row
	int						nz;					// the number of nonzero elements
	TT*						values;				// the values of the nonzero elements
	int*					row_ptr;			// the locations in the val vector that start a row
	int*					column_index;		// the column indices of the elements in the val vector

	int						value_ix;
	int						prev_row;

public:	// Constructors and Destructor
	CSR_MATRIX(void)
		: values(0), column_index(0), row_ptr(0)
	{}

	CSR_MATRIX(const int& N_input, const int& nz_input)
		: values(0), column_index(0), row_ptr(0)
	{
		Initialize(N_input, nz_input);
	}

	CSR_MATRIX(const CSR_MATRIX<T>& matrix_input)
		: values(0), column_index(0), row_ptr(0)
	{
		Initialize(matrix_input);
	}

	~CSR_MATRIX(void)
	{
		DeleteMemory();
	}

	void DeleteMemory(void)
	{
		if (values != 0)
		{
			delete [] values;
		}

		if (row_ptr != 0)
		{
			delete [] row_ptr;
		}

		if (column_index != 0)
		{
			delete [] column_index;
		}
	}

public: // Initialization Functions
	void Initialize(const int& N_input, const int& nz_input)
	{
		DeleteMemory();

		N = N_input;
		nz = nz_input;

		values = new TT [nz];
		row_ptr = new int [N+1];
		column_index = new int [nz];

		value_ix = 0;
		prev_row = -1;

		row_ptr[N_input] = nz_input;
	}
	
	void Initialize(const CSR_MATRIX<TT>& matrix_input)
	{
		Initialize(matrix_input.N, matrix_input.nz);
	
		for (int i = 0; i < nz; i++)
		{
			values[i] = matrix_input.values[i];
		}

		for (int i = 0; i < N + 1; i++)
		{
			row_ptr[i] = matrix_input.row_ptr[i];
		}

		for (int i = 0; i < nz; i++)
		{
			column_index[i] = matrix_input.column_index[i];
		}
		
		value_ix = matrix_input.value_ix;
		prev_row = matrix_input.prev_row;
	}

public: // Operator Overloading
	void operator *= (const TT& s)
	{
		for (int i = 0; i < nz; i++)
		{
			values[i] *= s;
		}
	}
	
	inline TT& operator()(const int& row_input, const int& column_input) const
	{
		bool is_nonzero(false);
		for (int vix = row_ptr[row_input]; vix < row_ptr[row_input + 1]; vix++)
		{
			if (column_index[vix] == column_input)
			{
				return values[vix];
				is_nonzero = true;
			}
		}

		if (is_nonzero == false)
		{
			TT zero(0);
			TT& temp = zero;
			return temp;
		}
	}

	CSR_MATRIX<TT>& operator=(const CSR_MATRIX<TT>& matrix_input)
	{
		Initialize(matrix_input);

		return (*this);
	}

public: // Member Functions
	void AssignValue(const int& row_input, const int& column_input, const TT& values_input)
	{
		values[value_ix] = values_input;

		if (row_input != prev_row)
		{
			row_ptr[row_input] = value_ix;
			prev_row = row_input;
		}

		column_index[value_ix] = column_input;

		value_ix++;
	}

	inline void Multiply(const VECTOR_ND<TT>& x, VECTOR_ND<TT>& b) const
	{
		assert(N == x.num_dimension);
		assert(x.num_dimension == b.num_dimension);

		T *bval(b.values), *xval(x.values);

		for (int row = 0; row < N; row++)
		{
			T v = 0;
			for (int vix = row_ptr[row]; vix < row_ptr[row+1]; vix++)
			{
				v += values[vix]*xval[column_index[vix]];
			}

			bval[row] = v;
         }
	}

	inline void ComputeResidual(const VECTOR_ND<TT>& x, const VECTOR_ND<TT>& b, VECTOR_ND<TT>& residual) const
	{
		assert(N == x.num_dimension);
		assert(x.num_dimension == b.num_dimension);
		assert(residual.num_dimension == N);

		TT *bval(b.values), *xval(x.values), *rval(residual.values);

		for (int row = 0; row < N; row++)
		{
			// Compute A*x
			TT v = 0;
			for (int vix = row_ptr[row]; vix < row_ptr[row+1]; vix++)
			{
				v += values[vix]*xval[column_index[vix]];
			}

			// residual = b - A*x
			rval[row] = bval[row] - v;
		}
	}

	TT* GetValue(const int& row, const int& column)
	{
		for (int vix = row_ptr[row]; vix < row_ptr[row+1]; vix++)
		{
			if (column_index[vix] == column)
			{
				return& values[vix];
			}
		}
	}
};

inline std::ostream& operator<<(std::ostream& output, const CSR_MATRIX<T>& A)
{
	output << "--Matrix Information in CSR Form--" << endl;
	
	output << "nz = [ ";
	for (int i = 0; i < A.nz; i++)
	{
		output << A.values[i] << " ";
	}
	output << "]" << endl;
	
	output << "ci = [ ";
	for (int i = 0; i < A.nz; i++)
	{
		output << A.column_index[i] << " ";
	}
	output << "]" << endl;

	output << "rp = [ ";
	for (int i = 0; i <= A.N; i++)
	{
		output << A.row_ptr[i] << " ";
	}
	output << "]" << endl;

	return output;
}

static CSR_MATRIX<T> CholeskyDecomposition(const CSR_MATRIX<T>& A)
{
	const int N = A.N;
	const int nz = A.nz;

	const int num_of_nz = N*(N + 1)/2;
	
	CSR_MATRIX<T> L;
	L.Initialize(N, num_of_nz);
		
	// Set up the column indices and row ptr since the given decomposition consist of lower triangular matrix
	for (int i = 0; i < N + 1; i++)
	{
		L.row_ptr[i] = i*(i + 1)/2;
	}

	for (int i = 0; i < N; i++)
	{
		for (int vix = L.row_ptr[i]; vix < L.row_ptr[i + 1]; vix++)
		{
			L.column_index[vix] = vix - L.row_ptr[i];
		}
	}

	for (int i = 0; i < N; i++)
	{
		for (int j = 0; j <= i; j++)
		{
			if (j == i)
			{
				T sum(0);
				if (j == 0)
				{
					L(i, j) = sqrt(A(i, j));
				}
				else
				{
					for (int k = 0; k <= j - 1; k++)
					{
						sum += POW2(L(i, k));
					}
					L(i, j) = sqrt(A(i, j) - sum);
				}
			}
			else
			{
				if (j == 0)
				{
					L(i, j) = (T)1/L(j, j)*(A(i, j));
				}
				else
				{
					T sum(0);
					for (int k = 0; k <= j - 1; k++)
					{
						sum += L(i, k)*(L(j, k));
					}
					L(i, j) = (T)1/L(j, j)*(A(i, j) - sum); 
				}
			}
		}
	}

	return L;
}

//template<class TT>
//static CSR_MATRIX<TT> IncompleteCholeskyDecomposition(const int& i_res_input, const int& j_res_input, const CSR_MATRIX<TT>& A)
//{
//	const int N = A.N;
//	const int nz = A.nz;
//
//	const int num_of_nz = (nz - N)/2 + N;
//	
//	FIELD_STRUCTURE_2D<int> index_field;
//	index_field.Initialize(i_res_input, j_res_input, 0, 0, 0, 0, 1, 1, 3);
//
//	const int i_start(index_field.i_start), i_end(index_field.i_end), j_start(index_field.j_start), j_end(index_field.j_end);
//	
//	int start_ix(0);
//	
//	int i, j;
//	LOOPS_2D(i, j, i_start, j_start, i_end, j_end)
//	{
//		index_field(i, j) = -1;
//	}
//
//	LOOPS_2D(i, j, index_field.i_start_g, index_field.j_start_g, index_field.i_end_g, index_field.j_end_g)
//	{
//		if (i < i_start || i > i_end || j < j_start || j > j_end)
//		{
//			index_field(i, j) = -1;
//		}
//	}
//
//	LOOPS_2D(i, j, i_start, j_start, i_end, j_end)
//	{
//		index_field(i, j) = start_ix++;
//	}
//	
//	CSR_MATRIX<TT> L;
//	L.Initialize(N, num_of_nz);
//
//	int number(0);
//	LOOPS_2D(i, j, i_start, j_start, i_end, j_end)
//	{
//		T sum(0), coef(0);
//		
//		if (index_field(i, j - 1) > -1)
//		{
//			if (index_field(i, j) == index_field.grid.i_res)
//			{
//				coef = (T)1/L.values[L.row_ptr[index_field(i, j) - index_field.grid.i_res]]*(A(index_field(i, j - 1), index_field(i, j)));
//			}
//			else if (((index_field(i, j) > index_field.grid.i_res) && (index_field(i, j) < 2*index_field.grid.i_res)) || (index_field(i, j) % index_field.grid.i_res == 0))
//			{
//				coef = (T)1/L.values[L.row_ptr[index_field(i, j) - index_field.grid.i_res] + 1]*(A(index_field(i, j - 1), index_field(i, j)));
//			}
//			else 
//			{
//				coef = (T)1/L.values[L.row_ptr[index_field(i, j) - index_field.grid.i_res] + 2]*(A(index_field(i, j - 1), index_field(i, j)));
//			}
//			
//			L.AssignValue(index_field(i, j), index_field(i, j - 1), coef);
//			number += 1;
//		}
//		
//		if (index_field(i - 1, j) > -1)
//		{
//			if (index_field(i, j) == 1)
//			{
//				coef = (T)1/L.values[L.row_ptr[index_field(i, j) - 1]]*(A(index_field(i - 1, j), index_field(i, j)));
//			}
//			else if (index_field(i, j) < index_field.grid.i_res)
//			{
//				coef = (T)1/L.values[L.row_ptr[index_field(i, j) - 1] + 1]*(A(index_field(i - 1, j), index_field(i, j)));
//			}
//			else
//			{
//				coef = (T)1/L.values[L.row_ptr[index_field(i, j)] - 1]*(A(index_field(i - 1, j), index_field(i, j)));
//			}
//			
//			L.AssignValue(index_field(i, j), index_field(i - 1, j), coef);
//			number += 1;
//		}
//		
//		if (index_field(i, j) == 0)
//		{
//			coef = sqrt(A(index_field(i, j), index_field(i, j)));
//		}
//		else
//		{
//			sum = 0;
//			
//			for (int k = L.row_ptr[index_field(i, j)]; k < number; k++)
//			{
//				sum += POW2(L.values[k]);
//			}
//			coef = sqrt(A(index_field(i, j), index_field(i, j)) - sum);
//		}
//		L.AssignValue(index_field(i, j), index_field(i, j), coef);
//			
//		number += 1;
//	}	
//		
//	return L;
//}

template<class TT>
static void IncompleteCholeskyDecomposition(const int& i_res_input, const int& j_res_input, const CSR_MATRIX<TT>& A, CSR_MATRIX<TT>& L)
{
	const int N = A.N;
	const int nz = A.nz;

	FIELD_STRUCTURE_2D<int> index_field;
	index_field.Initialize(i_res_input, j_res_input, 0, 0, 0, 0, 1, 1, 3);

	const int i_start(index_field.i_start), i_end(index_field.i_end), j_start(index_field.j_start), j_end(index_field.j_end);
	
	int start_ix(0);
	
	int i, j;
	LOOPS_2D(i, j, i_start, j_start, i_end, j_end)
	{
		index_field(i, j) = -1;
	}

	LOOPS_2D(i, j, index_field.i_start_g, index_field.j_start_g, index_field.i_end_g, index_field.j_end_g)
	{
		if (i < i_start || i > i_end || j < j_start || j > j_end)
		{
			index_field(i, j) = -1;
		}
	}

	LOOPS_2D(i, j, i_start, j_start, i_end, j_end)
	{
		index_field(i, j) = start_ix++;
	}
	
	int number(0);
	LOOPS_2D(i, j, i_start, j_start, i_end, j_end)
	{
		TT sum((TT)0), coef((TT)0);
		
		if (index_field(i, j - 1) > -1)
		{
			if (index_field(i, j) == index_field.grid.i_res)
			{
				coef = 1/L.values[L.row_ptr[index_field(i, j) - index_field.grid.i_res]]*(A(index_field(i, j - 1), index_field(i, j)));
			}
			else if (((index_field(i, j) > index_field.grid.i_res) && (index_field(i, j) < 2*index_field.grid.i_res)) || (index_field(i, j) % index_field.grid.i_res == 0))
			{
				coef = 1/L.values[L.row_ptr[index_field(i, j) - index_field.grid.i_res] + 1]*(A(index_field(i, j - 1), index_field(i, j)));
			}
			else 
			{
				coef = 1/L.values[L.row_ptr[index_field(i, j) - index_field.grid.i_res] + 2]*(A(index_field(i, j - 1), index_field(i, j)));
			}
			
			L.AssignValue(index_field(i, j), index_field(i, j - 1), coef);
			number += 1;
		}
		
		if (index_field(i - 1, j) > -1)
		{
			if (index_field(i, j) == 1)
			{
				coef = 1/L.values[L.row_ptr[index_field(i, j) - 1]]*(A(index_field(i - 1, j), index_field(i, j)));
			}
			else if (index_field(i, j) < index_field.grid.i_res)
			{
				coef = 1/L.values[L.row_ptr[index_field(i, j) - 1] + 1]*(A(index_field(i - 1, j), index_field(i, j)));
			}
			else
			{
				coef = 1/L.values[L.row_ptr[index_field(i, j)] - 1]*(A(index_field(i - 1, j), index_field(i, j)));
			}
			
			L.AssignValue(index_field(i, j), index_field(i - 1, j), coef);
			number += 1;
		}
		
		if (index_field(i, j) == 0)
		{
			coef = sqrt(A(index_field(i, j), index_field(i, j)));
		}
		else
		{
			sum = (TT)0;
			
			for (int k = L.row_ptr[index_field(i, j)]; k < number; k++)
			{
				sum += POW2(L.values[k]);
			}
			coef = sqrt(A(index_field(i, j), index_field(i, j)) - sum);
		}
		L.AssignValue(index_field(i, j), index_field(i, j), coef);
			
		number += 1;
	}	
}

template<class TT>
static CSR_MATRIX<TT> DiagonalPreconditioner(const CSR_MATRIX<TT>& A)
{
	const int N = A.N;
	
	CSR_MATRIX<TT> D;
	D.Initialize(N, N);

	for (int i = 0; i < N + 1; i++)
	{
		D.row_ptr[i] = i;
	}

	for (int i = 0; i < N; i++)
	{
		for (int vix = D.row_ptr[i]; vix < D.row_ptr[i + 1]; vix++)
		{
			D.column_index[vix] = vix;
		}
	}

	for (int i = 0; i < N; i++)
	{
		D(i, i) = A(i, i);
	}

	return D;
}