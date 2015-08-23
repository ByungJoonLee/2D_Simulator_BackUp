#include "CG_METHOD.h"

void CG_METHOD::Solve(const CSR_MATRIX<T>& A, VECTOR_ND<T>& x, const VECTOR_ND<T>& b, const FIELD_STRUCTURE_2D<int>& bc)
{
	CGMethod(A, x, b);
}

void CG_METHOD::Solve(const CSR_MATRIX<T>& A, VECTOR_ND<T>& x, const VECTOR_ND<T>& b, const FIELD_STRUCTURE_2D<int>& bc, const int& thread_id)
{
	CGMethod(A, x, b, thread_id);
}

void CG_METHOD::CGMethod(const CSR_MATRIX<T>& A, VECTOR_ND<T>& x, const VECTOR_ND<T>& b)
{
	const int N(x.num_dimension);

	res.Initialize(N);
	
	p.Initialize(N);
	
	Ap.Initialize(N);
	
	num_iteration = 0;

	T *rval(res.values), *pval(p.values), *Apval(Ap.values), *xval(x.values);

	T alpha, res_old, res_new;

	A.ComputeResidual(x, b, res);

	for (int i = 0; i < N; i++)
	{
		p.values[i] = res.values[i];	
	}

	DotProduct(res, p, res_old);

	num_iteration = 0;
	
	while (num_iteration < max_iteration)
	{
		A.Multiply(p, Ap);
		
		alpha = res_old/DotProduct(p, Ap);

		for (int i = 0; i < N; i++)
		{
			xval[i] += alpha*pval[i];
			rval[i] -= alpha*Apval[i];
		}
				
		res_new = DotProduct(res, res);
		
		if (res_new < sqr_tolerance)
		{
			cout << "Converge!!" << endl;
			cout << "Iteration Number : " << num_iteration << endl;
			break;
		}

		for (int i = 0; i < N; i++)
		{
			const T k = res_new/res_old;

			pval[i] = res.values[i] + k*p.values[i];
		}

		res_old = res_new;
		
		num_iteration++;
	}

	residual = sqrt(res_new);
	cout << "Residual : " << residual << endl;
}

void CG_METHOD::CGMethod(const CSR_MATRIX<T>& A, VECTOR_ND<T>& x, const VECTOR_ND<T>& b, const int& thread_id)
{
	BEGIN_HEAD_THREAD_WORK
	{
		multithreading->SplitDomainIndex1D(0, A.N);
	}
	END_HEAD_THREAD_WORK

	const int N(x.num_dimension), start_ix(multithreading->start_ix_1D[thread_id]), end_ix(multithreading->end_ix_1D[thread_id]);

	BEGIN_HEAD_THREAD_WORK
	{
		res.Initialize(N);
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

	T *rval(res.values), *pval(p.values), *Apval(Ap.values), *xval(x.values);

	T alpha, res_old, res_new;

	A.ComputeResidual(x, b, res, thread_id);

	for (int i = start_ix; i <= end_ix; i++)
	{
		p.values[i] = res.values[i];
	}
	multithreading->Sync(thread_id);

	res_old = DotProduct(res, p, multithreading, thread_id);
	
	while(num_iteration < max_iteration)
	{
		A.Multiply(p, Ap, thread_id);

		alpha = res_old/ DotProduct(p, Ap, multithreading, thread_id);
		
		for (int i = start_ix; i <= end_ix; i++)
		{
			xval[i] += alpha*pval[i];
			rval[i] -= alpha*Apval[i];
		}
		multithreading->Sync(thread_id);

		res_new = DotProduct(res, res, multithreading, thread_id);

		if(res_new < sqr_tolerance) break;					// In L2 Norm

		for (int i = start_ix; i <= end_ix; i++)
		{
			const T k = res_new/res_old;

			pval[i] = res.values[i] + k*p.values[i];
		}
		multithreading->Sync(thread_id);

		res_old = res_new;

        BEGIN_HEAD_THREAD_WORK
		{
			num_iteration++;
		}
		END_HEAD_THREAD_WORK
	}

	multithreading->Sync(thread_id);

	BEGIN_HEAD_THREAD_WORK
	{
		residual = sqrt(res_new);
		cout << "Residual : " << residual << endl;
	}
	END_HEAD_THREAD_WORK;
}
