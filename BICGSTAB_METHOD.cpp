#include "BICGSTAB_METHOD.h"

void BICGSTAB_METHOD::Solve(const CSR_MATRIX<T>& A, VECTOR_ND<T>& x, const VECTOR_ND<T>& b, const FIELD_STRUCTURE_2D<int>& bc)
{
	BICGSTABMethod(A, x, b);
}

void BICGSTAB_METHOD::Solve(const CSR_MATRIX<T>& A, VECTOR_ND<T>& x, const VECTOR_ND<T>& b, const FIELD_STRUCTURE_2D<int>& bc, const int& thread_id)
{
	BICGSTABMethod(A, x, b, thread_id);
}

void BICGSTAB_METHOD::BICGSTABMethod(const CSR_MATRIX<T>& A, VECTOR_ND<T>& x, const VECTOR_ND<T>& b)
{
	const int N(x.num_dimension);

	res.Initialize(N);

	rtilde.Initialize(N);

	p.Initialize(N);

	Ap.Initialize(N);

	s.Initialize(N);

	As.Initialize(N);

	num_iteration = 0;

	T *rval(res.values), *rtval(rtilde.values), *pval(p.values), *Apval(Ap.values), *sval(s.values), *xval(x.values);

	T rho_1(1), rho_2(1), alpha(1), beta(1), omega(1);

	A.ComputeResidual(x, b, res);

	T norm_of_b = b.Norm2();
	T norm_of_res = res.Norm2();

	if (norm_of_b == (T)0)
	{
		norm_of_b = 1;
	}

	if (norm_of_res/norm_of_b <= tolerance)
	{
		tolerance = norm_of_res/norm_of_b;
		max_iteration = 0;
		
		return;
	}

	for (int i = 0; i < N; i++)
	{
		rtilde.values[i] = res.values[i];
	}

	while (num_iteration < max_iteration)
	{
		rho_1 = DotProduct(rtilde, res);

		if (rho_1 == 0)
		{
			tolerance = norm_of_res/norm_of_b;
			return;
		}

		beta = (rho_1/rho_2)*(alpha/omega);

		if (num_iteration == 0)
		{
			for (int i = 0; i < N; i++)
			{
				pval[i] = rval[i];
			}
		}
		else
		{
			for (int i = 0; i < N; i++)
			{
				pval[i] = rval[i] + beta*(pval[i] - omega*Ap[i]);
			}
		}

		A.Multiply(p, Ap);

		alpha = rho_1/DotProduct(rtilde, Ap);
		
		for (int i = 0; i < N; i++)
		{
			sval[i] = rval[i] - alpha*Ap[i];
		}

		A.Multiply(s, As);

		omega = DotProduct(As, s)/DotProduct(As, As);

		for (int i = 0; i < N; i++)
		{
			xval[i] = xval[i] + alpha*pval[i] + omega*sval[i];
			rval[i] = sval[i] - omega*As[i];
		}

		rho_2 = rho_1;

		if (res.Norm2()/norm_of_b < tolerance)
		{
			cout << "Converge!!" << endl;
			cout << "Iteration Number : " << num_iteration << endl;
			break;
		}
		
		if (omega == 0)
		{
			tolerance = res.Norm2()/norm_of_b;
			break;
		}	

		num_iteration++;
	}

	residual = res.Norm2();
	cout << "Residual : " << residual << endl;
}

void BICGSTAB_METHOD::BICGSTABMethod(const CSR_MATRIX<T>& A, VECTOR_ND<T>& x, const VECTOR_ND<T>& b, const int& thread_id)
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
	END_HEAD_THREAD_WORK;

	BEGIN_HEAD_THREAD_WORK
	{
		rtilde.Initialize(N);
	}
	END_HEAD_THREAD_WORK;

	BEGIN_HEAD_THREAD_WORK
	{
		p.Initialize(N);
	}
	END_HEAD_THREAD_WORK;

	BEGIN_HEAD_THREAD_WORK
	{
		Ap.Initialize(N);
	}
	END_HEAD_THREAD_WORK;

	BEGIN_HEAD_THREAD_WORK
	{
		s.Initialize(N);
	}
	END_HEAD_THREAD_WORK;

	BEGIN_HEAD_THREAD_WORK
	{
		As.Initialize(N);
	}
	END_HEAD_THREAD_WORK;

	BEGIN_HEAD_THREAD_WORK
	{
		num_iteration = 0;
	}
	END_HEAD_THREAD_WORK;

	T *rval(res.values), *rtval(rtilde.values), *pval(p.values), *Apval(Ap.values), *sval(s.values), *xval(x.values);

	T rho_1(1), rho_2(1), alpha(1), beta(1), omega(1);

	A.ComputeResidual(x, b, res, thread_id);

	T norm_of_b, norm_of_res; 
	
	BEGIN_HEAD_THREAD_WORK
	{
		norm_of_b = b.Norm2();
	}
	END_HEAD_THREAD_WORK;
	
	BEGIN_HEAD_THREAD_WORK
	{
		norm_of_res = res.Norm2();
	}
	END_HEAD_THREAD_WORK;

	if (norm_of_b == (T)0)
	{
		norm_of_b = 1;
	}

	if (norm_of_res/norm_of_b <= tolerance)
	{
		tolerance = norm_of_res/norm_of_b;
		max_iteration = 0;
		
		return;
	}

	for (int i = start_ix; i <= end_ix; i++)
	{
		rtilde.values[i] = res.values[i];
	}
	multithreading->Sync(thread_id);

	while (num_iteration < max_iteration)
	{
		rho_1 = DotProduct(rtilde, res, multithreading, thread_id);

		if (rho_1 == 0)
		{
			tolerance = norm_of_res/norm_of_b;
			return;
		}

		beta = (rho_1/rho_2)*(alpha/omega);

		if (num_iteration == 0)
		{
			for (int i = start_ix; i <= end_ix; i++)
			{
				pval[i] = rval[i];
			}
			multithreading->Sync(thread_id);
		}
		else
		{
			for (int i = start_ix; i < end_ix; i++)
			{
				pval[i] = rval[i] + beta*(pval[i] - omega*Ap[i]);
			}
			multithreading->Sync(thread_id);
		}

		A.Multiply(p, Ap, thread_id);

		alpha = rho_1/DotProduct(rtilde, Ap, multithreading, thread_id);
		
		for (int i = start_ix; i <= end_ix; i++)
		{
			sval[i] = rval[i] - alpha*Ap[i];
		}
		multithreading->Sync(thread_id);

		A.Multiply(s, As, thread_id);

		omega = DotProduct(As, s, multithreading, thread_id)/DotProduct(As, As, multithreading, thread_id);

		for (int i = start_ix; i <= end_ix; i++)
		{
			xval[i] = xval[i] + alpha*pval[i] + omega*sval[i];
			rval[i] = sval[i] - omega*As[i];
		}
		multithreading->Sync(thread_id);

		rho_2 = rho_1;

		if (res.Norm2() < tolerance)
		{
			cout << "--------------BICGSTAB Method--------------" << endl;
			cout << "Converge!!" << endl;
			cout << "Iteration Number : " << num_iteration << endl;
			cout << "--------------------------------------------" << endl;
			break;
		}
		if (num_iteration == max_iteration)
		{
			cout << "--------------BICGSTAB Method--------------" << endl;
			cout << "Not Converge T^T" << endl;
			cout << "Residual: " << residual << endl;
			cout << "--------------------------------------------" << endl;
		}

		if (omega == 0)
		{
			tolerance = res.Norm2()/norm_of_b;
			break;
		}	

		BEGIN_HEAD_THREAD_WORK
		{
			num_iteration++;
		}
		END_HEAD_THREAD_WORK;
	}

	BEGIN_HEAD_THREAD_WORK
	{
		residual = res.Norm2();
		cout << "Residual : " << residual << endl;
	}
	END_HEAD_THREAD_WORK;
}