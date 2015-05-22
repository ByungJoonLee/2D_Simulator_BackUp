#include "CG_METHOD.h"

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
	