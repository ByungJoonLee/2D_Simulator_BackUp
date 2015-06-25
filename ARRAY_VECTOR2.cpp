#include "ARRAY_VECTOR2.h"
#include "math.h"

//////////////////////////////////////////////////////////////////////
// Basic 3 float Vector Math Functions
//////////////////////////////////////////////////////////////////////
namespace ARRAY_VECTOR2
{
	template<typename T>
	void convex(T *a, T *b, T *c, T s)
	{
		c[0] = a[0]*s + b[0]*(1.0f-s);
		c[1] = a[1]*s + b[1]*(1.0f-s);
	}

	template<typename T>
	void normalize(T *n)
	{
		T det = sqrt(n[0]*n[0] + n[1]*n[1]);
		if(det == 0.0f)	
			det = 1.0f;

		ARRAY_VECTOR2::div<T>(n, det);
	}

	template<typename T>
	void det(T *a, T* det)
	{
		*det = sqrt(a[0]*a[0] + a[1]*a[1]);
	}

	template<typename T>
	T det(T *a)
	{
		return sqrt(a[0]*a[0] + a[1]*a[1]);
	}

	template<typename T>
	void add(T *a, T *b, T *c)			// a + b = c
	{
		c[0] = a[0] + b[0];			
		c[1] = a[1] + b[1];			
	}

	template<typename T>
	void sub(T *a, T *b, T *c)			// a - b = c
	{
		c[0] = a[0] - b[0];			
		c[1] = a[1] - b[1];			
	}

	template<typename T>
	void set(T *a, T b)						// set a as b (scalar)
	{
		for(int i = 0; i < 2; i ++)
		{
			a[i] = b;
		}		
	}

	template<typename T>
	void set(T *a, T *b)					// set a as b (vector)
	{
		for(int i = 0; i < 2; i ++)
		{
			a[i] = b[i];		
		}
	}

	void set(int *a, int *b)
	{
		for(int i = 0; i < 2; i ++)
		{
			a[i] = b[i];
		}
	}

	void set(int *a, int b)
	{
		for(int i = 0; i < 2; i ++)
		{
			a[i] = b;
		}
	}

	void set(int *a, int i, int j)
	{
		a[0] = i;
		a[1] = j;
	}

	template<typename T>
	void div(T *a, T b)						// div a by b
	{
		for(int i = 0; i < 2; i ++)
		{
			a[i] /= b;
		}
	}

	template<typename T>
	int maxIndex(T *a)
	{
		int ans = 0;
		for(int i = 0; i < 2; i ++)
			if(a[i] > a[ans])
				ans = i;		
		return ans;
	}


	template<typename T>
	int minIndex(T *a)
	{
		int ans = 0;
		for(int i = 0; i < 2; i ++)
			if(a[i] < a[ans])
				ans = i;		
		return ans;
	}

	template<typename T>
	void mul(T *a, T b)
	{
		for(int i = 0; i < 2; i ++)
		{
			a[i] *= b;
		}
	}

	template<typename T>
	void cross(T *a, T *b, T *c)		// a cross b = c
	{	
		T temp[3];
		temp[0] = 0;
		temp[1] = 0;
		temp[2] = a[0]*b[1] - a[1]*b[0];		
		c[0] = temp[0];
		c[1] = temp[1];
		c[2] = temp[2];
	}

	template<typename T>
	T dot(T *a, T *b)
	{
		return a[0]*b[0] + a[1]*b[1];
	}

	template<typename T>
	T dist(T *a, T *b)
	{
		T diff[2] = {a[0]-b[0], a[1]-b[1]};
		return sqrt(diff[0]*diff[0] + diff[1]*diff[1]);
	}

	template<typename T>
	void weightedAverage(T *v0, T *v1, T *result, T w0, T w1)
	{
		T det = w0 + w1;
		det = 1.0f/det;
		result[0] = w1*v0[0] + w0*v1[0];
		result[1] = w1*v0[1] + w0*v1[1];
		
		result[0] *= det;
		result[1] *= det;
	}
}

#ifdef USE_FLOAT_T		
template	void	ARRAY_VECTOR2::add(float *a, float *b, float *c);
template	void	ARRAY_VECTOR2::convex(float *a, float *b, float *c, float s);
template	void	ARRAY_VECTOR2::cross(float* a, float* b, float* c);
template	void	ARRAY_VECTOR2::det(float *a, float *det);
template	float	ARRAY_VECTOR2::det(float *a);
template	float	ARRAY_VECTOR2::dist(float *a, float *b);
template	void	ARRAY_VECTOR2::div(float *a, float b);
template	float	ARRAY_VECTOR2::dot(float *a, float *b);

template	void	ARRAY_VECTOR2::set(float *a, float b);
template	void	ARRAY_VECTOR2::set(float *a, float *b);

template	void	ARRAY_VECTOR2::sub(float *a, float *b, float *c);	
template	int		ARRAY_VECTOR2::maxIndex(float *a);
template	int		ARRAY_VECTOR2::minIndex(float *a);
template	void	ARRAY_VECTOR2::mul(float *a, float b);
template	void	ARRAY_VECTOR2::normalize(float *n);	
template	void	ARRAY_VECTOR2::weightedAverage(float *v0, float *v1, float *result, float w0, float w1);
#else
template	void	ARRAY_VECTOR2::add(double *a, double *b, double *c);
template	void	ARRAY_VECTOR2::convex(double *a, double *b, double *c, double s);
template	void	ARRAY_VECTOR2::cross(double* a, double* b, double* c);
template	void	ARRAY_VECTOR2::det(double *a, double *det);
template	double	ARRAY_VECTOR2::det(double *a);
template	double	ARRAY_VECTOR2::dist(double *a, double *b);
template	void	ARRAY_VECTOR2::div(double *a, double b);
template	double	ARRAY_VECTOR2::dot(double *a, double *b);

template	void	ARRAY_VECTOR2::set(double *a, double b);
template	void	ARRAY_VECTOR2::set(double *a, double *b);

template	void	ARRAY_VECTOR2::sub(double *a, double *b, double *c);	
template	int		ARRAY_VECTOR2::maxIndex(double *a);
template	int		ARRAY_VECTOR2::minIndex(double *a);
template	void	ARRAY_VECTOR2::mul(double *a, double b);
template	void	ARRAY_VECTOR2::normalize(double *n);	
template	void	ARRAY_VECTOR2::weightedAverage(double *v0, double *v1, double *result, double w0, double w1);
#endif


