#pragma once

namespace ARRAY_VECTOR2
{
	template<typename T>	void	add(T *a, T *b, T *c);
	template<typename T>	void	convex(T *a, T *b, T *c, T s);			// c = a*s + b*(1-s)
	template<typename T>	void	cross(T* a, T* b, T* c);
	template<typename T>	void	det(T *a, T *det);
	template<typename T>	T		det(T *a);
	template<typename T>	T		dist(T *a, T *b);		// Distance between two points
	template<typename T>	void	div(T *a, T b);
	template<typename T>	T		dot(T *a, T *b);
	template<typename T>	void	set(T *a, T b);
	template<typename T>	void	set(T *a, T *b);
							void	set(int   *a, int   *b);
							void	set(int   *a, int    b);
							void	set(int	  *a, int i, int j);
	template<typename T>	void	sub(T *a, T *b, T *c);	
	template<typename T>	int		maxIndex(T *a);
	template<typename T>	int		minIndex(T *a);
	template<typename T>	void	mul(T *a, T b);
	template<typename T>	void	normalize(T *n);	
	template<typename T>	void	weightedAverage(T *v0, T *v1, T *result, T w0, T w1);
}
