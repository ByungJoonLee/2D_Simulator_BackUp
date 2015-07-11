#pragma once

#include "COMMON_DEFINITION.h"

class MATRIX_2X2
{
public: 
	T			x[4];

public: // Costructor and Destructor
	MATRIX_2X2(void)
	{}

	MATRIX_2X2(const T x11,const T x21,const T x12,const T x22)
    {
        x[0]=x11;
		x[1]=x21;
		x[2]=x12;
		x[3]=x22;
    }

	~MATRIX_2X2(void)
	{}

public:
	T& operator()(const int& i,const int& j)
    {
		assert(i >= 1 && i <= 2);
		assert(j >= 1 && j <= 2);
		return x[i-1+2*(j-1)]; // TODO : check
	}

	MATRIX_2X2 operator*(const MATRIX_2X2& A) const
    {
		return MATRIX_2X2(x[0]*A.x[0]+x[2]*A.x[1], x[1]*A.x[0]+x[3]*A.x[1], x[0]*A.x[2]+x[2]*A.x[3], x[1]*A.x[2]+x[3]*A.x[3]);
	}

	MATRIX_2X2 Transposed() const
    {
		return MATRIX_2X2(x[0], x[2], x[1], x[3]);
	}

	MATRIX_2X2 Inversed()
    {
		T determinant = x[0]*x[3] - x[1]*x[2];
		
		assert(determinant!=0);
		T s = 1/determinant;

		return MATRIX_2X2(x[3]*s, -x[1]*s, -x[2]*s, x[0]*s);
	}

   MATRIX_2X2 operator * (const T a) const
   {
		return MATRIX_2X2(a*x[0],a*x[1],a*x[2],a*x[3]);
   }

   static MATRIX_2X2 Identity()
   {
		return MATRIX_2X2(1, 0, 0, 1);
   }

   void operator += (const MATRIX_2X2& A)
   {
	   for(int i=0;i<4;i++) x[i] += A.x[i];
   }

   void operator *= (const T& a)
   {
	   for(int i=0;i<4;i++) x[i] *= a;
   }

   void operator = (const MATRIX_2X2& A)
   {
	   for(int i=0;i<4;i++) x[i] = A.x[i];
   }

   MATRIX_2X2 operator + (const MATRIX_2X2& A)
   {
		return MATRIX_2X2(x[0] + A.x[0], x[1] + A.x[1], x[2] + A.x[2], x[3] + A.x[3]);
   }

   MATRIX_2X2 operator - (const MATRIX_2X2& A)
   {
		return MATRIX_2X2(x[0] - A.x[0], x[1] - A.x[1], x[2] - A.x[2], x[3] - A.x[3]);
   }
};

inline std::ostream& operator<<(std::ostream& output, const MATRIX_2X2& A)
{
	output << "|" << A.x[0] << " " << A.x[2] << "|" << endl;
	output << "|" << A.x[1] << " " << A.x[3] << "|" << endl;

	return output;
}

inline static VT operator * (const MATRIX_2X2& A, const VT& v)
{
	return VT(v.x*A.x[0]+v.y*A.x[2], v.x*A.x[1]+v.y*A.x[3]);
}

inline static MATRIX_2X2 operator * (const T s, const MATRIX_2X2& A)
{
	return MATRIX_2X2(s*A.x[0], s*A.x[1], s*A.x[2], s*A.x[3]);
}
