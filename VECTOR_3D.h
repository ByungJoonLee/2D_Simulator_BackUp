#pragma once

#include <iostream>

template<class TT>
class VECTOR_3D
{
public: // Essential member variables
	union
	{
		struct{TT x, y, z;};
		struct{TT i, j, k;};
		TT values[3];
	};	// Zero based indexing

public: // Constructors and Destructor
	inline VECTOR_3D(void)
		: x(TT()), y(TT()), z(TT())
	{}

	inline VECTOR_3D(const TT& x_input, const TT& y_input, const TT& z_input)
		: x(x_input), y(y_input), z(z_input)
	{}

	inline VECTOR_3D(const VECTOR_3D<TT>& vector_3d_input)
		: x(vector_3d_input.x), y(vector_3d_input.y), z(vector_3d_input.z)
	{}

	inline VECTOR_3D(const TT values_input[3])
		: x(values_input[0]), y(values_input[1]), z(values_input[2])
	{}

	~VECTOR_3D(void)
	{}
	
public: // Initialization Function
	void Initialize(const TT& x_input, const TT& y_input, const TT& z_input)
	{
		x = x_input;
		y = y_input;
		z = z_input;
	}

	void Initialize(const VECTOR_3D<TT>& vector_3d_input)
	{
		x = vector_3d_input.x;
		y = vector_3d_input.y;
		z = vector_3d_input.z;
	}

	void Initialize(const TT values_input[3])
	{
		x = values_input[0];
		y = values_input[1];
		z = values_input[2];
	}

public: // Operation overloading
	inline TT operator [] (const int& i)
	{
		return values[i];
	}

	inline void operator = (const VECTOR_3D& v)
	{
		x = v.x;
		y = v.y;
		z = v.z;
	}

	inline void operator += (const VECTOR_3D& v)
	{
		x += v.x;
		y += v.y;
		z += v.z;
	}

	inline void operator -= (const VECTOR_3D& v)
	{
		x -= v.x;
		y -= v.y;
		z -= v.z;
	}

	inline void operator *= (const TT& s)
	{
		x *= s;
		y *= s;
		z *= s;
	}

	inline void operator /= (const TT& s)
	{
		TT one_over_s = (TT)1/s;
		x *= one_over_s;
		y *= one_over_s;
		z *= one_over_s;
	}

	inline VECTOR_3D operator + (const VECTOR_3D& v) const
	{
		return VECTOR_3D(x + v.x, y + v.y, z + v.z);
	}

	inline VECTOR_3D operator + (const TT& a) const
	{
		return VECTOR_3D(x + a, y + a, z + a);
	}

	inline VECTOR_3D operator - (const VECTOR_3D& v) const
	{
		return VECTOR_3D(x - v.x, y - v.y, z - v.z);
	}

	inline VECTOR_3D operator - () const
	{
		return VECTOR_3D(-x, -y, -z);
	}

	inline VECTOR_3D operator * (const TT& s) const
	{
		return VECTOR_3D(x*s, y*s, z*s);
	}

	inline VECTOR_3D operator * (const VECTOR_3D<TT>& v_input)
	{
		return VECTOR_3D(x*v_input.x, y*v_input.y, z*v_input.z);
	}

	inline VECTOR_3D operator / (const TT& s) const
	{
		TT one_over_s = (TT)1/s;
	
		return VECTOR_3D(x*one_over_s, y*one_over_s, z*one_over_s);
	}

	inline VECTOR_3D operator / (const VECTOR_3D<TT>& v) const
	{
		return VECTOR_3D(x/v.x, y/v.y, z/v.z);
	}

public: // Member Functions
	inline TT Magnitude() const
	{
		return sqrt(x*x + y*y + z*z);
	}

	inline TT SqrMagnitude() const
	{
		return x*x + y*y + z*z;
	}

	inline void Normalize()
	{
		TT M = Magnitude();
		
		if(M != 0)
		{
			TT one_over_M = (TT)1/M;
			
			x *= one_over_M;
			y *= one_over_M;
			z *= one_over_M;
		}
	}

	inline VECTOR_3D<TT> Normalized() const
	{
		VECTOR_3D<TT> normalized_vector(x, y, z);
		
		normalized_vector.Normalize();
		
		return normalized_vector;
	}

	inline bool IsSqrMagnitudeSmallerThan(const TT& sqrmagnitude) const
	{
		if(sqrmagnitude > (x*x + y*y)) 
			return true;
		else
			return false;
	}

	// Be caution let normal vector be a unit vector
	inline void ScalingComponents(const VECTOR_3D<TT>& normal, const TT& normal_coef, const TT& tangential_coef)
	{
		const TT alpha = DotProduct(*this, normal);
		(*this) -= alpha*normal;
		(*this) = (normal_coef*alpha)*normal + tangential_coef*(*this);
	}
	
	inline void Assign(const TT& x_input, const TT& y_input, const TT& z_input)
	{
		x = x_input;
		y = y_input;
		z = z_input;
	}

	inline void Assign(const VECTOR_3D& v_input)
	{
		x = v_input.x;
		y = v_input.y;
		z = v_input.z;
	}

	inline void Assign(const TT values_input[3])
	{
		x = values_input[0];
		y = values_input[1];
		z = values_input[2];
	}

	inline void AssignZeroVector()
	{
		x = (TT)0;
		y = (TT)0;
		z = (TT)0;
	}

	inline void AssignDifference(const VECTOR_3D<TT>& v1, const VECTOR_3D<TT>& v2)
	{
		x = v1.x - v2.x;
		y = v1.y - v2.y;
		z = v1.z - v2.z;
	}

	inline void AssignDifferencePlusScaledDifference(const VECTOR_3D<TT>& p1, const VECTOR_3D<TT>& p2, const VECTOR_3D<TT>& v1, const VECTOR_3D<TT>& v2, const TT& dt)
	{
		x = (p1.x - p2.x) + (v1.x - v2.x)*dt;
		y = (p1.y - p2.y) + (v1.y - v2.y)*dt;
		z = (p1.z - p2.z) + (v1.z - v2.z)*dt;
	}

	inline void AssignScaledDifference(const TT& scalar, const VECTOR_3D<TT>& v1, const VECTOR_3D<TT>& v2)
	{
		x = scalar*(v1.x - v2.x);
		y = scalar*(v1.y - v2.y);
		z = scalar*(v1.z - v2.z);
	}

	inline void AddScaledDifference(const TT& scalar, const VECTOR_3D<TT>& v1, const VECTOR_3D<TT>& v2)
	{
		x += scalar*(v1.x - v2.x);
		y += scalar*(v1.y - v2.y);
		z += scalar*(v1.z - v2.z);
	}

	inline void AssignScaledVector(const TT& scalar, const VECTOR_3D<TT>& v1)
	{
		x = scalar*v1.x;
		y = scalar*v1.y;
		z = scalar*v1.z;
	}

	inline void AddSum(const VECTOR_3D<TT>& v1, const VECTOR_3D<TT>& v2)
	{
		x += (v1.x + v2.x);
		y += (v1.y + v2.y);
		z += (v1.z + v2.z);
	}

	inline void SubtractSum(const VECTOR_3D<TT>& v1, const VECTOR_3D<TT>& v2)
	{
		x -= (v1.x + v2.x);
		y -= (v1.y + v2.y);
		z -= (v1.z + v2.z);
	}
};

// Miscellaneous Free Operators and Functions
template<class TT> inline static TT FastDistance(const VECTOR_3D<TT>& v1, const VECTOR_3D<TT>& v2)
{
	TT sqr_magnitude = (v1.x - v2.x)*(v1.x - v2.x) + (v1.y - v2.y)*(v1.y - v2.y) + (v1.z - v2.z)*(v1.z - v2.z);
	return sqrt(sqr_magnitude);
}

template<class TT> 
inline static VECTOR_3D<TT> operator*(const TT& a, const VECTOR_3D<TT>& v)
{
	return VECTOR_3D<TT>(a*v.x, a*v.y, a*v.z);
}

template<class TT>
inline static VECTOR_3D<TT> operator*(const int& a, const VECTOR_3D<TT>& v)
{
	return VECTOR_3D<TT>(a*v.x, a*v.y, a*v.z);
}

template<class TT>
inline static VECTOR_3D<TT> operator*(const VECTOR_3D<TT>& v1, const VECTOR_3D<TT>& v2)
{
	return VECTOR_3D<TT>(v1.x*v2.x, v1.y*v2.y, v1.z*v2.z);
}

template<class TT>
inline static VECTOR_3D<TT> operator+(const TT& a, const VECTOR_3D<TT>& v)
{
	return VECTOR_3D<TT>(a + v.x, a + v.y, a + v.z);
}

template<class TT>
inline static VECTOR_3D<TT> operator/(const TT& s, const VECTOR_3D<TT>& v)
{
	TT one_over_s = (TT)1/s;
	return VECTOR_3D<TT>(one_over_s*v.x, one_over_s*v.y, one_over_s*v.z);
}

template<class TT>
inline static TT DotProduct(const VECTOR_3D<TT>& v1, const VECTOR_3D<TT>& v2)
{
	return v1.x*v2.x + v1.y*v2.y + v1.z*v2.z;
}

template<class TT>
inline static VECTOR_3D<TT> CrossProduct(const VECTOR_3D<TT>& v1, const VECTOR_3D<TT>& v2)
{
	return VECTOR_3D<TT>(v1.y*v2.z - v1.z*v2.y, v1.z*v2.x - v1.x*v2.z, v1.x*v2.y - v1.y*v2.x);
}

template<class TT>
inline static bool IsSqrDistanceSmallerThan(const VECTOR_3D<TT>& v1, const VECTOR_3D<TT>& v2, const TT& sqrmagnitude)
{
	const TT diff_x(v1.x - v2.x), diff_y(v1.y - v2.y), diff_z(v1.z - v2.z);
	if(sqrmagnitude > (diff_x*diff_x + diff_y*diff_y + diff_z*diff_z))
		return true;
	else 
		return false;
}

template<class TT>
inline std::ostream& operator << (std::ostream& output, const VECTOR_3D<TT>& v)
{
	return output << v.x << " " << v.y << " " << v.z;
}







