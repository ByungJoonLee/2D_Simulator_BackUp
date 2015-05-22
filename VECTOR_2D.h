#pragma once

#include <iostream>

template<class TT>
class VECTOR_2D
{
public: // Essential member variable
	union
	{
		struct{TT x, y;};
		struct{TT i, j;};
		TT values[2];
	}; // Zero based indexing

public: // Constructors and Destructor
	VECTOR_2D(void)
		: x(TT()), y(TT())
	{}

	VECTOR_2D(const TT& x_input, const TT& y_input)
		: x(TT()), y(TT())
	{
		Initialize(x_input, y_input);
	}

	VECTOR_2D(const VECTOR_2D<TT>& vector_2d_input)
		: x(TT()), y(TT())
	{
		Initialize(vector_2d_input);
	}

	VECTOR_2D(const TT value_input[2])
		: x(TT()), y(TT())
	{
		Initialize(value_input);
	}

	~VECTOR_2D(void)
	{}

public: // Initialization Functions
	void Initialize(const TT& x_input, const TT& y_input)
	{
		x = x_input;
		y = y_input;
	}

	void Initialize(const VECTOR_2D<TT>& vector_2d_input)
	{
		x = vector_2d_input.x;
		y = vector_2d_input.y;
	}

	void Initialize(const TT values_input[2])
	{
		x = values_input[0];
		y = values_input[1];
	}

public: // Operator Overloading
	inline void operator = (const VECTOR_2D<TT>& v)
	{
		x = v.x;
		y = v.y;
	}

	inline void operator += (const VECTOR_2D<TT>& v)
	{
		x += v.x;
		y += v.y;
	}

	inline void operator -= (const VECTOR_2D<TT>& v)
	{
		x -= v.x;
		y -= v.y;
	}

	inline void operator *= (const TT& s)					
	{
		x *= s;
		y *= s;
	}

	inline void operator /= (const TT& s)
	{
		TT one_over_s = (TT)1/s;
		x *= one_over_s;
		y *= one_over_s;
	}

	inline VECTOR_2D operator + (const VECTOR_2D<TT>& v) const
	{
		return VECTOR_2D(x + v.x, y + v.y);
	}

	inline VECTOR_2D operator - (const VECTOR_2D<TT>& v) const
	{
		return VECTOR_2D(x - v.x, y - v.y);
	}

	inline VECTOR_2D operator * (const TT& s) const
	{
		return VECTOR_2D(x*s, y*s);
	}

	inline VECTOR_2D operator / (const TT& s) const
	{
		TT one_over_s = (TT)1/s;

		return VECTOR_2D(x*one_over_s, y*one_over_s);
	}

	inline VECTOR_2D operator * (const VECTOR_2D<TT>& v_input) const
	{
		return VECTOR_2D(x*v_input.x, y*v_input.y);
	}

	inline VECTOR_2D operator / (const VECTOR_2D<TT>& v_input) const
	{
		return VECTOR_2D(x/v_input.x, y/v_input.y);
	}

public: // Member Functions
	inline TT Magnitude(void) const
	{
		return sqrt(x*x + y*y);
	}

	inline TT SqrMagniude(void) const
	{
		return x*x + y*y;
	}

	void MakeThisUnit(void)
	{
		TT Magn = Magnitude();
		TT one_over_Magn = (TT)1/Magn;

		if (Magn != 0)
		{
			x *= one_over_Magn;
			y *= one_over_Magn;
		}
		else
		{
			x = 0;
			y = 0;
		}
	}
};

// Miscellaneous Free Operators and Functions
template<class TT>
inline static VECTOR_2D<TT> operator + (const TT& s, const VECTOR_2D<TT>& v)
{
	return VECTOR_2D<TT>(s + v.x, s + v.y);
}

template<class TT>
inline static VECTOR_2D<TT> operator * (const TT& s, const VECTOR_2D<TT>& v)
{
	return VECTOR_2D<TT>(s*v.x, s*v.y);
}

template<class TT>
inline static VECTOR_2D<TT> operator * (const int& s, const VECTOR_2D<TT>& v)
{
	return VECTOR_2D<TT>(s*v.x, s*v.y);
}

template<class TT>
inline static VECTOR_2D<TT> operator * (const VECTOR_2D<TT>& v1, const VECTOR_2D<TT>& v2)
{
	return VECTOR_2D<TT>(v1.x*v2.x, v1.y*v2.y);
}

template<class TT>
inline static VECTOR_2D<TT> operator / (const TT& s, const VECTOR_2D<TT>& v)
{
	TT one_over_s = (TT)1/s;
	return VECTOR_2D<TT>(one_over_s*v.x, one_over_s*v.y);
}

template<class TT>
inline static TT DotProduct(const VECTOR_2D<TT>& v1, const VECTOR_2D<TT>& v2)
{
	return v1.x*v2.x + v1.y*v2.y;
}

// Display Function just in case
template<class TT>
inline std::ostream& operator << (std::ostream& output, const VECTOR_2D<TT>& v)
{
	return output << v.x << " " << v.y;
}









