#pragma once

#include "COMMON_DEFINITION.h"

template<class TT>
class ARRAY_1D
{
public: // Index Information
	union 
	{
		struct{int i_start, i_end;};
	};

public: // Data Array
	TT* values;

public: // Speedup Variables
	int i_res;
	
public: // Constructors and Destructor
	ARRAY_1D(void)
		: values(0)
	{}

	ARRAY_1D(const int& i_start_input, const int& i_res_input, const bool& initialize = false)
		: values(0)
	{
		Initialize(i_start_input, i_res_input, initialize);
	}

	~ARRAY_1D(void)
	{
		if (values != 0)
		{
			delete [] values;
		}
	}
	
public: // Initialization Functions
	void Initialize(const int& i_start_input, const int& i_res_input, const bool& initialize = false)
	{
		if (values != 0)
		{
			delete [] values;
		}

		i_start = i_start_input;
		
		i_res = i_res_input;

		i_end = i_start + i_res - 1;


		assert(i_res > 0);
		values = new TT[i_res];

		if (initialize == true)
		{
			AssignAllValues(TT());
		}
	}

public: // Operator Overloading
	TT& operator [](const int& ix) const
	{
		assert(ix >= 0 && ix <= i_res);
		return values[ix];
	}

	TT& operator ()(const int& ix) const
	{
		assert(ix >= 0 && ix <= i_res);
		return values[ix];
	}

	void operator *=(const T& constant) 
	{
		for (int i = 0; i < i_res; i++)
		{
			values[i] *= constant;
		}
	}

	void operator +=(const T& constant)
	{
		for (int i = 0; i < i_res; i++)
		{
			values[i] += constant;
		}
	}

	void operator -=(const T& constant)
	{
		for (int i = 0; i < i_res; i++)
		{
			values[i] -= constant;
		}
	}

	void operator /=(const T& constant)
	{
		T one_over_cons = (T)1/constant;
		for (int i = 0; i < i_res; i++)
		{
			values[i] *= one_over_cons;
		}
	}

public: // Indexing Functions
	inline const VI ClampedIndex(const int& i) const
	{
		return CLAMP(i, i_start, i_end);
	}

public: // Member Functions
	TT& Deviated(const int& base, const int& i_deviation) const
	{
		assert((base + i_deviation) >= 0 && (base + i_deviation) <= i_res);

		return values[base + i_deviation];
	}

	void AssignAllValues(const TT& constant)
	{
		for (int i = 0; i < i_res; i++)
		{
			values[i] = constant;
		}
	}

	void AssignRegionalValues(const TT& constant, const int& i_start_input, const int& i_end_input)
	{
		for (int i = i_start_input; i <= i_end_input; i++)
		{
			values[i] = constant;
		}
	}
};






