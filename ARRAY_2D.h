#pragma once

#include "COMMON_DEFINITION.h"

template<class TT>
class ARRAY_2D
{
public: // Index Information
	union 
	{
		struct{int i_start, j_start, i_end, j_end;};
		struct{int ix_start[2], ix_end[2];};
	};

public: // Data Array
	TT* values;

public: // Speedup Variables
	int i_res, j_res;
	int ij_res;

public: // Constructors and Destructor
	ARRAY_2D(void)
		: values(0)
	{}

	ARRAY_2D(const int& i_start_input, const int& j_start_input, const int& i_res_input, const int& j_res_input, const bool& initialize = false)
		: values(0)
	{
		Initialize(i_start_input, j_start_input, i_res_input, j_res_input, initialize);
	}

	ARRAY_2D(const VI& start_ix_input, const VI& res_ix_input, const bool& initialize = false)
		: values(0)
	{
		Initialize(start_ix_input, res_ix_input, initialize);
	}

	~ARRAY_2D(void)
	{
		if (values != 0)
		{
			delete [] values;
		}
	}
	
public: // Initialization Functions
	void Initialize(const int& i_start_input, const int& j_start_input, const int& i_res_input, const int& j_res_input, const bool& initialize = false)
	{
		if (values != 0)
		{
			delete [] values;
		}

		i_start = i_start_input;
		j_start = j_start_input;

		i_res = i_res_input;
		j_res = j_res_input;

		i_end = i_start + i_res - 1;
		j_end = j_start + j_res - 1;

		ij_res = i_res*j_res;

		assert(i_res > 0 && j_res > 0);
		values = new TT[ij_res];

		if (initialize == true)
		{
			AssignAllValues(TT());
		}
	}

	void Initialize(const VI& start_ix_input, const VI& res_ix_input, const bool& Initialize = false)
	{
		Initialize(start_ix_input.i, start_ix_input.j, res_ix_input.i, res_ix_input.j, Initialize);
	}

public: // Operator Overloading
	TT& operator [](const int& ix) const
	{
		assert(ix >= 0 && ix <= ij_res);
		return values[ix];
	}

	TT& operator ()(const int& ix) const
	{
		assert(ix >= 0 && ix <= ij_res);
		return values[ix];
	}

	TT& operator ()(const int& i, const int& j) const
	{
		assert(i >= i_start && i <= i_end);
		assert(j >= j_start && j <= j_end);

		return values[Index1D(i,j)];
	}

	TT& operator ()(const VI& ix) const
	{
		return (*this)(ix.i, ix.j);
	}

	void operator *=(const T& constant) 
	{
		for (int i = 0; i < ij_res; i++)
		{
			values[i] *= constant;
		}
	}

	void operator +=(const T& constant)
	{
		for (int i = 0; i < ij_res; i++)
		{
			values[i] += constant;
		}
	}

	void operator -=(const T& constant)
	{
		for (int i = 0; i < ij_res; i++)
		{
			values[i] -= constant;
		}
	}

	void operator /=(const T& constant)
	{
		T one_over_cons = (T)1/constant;
		for (int i = 0; i < ij_res; i++)
		{
			values[i] *= one_over_cons;
		}
	}

public: // Indexing Functions
	const int Index1D(const int& i, const int& j) const
	{
		assert(i >= i_start && i <= i_end);
		assert(j >= j_start && j <= j_end);

		return (i - i_start) + (j - j_start)*i_res;
	}

	const int Index1D(const VI& ix) const
	{
		assert(ix.i >= i_start && ix.i <= i_end);
		assert(ix.j >= j_start && ix.j <= j_end);
		
		return (ix.i - i_start) + (ix.j - j_start)*i_res;
	}

	inline const VI ClampedIndex(const int& i, const int& j) const
	{
		return VI(CLAMP(i, i_start, i_end), CLAMP(j, j_start, j_end));
	}

	inline const VI ClampedIndex(const VI& ix) const
	{
		return ClampedIndex(ix.i, ix.j);
	}

public: // Member Functions
	TT& DeviatedX(const int& base, const int& i_deviation) const
	{
		assert((base + i_deviation) >= 0 && (base + i_deviation) <= ij_res);

		return values[base + i_deviation];
	}

	TT& DeviatedY(const int& base, const int& j_deviation) const
	{
		assert((base + j_deviation) >= 0 && (base + j_deviation) <= ij_res);

		return values[base + j_deviation];
	}

	void AssignAllValues(const TT& constant)
	{
		for (int i = 0; i < ij_res; i++)
		{
			values[i] = constant;
		}
	}

	void AssignRegionalValues(const TT& constant, const int& i_start, const int& j_start, const int& i_end, const int& j_end)
	{
		for (int i = i_start; i <= i_end; i++)
		{
			for (int j = j_start; j <= j_end; j++)
			{
				values[Index1D(i,j)] = constant;
			}
		}
	}
};






