#pragma once

#include "COMMON_DEFINITION.h"

template<class TT>
class ARRAY
{
public: // Essential Data
	int length;
	TT* values;

public: // Constructors and Destructor
	ARRAY(void)
		: length(0), values(0)
	{}

	ARRAY(const int& length_input)
		: length(0), values(0)
	{
		Initialize(length_input);
	}

	ARRAY(const int& length_input, const TT& values_input)
		: length(0), values(0)
	{
		Initialize(length_input, values_input);
	}

	ARRAY(const ARRAY<TT>& array_input)
		: length(0), values(0)
	{
		Initialize(array_input);
	}

	~ARRAY(void)
	{
		if (values != 0)
		{
			delete [] values;
		}
		
		length = 0;
	}

public: // Initialization Functions
	inline void Initialize(const int& length_input)
	{
		length = length_input;

		if (values != 0)
		{
			delete [] values;
			values = 0;
		}

		if (length > 0)
		{
			values = new TT[length];
		}
	}

	inline void Initialize(const int& length_input, const TT& values_input)
	{
		length = length_input;

		if (values != 0)
		{
			delete [] values;
			values = 0;
		}

		if (length > 0)
		{
			values = new TT[length];
			AssignAllValues(values_input);
		}
	}

	inline void Initialize(const ARRAY<TT>& array_input)
	{
		Initialize(array_input.length);

		CopyFrom(array_input);
	}

public: // Operator Overloading
	inline TT& operator [] (const int& i) const
	{
		return values[i];
	}

public: // Member Functions
	void AssignAllValues(const TT& values_input)
	{
		for (int i = 0; i < length; i++)
		{
			values[i] = values_input;
		}
	}

	void AssignValues(const int& start_ix, const int& end_ix, const TT& value_input)
	{
		for (int i = start_ix; i <= end_ix; i++)
		{
			values[i] = values_input;
		}
	}

	void CopyFrom(const ARRAY<TT>& array_input)
	{
		assert(length == array_input.length);

		TT* from_val = array_input.values;

		for (int i = 0; i < length; i++)
		{
			values[i] = from_val[i];
		}
	}
};



