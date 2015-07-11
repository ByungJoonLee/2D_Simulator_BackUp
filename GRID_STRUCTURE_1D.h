#pragma once

#include "COMMON_DEFINITION.h"
#include "ARRAY.h"
#include "SCRIPT_READER.h"

class GRID_STRUCTURE_1D
{
public: // Essential Data
	int i_res;
	int i_start, i_end;
	T	x_min, x_max;
	T	dx;
	T	two_dx;
	T	dx2;
	T	one_over_dx;
	T	one_over_2dx;
	T	one_over_dx2;

public: // Constructors and Destructor
	GRID_STRUCTURE_1D(void)
	{}

	GRID_STRUCTURE_1D(const int& i_res_input, const int& i_start_input, const T& x_min_input, const T& x_max_input)
	{
		Initialize(i_res_input, i_start_input, x_min_input, x_max_input);
	}

	GRID_STRUCTURE_1D(const GRID_STRUCTURE_1D& grid_input)
	{
		Initialize(grid_input);
	}

	~GRID_STRUCTURE_1D(void)
	{}

public: // Initialization Functions
	void Initialize(const int& i_res_input, const int& i_start_input, const T& x_min_input, const T& x_max_input)
	{
		i_res = i_res_input;
		i_start = i_start_input;
		i_end = i_start + i_res - 1;
		
		x_min = x_min_input;
		x_max = x_max_input;

		T diff_x = x_max - x_min;

		dx = diff_x/(i_res - 1);
		
		two_dx = 2*dx;
		
		dx2 = dx*dx;

		one_over_dx = (T)1/dx;

		one_over_2dx = (T)1/two_dx;

		one_over_dx2 = (T)1/dx2;
	}

	void Initialize(const GRID_STRUCTURE_1D& grid_input)
	{
		Initialize(grid_input.i_res, grid_input.i_start, grid_input.x_min, grid_input.x_max);
	}

	void InitializeFromBlock(const SCRIPT_BLOCK& block)
	{
		int start = block.GetInteger("start_indices");
		int res = block.GetInteger("base_grid_resolution");
		T min = block.GetFloat("base_grid_min");
		T max = block.GetFloat("base_grid_max");

		Initialize(res, start, min, max);
	}

public: // Operator Overloading
	void operator = (const GRID_STRUCTURE_1D& grid_input)
	{
		Initialize(grid_input);
	}

public: // Index Correction Functions
	inline const int ClampedIndex(const int& i) const
	{
		return CLAMP(i, i_start, i_end);
	}

public: // Member Functions
	inline const T CellCenter(const int& i) const
	{
		return x_min + ((T)0.5 + (T)(i - i_start))*dx;
	}

	inline const T GridPoint(const int& i) const
	{
		return x_min + (i - i_start)*dx;
	}

	inline void Center(const int& i, T& position) const
	{
		position = x_min + ((T)0.5 + (i - i_start))*dx;
	}

	// Note that this one returns index of the cell containing this position
	inline const int Cell(const T& position) const
	{
		return (int)((position - x_min)*one_over_dx);
	}

	inline const void Cell(const T& position, int& i) const
	{
		i = (int)((position - x_min)*one_over_dx);

		assert(i >= i_start && i <= i_end);
	}

	inline const int ClampedCell(const T& position) const
	{
		int index = (int)((position - x_min)*one_over_dx);

		index = CLAMP(index, i_start, i_end);
	
		return index;
	}

	inline const void ClampedCell(const T& position, int& i) const
	{
		i = (int)((position - x_min)*one_over_dx);

		i = CLAMP(i, i_start, i_end);
	}

	inline bool Inside(const T& position) const
	{
		if (position <= x_min)
		{
			return false;
		}
		else if (position >= x_max)
		{
			return false;
		}
		else
		{
			return true;
		}
	}

	inline bool Inside(const T& position, const T& width) const
	{
		if (position <= x_min + width)
		{
			return false;
		}
		else if (position >= x_max - width)
		{
			return false;
		}
		else
		{
			return true;
		}
	}

	inline bool Inside(const int& i) const
	{
		if (i < i_start)
		{
			return false;
		}
		else if (i > i_end)
		{
			return false;
		}
		else
		{
			return true;
		}
	}

	GRID_STRUCTURE_1D Enlarged(const int& width) const
	{
		return GRID_STRUCTURE_1D(i_res + 2*width, i_start - width, x_min - (T)width*dx, x_max + (T)width*dx);
	}

	void Enlarge(const int& width) 
	{
		Initialize(i_res + 2*width, i_start - width, x_min - (T)width*dx, x_max + (T)width*dx);
	}

	void Translate(const T& variation)
	{
		x_min += variation;
		x_max += variation;
	}

public: // For Multithreading
	void SplitInXDirection(const int& num_threads, ARRAY<GRID_STRUCTURE_1D>& partial_grids);
};

// Osteam object overloading
inline std::ostream& operator<<(std::ostream& output, const GRID_STRUCTURE_1D& grid)
{
	output << "GRID_STRUCTURE_1D" << endl;
	output << "- Resolution = " << grid.i_res << endl;
	output << "- Index range =(" << grid.i_start << ") to (" << grid.i_end  << ") " << endl;
	output << "- Range = (" << grid.x_min << ") to (" << grid.x_max << ") " << endl;
	output << "- dx = " << grid.dx << endl;

	return output;
}


