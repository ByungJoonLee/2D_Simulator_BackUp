#pragma once

#include "COMMON_DEFINITION.h"
#include "ARRAY.h"
#include "SCRIPT_READER.h"

class GRID_STRUCTURE_2D
{
public: // Essential Data
	union // Grid Resolution								
	{
		struct{int i_res, j_res;};
		int res[2];
	};
	
	union // Start Indices
	{
		struct{int i_start, j_start;};
		int start_ix[2];
	};

	union // End Indices
	{
		struct{int i_end, j_end;};
		int end_ix[2];
	};

	union // Grid Domain (Excluding the ghost cells)
	{
		struct{T x_min, y_min, x_max, y_max;};
		struct{T min[2], max[2];};
	};

	union // Grid Spacing
	{
		struct{T dx, dy;};
		T dxdy[2];
	};

	union // Grid Spacing - 2dx, 2dy
	{
		struct{T two_dx, two_dy;};
		T two_dxdy[2];
	};

	union // Grid Spacing - dx^2, dy^2
	{
		struct{T dx2, dy2;};
		T dx2dy2[2];
	};

	union // Inverse of dx, dy
	{
		struct{T one_over_dx, one_over_dy;};
		T one_over_dxdy[2];
	};

	union // Inverse of 2dx, 2dy
	{
		struct{T one_over_2dx, one_over_2dy;};
		T one_over_2dx2dy[2];
	};

	union // Inverse of dx2, dy2
	{
		struct{T one_over_dx2, one_over_dy2;};
		T one_over_dx2dy2[2];
	};

public: // Constructors and Destructor
	GRID_STRUCTURE_2D(void)
	{}

	GRID_STRUCTURE_2D(const int& i_res_input, const int& j_res_input, const int& i_start_input, const int& j_start_input, const T& x_min_input, const T& y_min_input, const T& x_max_input, const T& y_max_input)
	{
		Initialize(i_res_input, j_res_input, i_start_input, j_start_input, x_min_input, y_min_input, x_max_input, y_max_input);
	}

	GRID_STRUCTURE_2D(const GRID_STRUCTURE_2D& grid_input)
	{
		Initialize(grid_input);
	}

	~GRID_STRUCTURE_2D(void)
	{}

public: // Initialization Functions
	void Initialize(const int& i_res_input, const int& j_res_input, const int& i_start_input, const int& j_start_input, const T& x_min_input, const T& y_min_input, const T& x_max_input, const T& y_max_input)
	{
		i_res = i_res_input;
		j_res = j_res_input;
		
		i_start = i_start_input;
		j_start = j_start_input;

		i_end = i_start + i_res - 1;
		j_end = j_start + j_res - 1;

		x_min = x_min_input;
		y_min = y_min_input;
		x_max = x_max_input;
		y_max = y_max_input;

		T diff_x = x_max - x_min;
		T diff_y = y_max - y_min;

		dx = diff_x/(i_res - 1);
		dy = diff_y/(j_res - 1);

		two_dx = 2*dx;
		two_dy = 2*dy;

		dx2 = dx*dx;
		dy2 = dy*dy;

		one_over_dx = (T)1/dx;
		one_over_dy = (T)1/dy;

		one_over_2dx = (T)1/two_dx;
		one_over_2dy = (T)1/two_dy;

		one_over_dx2 = (T)1/dx2;
		one_over_dy2 = (T)1/dy2;
	}

	void Initialize(const GRID_STRUCTURE_2D& grid_input)
	{
		Initialize(grid_input.i_res, grid_input.j_res, grid_input.i_start, grid_input.j_start, grid_input.x_min, grid_input.y_min, grid_input.x_max, grid_input.y_max);
	}

	void InitializeFromBlock(const SCRIPT_BLOCK& block)
	{
		VI start = block.GetInt2("start_indices");
		VI res = block.GetInt2("base_grid_resolution");
		VT min = block.GetVector2("base_grid_min");
		VT max = block.GetVector2("base_grid_max");

		Initialize(res.i, res.j, start.i, start.j, min.i, min.j, max.i, max.j);
	}

public: // Operator Overloading
	void operator = (const GRID_STRUCTURE_2D& grid_input)
	{
		Initialize(grid_input);
	}

public: // Index Correction Functions
	inline const VI ClampedIndex(const int& i, const int& j) const
	{
		return VI(CLAMP(i, i_start, i_end), CLAMP(j, j_start, j_end));
	}

	inline const VI ClampedIndex(const VI& ix) const
	{
		return VI(CLAMP(ix.i, i_start, i_end), CLAMP(ix.j, j_start, j_end));
	}

	inline const int Index1D(const int& i, const int& j) const
	{
		return (i - i_start) + (j - j_start)*i_res;
	}

	inline const int ClampedIndexI(const int& i) const
	{
		return CLAMP(i, i_start, i_end);
	}

	inline const int ClampedIndexJ(const int& j) const
	{
		return CLAMP(j, j_start, j_end);
	}

public: // Member Functions
	inline const VT CellCenter(const int& i, const int& j) const
	{
		return VT(x_min + ((T)0.5 + (T)(i - i_start))*dx, y_min + ((T)0.5 + (T)(j - j_start))*dy);
	}

	inline const VT CellCenter(const VI& ix) const
	{
		return VT(x_min + ((T)0.5 + (T)(ix.i - i_start))*dx, y_min + ((T)0.5 + (T)(ix.j - j_start))*dy);
	}

	inline const VT GridPoint(const int& i, const int& j) const
	{
		return VT(x_min + (i - i_start)*dx, y_min + (j - j_start)*dy);
	}

	inline const VT GridPoint(const VI& ix) const
	{
		return VT(x_min + (ix.i - i_start)*dx, y_min + (ix.j - j_start)*dy);
	}

	inline void Center(const int& i, const int& j, VT& position) const
	{
		position.x = x_min + ((T)0.5 + (i - i_start))*dx;
		position.y = y_min + ((T)0.5 + (j - j_start))*dy;
	}

	inline void Center(const VI& ix, VT& position) const
	{
		position.x = x_min + ((T)0.5 + (ix.i - i_start))*dx;
		position.y = y_min + ((T)0.5 + (ix.j - j_start))*dy;
	}

	// Note that this one returns index of the cell containing this position
	inline const VI Cell(const VI& position) const
	{
		return VI((int)((position.x - x_min)*one_over_dx), (int)((position.y - y_min)*one_over_dy));
	}

	inline const void Cell(const VI& position, int& i, int& j) const
	{
		i = (int)((position.x - x_min)*one_over_dx);
		j = (int)((position.y - y_min)*one_over_dy);

		assert(i >= i_start && i <= i_end);
		assert(j >= j_start && j <= j_end);
	}

	inline const VI ClampedCell(const VT& position) const
	{
		VI index((int)((position.x - x_min)*one_over_dx), (int)((position.y - y_min)*one_over_dy));

		index.i = CLAMP(index.i, i_start, i_end);
		index.j = CLAMP(index.j, j_start, j_end);
		
		return index;
	}

	inline const VI ClampedCell(const VT& position, VI& index) const
	{
		index.i = (int)((position.x - x_min)*one_over_dx);
		index.j = (int)((position.y - y_min)*one_over_dy);

		index.i = CLAMP(index.i, i_start, i_end);
		index.j = CLAMP(index.j, j_start, j_end);

		return index;
	}

	inline const void ClampedCell(const VT& position, int& i, int& j) const
	{
		i = (int)((position.x - x_min)*one_over_dx);
		j = (int)((position.y - y_min)*one_over_dy);

		i = CLAMP(i, i_start, i_end);
		j = CLAMP(j, j_start, j_end);
	}

	inline VI LeftBottomCell(const VT& position) const
	{
		return VI((int)((position.x - x_min)*one_over_dx - (T)0.5), (int)((position.y - y_min)*one_over_dy - (T)0.5));
	}

	inline void LeftBottomCell(const VT& position, int& i, int& j) const
	{
		i = (int)((position.x - x_min)*one_over_dx - (T)0.5);
		j = (int)((position.y - y_min)*one_over_dy - (T)0.5);
	}

	inline void LeftBottomCell(const VT& position, VI& ix) const
	{
		LeftBottomCell(position, ix.i, ix.j);
	}

	inline bool Inside(const VT& position) const
	{
		if (position.x <= x_min)
		{
			return false;
		}
		else if (position.x >= x_max)
		{
			return false;
		}
		else if (position.y <= y_min)
		{
			return false;
		}
		else if (position.y >= y_max)
		{
			return false;
		}
		else
		{
			return true;
		}
	}

	inline bool Inside(const VT& position, const T& width) const
	{
		if (position.x <= x_min + width)
		{
			return false;
		}
		else if (position.x >= x_max - width)
		{
			return false;
		}
		else if (position.y <= y_min + width)
		{
			return false;
		}
		else if (position.y >= y_max - width)
		{
			return false;
		}
		else
		{
			return true;
		}
	}

	inline bool Inside(const int& i, const int& j) const
	{
		if (i < i_start)
		{
			return false;
		}
		else if (i > i_end)
		{
			return false;
		}
		else if (j < j_start)
		{
			return false;
		}
		else if (j > j_end)
		{
			return false;
		}
		else
		{
			return true;
		}
	}

	inline bool Inside(const VI& ix) const
	{
		if (ix.i < i_start)
		{
			return false;
		}
		else if (ix.i > i_end)
		{
			return false;
		}
		else if (ix.j < j_start)
		{
			return false;
		}
		else if (ix.j > j_end)
		{
			return false;
		}
		else
		{
			return true;
		}
	}

	inline bool Inside(const VI& ix, const int& inner_width) const
	{
		if (ix.i < i_start + inner_width)
		{
			return false;
		}
		else if (ix.i > i_end - inner_width)
		{
			return false;
		}
		else if (ix.j < j_start + inner_width)
		{
			return false;
		}
		else if (ix.j > j_end - inner_width)
		{
			return false;
		}
		else
		{
			return true;
		}
	}

	inline bool InsideI(const int& i) const
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

	inline bool InsideJ(const int& j) const
	{
		if (j < j_start)
		{
			return false;
		}
		else if (j > j_end)
		{
			return false;
		}
		else
		{
			return true;
		}
	}

	GRID_STRUCTURE_2D Enlarged(const int& width) const
	{
		return GRID_STRUCTURE_2D(i_res + 2*width, j_res + 2*width, i_start - width, j_start - width, x_min - (T)width*dx, y_min - (T)width*dy, x_max + (T)width*dx, y_max + (T)width*dy);
	}

	void Enlarge(const int& width) 
	{
		Initialize(i_res + 2*width, j_res + 2*width, i_start - width, j_start - width, x_min - (T)width*dx, y_min - (T)width*dy, x_max + (T)width*dx, y_max + (T)width*dy);
	}

	void Translate(const VT& variation)
	{
		x_min += variation.x;
		y_min += variation.y;

		x_max += variation.x;
		y_max += variation.y;
	}

public: // For Multithreading
	void SplitInYDirection(const int& num_threads, ARRAY<GRID_STRUCTURE_2D>& partial_grids);
};

// Osteam object overloading
inline std::ostream& operator<<(std::ostream& output, const GRID_STRUCTURE_2D& grid)
{
	output << "GRID_STRUCTURE_2D" << endl;
	output << "- Resolution = (" << grid.i_res << "," << grid.j_res << ")" << endl;
	output << "- Index range =(" << grid.i_start << "," << grid.j_start << ") to (" << grid.i_end << "," << grid.j_end << ") " << endl;
	output << "- Range = (" << grid.x_min << "," << grid.y_min << ") to (" << grid.x_max << "," << grid.y_max << ") " << endl;
	output << "- dx = (" << grid.dx << "," << grid.dy << ")" << endl;

	return output;
}


