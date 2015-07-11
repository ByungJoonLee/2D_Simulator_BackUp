#pragma once

#include "GRID_STRUCTURE_1D.h"
#include "ARRAY_1D.h"
#include "MULTITHREADING.h"

template<class TT>
class FIELD_STRUCTURE_1D
{
public: // Essential Data
	GRID_STRUCTURE_1D			grid;
	GRID_STRUCTURE_1D			grid_ghost;
	ARRAY<GRID_STRUCTURE_1D>	partial_grids;
	ARRAY<GRID_STRUCTURE_1D>	partial_grids_ghost;
	ARRAY_1D<TT>				array_for_this;

public: // Properties
	int						ghost_width;

public: // Speedup Variables
	int						i_start, i_end, i_start_g, i_end_g;
	int						i_res_g;
	T						x_min;
	T						dx;
	T						one_over_dx;
	T						one_over_2dx;
	T						one_over_dx2;
	TT*						values;

public: // Multithreading
	MULTITHREADING*			multithreading;

public: // Constructors and Destructor
	FIELD_STRUCTURE_1D(void)
		: values(0), multithreading(0) 
	{}

	FIELD_STRUCTURE_1D(const int& i_res_input, const int& i_start_input, const T& x_min_input, const T& x_max_input, const int& ghost_width_input = 0, MULTITHREADING* multithreading_input = 0)
		: partial_grids(0), values(0)
	{
		Initialize(i_res_input, i_start_input, x_min_input, x_max_input, ghost_width_input, multithreading_input);
	}
			
	~FIELD_STRUCTURE_1D(void)
	{}

public: // Initialization Functions
	void Initialize(const int& i_res_input, const int& i_start_input, const T& x_min_input, const T& x_max_input, const int& ghost_width_input = 0, MULTITHREADING* multithreading_input = 0)
	{
		grid.Initialize(i_res_input, i_start_input, x_min_input, x_max_input);

		ghost_width = ghost_width_input;

		grid_ghost.Initialize(grid.Enlarged(ghost_width));

		i_start = grid.i_start;
		i_end = grid.i_end;
		i_start_g = grid_ghost.i_start;
		i_end_g = grid_ghost.i_end;

		x_min = grid.x_min;

		dx = grid.dx;

		one_over_dx = grid.one_over_dx;

		one_over_2dx = grid.one_over_2dx;
		
		one_over_dx2 = grid.one_over_dx2;

		array_for_this.Initialize(grid.i_start - ghost_width, grid.i_res + 2*ghost_width, true);

		values = array_for_this.values;

		i_res_g = array_for_this.i_res;

		if (multithreading_input != 0)
		{
			multithreading = multithreading_input;
			grid.SplitInXDirection(multithreading->num_threads, partial_grids);
			grid_ghost.SplitInXDirection(multithreading->num_threads, partial_grids_ghost);
		}
	}

	void Initialize(const GRID_STRUCTURE_1D& grid_input, const int& ghost_width_input = 0, MULTITHREADING* multithreading_input = 0)
	{
		Initialize(grid_input.i_res, grid_input.i_start, grid_input.x_min, grid_input.x_max, ghost_width_input, multithreading_input);
	}

public: // Operator Overloading
	inline TT& operator [](const int& ix) const
	{
		return array_for_this(ix);
	}

	inline TT& operator ()(const int& ix) const
	{
		return array_for_this(ix);
	}

	inline TT operator ()(const T& position) const
	{
		return LinearInterpolation(position);
	}

public: // Indexing Functions
	inline int ClampI(const int& i) const
	{
		if (i < i_start_g)
		{
			return i_start_g;
		}
		else if (i > i_end_g)
		{
			return i_end_g;
		}

		return i;
	}

public: // Member Functions
	inline T CellCenter(const int& i) const
	{
		return grid.CellCenter(i);
	}

	inline T GridPoint(const int& i) const
	{
		return grid.GridPoint(i);
	}

	inline TT LinearInterpolation(const T& pos) const
	{
		const T &posx(pos);

		const int i0 = (int)((posx - x_min)*one_over_dx);
		
		const T a = (posx - (x_min + (T)(i0 - i_start)*dx))*one_over_dx;
		
		const int ci0(ClampI(i0));
		const int ci1(ClampI(i0+1));

		const TT v0(array_for_this(ci0)), v1(array_for_this(ci1));

		return (1-a)*v0 + a*v1;
	}

	void AssignAllValue(const TT& value)
	{
		for (int i = 0; i < i_res_g; i++)
		{
			array_for_this[i] = value;
		}
	}

	void AssignAllValueGhost(const TT& value)
	{
		for (int i = 0; i < array_for_this.i_res; i++)
		{
			array_for_this.values[i] = value;
		}
	}

	template<class TTT>
	void AssignAllValue(const ARRAY_1D<TTT>& arr, const TTT& value)
	{
		for (int i = 0; i < i_res_g; i++)
		{
			arr[i] = value;
		}
	}

	void AddAllValue(const TT& value)
	{
		for (int i = 0; i < i_res_g; i++)
		{
			array_for_this[i] += value;
		}
	}

public: // Speedup Functions
	inline const int ClampedArrayIndex(const int& i) const
	{
		return CLAMP(i, i_start_g, i_end_g);
	}

	void Translate(const T& deviation) 
	{
		grid.Translate(deviation);
		grid_ghost.Translate(deviation);
		
		x_min += deviation;
	}

public: // Functions filling the ghost cell 
	void FillGhostCellsFrom(ARRAY_1D<TT>& phi_real, const bool& copy_real_data, const int& thread_id)
	{
		// Fill real region
		if (copy_real_data)
		{
			const int i_start(partial_grids[thread_id].i_start), i_end(partial_grids[thread_id].i_end);

			for (int i = i_start; i <= i_end; i++)
			{
				array_for_this(i) = phi_real(i);
			}
		}

		if (ghost_width == 0)
		{
			multithreading->Sync(thread_id);
			return;
		}

		// Fill Ghost Region
		const int i_start(grid.i_start), i_end(grid.i_end);
		const int num_threads(partial_grids.length);
		
		// Face left
		if (0 % num_threads == thread_id)
		{
			for (int i = i_start - ghost_width; i <= i_start - 1; i++)
			{
				array_for_this(i) = phi_real(i_start + 1);
			}
		}
		
		// Face right
		if (1 % num_threads == thread_id)
		{
			for (int i = i_end + 1; i <= i_end + ghost_width; i++)
			{
				array_for_this(i) = phi_real(i_end - 1);
			}
		}

		multithreading->Sync(thread_id);
	}

	void FillGhostCellsFrom(const int& thread_id, ARRAY_1D<TT>& phi_real, const bool& copy_real_data)
	{
		// Fill real region
		if (copy_real_data)
		{
			const int i_start(partial_grids[thread_id].i_start), i_end(partial_grids[thread_id].i_end);

			for (int i = i_start; i <= i_end; i++)
			{
				array_for_this(i) = phi_real(i);
			}
		}

		if (ghost_width == 0)
		{
			multithreading->Sync(thread_id);
			return;
		}

		// Fill Ghost Region
		const int i_start(grid.i_start), i_end(grid.i_end);
		const int num_threads(partial_grids.length);
		
		// Face left
		if (0 % num_threads == thread_id)
		{
			for (int i = i_start - ghost_width; i <= i_start - 1; i++)
			{
				array_for_this(i) = phi_real(i_start + 1);
			}
		}
		
		// Face right
		if (1 % num_threads == thread_id)
		{
			for (int i = i_end + 1; i <= i_end + ghost_width; i++)
			{
				array_for_this(i) = phi_real(i_end - 1);
			}
		}

		multithreading->Sync(thread_id);
	}

	void FillGhostCellsFrom(ARRAY_1D<TT>& phi_real, const bool& copy_real_data)
	{
		// Fill real region
		if (copy_real_data)
		{
			const int i_start(grid.i_start), i_end(grid.i_end);

			for (int i = i_start; i <= i_end; i++)
			{
				array_for_this(i) = phi_real(i);
			}
		}

		// Fill Ghost Region
		const int i_start(grid.i_start), i_end(grid.i_end);

		// Face left
		for (int i = i_start - ghost_width; i <= i_start - 1; i++)
		{
			array_for_this(i) = phi_real(i_start + 1);
		}

		// Face right
		for (int i = i_end + 1; i <= i_end + ghost_width; i++)
		{
			array_for_this(i) = phi_real(i_end - 1);
		}
	}

	void FillGhostCellsContinuousDerivativesFrom(ARRAY_1D<TT>& phi_real, const bool& copy_real_data)
	{
		// Fill real region
		if (copy_real_data)
		{
			const int i_start(grid.i_start), i_end(grid.i_end);

			for (int i = i_start; i <= i_end; i++)
			{
				array_for_this(i) = phi_real(i);
			}
		}

		// Fill Ghost Region
		const int i_start(grid.i_start), i_end(grid.i_end);

		// Face left
		for (int i = i_start - 1; i >= i_start - ghost_width; i--)
		{
			array_for_this(i) = (T)2*phi_real(i + 1) - phi_real(i + 2);
		}

		// Face right
		for (int i = i_end + 1; i <= i_end + ghost_width; i++)
		{
			array_for_this(i) = (T)2*phi_real(i - 1) - phi_real(i - 2);
		}
	}

	void FillGhostCellsContinuousDerivativesFrom(ARRAY_1D<TT>& phi_real, const bool& copy_real_data, const int& thread_id)
	{
		// Fill real region
		if (copy_real_data)
		{
			const int i_start(partial_grids[thread_id].i_start), i_end(partial_grids[thread_id].i_end);

			for (int i = i_start; i <= i_end; i++)
			{
				array_for_this(i) = phi_real(i);
			}
		}

		if (ghost_width == 0)
		{
			multithreading->Sync(thread_id);
			return;
		}

		// Fill Ghost Region
		const int i_start(grid.i_start), i_end(grid.i_end);
		const int num_threads(partial_grids.length);
		
		// Face left
		if (0 % num_threads == thread_id)
		{
			for (int i = i_start - ghost_width; i <= i_start - 1; i++)
			{
				array_for_this(i) = (T)2*phi_real(i + 1) - phi_real(i + 2);
			}
		}
		
		// Face right
		if (1 % num_threads == thread_id)
		{
			for (int i = i_end + 1; i <= i_end + ghost_width; i++)
			{
				array_for_this(i) = (T)2*phi_real(i - 1) - phi_real(i - 2);
			}
		}
		
		multithreading->Sync(thread_id);
	}
};






