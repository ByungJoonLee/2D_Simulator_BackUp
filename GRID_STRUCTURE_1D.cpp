#include "GRID_STRUCTURE_1D.h"

void GRID_STRUCTURE_1D::SplitInXDirection(const int& num_threads, ARRAY<GRID_STRUCTURE_1D>& partial_grids)
	{
		partial_grids.Initialize(num_threads);

		const int quotient = i_res / num_threads;
		const int remainder = i_res % num_threads;

		int i_start_p = i_start;
		T x_min_p = x_min;

		for (int i = 0; i < num_threads; i++)
		{
			const int i_res_p = (i < remainder) ? (quotient + 1) : (quotient);
			const T x_max_p = dx*(T)i_res_p + x_min_p;
			partial_grids[i].Initialize(i_res, i_start_p, x_min_p, x_max_p);
			i_start_p += i_res_p;
			x_min_p = x_max_p;
		}
	}