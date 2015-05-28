#include "GRID_STRUCTURE_2D.h"

void GRID_STRUCTURE_2D::SplitInYDirection(const int& num_threads, ARRAY<GRID_STRUCTURE_2D>& partial_grids)
	{
		partial_grids.Initialize(num_threads);

		const int quotient = j_res / num_threads;
		const int remainder = j_res % num_threads;

		int j_start_p = j_start;
		T y_min_p = y_min;

		for (int j = 0; j < num_threads; j++)
		{
			const int j_res_p = (j < remainder) ? (quotient + 1) : (quotient);
			const T y_max_p = dy*(T)j_res_p + y_min_p;
			partial_grids[j].Initialize(i_res, j_res_p, i_start, j_start_p, x_min, y_min_p, x_max, y_max_p);
			j_start_p += j_res_p;
			y_min_p = y_max_p;
		}
	}