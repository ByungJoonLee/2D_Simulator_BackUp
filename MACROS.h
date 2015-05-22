#pragma once

enum  POISSON_SOLVER_TYPE                                   {NO_SOLVER, CG, PCG};

#define PI													(T)3.141592
#define BC_FULL												0
#define BC_DIR												-1
#define BC_OBJ												-2
#define BC_NULL												-3
#define BC_NEUM												-4

#define CLAMP(v, min, max)									((v) < (min) ? (min) : ((v) > (max) ? (max) : (v)))

#define POW2(a)												((a)*(a))
#define POW3(a)												((a)*(a)*(a))

#define MAX(a, b)											((a) > (b) ? (a) : (b))
#define MIN(a, b)											((a) > (b) ? (b) : (a))

#define INCREASING_SORT2(a, b, a1, a2)						if (a <= b){a1 = a; a2 = b;}                 \
															else{a1 = b; a2 = a;}                        

#define DELETE_POINTER(pointer)								if (pointer != 0) {delete pointer; pointer = 0;}
#define DELETE_ARRAY(pointer)								if (pointer != 0) {delete [] pointer; pointer = 0;}

#define LOOPS_2D(i, j, i_start, j_start, i_end, j_end)      for((j) = (j_start) ; (j) <= (j_end); ++(j)) for((i) = (i_start); (i) <= (i_end); ++(i))

#define GRID_ITERATION_2D(grid_2d_input)					for(int j_start = (grid_2d_input).j_start, j_end = (grid_2d_input).j_end, i_start = (grid_2d_input).i_start, i_end = (grid_2d_input).i_end, i, j = j_start; j <= j_end; ++j) for(i = i_start; i <= i_end; ++i)

#define BEGIN_GRID_ITERATION_2D(grid_2d_input)				{GRID_STRUCTURE_2D& grid_2d_itr(grid_2d_input);																								\
															 int i, j;																																	\
															 const int j_start = grid_2d_itr.j_start, j_end = grid_2d_itr.j_end, i_start = grid_2d_itr.i_start, i_end = grid_2d_itr.i_end;				\
															 for (int j = j_start; j <= j_end; ++j) for (int i = i_start; i <= i_end; ++i) 
#define END_GRID_ITERATION_2D								 multithreading->Sync(thread_id);}


															