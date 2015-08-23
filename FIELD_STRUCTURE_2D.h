#pragma once

#include "GRID_STRUCTURE_2D.h"
#include "ARRAY_2D.h"
#include "MULTITHREADING.h"

template<class TT>
class FIELD_STRUCTURE_2D
{
public: // Essential Data
	GRID_STRUCTURE_2D			grid;
	GRID_STRUCTURE_2D			grid_ghost;
	ARRAY<GRID_STRUCTURE_2D>	partial_grids;
	ARRAY<GRID_STRUCTURE_2D>	partial_grids_ghost;
	ARRAY_2D<TT>				array_for_this;

public: // Properties
	int						ghost_width;

public: // Speedup Variables
	int						i_start, i_end, i_start_g, i_end_g;
	int						j_start, j_end, j_start_g, j_end_g;
	int						i_res_g, ij_res_g;
	T						x_min, y_min;
	T						dx, dy;
	T						one_over_dx, one_over_dy;
	T						one_over_2dx, one_over_2dy;
	T						one_over_dx2, one_over_dy2;
	TT*						values;

public: // Choosing the type
	bool					is_scalar;
	bool					is_vector;

	// Choosing the component
	bool					is_x_component;
	bool					is_y_component;

public: // Multithreading
	MULTITHREADING*			multithreading;

public: // Constructors and Destructor
	FIELD_STRUCTURE_2D(void)
		: values(0), multithreading(0), is_scalar(true), is_vector(false), is_x_component(false), is_y_component(false)
	{}

	FIELD_STRUCTURE_2D(const int& i_res_input, const int& j_res_input, const int& i_start_input, const int& j_start_input, const T& x_min_input, const T& y_min_input, const T& x_max_input, const T& y_max_input, const int& ghost_width_input = 0, const bool& is_scalar_input = true, const bool& is_vector_input = false, MULTITHREADING* multithreading_input = 0)
		: partial_grids(0), is_scalar(true), is_vector(false), is_x_component(false), is_y_component(false), values(0)
	{
		Initialize(i_res_input, j_res_input, i_start_input, j_start_input, x_min_input, y_min_input, x_max_input, y_max_input, ghost_width_input, is_scalar_input, is_vector_input, multithreading_input);
	}
			
	~FIELD_STRUCTURE_2D(void)
	{}

public: // Initialization Functions
	void Initialize(const int& i_res_input, const int& j_res_input, const int& i_start_input, const int& j_start_input, const T& x_min_input, const T& y_min_input, const T& x_max_input, const T& y_max_input, const int& ghost_width_input = 0, const bool& is_scalar_input = true, const bool& is_vector_input = false, MULTITHREADING* multithreading_input = 0)
	{
		grid.Initialize(i_res_input, j_res_input, i_start_input, j_start_input, x_min_input, y_min_input, x_max_input, y_max_input);

		ghost_width = ghost_width_input;

		grid_ghost.Initialize(grid.Enlarged(ghost_width));

		i_start = grid.i_start;
		j_start = grid.j_start;

		i_end = grid.i_end;
		j_end = grid.j_end;

		i_start_g = grid_ghost.i_start;
		j_start_g = grid_ghost.j_start;

		i_end_g = grid_ghost.i_end;
		j_end_g = grid_ghost.j_end;

		x_min = grid.x_min;
		y_min = grid.y_min;

		dx = grid.dx;
		dy = grid.dy;

		one_over_dx = grid.one_over_dx;
		one_over_dy = grid.one_over_dy;

		one_over_2dx = grid.one_over_2dx;
		one_over_2dy = grid.one_over_2dy;
		
		one_over_dx2 = grid.one_over_dx2;
		one_over_dy2 = grid.one_over_dy2;

		array_for_this.Initialize(grid.i_start - ghost_width, grid.j_start - ghost_width, grid.i_res + 2*ghost_width, grid.j_res + 2*ghost_width, true);

		values = array_for_this.values;

		i_res_g = array_for_this.i_res;
		ij_res_g = array_for_this.ij_res;

		is_scalar = is_scalar_input;
		is_vector = is_vector_input;

		if (multithreading_input != 0)
		{
			multithreading = multithreading_input;
			grid.SplitInYDirection(multithreading->num_threads, partial_grids);
			grid_ghost.SplitInYDirection(multithreading->num_threads, partial_grids_ghost);
		}
	}

	void Initialize(const VI& ij_res, const VI& ij_start, const VT& min, const VT& max, const int& ghost_width_input = 0, const bool& is_scalar_input = true, const bool& is_vector_input = false, MULTITHREADING* multithreading_input = 0)
	{
		Initialize(ij_res.i, ij_res.j, ij_start.i, ij_start.j, min.i, min.j, max.i, max.j, ghost_width_input, is_scalar_input, is_vector_input, multithreading_input);
	}

	void Initialize(const GRID_STRUCTURE_2D& grid_input, const int& ghost_width_input = 0, const bool& is_scalar_input = true, const bool& is_vector_input = false, MULTITHREADING* multithreading_input = 0)
	{
		Initialize(grid_input.i_res, grid_input.j_res, grid_input.i_start, grid_input.j_start, grid_input.x_min, grid_input.y_min, grid_input.x_max, grid_input.y_max, ghost_width_input, is_scalar_input, is_vector_input, multithreading_input);
	}

	void Initialize(const GRID_STRUCTURE_2D& grid_input, const int& ghost_width_input, MULTITHREADING* multithreading_input)
	{
		Initialize(grid_input.i_res, grid_input.j_res, grid_input.i_start, grid_input.j_start, grid_input.x_min, grid_input.y_min, grid_input.x_max, grid_input.y_max, ghost_width_input, true, false, multithreading_input);
	}

public: // Operator Overloading
	inline TT& operator [](const int& ix) const
	{
		return array_for_this(ix);
	}

	inline TT& operator ()(const VI& ix) const
	{
		return array_for_this(ix);
	}

	inline TT& operator ()(const int& ix) const
	{
		return array_for_this(ix);
	}

	inline TT& operator ()(const int& i, const int& j) const
	{
		return array_for_this(i, j);
	}

	inline TT operator ()(const VT& position) const
	{
		return BiLinearInterpolation(position);
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

	inline int ClampJ(const int& j) const
	{
		if (j < j_start_g)
		{
			return j_start_g;
		}
		else if (j > j_end_g)
		{
			return j_end_g;
		}

		return j;
	}

	inline const int Index1D(const int& i, const int& j) const
	{
		return array_for_this.Index1D(i, j);
	}

	inline const int Index1D(const VI& ix) const
	{
		return array_for_this.Index1D(ix);
	}

public: // Member Functions
	inline VT CellCenter(const int& i, const int& j) const
	{
		return grid.CellCenter(i, j);
	}

	inline VT GridPoint(const int& i, const int& j) const
	{
		return grid.GridPoint(i, j);
	}

	inline TT BilinearInterpolation(const VT& pos) const
	{
		const T &posx(pos.x), &posy(pos.y);

		const int i0 = (int)((posx - x_min)*one_over_dx);
		const int j0 = (int)((posy - y_min)*one_over_dy);

		const T a = (posx - (x_min + (T)(i0 - i_start)*dx))*one_over_dx;
		const T b = (posy - (y_min + (T)(j0 - j_start)*dy))*one_over_dy;

		const int ci0(ClampI(i0)), cj0(ClampJ(j0));
		const int ci1(ClampI(i0+1)), cj1(ClampJ(j0+1));

		const TT v00(array_for_this(ci0,cj0)), v10(array_for_this(ci1,cj0)), v01(array_for_this(ci0,cj1)), v11(array_for_this(ci1,cj1));

		const T mamb = ((T)1 - a)*((T)1 - b);

		return mamb*v00 + a*((T)1 - b)*v10 + b*((T)1 - a)*v01 + b*a*v11;
	}

	inline void HermiteCubicInterpolation(const VT& pos, const FIELD_STRUCTURE_2D<VT>& gradient_input, T& hermite_cubic, T& hermite_cubic_x, T& hermite_cubic_y) const
	{
		const T &posx(pos.x), &posy(pos.y);

		const int i0 = (int)((posx - x_min)*one_over_dx);
		const int j0 = (int)((posy - y_min)*one_over_dy);

		// Values for derivation
		T center, right, up, rightup;
		T center_x, right_x, up_x, rightup_x;
		T center_y, right_y, up_y, rightup_y;
		T center_xy, right_xy, up_xy, rightup_xy;

		center = array_for_this(i0, j0);
		right = array_for_this(i0 + 1, j0);
		up = array_for_this(i0, j0 + 1);
		rightup = array_for_this(i0 + 1, j0 + 1);

		center_x = gradient_input(i0, j0).x;
		right_x = gradient_input(i0 + 1, j0).x;
		up_x = gradient_input(i0, j0 + 1).x;
		rightup_x = gradient_input(i0 + 1, j0 + 1).x;

		center_y = gradient_input(i0, j0).y;
		right_y = gradient_input(i0 + 1, j0).y;
		up_y = gradient_input(i0, j0 + 1).y;
		rightup_y = gradient_input(i0 + 1, j0 + 1).y;

		center_xy = (gradient_input(i0 + 1, j0).y - gradient_input(i0 - 1, j0).y)*one_over_2dx;
		right_xy = (gradient_input(i0 + 2, j0).y - gradient_input(i0, j0).y)*one_over_2dx;
		up_xy = (gradient_input(i0 + 1, j0 + 1).y - gradient_input(i0 - 1, j0 + 1).y)*one_over_2dx;
		rightup_xy = (gradient_input(i0 + 2, j0 + 1).y - gradient_input(i0, j0 + 1).y)*one_over_2dx;

		// 1D cubic polynomial functions
		T zeta, eta, f_zeta, f_zeta_c, f_eta, f_eta_c, g_zeta, g_zeta_c, g_eta, g_eta_c;
		
		zeta = (posx - (x_min + i0*dx))*grid.one_over_dx;
		eta = (posy - (y_min + j0*dx))*grid.one_over_dy;

		f_zeta = (T)1 - (T)3*POW2(zeta) + (T)2*POW3(zeta);
		f_zeta_c = (T)1 - (T)3*POW2((T)1 - zeta) + (T)2*POW3((T)1 - zeta);
		f_eta = (T)1 - (T)3*POW2(eta) + (T)2*POW3(eta);
		f_eta_c = (T)1 - (T)3*POW2((T)1 - eta) + (T)2*POW3((T)1 - eta);
		g_zeta = zeta*POW2((T)1 - zeta);
		g_zeta_c = ((T)1 - zeta)*POW2(zeta);
		g_eta = eta*POW2((T)1 - eta);
		g_eta_c = ((T)1 - eta)*POW2(eta);

		T f_zeta_d, f_zeta_c_d, g_zeta_d, g_zeta_c_d, f_eta_d, f_eta_c_d, g_eta_d, g_eta_c_d;
		 
		f_zeta_d = (T)-6*zeta + (T)6*POW2(zeta);
		f_zeta_c_d = (T)-6*((T)1 - zeta) + (T)6*POW2((T)1 - zeta);
		g_zeta_d = (T)1 - (T)4*zeta + (T)3*POW2(zeta);
		g_zeta_c_d = (T)1 - (T)4*((T)1 - zeta) + (T)3*POW2((T)1 - zeta); 
		f_eta_d = (T)-6*eta + (T)6*POW2(eta);
		f_eta_c_d = (T)-6*((T)1 - eta) + (T)6*POW2((T)1 - eta);
		g_eta_d = (T)1 - (T)4*eta + (T)3*POW2(eta);
		g_eta_c_d = (T)1 - (T)4*((T)1 - eta) + (T)3*POW2((T)1 - eta); 

		hermite_cubic = (T)0;
		hermite_cubic_x = (T)0;
		hermite_cubic_y = (T)0;
		
		hermite_cubic += center*f_zeta*f_eta + right*f_zeta_c*f_eta + up*f_zeta*f_eta_c + rightup*f_zeta_c*f_eta_c;
		hermite_cubic += grid.dx*(center_x*g_zeta*f_eta + center_y*f_zeta*g_eta - right_x*g_zeta_c*f_eta + right_y*f_zeta_c*g_eta + up_x*g_zeta*f_eta_c - up_y*f_zeta*g_eta_c - rightup_x*g_zeta_c*f_eta_c - rightup_y*f_zeta_c*g_eta_c);
		hermite_cubic += grid.dx2*(center_xy*g_zeta*g_eta - right_xy*g_zeta_c*g_eta - up_xy*g_zeta*g_eta_c + rightup_xy*g_zeta_c*g_eta_c);
		
		hermite_cubic_x += grid.one_over_dx*(center*f_zeta_d*f_eta - right*f_zeta_c_d*f_eta + up*f_zeta_d*f_eta_c - rightup*f_zeta_c_d*f_eta_c);
		hermite_cubic_x += center_x*g_zeta_d*f_eta + center_y*f_zeta_d*g_eta + right_x*g_zeta_c_d*f_eta - right_y*f_zeta_c_d*g_eta + up_x*g_zeta_d*f_eta_c - up_y*f_zeta_d*g_eta_c + rightup_x*g_zeta_c_d*f_eta_c + rightup_y*f_zeta_c_d*g_eta_c;
		hermite_cubic_x += grid.dx*(center_xy*g_zeta_d*g_eta + right_xy*g_zeta_c_d*g_eta - up_xy*g_zeta_d*g_eta_c - rightup_xy*g_zeta_c_d*g_eta_c);

		hermite_cubic_y += grid.one_over_dy*(center*f_zeta*f_eta_d + right*f_zeta_c*f_eta_d - up*f_zeta*f_eta_c_d - rightup*f_zeta_c*f_eta_c_d);
		hermite_cubic_y += center_x*g_zeta*f_eta_d + center_y*f_zeta*g_eta_d - right_x*g_zeta_c*f_eta_d + right_y*f_zeta_c*g_eta_d - up_x*g_zeta*f_eta_c_d + up_y*f_zeta*g_eta_c_d + rightup_x*g_zeta_c*f_eta_c_d + rightup_y*f_zeta_c*g_eta_c_d;
		hermite_cubic_y += grid.dy*(center_xy*g_zeta*g_eta_d - right_xy*g_zeta_c*g_eta_d + up_xy*g_zeta*g_eta_c_d - rightup_xy*g_zeta_c*g_eta_c_d);
	}

	void AssignAllValue(const TT& value)
	{
		for (int i = 0; i < ij_res_g; i++)
		{
			array_for_this[i] = value;
		}
	}

	void AssignAllValueGhost(const TT& value)
	{
		for (int i = 0; i < array_for_this.ij_res; i++)
		{
			array_for_this.values[i] = value;
		}
	}

	template<class TTT>
	void AssignAllValue(const ARRAY_2D<TTT>& arr, const TTT& value)
	{
		for (int i = 0; i < ij_res_g; i++)
		{
			arr[i] = value;
		}
	}

	void AddAllValue(const TT& value)
	{
		for (int i = 0; i < ij_res_g; i++)
		{
			array_for_this[i] += value;
		}
	}
	
public: // Speedup Functions
	inline const TT ArrayValue(const int& i, const int& j) const
	{
		return *(value + (i - i_start_g) + (j - j_start_g)*i_res_g);
	}

	inline const VI ClampedArrayIndex(const int& i, const int& j) const
	{
		return VI(CLAMP(i, i_start_g, i_end_g), CLAMP(j, j_start_g, j_end_g));
	}

	inline const VI ClampedArrayIndex(const VI& ix) const
	{
		return VI(CLAMP(ix.i, i_start_g, i_end_g), CLAMP(ix.j, j_start_g, j_end_g));
	}

	void Translate(const VT& deviation) 
	{
		grid.Translate(deviation);
		grid_ghost.Translate(deviation);
		
		x_min += deviation.x;
		y_min += deviation.y;
	}

public: // Functions filling the ghost cell 
	void FillGhostCellsFrom(ARRAY_2D<TT>& phi_real, const bool& copy_real_data, const int& thread_id)
	{
		// Fill real region
		if (copy_real_data)
		{
			const int i_start(partial_grids[thread_id].i_start), j_start(partial_grids[thread_id].j_start), i_end(partial_grids[thread_id].i_end), j_end(partial_grids[thread_id].j_end);

			for (int j = j_start; j <= j_end; j++)
			{
				for (int i = i_start; i <= i_end; i++)
				{
					array_for_this(i, j) = phi_real(i, j);
				}
			}
		}

		if (ghost_width == 0)
		{
			multithreading->Sync(thread_id);
			return;
		}

		// Fill Ghost Region
		const int i_start(grid.i_start), j_start(grid.j_start), i_end(grid.i_end), j_end(grid.j_end);
		const int num_threads(partial_grids.length);
		
		// Face left
		if (0 % num_threads == thread_id)
		{
			for (int j = j_start; j <= j_end; j++)
			{
				for (int i = i_start - ghost_width; i <= i_start - 1; i++)
				{
					array_for_this(i, j) = phi_real(i_start + 1, j);
				}
			}
		}
		
		// Face right
		if (1 % num_threads == thread_id)
		{
			for (int j = j_start; j <= j_end; j++)
			{
				for (int i = i_end + 1; i <= i_end + ghost_width; i++)
				{
					array_for_this(i, j) = phi_real(i_end - 1 , j);
				}
			}
		}
		
		// Face top
		if (2 % num_threads == thread_id)
		{
			for (int i = i_start; i <= i_end; i++)
			{
				for (int j = j_end + 1; j <= j_end + ghost_width; j++)
				{
					array_for_this(i, j) = phi_real(i, j_end - 1);
				}
			}
		}
		
		// Face bottom
		if (3 % num_threads == thread_id)
		{
			for (int i = i_start; i <= i_end; i++)
			{
				for (int j = j_start - ghost_width; j <= j_start - 1; j++)
				{
					array_for_this(i, j) = phi_real(i, j_start + 1);
				}
			}
		}

		multithreading->Sync(thread_id);
	}

	void FillGhostCellsFrom(const int& thread_id, ARRAY_2D<TT>& phi_real, const bool& copy_real_data)
	{
		// Fill real region
		if (copy_real_data)
		{
			const int i_start(partial_grids[thread_id].i_start), j_start(partial_grids[thread_id].j_start), i_end(partial_grids[thread_id].i_end), j_end(partial_grids[thread_id].j_end);

			for (int j = j_start; j <= j_end; j++)
			{
				for (int i = i_start; i <= i_end; i++)
				{
					array_for_this(i, j) = phi_real(i, j);
				}
			}
		}

		if (ghost_width == 0)
		{
			multithreading->Sync(thread_id);
			return;
		}

		// Fill Ghost Region
		const int i_start(grid.i_start), j_start(grid.j_start), i_end(grid.i_end), j_end(grid.j_end);
		const int num_threads(partial_grids.length);
		
		// Face left
		if (0 % num_threads == thread_id)
		{
			for (int j = j_start; j <= j_end; j++)
			{
				for (int i = i_start - ghost_width; i <= i_start - 1; i++)
				{
					array_for_this(i, j) = phi_real(i_start + 1, j);
				}
			}
		}
		
		// Face right
		if (1 % num_threads == thread_id)
		{
			for (int j = j_start; j <= j_end; j++)
			{
				for (int i = i_end + 1; i <= i_end + ghost_width; i++)
				{
					array_for_this(i, j) = phi_real(i_end - 1 , j);
				}
			}
		}
		
		// Face top
		if (2 % num_threads == thread_id)
		{
			for (int i = i_start; i <= i_end; i++)
			{
				for (int j = j_end + 1; j <= j_end + ghost_width; j++)
				{
					array_for_this(i, j) = phi_real(i, j_end - 1);
				}
			}
		}
		
		// Face bottom
		if (3 % num_threads == thread_id)
		{
			for (int i = i_start; i <= i_end; i++)
			{
				for (int j = j_start - ghost_width; j <= j_start - 1; j++)
				{
					array_for_this(i, j) = phi_real(i, j_start + 1);
				}
			}
		}

		multithreading->Sync(thread_id);
	}

	void FillGhostCellsFrom(ARRAY_2D<TT>& phi_real, const bool& copy_real_data)
	{
		// Fill real region
		if (copy_real_data)
		{
			const int i_start(grid.i_start), j_start(grid.j_start), i_end(grid.i_end), j_end(grid.j_end);

			for (int j = j_start; j <= j_end; j++)
			{
				for (int i = i_start; i <= i_end; i++)
				{
					array_for_this(i, j) = phi_real(i, j);
				}
			}
		}

		// Fill Ghost Region
		const int i_start(grid.i_start), j_start(grid.j_start), i_end(grid.i_end), j_end(grid.j_end);

		// Face left
		for (int j = j_start; j <= j_end; j++)
		{
			for (int i = i_start - ghost_width; i <= i_start - 1; i++)
			{
				array_for_this(i, j) = phi_real(i_start + 1, j);
			}
		}

		// Face right
		for (int j = j_start; j <= j_end; j++)
		{
			for (int i = i_end + 1; i <= i_end + ghost_width; i++)
			{
				array_for_this(i, j) = phi_real(i_end - 1 , j);
			}
		}

		// Face top
		for (int i = i_start; i <= i_end; i++)
		{
			for (int j = j_end + 1; j <= j_end + ghost_width; j++)
			{
				array_for_this(i, j) = phi_real(i, j_end - 1);
			}
		}

		// Face bottom
		for (int i = i_start; i <= i_end; i++)
		{
			for (int j = j_start - ghost_width; j <= j_start - 1; j++)
			{
				array_for_this(i, j) = phi_real(i, j_start + 1);
			}
		}
	}

	void FillGhostCellsContinuousDerivativesFrom(ARRAY_2D<TT>& phi_real, const bool& copy_real_data)
	{
		// Fill real region
		if (copy_real_data)
		{
			const int i_start(grid.i_start), j_start(grid.j_start), i_end(grid.i_end), j_end(grid.j_end);

			for (int j = j_start; j <= j_end; j++)
			{
				for (int i = i_start; i <= i_end; i++)
				{
					array_for_this(i, j) = phi_real(i, j);
				}
			}
		}

		// Fill Ghost Region
		const int i_start(grid.i_start), j_start(grid.j_start), i_end(grid.i_end), j_end(grid.j_end);

		// Face left
		for (int j = j_start; j <= j_end; j++)
		{
			for (int i = i_start - 1; i >= i_start - ghost_width; i--)
			{
				array_for_this(i, j) = (T)2*phi_real(i + 1, j) - phi_real(i + 2, j);
			}
		}

		// Face right
		for (int j = j_start; j <= j_end; j++)
		{
			for (int i = i_end + 1; i <= i_end + ghost_width; i++)
			{
				array_for_this(i, j) = (T)2*phi_real(i - 1 , j) - phi_real(i - 2, j);
			}
		}

		// Face top
		for (int i = i_start; i <= i_end; i++)
		{
			for (int j = j_end + 1; j <= j_end + ghost_width; j++)
			{
				array_for_this(i, j) = (T)2*phi_real(i, j - 1) - phi_real(i, j - 2);
			}
		}

		// Face bottom
		for (int i = i_start; i <= i_end; i++)
		{
			for (int j = j_start - 1; j >= j_start - ghost_width; j--)
			{
				array_for_this(i, j) = (T)2*phi_real(i, j + 1) - phi_real(i, j + 2);
			}
		}
        
	}

	void FillGhostCellsContinuousDerivativesFrom(ARRAY_2D<TT>& phi_real, const bool& copy_real_data, const int& thread_id)
	{
		// Fill real region
		if (copy_real_data)
		{
			const int i_start(partial_grids[thread_id].i_start), j_start(partial_grids[thread_id].j_start), i_end(partial_grids[thread_id].i_end), j_end(partial_grids[thread_id].j_end);

			for (int j = j_start; j <= j_end; j++)
			{
				for (int i = i_start; i <= i_end; i++)
				{
					array_for_this(i, j) = phi_real(i, j);
				}
			}
		}

		if (ghost_width == 0)
		{
			multithreading->Sync(thread_id);
			return;
		}

		// Fill Ghost Region
		const int i_start(grid.i_start), j_start(grid.j_start), i_end(grid.i_end), j_end(grid.j_end);
		const int num_threads(partial_grids.length);
		
		// Face left
		if (0 % num_threads == thread_id)
		{
			for (int j = j_start; j <= j_end; j++)
			{
				for (int i = i_start - ghost_width; i <= i_start - 1; i++)
				{
					array_for_this(i, j) = (T)2*phi_real(i + 1, j) - phi_real(i + 2, j);
				}
			}
		}
		
		// Face right
		if (1 % num_threads == thread_id)
		{
			for (int j = j_start; j <= j_end; j++)
			{
				for (int i = i_end + 1; i <= i_end + ghost_width; i++)
				{
					array_for_this(i, j) = (T)2*phi_real(i - 1 , j) - phi_real(i - 2, j);
				}
			}
		}
		
		// Face top
		if (2 % num_threads == thread_id)
		{
			for (int i = i_start; i <= i_end; i++)
			{
				for (int j = j_end + 1; j <= j_end + ghost_width; j++)
				{
					array_for_this(i, j) = (T)2*phi_real(i, j - 1) - phi_real(i, j - 2);
				}
			}
		}
		
		// Face bottom
		if (3 % num_threads == thread_id)
		{
			for (int i = i_start; i <= i_end; i++)
			{
				for (int j = j_start - ghost_width; j <= j_start - 1; j++)
				{
					array_for_this(i, j) = (T)2*phi_real(i, j + 1) - phi_real(i, j + 2);
				}
			}
		}

		multithreading->Sync(thread_id);
	}
};






