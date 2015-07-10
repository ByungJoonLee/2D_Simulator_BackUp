#pragma once

#include "COMMON_DEFINITION.h"
#include "FIELD_STRUCTURE_1D.h"
#include "DYNAMIC_ARRAY.h"
#include "LEVELSET_OBJECT.h"

class LEVELSET_1D;

class LEVELSET_1D : public LEVELSET_OBJECT
{
public: // Essential Data
	FIELD_STRUCTURE_1D<T>			signed_distance_field;
	FIELD_STRUCTURE_1D<T>			scalar_field_ghost;
	
	// For test
	FIELD_STRUCTURE_1D<T>			exact_solution;

	GRID_STRUCTURE_1D&				grid;
	ARRAY<GRID_STRUCTURE_1D>&		partial_grids;
	ARRAY<GRID_STRUCTURE_1D>&		partial_grids_ghost;
	int&							ghost_width;
	ARRAY_1D<T>						&phi, &arr, phi_true;

	FIELD_STRUCTURE_1D<T>			normal;
	FIELD_STRUCTURE_1D<T>			curvature;
	
	// For curvature calculation
	FIELD_STRUCTURE_1D<T>           phi_x, phi_xx;
	
	int								sweep_direction;

public: // The way of calculating curvature vector
	bool							curvature_by_normal_vector;
	bool							curvature_by_levelset;

public: // Multithreading
	MULTITHREADING*					multithreading;

public: // Constructors and Destructor
	LEVELSET_1D(void)
		: grid(signed_distance_field.grid), partial_grids(signed_distance_field.partial_grids), partial_grids_ghost(signed_distance_field.partial_grids_ghost), phi(signed_distance_field.array_for_this), arr(signed_distance_field.array_for_this), phi_true(signed_distance_field.array_for_this), sweep_direction(0), ghost_width(signed_distance_field.ghost_width)
	{}
	
	~LEVELSET_1D(void)
	{}

public: // Initialization Functions
	void Initialize(const GRID_STRUCTURE_1D& grid_input, const int& ghost_width_input, MULTITHREADING* multithreading_input = 0)
	{
		assert(ghost_width_input >= 2);

		multithreading = multithreading_input;
		signed_distance_field.Initialize(grid_input, ghost_width_input, multithreading_input);
		scalar_field_ghost.Initialize(grid_input, ghost_width_input, multithreading_input);
		exact_solution.Initialize(grid_input, ghost_width_input, multithreading_input);

		normal.Initialize(grid_input, 1, multithreading_input);
		curvature.Initialize(grid_input, 0, multithreading_input);

		phi_x.Initialize(signed_distance_field.grid, 2, multithreading_input);
		phi_xx.Initialize(signed_distance_field.grid, 2, multithreading_input);
		
		phi.AssignAllValues(grid.dx);
		phi_true.Initialize(phi.i_start, phi.i_res, true);

		sweep_direction = 0;

		curvature_by_normal_vector = false;
		curvature_by_levelset = false;
	}

	void Initialize(const int& i_start_input, const int& i_res_input, const T& x_min_input, const T& x_max_input, const int& ghost_width_input, MULTITHREADING* multithreading_input = 0)
	{
		assert(ghost_width_input >= 2);

		multithreading = multithreading_input;
		
		signed_distance_field.Initialize(i_res_input, i_start_input, x_min_input, x_max_input, ghost_width_input, multithreading_input);
		scalar_field_ghost.Initialize(i_res_input, i_start_input, x_min_input, x_max_input, ghost_width_input, multithreading_input);
		exact_solution.Initialize(i_res_input, i_start_input, x_min_input, x_max_input, ghost_width_input, multithreading_input);

		normal.Initialize(i_res_input + 2, i_start_input - 1, x_min_input - grid.dx, x_max_input + grid.dx, 0, multithreading_input);
		curvature.Initialize(i_res_input, i_start_input, x_min_input, x_max_input, 0, multithreading_input);

		phi_x.Initialize(signed_distance_field.grid, 2, multithreading_input);
		phi_xx.Initialize(signed_distance_field.grid, 2, multithreading_input);

		phi.AssignAllValues(grid.dx);
		phi_true.Initialize(phi.i_start, phi.i_res, true);

		sweep_direction = 0;

		curvature_by_normal_vector = false;
		curvature_by_levelset = false;
	}

public: // Operator Overloading
	inline T& operator ()(const int& i) const
	{
		return phi(i);
	}

	inline T operator ()(const T& position) const
	{
		return LinearInterpolation(position);
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

	inline T LinearInterpolation(const T& position) const
	{
		return signed_distance_field.LinearInterpolation(position);
	}

	void Translate(const T& deviation)
	{
		signed_distance_field.Translate(deviation);
		normal.Translate(deviation);
		curvature.Translate(deviation);
	}

	inline bool InsideOBB(const T& position) const
	{
		return grid.Inside(position);
	}

	inline const T SignedDistance(const T& position) const
	{
		return LinearInterpolation(position);
	}

	inline const T SignedDistance(const int& i) const
	{
		return signed_distance_field(i);
	}
	
	void AssignAllValuesLevelset(const T& value)
	{
		for (int i = grid.i_start; i <= grid.i_end; i++)
		{
			phi(i) = value;
		}
	}

	void CopyAllValuesFrom(const LEVELSET_1D& levelset_input)
	{
		for (int i = grid.i_start; i <= grid.i_end; i++)
		{
			phi(i) = levelset_input.phi(i);
		}
	}

	void FillGhostCellsFrom(ARRAY_1D<T>& phi_real, const bool& copy_real_data)
	{
		signed_distance_field.FillGhostCellsFrom(phi_real, copy_real_data);
	}

	void FillGhostCellsFrom(const int& thread_id, ARRAY_1D<T>& phi_real, const bool& copy_real_data)
	{
		signed_distance_field.FillGhostCellsFrom(phi_real, copy_real_data, thread_id);
	}

	void FillGhostCellsFromPointer(ARRAY_1D<T>* phi_real, const bool& copy_real_data)
	{
		signed_distance_field.FillGhostCellsFrom(*phi_real, copy_real_data);
	}
	
	void FillGhostCellsFromPointerThreaded(const int& thread_id, ARRAY_1D<T>* phi_real, const bool& copy_real_data)
	{
		signed_distance_field.FillGhostCellsFrom(thread_id, *phi_real, copy_real_data);
	}

	void FillGhostCellsFromThreaded(ARRAY_1D<T>* phi_real, const bool& copy_real_data)
	{
		multithreading->RunThreads(&LEVELSET_1D::FillGhostCellsFromPointerThreaded, this, phi_real, copy_real_data);
	}

	void FillGhostCellsContinuousDerivativesFrom(ARRAY_1D<T>& phi_real, const bool& copy_real_data)
	{
		signed_distance_field.FillGhostCellsContinuousDerivativesFrom(phi_real, copy_real_data);
	}

	void FillGhostCellsContinuousDerivativesFrom(ARRAY_1D<T>& phi_real, const bool& copy_real_data, const int& thread_id)
	{
		signed_distance_field.FillGhostCellsContinuousDerivativesFrom(phi_real, copy_real_data, thread_id);
	}

	void FillGhostCellsContinuousDerivativesFromPointer(ARRAY_1D<T>* phi_real, const bool& copy_real_data)
	{
		signed_distance_field.FillGhostCellsContinuousDerivativesFrom(*phi_real, copy_real_data);
	}

	void ComputeNormals()
	{
		for (int i = grid.i_start; i <= grid.i_end; i++)
		{
			T& nor(normal(i));

			T nor_x = (phi(i+1) - phi(i-1))*signed_distance_field.one_over_2dx;
			T mag_of_normal = abs(nor_x);

			if (mag_of_normal != 0)
			{
				nor = nor_x;
			}
			else
			{
				nor = (phi(i+1) - phi(i))*signed_distance_field.one_over_dx;
			}

			nor = nor/abs(nor);			
		}
	}

	void ComputeNormals(const int& thread_id)
	{
		signed_distance_field.FillGhostCellsFrom(phi, false, thread_id);

		BEGIN_GRID_ITERATION_1D(normal.partial_grids[thread_id])
		{
			T& nor(normal(i));

			T nor_x = (phi(i+1) - phi(i-1))*signed_distance_field.one_over_2dx;
			T mag_of_normal = abs(nor_x);

			if (mag_of_normal != 0)
			{
				nor = nor_x;
			}
			else
			{
				nor = (phi(i+1) - phi(i))*signed_distance_field.one_over_dx;
			}
			nor = nor/abs(nor);
		}	
		END_GRID_ITERATION_1D;
	}

	void ComputeCurvatures(void)
	{
		T tolerance = (T)1/grid.dx;
		
		// Curvature by the definition of divergence of normal
		if (curvature_by_normal_vector == true)
		{
			normal.FillGhostCellsFrom(normal.array_for_this, true);

			GRID_ITERATION_1D(grid)
			{
				T curv = (normal(i+1) - normal(i-1))*signed_distance_field.one_over_2dx;
								
				if (abs(curv) <= tolerance)
				{
					curvature(i) = -curv;
				}
				else
				{
					curvature(i) = -tolerance;
				}
			}
		}
		// Curvature calculated by levelset
		//if (curvature_by_levelset == true)
		//{
		//	GRID_ITERATION_1D(signed_distance_field.grid)
		//		{
		//			phi_x(i) = (phi(i + 1) - phi(i - 1))*signed_distance_field.one_over_2dx;
		//			phi_xx(i) = (phi(i + 1) - 2*phi(i) + phi(i - 1))*signed_distance_field.one_over_dx2;
		//		}
		//	
		//		T px, pxx;
		//	
		//		for (int i = grid.i_start; i <= grid.i_end; i++)
		//		{
		//			px = phi_x(i);
		//			pxx = phi_xx(i);
		//									
		//			T magnitude = abs(px);
		//			T deno = POW3(magnitude);
		//				
		//			T curv(0);
	
		//			if (deno != 0)
		//			{
		//				T one_over_deno = (T)1/deno;
		//				curv = -(POW2(px)*pyy - (T)2*px*py*pxy + POW2(py)*pxx)*one_over_deno;
		//			}
		//			else
		//				{
		//					px = (phi(i + 1, j) - phi(i, j))*grid.one_over_dx, py = (phi(i, j + 1) - phi(i, j))*grid.one_over_dy;
		//					magnitude = sqrt(POW2(px) + POW2(py));
		//					deno = POW3(magnitude);
		//					if (deno == 0)
		//					{
		//						cout << "Denominator cannot be zero!!" << endl;
		//						exit(0);
		//					}
		//					T one_over_deno = (T)1/deno;
		//					curv = -(POW2(px)*pyy - (T)2*px*py*pxy + POW2(py)*pxx)*one_over_deno;
		//				}
		//
		//				if (abs(curv) <= tolerance)
		//				{
		//					curvature(i, j) = curv;
		//				}
		//				else if (curv > tolerance)
		//				{
		//					curvature(i, j) = tolerance;
		//				}
		//				else if (curv < -tolerance)
		//				{
		//					curvature(i, j) = -tolerance;
		//				}
	
		//				//curvature(i, j) = curv;
		//			}
		//		}
		//	}
		//}
	}
	
	void ComputeCurvaturesThread(const int& thread_id)
	{
		T tolerance = (T)1/grid.dx;
		
		// Curvature by the definition of divergence of normal
		if (curvature_by_normal_vector == true)
		{
			normal.FillGhostCellsFrom(normal.array_for_this, true, thread_id);

			BEGIN_GRID_ITERATION_1D(curvature.partial_grids[thread_id])
			{
				T curv = (normal(i+1) - normal(i-1))*signed_distance_field.one_over_2dx;
								
				if (abs(curv) <= tolerance)
				{
					curvature(i) = -curv;
				}
				else
				{
					curvature(i) = -tolerance;
				}
			}
			END_GRID_ITERATION_2D;
		}
		
		signed_distance_field.FillGhostCellsContinuousDerivativesFrom(signed_distance_field.array_for_this, true, thread_id);

		// Curvature calculated by levelset
		//if (curvature_by_levelset == true)
		//{
		//	BEGIN_GRID_ITERATION_2D(signed_distance_field.partial_grids[thread_id])
		//	{
		//		phi_x(i, j) = (phi(i + 1, j) - phi(i - 1, j))*signed_distance_field.one_over_2dx;
		//		phi_y(i, j) = (phi(i, j + 1) - phi(i, j - 1))*signed_distance_field.one_over_2dy;
		//		phi_xx(i, j) = (phi(i + 1, j) - 2*phi(i, j) + phi(i - 1, j))*signed_distance_field.one_over_dx2;
		//		phi_yy(i, j) = (phi(i, j + 1) - 2*phi(i, j) + phi(i, j - 1))*signed_distance_field.one_over_dy2;
		//		phi_xy(i, j) = (phi(i + 1, j + 1) - phi(i + 1, j - 1) - phi(i - 1, j + 1) + phi(i - 1, j - 1))*signed_distance_field.one_over_2dx*signed_distance_field.one_over_2dy;		
		//	}
		//	END_GRID_ITERATION_2D;

		//	T px, py, pxx, pyy, pxy;
		//	
		//	BEGIN_GRID_ITERATION_2D(signed_distance_field.partial_grids[thread_id])
		//	{
		//		px = phi_x(i, j);
		//		py = phi_y(i, j);
		//		pxx = phi_xx(i, j);
		//		pyy = phi_yy(i, j);
		//		pxy = phi_xy(i, j);
		//				
		//		T magnitude = sqrt(POW2(px) + POW2(py));
		//		T deno = POW3(magnitude);
		//				
		//		T curv(0);
	
		//		if (deno != 0)
		//		{
		//			T one_over_deno = (T)1/deno;
		//			curv = -(POW2(px)*pyy - (T)2*px*py*pxy + POW2(py)*pxx)*one_over_deno;
		//		}
		//		else
		//		{
		//			px = (phi(i + 1, j) - phi(i, j))*grid.one_over_dx, py = (phi(i, j + 1) - phi(i, j))*grid.one_over_dy;
		//			magnitude = sqrt(POW2(px) + POW2(py));
		//			deno = POW3(magnitude);
		//			if (deno == 0)
		//			{
		//				cout << "Denominator cannot be zero!!" << endl;
		//				exit(0);
		//			}
		//			T one_over_deno = (T)1/deno;
		//			curv = -(POW2(px)*pyy - (T)2*px*py*pxy + POW2(py)*pxx)*one_over_deno;
		//		}
		//
		//		if (abs(curv) <= tolerance)
		//		{
		//			curvature(i, j) = curv;
		//		}
		//		else if (curv > tolerance)
		//		{
		//			curvature(i, j) = tolerance;
		//		}
		//		else if (curv < -tolerance)
		//		{
		//			curvature(i, j) = -tolerance;
		//		}
	
		//		//curvature(i, j) = curv;
		//	}
		//	END_GRID_ITERATION_2D;
		//}
	}

	void ComputeCurvaturesThreaded()
	{
		multithreading->RunThreads(&LEVELSET_1D::ComputeCurvaturesThread, this);
	}

	inline const T Normal(const T& position) const
	{
		return normal.LinearInterpolation(position);
	}

	inline void UnitNormal(const T& position, T& normal_output) const
	{
		normal_output = normal.LinearInterpolation(position);
		normal_output = normal_output/abs(normal_output);
	}
	
	inline T Curvature(const T& position) const
	{
		return curvature.LinearInterpolation(position); 
	}

public: // Error Estimate Functions (For the Test)
	T ErrorInf(void)
	{
		T err_inf(0);

		for (int i = grid.i_start; i <= grid.i_end; i++)
		{
			T temp = abs(phi(i) - phi_true(i));
			err_inf = MAX(err_inf, temp);
		}

		return err_inf;
	}

	T Error2Norm(void)
	{
		T err_2norm(0);

		for (int i = grid.i_start; i <= grid.i_end; i++)
		{
			err_2norm += POW2(abs(phi(i) - phi_true(i)));
		}

		err_2norm = grid.dx*sqrt(err_2norm);

		return err_2norm;
	}

	T Error1Norm(void)
	{
		T err_1norm(0);

		for (int i = grid.i_start; i <= grid.i_end; i++)
		{
			err_1norm += abs(phi(i) - phi_true(i));
		}

		err_1norm = grid.dx*grid.dx*err_1norm;

		return err_1norm;
	}
};

