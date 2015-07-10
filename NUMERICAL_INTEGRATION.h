#pragma once

#include "MATRIX_2X2.h"
#include "LEVELSET_2D.h"
#include "WORLD_DISCRETIZATION_2D.h"

class NUMERICAL_INTEGRATION
{
public: // Grids
	GRID_STRUCTURE_2D				base_grid;
	GRID_STRUCTURE_2D				ghost_grid;

	ARRAY<GRID_STRUCTURE_2D>		partial_base_grids;
	ARRAY<GRID_STRUCTURE_2D>		partial_ghost_grids;

public: 
	LEVELSET_2D*					object_levelset;

	WORLD_DISCRETIZATION_2D*		world_discretization;

	FIELD_STRUCTURE_2D<T>			integral_value_array;
	T								integral_value;
	
public: // For error check
	T								true_integral_value;
	T								relative_error;

public: // For Debug
	DYNAMIC_ARRAY<VECTOR_2D<T>>		inserted_point;

public: // Multithreading
	MULTITHREADING*					multithreading;

public: // Constructor and Destructor
	NUMERICAL_INTEGRATION(void)
		: object_levelset(0), multithreading(0), integral_value((T)0)
	{}

	~NUMERICAL_INTEGRATION(void)
	{
		DELETE_POINTER(object_levelset);
	}

public: // Initialization Functions
	void InitializeFromScriptBlock(const SCRIPT_BLOCK& script_block, MULTITHREADING* multithreading_input)
	{
		// Initialize grids from outside
		SCRIPT_BLOCK grid_sb = script_block.FindBlock("GRID_STRUCTURE_2D");

		base_grid.InitializeFromBlock(grid_sb);

		InitializeFromScriptBlock(base_grid, script_block, multithreading_input);
	}

	void InitializeFromScriptBlock(const GRID_STRUCTURE_2D& world_grid, const SCRIPT_BLOCK& numerical_integration_block, MULTITHREADING* multithreading_input)
	{
		// Simulation properties from script
		int ghost_width = numerical_integration_block.GetInteger("ghost_width", 3);
		
		// Grid
		base_grid.Initialize(world_grid.i_res - 1, world_grid.j_res - 1, world_grid.i_start, world_grid.j_start, world_grid.x_min + (T)0.5*world_grid.dx, world_grid.y_min + (T)0.5*world_grid.dy, world_grid.x_max + (T)0.5*world_grid.dx, world_grid.y_max + (T)0.5*world_grid.dy);
		ghost_grid.Initialize(base_grid.Enlarged(ghost_width));

		// Multithreading
		multithreading = multithreading_input;
		base_grid.SplitInYDirection(multithreading->num_threads, partial_base_grids);
		ghost_grid.SplitInYDirection(multithreading->num_threads, partial_ghost_grids);

		// Integral Value of each cell
		integral_value_array.Initialize(base_grid.i_res - 1, base_grid.j_res - 1, base_grid.i_start, base_grid.j_start, base_grid.x_min + (T)0.5*base_grid.dx, base_grid.y_min + (T)0.5*base_grid.dy, base_grid.x_max - (T)0.5*base_grid.dx, base_grid.y_max - (T)0.5*base_grid.dy, 2, true, false, multithreading_input);
		
		// For debug
		int ij_res = base_grid.i_res*base_grid.j_res;
		inserted_point.Initialize(ij_res);
		
		// Object Levelset
		DELETE_POINTER(object_levelset);
		object_levelset = new LEVELSET_2D();
		object_levelset->Initialize(world_grid, 2, multithreading);
		
		// Initialize object levelset
		object_levelset->AssignAllValuesLevelset(world_grid.dx*(T)3);
		object_levelset->FillGhostCellsFromPointer(&(object_levelset->phi), false);
		object_levelset->curvature_by_normal_vector = numerical_integration_block.GetBoolean("curvature_by_normal_vector", false);
		object_levelset->curvature_by_levelset = numerical_integration_block.GetBoolean("curvature_by_levelset", false);
	}

public: // Calculation of Numerical Integration
	void IntegrandFunctionValue(const VT& input_value, T& integrand_value, const int& thread_id)
	{
		integrand_value = (T)1;
	}

	void SurfaceIntegralThreaded(const int& thread_id)
	{
		BEGIN_HEAD_THREAD_WORK
		{
			integral_value = (T)0;
		}
		END_HEAD_THREAD_WORK;

		BEGIN_GRID_ITERATION_2D(partial_base_grids[thread_id])
		{
			// Levelset Value on the cell
			T phi[4];

			phi[0] = (*object_levelset)(i    , j    );
			phi[1] = (*object_levelset)(i + 1, j    );
			phi[2] = (*object_levelset)(i + 1, j + 1);
			phi[3] = (*object_levelset)(i    , j + 1);

			T epsilon = 1e-20;
			
			for (int k = 0; k < 4; k++)
			{
				if (abs(phi[k]) < (double)1e-20 && phi[k] > (T)0)
				{
					phi[k] = epsilon;
				}
				else if (abs(phi[k]) < (double)1e-20 && phi[k] <= (T)0)
				{
					phi[k] = -epsilon;
				}
			}

			// Speed up variables
			T x_min = (*object_levelset).signed_distance_field.x_min, y_min = (*object_levelset).signed_distance_field.y_min;
			T dx = (*object_levelset).signed_distance_field.dx, dy = (*object_levelset).signed_distance_field.dy;
			T x0 = x_min + i*dx, y0 = y_min + j*dy, x1 = x_min + (i + 1)*dx, y1 = y_min + (j + 1)*dy;

			VT v0 = VT(x0, y0), v1 = VT(x1, y0), v2 = VT(x1, y1), v3 = VT(x0, y1);
			if (phi[0] > (T)0 && phi[1] > (T)0 && phi[2] > (T)0 && phi[3] > (T)0)
			{
				integral_value_array(i, j) = (T)0;
			}
			else if (phi[0] < (T)0 && phi[1] > (T)0 && phi[2] > (T)0 && phi[3] > (T)0)				// only one point inside of levelset
			{
				VT v01 = v1*phi[0]/(phi[0] - phi[1]) + v0*phi[1]/(phi[1] - phi[0]);
				VT v02 = v2*phi[0]/(phi[0] - phi[2]) + v0*phi[2]/(phi[2] - phi[0]);
				VT v03 = v3*phi[0]/(phi[0] - phi[3]) + v0*phi[3]/(phi[3] - phi[0]);

				T f1(0), f2(0), f3(0);

				IntegrandFunctionValue(v01, f1, thread_id);
				IntegrandFunctionValue(v02, f2, thread_id);
				IntegrandFunctionValue(v03, f3, thread_id);

				VT l1 = v02 - v01, l2 = v03 - v02;
				
				T mag_l1 = l1.Magnitude(), mag_l2 = l2.Magnitude();

				integral_value_array(i, j) = mag_l1*(f1 + f2)*(T)0.5 + mag_l2*(f2 + f3)*(T)0.5;
			}
			else if (phi[0] > (T)0 && phi[1] < (T)0 && phi[2] > (T)0 && phi[3] > (T)0)
			{
				VT v10 = v0*phi[1]/(phi[1] - phi[0]) + v1*phi[0]/(phi[0] - phi[1]);
				VT v13 = v3*phi[1]/(phi[1] - phi[3]) + v1*phi[3]/(phi[3] - phi[1]);
				VT v12 = v2*phi[1]/(phi[1] - phi[2]) + v1*phi[2]/(phi[2] - phi[1]);

				T f1(0), f2(0), f3(0);

				IntegrandFunctionValue(v10, f1, thread_id);
				IntegrandFunctionValue(v13, f2, thread_id);
				IntegrandFunctionValue(v12, f3, thread_id);

				VT l1 = v10 - v13, l2 = v13 - v12;

				T mag_l1 = l1.Magnitude(), mag_l2 = l2.Magnitude();

				integral_value_array(i, j) = mag_l1*(f1 + f2)*(T)0.5 + mag_l2*(f2 + f3)*(T)0.5;
			}
			else if (phi[0] > (T)0 && phi[1] > (T)0 && phi[2] < (T)0 && phi[3] > (T)0)
			{
				VT v21 = v1*phi[2]/(phi[2] - phi[1]) + v2*phi[1]/(phi[1] - phi[2]);
				VT v20 = v0*phi[2]/(phi[2] - phi[0]) + v2*phi[0]/(phi[0] - phi[2]);
				VT v23 = v3*phi[2]/(phi[2] - phi[3]) + v2*phi[3]/(phi[3] - phi[2]);

				T f1(0), f2(0), f3(0);

				IntegrandFunctionValue(v21, f1, thread_id);
				IntegrandFunctionValue(v20, f2, thread_id);
				IntegrandFunctionValue(v23, f3, thread_id);

				VT l1 = v21 - v20, l2 = v23 - v20;

				T mag_l1 = l1.Magnitude(), mag_l2 = l2.Magnitude();

				integral_value_array(i, j) = mag_l1*(f1 + f2)*(T)0.5 + mag_l2*(f2 + f3)*(T)0.5;
			}
			else if (phi[0] > (T)0 && phi[1] > (T)0 && phi[2] > (T)0 && phi[3] < (T)0)
			{
				VT v30 = v0*phi[3]/(phi[3] - phi[0]) + v3*phi[0]/(phi[0] - phi[3]);
				VT v31 = v1*phi[3]/(phi[3] - phi[1]) + v3*phi[1]/(phi[1] - phi[3]);
				VT v32 = v2*phi[3]/(phi[3] - phi[2]) + v3*phi[2]/(phi[2] - phi[3]);

				T f1(0), f2(0), f3(0);

				IntegrandFunctionValue(v30, f1, thread_id);
				IntegrandFunctionValue(v31, f2, thread_id);
				IntegrandFunctionValue(v32, f3, thread_id);

				VT l1 = v30 - v31, l2 = v32 - v31;

				T mag_l1 = l1.Magnitude(), mag_l2 = l2.Magnitude();

				integral_value_array(i, j) = mag_l1*(f1 + f2)*(T)0.5 + mag_l2*(f2 + f3)*(T)0.5;
			}
			else if (phi[0] < (T)0 && phi[1] < (T)0 && phi[2] > (T)0 && phi[3] > (T)0)				// Two points inside levelset	
			{
				VT v03 = v3*phi[0]/(phi[0] - phi[3]) + v0*phi[3]/(phi[3] - phi[0]);
				VT v02 = v2*phi[0]/(phi[0] - phi[2]) + v0*phi[2]/(phi[2] - phi[0]);
				VT v12 = v2*phi[1]/(phi[1] - phi[2]) + v1*phi[2]/(phi[2] - phi[1]);

				T f1(0), f2(0), f3(0);

				IntegrandFunctionValue(v03, f1, thread_id);
				IntegrandFunctionValue(v02, f2, thread_id);
				IntegrandFunctionValue(v12, f3, thread_id);

				VT l1 = v03 - v02, l2 = v02 - v12;

				T mag_l1 = l1.Magnitude(), mag_l2 = l2.Magnitude();

				integral_value_array(i, j) = mag_l1*(f1 + f2)*(T)0.5 + mag_l2*(f2 + f3)*(T)0.5;
			}
			else if (phi[0] > (T)0 && phi[1] < (T)0 && phi[2] < (T)0 && phi[3] > (T)0)
			{
				VT v23 = v3*phi[2]/(phi[2] - phi[3]) + v2*phi[3]/(phi[3] - phi[2]);
				VT v20 = v0*phi[2]/(phi[2] - phi[0]) + v2*phi[0]/(phi[0] - phi[2]);
				VT v01 = v1*phi[0]/(phi[0] - phi[1]) + v0*phi[1]/(phi[1] - phi[0]);

				T f1(0), f2(0), f3(0);

				IntegrandFunctionValue(v23, f1, thread_id);
				IntegrandFunctionValue(v20, f2, thread_id);
				IntegrandFunctionValue(v01, f3, thread_id);

				VT l1 = v23 - v20, l2 = v20 - v01;

				T mag_l1 = l1.Magnitude(), mag_l2 = l2.Magnitude();

				integral_value_array(i, j) = mag_l1*(f1 + f2)*(T)0.5 + mag_l2*(f2 + f3)*(T)0.5;
			}
			else if (phi[0] > (T)0 && phi[1] > (T)0 && phi[2] < (T)0 && phi[3] < (T)0)
			{
				VT v21 = v1*phi[2]/(phi[2] - phi[1]) + v2*phi[1]/(phi[1] - phi[2]);
				VT v20 = v0*phi[2]/(phi[2] - phi[0]) + v2*phi[0]/(phi[0] - phi[2]);
				VT v30 = v0*phi[3]/(phi[3] - phi[0]) + v3*phi[0]/(phi[0] - phi[3]);

				T f1(0), f2(0), f3(0);

				IntegrandFunctionValue(v21, f1, thread_id);
				IntegrandFunctionValue(v20, f2, thread_id);
				IntegrandFunctionValue(v30, f3, thread_id);

				VT l1 = v21 - v20, l2 = v20 - v30;

				T mag_l1 = l1.Magnitude(), mag_l2 = l2.Magnitude();

				integral_value_array(i, j) = mag_l1*(f1 + f2)*(T)0.5 + mag_l2*(f2 + f3)*(T)0.5;
			}
			else if (phi[0] < (T)0 && phi[1] > (T)0 && phi[2] > (T)0 && phi[3] < (T)0)
			{
				VT v23 = v3*phi[2]/(phi[2] - phi[3]) + v2*phi[3]/(phi[3] - phi[2]);
				VT v20 = v0*phi[2]/(phi[2] - phi[0]) + v2*phi[0]/(phi[0] - phi[2]);
				VT v01 = v1*phi[0]/(phi[0] - phi[1]) + v0*phi[1]/(phi[1] - phi[0]);

				T f1(0), f2(0), f3(0);

				IntegrandFunctionValue(v23, f1, thread_id);
				IntegrandFunctionValue(v20, f2, thread_id);
				IntegrandFunctionValue(v01, f3, thread_id);

				VT l1 = v23 - v20, l2 = v20 - v01;

				T mag_l1 = l1.Magnitude(), mag_l2 = l2.Magnitude();

				integral_value_array(i, j) = mag_l1*(f1 + f2)*(T)0.5 + mag_l2*(f2 + f3)*(T)0.5;
			}
			else if (phi[0] > (T)0 && phi[1] < (T)0 && phi[2] < (T)0 && phi[3] < (T)0)			// Three points inside levelset
			{
				VT v01 = v1*phi[0]/(phi[0] - phi[1]) + v0*phi[1]/(phi[1] - phi[0]);
				VT v02 = v2*phi[0]/(phi[0] - phi[2]) + v0*phi[2]/(phi[2] - phi[0]);
				VT v03 = v3*phi[0]/(phi[0] - phi[3]) + v0*phi[3]/(phi[3] - phi[0]);

				T f1(0), f2(0), f3(0);

				IntegrandFunctionValue(v01, f1, thread_id);
				IntegrandFunctionValue(v02, f2, thread_id);
				IntegrandFunctionValue(v03, f3, thread_id);

				VT l1 = v02 - v01, l2 = v03 - v02;

				T mag_l1 = l1.Magnitude(), mag_l2 = l2.Magnitude();

				integral_value_array(i, j) = mag_l1*(f1 + f2)*(T)0.5 + mag_l2*(f2 + f3)*(T)0.5;
			}
			else if (phi[0] < (T)0 && phi[1] > (T)0 && phi[2] < (T)0 && phi[3] < (T)0)
			{
				VT v10 = v0*phi[1]/(phi[1] - phi[0]) + v1*phi[0]/(phi[0] - phi[1]);
				VT v13 = v3*phi[1]/(phi[1] - phi[3]) + v1*phi[3]/(phi[3] - phi[1]);
				VT v12 = v2*phi[1]/(phi[1] - phi[2]) + v1*phi[2]/(phi[2] - phi[1]);

				T f1(0), f2(0), f3(0);

				IntegrandFunctionValue(v10, f1, thread_id);
				IntegrandFunctionValue(v13, f2, thread_id);
				IntegrandFunctionValue(v12, f3, thread_id);

				VT l1 = v10 - v13, l2 = v13 - v12;

				T mag_l1 = l1.Magnitude(), mag_l2 = l2.Magnitude();

				integral_value_array(i, j) = mag_l1*(f1 + f2)*(T)0.5 + mag_l2*(f2 + f3)*(T)0.5;
			}
			else if (phi[0] < (T)0 && phi[1] < (T)0 && phi[2] > (T)0 && phi[3] < (T)0)
			{
				VT v21 = v1*phi[2]/(phi[2] - phi[1]) + v2*phi[1]/(phi[1] - phi[2]);
				VT v20 = v0*phi[2]/(phi[2] - phi[0]) + v2*phi[0]/(phi[0] - phi[2]);
				VT v23 = v3*phi[2]/(phi[2] - phi[3]) + v2*phi[3]/(phi[3] - phi[2]);

				T f1(0), f2(0), f3(0);

				IntegrandFunctionValue(v21, f1, thread_id);
				IntegrandFunctionValue(v20, f2, thread_id);
				IntegrandFunctionValue(v23, f3, thread_id);

				VT l1 = v21 - v20, l2 = v23 - v20;

				T mag_l1 = l1.Magnitude(), mag_l2 = l2.Magnitude();

				integral_value_array(i, j) = mag_l1*(f1 + f2)*(T)0.5 + mag_l2*(f2 + f3)*(T)0.5;
			}
			else if (phi[0] < (T)0 && phi[1] < (T)0 && phi[2] < (T)0 && phi[3] > (T)0)
			{
				VT v30 = v0*phi[3]/(phi[3] - phi[0]) + v3*phi[0]/(phi[0] - phi[3]);
				VT v31 = v1*phi[3]/(phi[3] - phi[1]) + v3*phi[1]/(phi[1] - phi[3]);
				VT v32 = v2*phi[3]/(phi[3] - phi[2]) + v3*phi[2]/(phi[2] - phi[3]);

				T f1(0), f2(0), f3(0);

				IntegrandFunctionValue(v30, f1, thread_id);
				IntegrandFunctionValue(v31, f2, thread_id);
				IntegrandFunctionValue(v32, f3, thread_id);

				VT l1 = v30 - v31, l2 = v32 - v31;

				T mag_l1 = l1.Magnitude(), mag_l2 = l2.Magnitude();

				integral_value_array(i, j) = mag_l1*(f1 + f2)*(T)0.5 + mag_l2*(f2 + f3)*(T)0.5;
			}
			else
			{
				integral_value_array(i, j) = 0;
			} 
			
			integral_value += integral_value_array(i, j);
		}
		END_GRID_ITERATION_2D;	
	}

	void SurfaceIntegralThreadedBJ(const int& thread_id)
	{
		BEGIN_HEAD_THREAD_WORK
		{
			integral_value = (T)0;
		}
		END_HEAD_THREAD_WORK;

		BEGIN_GRID_ITERATION_2D(partial_base_grids[thread_id])
		{
			// Levelset Value on the cell
			T phi[4];

			phi[0] = (*object_levelset)(i    , j    );
			phi[1] = (*object_levelset)(i + 1, j    );
			phi[2] = (*object_levelset)(i + 1, j + 1);
			phi[3] = (*object_levelset)(i    , j + 1);

			// Levelset Normal Vector on the cell
			VT normal_phi[4];

			normal_phi[0] = (*object_levelset).normal(i    , j    );
			normal_phi[1] = (*object_levelset).normal(i + 1, j    );
			normal_phi[2] = (*object_levelset).normal(i + 1, j + 1);
			normal_phi[3] = (*object_levelset).normal(i    , j + 1);

			T epsilon = 1e-20;
			
			for (int k = 0; k < 4; k++)
			{
				if (abs(phi[k]) < (double)1e-20 && phi[k] > (T)0)
				{
					phi[k] = epsilon;
				}
				else if (abs(phi[k]) < (double)1e-20 && phi[k] <= (T)0)
				{
					phi[k] = -epsilon;
				}
			}

			// Speed up variables
			T x_min = (*object_levelset).signed_distance_field.x_min, y_min = (*object_levelset).signed_distance_field.y_min;
			T dx = (*object_levelset).signed_distance_field.dx, dy = (*object_levelset).signed_distance_field.dy;
			T x0 = x_min + i*dx, y0 = y_min + j*dy, x1 = x_min + (i + 1)*dx, y1 = y_min + (j + 1)*dy;

			VT v0 = VT(x0, y0), v1 = VT(x1, y0), v2 = VT(x1, y1), v3 = VT(x0, y1);
			
			if (phi[0] > (T)0 && phi[1] > (T)0 && phi[2] > (T)0 && phi[3] > (T)0)
			{
				integral_value_array(i, j) = (T)0;
			}
			else if (phi[0] < (T)0 && phi[1] > (T)0 && phi[2] > (T)0 && phi[3] > (T)0)				// only one point inside of levelset
			{
				VT v01 = v1*phi[0]/(phi[0] - phi[1]) + v0*phi[1]/(phi[1] - phi[0]);
				VT v03 = v3*phi[0]/(phi[0] - phi[3]) + v0*phi[3]/(phi[3] - phi[0]);

				/*VT n01 = normal_phi[1]*abs(phi[0])/(abs(phi[0]) + abs(phi[1])) + normal_phi[0]*abs(phi[1])/(abs(phi[0]) + abs(phi[1]));
				VT n03 = normal_phi[3]*abs(phi[0])/(abs(phi[0]) + abs(phi[3])) + normal_phi[0]*abs(phi[3])/(abs(phi[0]) + abs(phi[3]));*/
				VT n01 = normal_phi[1]*v0.x/(v0.x - v1.x) + normal_phi[0]*v1.x/(v1.x - v0.x);
				VT n03 = normal_phi[3]*v0.y/(v0.y - v3.y) + normal_phi[0]*v3.y/(v3.y - v0.y);
				
				MATRIX_2X2 NN(n01.y, n03.y, -n01.x, -n03.x);
				VT RHS(n01.y*v01.x - n01.x*v01.y, n03.y*v03.x - n03.x*v03.y);

				VT v_intersection = NN.Inversed()*RHS;
				T phi_intersection = (*object_levelset)(v_intersection);

				VT v02;

				if (abs(NN.Determinant()) < (T)1e-20)
				{
					v02 = v2*phi[0]/(phi[0] - phi[2]) + v0*phi[2]/(phi[2] - phi[0]);
				}
				else
				{
					v02 = v_intersection*phi[0]/(phi[0] - phi_intersection) + v0*phi_intersection/(phi_intersection - phi[0]);
				}
				
				inserted_point.Push(v02);
				
				T f1(0), f2(0), f3(0);

				IntegrandFunctionValue(v01, f1, thread_id);
				IntegrandFunctionValue(v02, f2, thread_id);
				IntegrandFunctionValue(v03, f3, thread_id);

				VT l1 = v02 - v01, l2 = v03 - v02;
				
				T mag_l1 = l1.Magnitude(), mag_l2 = l2.Magnitude();

				integral_value_array(i, j) = mag_l1*(f1 + f2)*(T)0.5 + mag_l2*(f2 + f3)*(T)0.5;
			}
			else if (phi[0] > (T)0 && phi[1] < (T)0 && phi[2] > (T)0 && phi[3] > (T)0)
			{
				VT v10 = v0*phi[1]/(phi[1] - phi[0]) + v1*phi[0]/(phi[0] - phi[1]);
				VT v12 = v2*phi[1]/(phi[1] - phi[2]) + v1*phi[2]/(phi[2] - phi[1]);

				/*VT n10 = normal_phi[1]*abs(phi[0])/(abs(phi[0]) + abs(phi[1])) + normal_phi[0]*abs(phi[1])/(abs(phi[0]) + abs(phi[1]));
				VT n12 = normal_phi[2]*abs(phi[1])/(abs(phi[1]) + abs(phi[2])) + normal_phi[1]*abs(phi[2])/(abs(phi[1]) + abs(phi[2]));*/
				VT n10 = normal_phi[1]*v0.x/(v0.x - v1.x) + normal_phi[0]*v1.x/(v1.x - v0.x);
				VT n12 = normal_phi[2]*v1.y/(v1.y - v2.y) + normal_phi[1]*v2.y/(v2.y - v1.y);

				MATRIX_2X2 NN(n10.y, n12.y, -n10.x, -n12.x);
				VT RHS(n10.y*v10.x - n10.x*v10.y, n12.y*v12.x - n12.x*v12.y);

				VT v_intersection = NN.Inversed()*RHS;
				T phi_intersection = (*object_levelset)(v_intersection);

				VT v13;

				if (abs(NN.Determinant()) < (T)1e-20)
				{
					v13 = v3*phi[1]/(phi[1] - phi[3]) + v1*phi[3]/(phi[3] - phi[1]);
				}
				else
				{
					v13 = v_intersection*phi[1]/(phi[1] - phi_intersection) + v1*phi_intersection/(phi_intersection - phi[1]);
				}
				
				inserted_point.Push(v13);

				T f1(0), f2(0), f3(0);

				IntegrandFunctionValue(v10, f1, thread_id);
				IntegrandFunctionValue(v13, f2, thread_id);
				IntegrandFunctionValue(v12, f3, thread_id);

				VT l1 = v10 - v13, l2 = v13 - v12;

				T mag_l1 = l1.Magnitude(), mag_l2 = l2.Magnitude();

				integral_value_array(i, j) = mag_l1*(f1 + f2)*(T)0.5 + mag_l2*(f2 + f3)*(T)0.5;
			}
			else if (phi[0] > (T)0 && phi[1] > (T)0 && phi[2] < (T)0 && phi[3] > (T)0)
			{
				VT v21 = v1*phi[2]/(phi[2] - phi[1]) + v2*phi[1]/(phi[1] - phi[2]);
				VT v23 = v3*phi[2]/(phi[2] - phi[3]) + v2*phi[3]/(phi[3] - phi[2]);

				T f1(0), f2(0), f3(0);

				/*VT n23 = normal_phi[2]*abs(phi[3])/(abs(phi[2]) + abs(phi[3])) + normal_phi[3]*abs(phi[2])/(abs(phi[2]) + abs(phi[3]));
				VT n21 = normal_phi[2]*abs(phi[1])/(abs(phi[1]) + abs(phi[2])) + normal_phi[1]*abs(phi[2])/(abs(phi[1]) + abs(phi[2]));*/
				VT n23 = normal_phi[2]*v3.x/(v3.x - v2.x) + normal_phi[3]*v2.x/(v2.x - v3.x);
				VT n21 = normal_phi[2]*v1.y/(v1.y - v2.y) + normal_phi[1]*v2.y/(v2.y - v1.y);

				MATRIX_2X2 NN(n23.y, n21.y, -n23.x, -n21.x);
				VT RHS(n23.y*v23.x - n23.x*v23.y, n21.y*v21.x - n21.x*v21.y);

				VT v_intersection = NN.Inversed()*RHS;
				T phi_intersection = (*object_levelset)(v_intersection);

				VT v20;

				if (abs(NN.Determinant()) < (T)1e-20)
				{
					v20 = v0*phi[2]/(phi[2] - phi[0]) + v2*phi[0]/(phi[0] - phi[2]);
				}
				else
				{
					v20 = v_intersection*phi[2]/(phi[2] - phi_intersection) + v2*phi_intersection/(phi_intersection - phi[2]);
				}
				
				inserted_point.Push(v20);

				IntegrandFunctionValue(v21, f1, thread_id);
				IntegrandFunctionValue(v20, f2, thread_id);
				IntegrandFunctionValue(v23, f3, thread_id);

				VT l1 = v21 - v20, l2 = v23 - v20;

				T mag_l1 = l1.Magnitude(), mag_l2 = l2.Magnitude();

				integral_value_array(i, j) = mag_l1*(f1 + f2)*(T)0.5 + mag_l2*(f2 + f3)*(T)0.5;
			}
			else if (phi[0] > (T)0 && phi[1] > (T)0 && phi[2] > (T)0 && phi[3] < (T)0)
			{
				VT v30 = v0*phi[3]/(phi[3] - phi[0]) + v3*phi[0]/(phi[0] - phi[3]);
				VT v32 = v2*phi[3]/(phi[3] - phi[2]) + v3*phi[2]/(phi[2] - phi[3]);

				VT n32 = normal_phi[3]*v2.x/(v2.x - v3.x) + normal_phi[2]*v3.x/(v3.x - v2.x);
				VT n30 = normal_phi[3]*v0.y/(v0.y - v3.y) + normal_phi[0]*v3.y/(v3.y - v0.y);
				/*VT n32 = normal_phi[3]*abs(phi[2])/(abs(phi[2]) + abs(phi[3])) + normal_phi[2]*abs(phi[3])/(abs(phi[2]) + abs(phi[3]));
				VT n30 = normal_phi[3]*abs(phi[0])/(abs(phi[0]) + abs(phi[3])) + normal_phi[0]*abs(phi[3])/(abs(phi[0]) + abs(phi[3]));*/

				MATRIX_2X2 NN(n32.y, n30.y, -n32.x, -n30.x);
				VT RHS(n32.y*v32.x - n32.x*v32.y, n30.y*v30.x - n30.x*v30.y);

				VT v_intersection = NN.Inversed()*RHS;
				T phi_intersection = (*object_levelset)(v_intersection);

				VT v31;

				if (abs(NN.Determinant()) < (T)1e-20)
				{
					v31 = v1*phi[3]/(phi[3] - phi[1]) + v3*phi[1]/(phi[1] - phi[3]);
				}
				else
				{
					v31 = v_intersection*phi[3]/(phi[3] - phi_intersection) + v3*phi_intersection/(phi_intersection - phi[3]);
				}
				
				inserted_point.Push(v31);

				T f1(0), f2(0), f3(0);

				IntegrandFunctionValue(v30, f1, thread_id);
				IntegrandFunctionValue(v31, f2, thread_id);
				IntegrandFunctionValue(v32, f3, thread_id);

				VT l1 = v30 - v31, l2 = v32 - v31;

				T mag_l1 = l1.Magnitude(), mag_l2 = l2.Magnitude();

				integral_value_array(i, j) = mag_l1*(f1 + f2)*(T)0.5 + mag_l2*(f2 + f3)*(T)0.5;
			}
			else if (phi[0] < (T)0 && phi[1] < (T)0 && phi[2] > (T)0 && phi[3] > (T)0)				// Two points inside levelset	
			{
				VT v03 = v3*phi[0]/(phi[0] - phi[3]) + v0*phi[3]/(phi[3] - phi[0]);
				VT v02 = v2*phi[0]/(phi[0] - phi[2]) + v0*phi[2]/(phi[2] - phi[0]);
				VT v12 = v2*phi[1]/(phi[1] - phi[2]) + v1*phi[2]/(phi[2] - phi[1]);

				/*VT n03 = normal_phi[0]*abs(phi[3])/(abs(phi[0]) + abs(phi[3])) + normal_phi[3]*abs(phi[0])/(abs(phi[0]) + abs(phi[3]));
				VT n12 = normal_phi[1]*abs(phi[2])/(abs(phi[1]) + abs(phi[2])) + normal_phi[2]*abs(phi[1])/(abs(phi[1]) + abs(phi[2]));

				MATRIX_2X2 NN(n03.y, n12.y, -n03.x, -n12.x);
				VT RHS(n03.y*v03.x - n03.x*v03.y, n12.y*v12.x - n12.x*v12.y);

				VT v_intersection = NN.Inversed()*RHS;
				T phi_intersection = (*object_levelset)(v_intersection);

				VT v02;

				if (abs(NN.Determinant()) < (T)1e-20)
				{
					v02 = v2*phi[0]/(phi[0] - phi[2]) + v0*phi[2]/(phi[2] - phi[0]);
				}
				else
				{
					v02 = v_intersection*phi[0]/(phi[0] - phi_intersection) + v0*phi_intersection/(phi_intersection - phi[0]);
				}*/

				inserted_point.Push(v02);

				T f1(0), f2(0), f3(0);
				
				IntegrandFunctionValue(v03, f1, thread_id);
				IntegrandFunctionValue(v02, f2, thread_id);
				IntegrandFunctionValue(v12, f3, thread_id);
				
				VT l1 = v03 - v02, l2 = v02 - v12;

				T mag_l1 = l1.Magnitude(), mag_l2 = l2.Magnitude();

				integral_value_array(i, j) = mag_l1*(f1 + f2)*(T)0.5 + mag_l2*(f2 + f3)*(T)0.5;
			}
			else if (phi[0] > (T)0 && phi[1] < (T)0 && phi[2] < (T)0 && phi[3] > (T)0)
			{
				VT v23 = v3*phi[2]/(phi[2] - phi[3]) + v2*phi[3]/(phi[3] - phi[2]);
				VT v20 = v0*phi[2]/(phi[2] - phi[0]) + v2*phi[0]/(phi[0] - phi[2]);
				VT v01 = v1*phi[0]/(phi[0] - phi[1]) + v0*phi[1]/(phi[1] - phi[0]);

				/*VT n23 = normal_phi[2]*abs(phi[3])/(abs(phi[2]) + abs(phi[3])) + normal_phi[3]*abs(phi[2])/(abs(phi[2]) + abs(phi[3]));
				VT n01 = normal_phi[0]*abs(phi[1])/(abs(phi[0]) + abs(phi[1])) + normal_phi[1]*abs(phi[0])/(abs(phi[0]) + abs(phi[1]));

				MATRIX_2X2 NN(n23.y, n01.y, -n23.x, -n01.x);
				VT RHS(n23.y*v23.x - n23.x*v23.y, n01.y*v01.x - n01.x*v01.y);

				VT v_intersection = NN.Inversed()*RHS;
				T phi_intersection = (*object_levelset)(v_intersection);

				VT v20;

				if (abs(NN.Determinant()) < (T)1e-20)
				{
					v20 = v2*phi[0]/(phi[0] - phi[2]) + v0*phi[2]/(phi[2] - phi[0]);
				}
				else
				{
					v20 = v_intersection*phi[0]/(phi[0] - phi_intersection) + v0*phi_intersection/(phi_intersection - phi[0]);
				}*/

				inserted_point.Push(v20);

				T f1(0), f2(0), f3(0);

				IntegrandFunctionValue(v23, f1, thread_id);
				IntegrandFunctionValue(v20, f2, thread_id);
				IntegrandFunctionValue(v01, f3, thread_id);

				VT l1 = v23 - v20, l2 = v20 - v01;

				T mag_l1 = l1.Magnitude(), mag_l2 = l2.Magnitude();

				integral_value_array(i, j) = mag_l1*(f1 + f2)*(T)0.5 + mag_l2*(f2 + f3)*(T)0.5;
			}
			else if (phi[0] > (T)0 && phi[1] > (T)0 && phi[2] < (T)0 && phi[3] < (T)0)
			{
				VT v21 = v1*phi[2]/(phi[2] - phi[1]) + v2*phi[1]/(phi[1] - phi[2]);
				VT v20 = v0*phi[2]/(phi[2] - phi[0]) + v2*phi[0]/(phi[0] - phi[2]);
				VT v30 = v0*phi[3]/(phi[3] - phi[0]) + v3*phi[0]/(phi[0] - phi[3]);

				//inserted_point.Push(v20);

				T f1(0), f2(0), f3(0);

				IntegrandFunctionValue(v21, f1, thread_id);
				IntegrandFunctionValue(v20, f2, thread_id);
				IntegrandFunctionValue(v30, f3, thread_id);

				VT l1 = v21 - v20, l2 = v20 - v30;

				T mag_l1 = l1.Magnitude(), mag_l2 = l2.Magnitude();

				integral_value_array(i, j) = mag_l1*(f1 + f2)*(T)0.5 + mag_l2*(f2 + f3)*(T)0.5;
			}
			else if (phi[0] < (T)0 && phi[1] > (T)0 && phi[2] > (T)0 && phi[3] < (T)0)
			{
				VT v23 = v3*phi[2]/(phi[2] - phi[3]) + v2*phi[3]/(phi[3] - phi[2]);
				VT v20 = v0*phi[2]/(phi[2] - phi[0]) + v2*phi[0]/(phi[0] - phi[2]);
				VT v01 = v1*phi[0]/(phi[0] - phi[1]) + v0*phi[1]/(phi[1] - phi[0]);

				//inserted_point.Push(v20);

				T f1(0), f2(0), f3(0);

				IntegrandFunctionValue(v23, f1, thread_id);
				IntegrandFunctionValue(v20, f2, thread_id);
				IntegrandFunctionValue(v01, f3, thread_id);

				VT l1 = v23 - v20, l2 = v20 - v01;

				T mag_l1 = l1.Magnitude(), mag_l2 = l2.Magnitude();

				integral_value_array(i, j) = mag_l1*(f1 + f2)*(T)0.5 + mag_l2*(f2 + f3)*(T)0.5;
			}
			else if (phi[0] > (T)0 && phi[1] < (T)0 && phi[2] < (T)0 && phi[3] < (T)0)			// Three points inside levelset
			{
				VT v01 = v1*phi[0]/(phi[0] - phi[1]) + v0*phi[1]/(phi[1] - phi[0]);
				VT v02 = v2*phi[0]/(phi[0] - phi[2]) + v0*phi[2]/(phi[2] - phi[0]);
				VT v03 = v3*phi[0]/(phi[0] - phi[3]) + v0*phi[3]/(phi[3] - phi[0]);

				/*VT n01 = normal_phi[1]*abs(phi[0])/(abs(phi[0]) + abs(phi[1])) + normal_phi[0]*abs(phi[1])/(abs(phi[0]) + abs(phi[1]));
				VT n03 = normal_phi[3]*abs(phi[0])/(abs(phi[0]) + abs(phi[3])) + normal_phi[0]*abs(phi[3])/(abs(phi[0]) + abs(phi[3]));

				MATRIX_2X2 NN(n01.y, n03.y, -n01.x, -n03.x);
				VT RHS(n01.y*v01.x - n01.x*v01.y, n03.y*v03.x - n03.x*v03.y);

				VT v_intersection = NN.Inversed()*RHS;
				T phi_intersection = (*object_levelset)(v_intersection);

				VT v02;

				if (abs(NN.Determinant()) < (T)1e-20)
				{
					v02 = v2*phi[0]/(phi[0] - phi[2]) + v0*phi[2]/(phi[2] - phi[0]);
				}
				else
				{
					if (phi_intersection > (T)0)
					{
						v02 = v_intersection*phi[2]/(phi[2] - phi_intersection) + v2*phi_intersection/(phi_intersection - phi[2]);
					}
					else
					{
						v02 = v_intersection*phi[0]/(phi[0] - phi_intersection) + v0*phi_intersection/(phi_intersection - phi[0]);
					}
				}
				
				inserted_point.Push(v02);*/

				T f1(0), f2(0), f3(0);

				IntegrandFunctionValue(v01, f1, thread_id);
				IntegrandFunctionValue(v02, f2, thread_id);
				IntegrandFunctionValue(v03, f3, thread_id);

				VT l1 = v02 - v01, l2 = v03 - v02;
				
				T mag_l1 = l1.Magnitude(), mag_l2 = l2.Magnitude();

				integral_value_array(i, j) = mag_l1*(f1 + f2)*(T)0.5 + mag_l2*(f2 + f3)*(T)0.5;
			}
			else if (phi[0] < (T)0 && phi[1] > (T)0 && phi[2] < (T)0 && phi[3] < (T)0)
			{
				VT v10 = v0*phi[1]/(phi[1] - phi[0]) + v1*phi[0]/(phi[0] - phi[1]);
				VT v13 = v3*phi[1]/(phi[1] - phi[3]) + v1*phi[3]/(phi[3] - phi[1]);
				VT v12 = v2*phi[1]/(phi[1] - phi[2]) + v1*phi[2]/(phi[2] - phi[1]);

				/*VT n10 = normal_phi[1]*abs(phi[0])/(abs(phi[0]) + abs(phi[1])) + normal_phi[0]*abs(phi[1])/(abs(phi[0]) + abs(phi[1]));
				VT n12 = normal_phi[1]*abs(phi[2])/(abs(phi[1]) + abs(phi[2])) + normal_phi[2]*abs(phi[1])/(abs(phi[1]) + abs(phi[2]));

				MATRIX_2X2 NN(n10.y, n12.y, -n10.x, -n12.x);
				VT RHS(n10.y*v10.x - n10.x*v10.y, n12.y*v12.x - n12.x*v12.y);

				VT v_intersection = NN.Inversed()*RHS;
				T phi_intersection = (*object_levelset)(v_intersection);

				VT v13;

				if (abs(NN.Determinant()) < (T)1e-20)
				{
					v13 = v3*phi[1]/(phi[1] - phi[3]) + v1*phi[3]/(phi[3] - phi[1]);
				}
				else
				{
					if (phi_intersection > (T)0)
					{
						v13 = v3*phi_intersection/(phi_intersection - phi[3]) + v_intersection*phi[3]/(phi[3] - phi_intersection);
					}
					else
					{
						v13 = v_intersection*phi[1]/(phi[1] - phi_intersection) + v1*phi_intersection/(phi_intersection - phi[1]);
					}
				}
				
				inserted_point.Push(v13);*/

				T f1(0), f2(0), f3(0);

				IntegrandFunctionValue(v10, f1, thread_id);
				IntegrandFunctionValue(v13, f2, thread_id);
				IntegrandFunctionValue(v12, f3, thread_id);

				VT l1 = v10 - v13, l2 = v13 - v12;

				T mag_l1 = l1.Magnitude(), mag_l2 = l2.Magnitude();

				integral_value_array(i, j) = mag_l1*(f1 + f2)*(T)0.5 + mag_l2*(f2 + f3)*(T)0.5;
			}
			else if (phi[0] < (T)0 && phi[1] < (T)0 && phi[2] > (T)0 && phi[3] < (T)0)
			{
				VT v21 = v1*phi[2]/(phi[2] - phi[1]) + v2*phi[1]/(phi[1] - phi[2]);
				VT v20 = v0*phi[2]/(phi[2] - phi[0]) + v2*phi[0]/(phi[0] - phi[2]);
				VT v23 = v3*phi[2]/(phi[2] - phi[3]) + v2*phi[3]/(phi[3] - phi[2]);

				//inserted_point.Push(v20);

				/*VT n21 = normal_phi[2]*abs(phi[1])/(abs(phi[1]) + abs(phi[2])) + normal_phi[1]*abs(phi[2])/(abs(phi[1]) + abs(phi[2]));
				VT n23 = normal_phi[2]*abs(phi[3])/(abs(phi[2]) + abs(phi[3])) + normal_phi[3]*abs(phi[2])/(abs(phi[2]) + abs(phi[3]));

				MATRIX_2X2 NN(n21.y, n23.y, -n21.x, -n23.x);
				VT RHS(n21.y*v21.x - n21.x*v21.y, n23.y*v23.x - n23.x*v23.y);

				VT v_intersection = NN.Inversed()*RHS;
				T phi_intersection = (*object_levelset)(v_intersection);

				VT v20;

				if (abs(NN.Determinant()) < (T)1e-20)
				{
					v20 = v0*phi[2]/(phi[2] - phi[0]) + v2*phi[0]/(phi[0] - phi[2]);
				}
				else
				{
					if (phi_intersection > (T)0)
					{
						v20 = v0*phi_intersection/(phi_intersection - phi[0]) + v_intersection*phi[0]/(phi[0] - phi_intersection);
					}
					else
					{
						v20 = v_intersection*phi[2]/(phi[2] - phi_intersection) + v2*phi_intersection/(phi_intersection - phi[2]);
					}
				}
				
				inserted_point.Push(v20);*/

				T f1(0), f2(0), f3(0);

				IntegrandFunctionValue(v21, f1, thread_id);
				IntegrandFunctionValue(v20, f2, thread_id);
				IntegrandFunctionValue(v23, f3, thread_id);

				VT l1 = v21 - v20, l2 = v23 - v20;

				T mag_l1 = l1.Magnitude(), mag_l2 = l2.Magnitude();

				integral_value_array(i, j) = mag_l1*(f1 + f2)*(T)0.5 + mag_l2*(f2 + f3)*(T)0.5;
			}
			else if (phi[0] < (T)0 && phi[1] < (T)0 && phi[2] < (T)0 && phi[3] > (T)0)
			{
				VT v30 = v0*phi[3]/(phi[3] - phi[0]) + v3*phi[0]/(phi[0] - phi[3]);
				VT v31 = v1*phi[3]/(phi[3] - phi[1]) + v3*phi[1]/(phi[1] - phi[3]);
				VT v32 = v2*phi[3]/(phi[3] - phi[2]) + v3*phi[2]/(phi[2] - phi[3]);

				/*VT n30 = normal_phi[3]*abs(phi[0])/(abs(phi[0]) + abs(phi[3])) + normal_phi[0]*abs(phi[3])/(abs(phi[0]) + abs(phi[3]));
				VT n32 = normal_phi[3]*abs(phi[2])/(abs(phi[2]) + abs(phi[3])) + normal_phi[2]*abs(phi[3])/(abs(phi[2]) + abs(phi[3]));

				MATRIX_2X2 NN(n30.y, n32.y, -n30.x, -n32.x);
				VT RHS(n30.y*v30.x - n30.x*v30.y, n32.y*v32.x - n32.x*v32.y);

				VT v_intersection = NN.Inversed()*RHS;
				T phi_intersection = (*object_levelset)(v_intersection);

				VT v31;

				if (abs(NN.Determinant()) < (T)1e-20)
				{
					v31 = v1*phi[3]/(phi[3] - phi[1]) + v3*phi[1]/(phi[1] - phi[3]);
				}
				else
				{
					if (phi_intersection > (T)0)
					{
						v31 = v1*phi_intersection/(phi_intersection - phi[1]) + v_intersection*phi[1]/(phi[1] - phi_intersection);
					}
					else
					{
						v31 = v_intersection*phi[3]/(phi[3] - phi_intersection) + v3*phi_intersection/(phi_intersection - phi[3]);
					}
				}

				inserted_point.Push(v31);*/

				T f1(0), f2(0), f3(0);

				IntegrandFunctionValue(v30, f1, thread_id);
				IntegrandFunctionValue(v31, f2, thread_id);
				IntegrandFunctionValue(v32, f3, thread_id);

				VT l1 = v30 - v31, l2 = v32 - v31;

				T mag_l1 = l1.Magnitude(), mag_l2 = l2.Magnitude();

				integral_value_array(i, j) = mag_l1*(f1 + f2)*(T)0.5 + mag_l2*(f2 + f3)*(T)0.5;
			}
			else
			{
				integral_value_array(i, j) = 0;
			} 
			
			integral_value += integral_value_array(i, j);
		}
		END_GRID_ITERATION_2D;	
	}

	void CalculateRelativeError()
	{
		relative_error = abs((integral_value - true_integral_value)/true_integral_value);
	}	
};