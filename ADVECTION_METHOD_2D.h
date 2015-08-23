#pragma once

#include "LEVELSET_2D.h"
#include "MATRIX_2X2.h"

template<class TT>
class ADVECTION_METHOD_2D
{
public: 
	// Semi-Legrangian Method
	static void SL1stOrder(FIELD_STRUCTURE_2D<TT>& rho, FIELD_STRUCTURE_2D<TT>& rho_ghost, const FIELD_STRUCTURE_2D<VT>& velocity, const T& dt);
	static void SL1stOrder(LEVELSET_2D& rho, FIELD_STRUCTURE_2D<TT>& rho_ghost, const FIELD_STRUCTURE_2D<VT>& velocity, const T& dt);
	
	// Gradient Augmented_Method
	static void Gradient_Augmented_Method(LEVELSET_2D& rho, FIELD_STRUCTURE_2D<TT>& rho_ghost, const FIELD_STRUCTURE_2D<VT>& velocity, const T& dt)
	{
		const ARRAY_2D<VT>& velocity_array(velocity.array_for_this);
		ARRAY_2D<TT>& rho_array(rho.arr);

		rho_ghost.FillGhostCellsContinuousDerivativesFrom(rho_array, true);

		VT characteristic_point;
		
		T hermite_cubic, hermite_cubic_x, hermite_cubic_y;

		// For temproary value for phi_x, phi_y
		FIELD_STRUCTURE_2D<VT> gradient_new;
		gradient_new.Initialize(rho.grid, 2);
		
		// Boundary Condition - Inflow homogenious Neumann Boundary Condition
		int i_start = rho.grid.i_start, i_end = rho.grid.i_end;
		int j_start = rho.grid.j_start, j_end = rho.grid.j_end;

		// For i_start and i_end
		for (int j = rho.grid.j_start + 1; j <= rho.grid.j_end - 1; j++)
		{
			rho.gradient(i_start, j).x = 0;
			rho.gradient(i_end, j).x = 0;
		}

		for (int j = rho.grid.j_start + 1; j <= rho.grid.j_end - 1; j++)
		{
			// Define the deformation matrix
			T u_x = (velocity_array(i_end + 1, j).x - velocity_array(i_end - 1, j).x)*rho.grid.one_over_2dx;
			T u_y = (velocity_array(i_end, j + 1).x - velocity_array(i_end, j - 1).x)*rho.grid.one_over_2dy;
			T v_x = (velocity_array(i_end + 1, j).y - velocity_array(i_end - 1, j).y)*rho.grid.one_over_2dx;
			T v_y = (velocity_array(i_end, j + 1).y - velocity_array(i_end, j - 1).y)*rho.grid.one_over_2dy;

			MATRIX_2X2 deformation_matrix(u_x, u_y, v_x, v_y);

			T y_coordinate = rho.grid.y_min + j*rho.grid.dy;

			if (y_coordinate > 50)
			{
				VT step_1 = rho.GridPoint(i_end, j) - velocity_array(i_end, j)*dt;

				MATRIX_2X2 step_1_m = step_1_m.Identity() - dt*deformation_matrix;

				VT step_2 = rho.GridPoint(i_end, j) - dt*((T)0.25*velocity_array(i_end, j) + (T)0.25*VT(PI/314*((T)50 - step_1.y), PI/314*(step_1.x - (T)50)));
								
				MATRIX_2X2 step_2_m = step_2_m.Identity() - dt*((T)1/4*deformation_matrix + (T)1/4*step_1_m*deformation_matrix);

				characteristic_point = rho.GridPoint(i_end, j) - dt*((T)1/6*velocity_array(i_end, j) + (T)1/6*VT(PI/314*((T)50 - step_1.y), PI/314*(step_1.x - (T)50)) + (T)2/3*VT(PI/314*((T)50 - step_2.y), PI/314*(step_2.x - (T)50)));

				MATRIX_2X2 characteristic_matrix = characteristic_matrix.Identity() - dt*((T)1/6*deformation_matrix + (T)1/6*step_1_m*deformation_matrix + (T)2/3*step_2_m*deformation_matrix);

				rho_ghost.HermiteCubicInterpolation(characteristic_point, rho.gradient, hermite_cubic, hermite_cubic_x, hermite_cubic_y);

				rho_array(i_end, j) = hermite_cubic;

				gradient_new(i_end, j) = characteristic_matrix*VT(hermite_cubic_x, hermite_cubic_y);
			}

			if (y_coordinate < 50)
			{
				VT step_1 = rho.GridPoint(i_start, j) - velocity_array(i_start, j)*dt;

				MATRIX_2X2 step_1_m = step_1_m.Identity() - dt*deformation_matrix;

				VT step_2 = rho.GridPoint(i_start, j) - dt*((T)0.25*velocity_array(i_start, j) - (T)0.25*VT(PI/314*((T)50 - step_1.y), PI/314*(step_1.x - (T)50)));

				MATRIX_2X2 step_2_m = step_2_m.Identity() - dt*((T)1/4*deformation_matrix + (T)1/4*step_1_m*deformation_matrix);

				characteristic_point = rho.GridPoint(i_start, j) - dt*((T)1/6*velocity_array(i_start, j) + (T)1/6*VT(PI/314*((T)50 - step_1.y), PI/314*(step_1.x - (T)50)) + (T)2/3*VT(PI/314*((T)50 - step_2.y), PI/314*(step_2.x - (T)50)));

				MATRIX_2X2 characteristic_matrix = characteristic_matrix.Identity() - dt*((T)1/6*deformation_matrix + (T)1/6*step_1_m*deformation_matrix + (T)2/3*step_2_m*deformation_matrix);

				rho_ghost.HermiteCubicInterpolation(characteristic_point, rho.gradient, hermite_cubic, hermite_cubic_x, hermite_cubic_y);

				rho_array(i_start, j) = hermite_cubic;

				gradient_new(i_start, j) = characteristic_matrix*VT(hermite_cubic_x, hermite_cubic_y);
			}
		}

		for (int j = rho.grid.j_start + 1; j <= rho.grid.j_end - 1; j++)
		{
			rho.gradient(i_start, j) = gradient_new(i_start, j);
			rho.gradient(i_end, j) = gradient_new(i_end, j);
		}

		// For j_start and j_end
		for (int i = rho.grid.i_start + 1; i <= rho.grid.i_end - 1; i++)
		{
			rho.gradient(i, j_start).y = 0;
			rho.gradient(i, j_end).y = 0;
		}

		for (int i = rho.grid.i_start + 1; i <= rho.grid.i_end - 1; i++)
		{
			// Define the deformation matrix
			T u_x = (velocity_array(i + 1, j_end).x - velocity_array(i - 1, j_end).x)*rho.grid.one_over_2dx;
			T u_y = (velocity_array(i, j_end + 1).x - velocity_array(i, j_end - 1).x)*rho.grid.one_over_2dy;
			T v_x = (velocity_array(i + 1, j_end).y - velocity_array(i - 1, j_end).y)*rho.grid.one_over_2dx;
			T v_y = (velocity_array(i, j_end + 1).y - velocity_array(i, j_end - 1).y)*rho.grid.one_over_2dy;

			MATRIX_2X2 deformation_matrix(u_x, u_y, v_x, v_y);

			T x_coordinate = rho.grid.x_min + i*rho.grid.dx;

			if (x_coordinate > 50)
			{
				VT step_1 = rho.GridPoint(i, j_start) - VT(velocity_array(i, j_start).x, 0)*dt;

				MATRIX_2X2 step_1_m = step_1_m.Identity() - dt*deformation_matrix;

				VT step_2 = rho.GridPoint(i, j_start) - dt*((T)0.25*velocity_array(i, j_start) + (T)0.25*VT(PI/314*((T)50 - step_1.y), 0));

				MATRIX_2X2 step_2_m = step_2_m.Identity() - dt*((T)1/4*deformation_matrix + (T)1/4*step_1_m*deformation_matrix);

				characteristic_point = rho.GridPoint(i, j_start) - dt*((T)1/6*velocity_array(i, j_start) + (T)1/6*VT(PI/314*((T)50 - step_1.y), 0) + (T)2/3*VT(PI/314*((T)50 - step_2.y), 0));

				MATRIX_2X2 characteristic_matrix = characteristic_matrix.Identity() - dt*((T)1/6*deformation_matrix + (T)1/6*step_1_m*deformation_matrix + (T)2/3*step_2_m*deformation_matrix);

				rho_ghost.HermiteCubicInterpolation(characteristic_point, rho.gradient, hermite_cubic, hermite_cubic_x, hermite_cubic_y);

				rho_array(i, j_start) = hermite_cubic;

				gradient_new(i, j_start) = characteristic_matrix*VT(hermite_cubic_x, hermite_cubic_y);
			}
			if (x_coordinate < 50)
			{
				VT step_1 = rho.GridPoint(i, j_end) - VT(velocity_array(i, j_end).x, 0)*dt;

				MATRIX_2X2 step_1_m = step_1_m.Identity() - dt*deformation_matrix;

				VT step_2 = rho.GridPoint(i, j_end) - dt*((T)0.25*velocity_array(i, j_end) + (T)0.25*VT(PI/314*((T)50 - step_1.y), 0));

				MATRIX_2X2 step_2_m = step_2_m.Identity() - dt*((T)1/4*deformation_matrix + (T)1/4*step_1_m*deformation_matrix);

				characteristic_point = rho.GridPoint(i, j_end) - dt*((T)1/6*velocity_array(i, j_end) + (T)1/6*VT(PI/314*((T)50 - step_1.y), 0) + (T)2/3*VT(PI/314*((T)50 - step_2.y), 0));

				MATRIX_2X2 characteristic_matrix = characteristic_matrix.Identity() - dt*((T)1/6*deformation_matrix + (T)1/6*step_1_m*deformation_matrix + (T)2/3*step_2_m*deformation_matrix);

				rho_ghost.HermiteCubicInterpolation(characteristic_point, rho.gradient, hermite_cubic, hermite_cubic_x, hermite_cubic_y);

				rho_array(i, j_end) = hermite_cubic;

				gradient_new(i, j_end) = characteristic_matrix*VT(hermite_cubic_x, hermite_cubic_y);
			}
		}

		for (int i = rho.grid.i_start + 1; i <= rho.grid.i_end - 1; i++)
		{
			rho.gradient(i, j_start) = gradient_new(i, j_start);
			rho.gradient(i, j_end) = gradient_new(i, j_end);
		}

		// 1st step - Solve the characteristic equation and assign levelset function
		for (int j = rho.grid.j_start + 1; j <= rho.grid.j_end - 1; j++)
		{
			for (int i = rho.grid.i_start + 1; i <= rho.grid.i_end - 1; i++)
			{
				// Define the deformation matrix
				T u_x = (velocity_array(i + 1, j).x - velocity_array(i - 1, j).x)*rho.grid.one_over_2dx;
				T u_y = (velocity_array(i, j + 1).x - velocity_array(i, j - 1).x)*rho.grid.one_over_2dy;
				T v_x = (velocity_array(i + 1, j).y - velocity_array(i - 1, j).y)*rho.grid.one_over_2dx;
				T v_y = (velocity_array(i, j + 1).y - velocity_array(i, j - 1).y)*rho.grid.one_over_2dy;
				
				MATRIX_2X2 deformation_matrix(u_x, u_y, v_x, v_y);
				
				// Third order RK3rd - Need to change this into the general case
				VT step_1 = rho.GridPoint(i, j) - velocity_array(i, j)*dt;
					
				MATRIX_2X2 step_1_m = step_1_m.Identity() - dt*deformation_matrix;
					
				VT step_2 = rho.GridPoint(i, j) - dt*((T)0.25*velocity_array(i, j) + (T)0.25*VT(PI/314*((T)50 - step_1.y),PI/314*(step_1.x - (T)50)));
				
				MATRIX_2X2 step_2_m = step_2_m.Identity() - dt*((T)1/4*deformation_matrix + (T)1/4*step_1_m*deformation_matrix);

				characteristic_point = rho.GridPoint(i, j) - dt*((T)1/6*velocity_array(i, j) + (T)1/6*VT(PI/314*((T)50 - step_1.y),PI/314*(step_1.x - (T)50)) + (T)2/3*VT(PI/314*((T)50 - step_2.y),PI/314*(step_2.x - (T)50)));
				
				MATRIX_2X2 characteristic_matrix = characteristic_matrix.Identity() - dt*((T)1/6*deformation_matrix + (T)1/6*step_1_m*deformation_matrix + (T)2/3*step_2_m*deformation_matrix);
				
				rho_ghost.HermiteCubicInterpolation(characteristic_point, rho.gradient, hermite_cubic, hermite_cubic_x, hermite_cubic_y);
				
				rho_array(i, j) = hermite_cubic;
				
				gradient_new(i, j) = characteristic_matrix*VT(hermite_cubic_x, hermite_cubic_y);
			}
		}
		
		// 2nd step - Interpolate the gradient vector
		for (int j = rho.grid.j_start + 1; j <= rho.grid.j_end - 1; j++)
		{
			for (int i = rho.grid.i_start + 1; i <= rho.grid.i_end - 1; i++)
			{
				rho.gradient(i, j) = gradient_new(i, j);
			}
		}	
	}

	// WENO 5th 
	static void WENO5th(FIELD_STRUCTURE_2D<TT>& rho, FIELD_STRUCTURE_2D<TT>& rho_ghost, const FIELD_STRUCTURE_2D<VT>& velocity, const T& dt, const T& epsilon)
	{
		ARRAY_2D<TT>& rho_array(rho.array_for_this);

		T one_over_dx(rho.one_over_dx), one_over_dy(rho.one_over_dy);

		// x-components
		for (int j = rho.grid.j_start; j <= rho.grid.j_end; j++)
		{
			for (int i = rho.grid.i_start; i <= rho.grid.i_end; i++)
			{
				TT rho_m_x_1 = one_over_dx*(rho_ghost(i - 2, j) - rho_ghost(i - 3, j));
				TT rho_m_x_2 = one_over_dx*(rho_ghost(i - 1, j) - rho_ghost(i - 2, j));
				TT rho_m_x_3 = one_over_dx*(rho_ghost(i, j) - rho_ghost(i - 1, j));
				TT rho_m_x_4 = one_over_dx*(rho_ghost(i + 1, j) - rho_ghost(i, j));
				TT rho_m_x_5 = one_over_dx*(rho_ghost(i + 2, j) - rho_ghost(i + 1, j));

				TT rho_p_x_1 = one_over_dx*(rho_ghost(i + 3, j) - rho_ghost(i + 2, j));
				TT rho_p_x_2 = one_over_dx*(rho_ghost(i + 2, j) - rho_ghost(i + 1, j));
				TT rho_p_x_3 = one_over_dx*(rho_ghost(i + 1, j) - rho_ghost(i , j));
				TT rho_p_x_4 = one_over_dx*(rho_ghost(i , j) - rho_ghost(i - 1, j));
				TT rho_p_x_5 = one_over_dx*(rho_ghost(i - 1, j) - rho_ghost(i - 2, j));

				// y-components
				TT rho_m_y_1 = one_over_dy*(rho_ghost(i, j - 2) - rho_ghost(i, j - 3));
				TT rho_m_y_2 = one_over_dy*(rho_ghost(i, j - 1) - rho_ghost(i, j - 2));
				TT rho_m_y_3 = one_over_dy*(rho_ghost(i, j) - rho_ghost(i, j - 1));
				TT rho_m_y_4 = one_over_dy*(rho_ghost(i, j + 1) - rho_ghost(i, j));
				TT rho_m_y_5 = one_over_dy*(rho_ghost(i, j + 2) - rho_ghost(i, j + 1));

				TT rho_p_y_1 = one_over_dy*(rho_ghost(i, j + 3) - rho_ghost(i, j + 2));
				TT rho_p_y_2 = one_over_dy*(rho_ghost(i, j + 2) - rho_ghost(i, j + 1));
				TT rho_p_y_3 = one_over_dy*(rho_ghost(i, j + 1) - rho_ghost(i, j));
				TT rho_p_y_4 = one_over_dy*(rho_ghost(i, j) - rho_ghost(i, j - 1));
				TT rho_p_y_5 = one_over_dy*(rho_ghost(i, j - 1) - rho_ghost(i, j - 2));
				
				// Smoothness
				// x-components
				TT s_m_x_1 = (T)13/12*(rho_m_x_1 - 2*rho_m_x_2 + rho_m_x_3)*(rho_m_x_1 - 2*rho_m_x_2 + rho_m_x_3) + (T)1/4*(rho_m_x_1 - 4*rho_m_x_2 + 3*rho_m_x_3)*(rho_m_x_1 - 4*rho_m_x_2 + 3*rho_m_x_3);
				TT s_m_x_2 = (T)13/12*(rho_m_x_2 - 2*rho_m_x_3 + rho_m_x_4)*(rho_m_x_2 - 2*rho_m_x_3 + rho_m_x_4) + (T)1/4*(rho_m_x_2 - rho_m_x_4)*(rho_m_x_2 - rho_m_x_4);
				TT s_m_x_3 = (T)13/12*(rho_m_x_3 - 2*rho_m_x_4 + rho_m_x_5)*(rho_m_x_3 - 2*rho_m_x_4 + rho_m_x_5) + (T)1/4*(3*rho_m_x_3 - 4*rho_m_x_4 + rho_m_x_5)*(3*rho_m_x_3 - 4*rho_m_x_4 + rho_m_x_5);

				TT s_p_x_1 = (T)13/12*(rho_p_x_1 - 2*rho_p_x_2 + rho_p_x_3)*(rho_p_x_1 - 2*rho_p_x_2 + rho_p_x_3) + (T)1/4*(rho_p_x_1 - 4*rho_p_x_2 + 3*rho_p_x_3)*(rho_p_x_1 - 4*rho_p_x_2 + 3*rho_p_x_3);
				TT s_p_x_2 = (T)13/12*(rho_p_x_2 - 2*rho_p_x_3 + rho_p_x_4)*(rho_p_x_2 - 2*rho_p_x_3 + rho_p_x_4) + (T)1/4*(rho_p_x_2 - rho_p_x_4)*(rho_p_x_2 - rho_p_x_4);
				TT s_p_x_3 = (T)13/12*(rho_p_x_3 - 2*rho_p_x_4 + rho_p_x_5)*(rho_p_x_3 - 2*rho_p_x_4 + rho_p_x_5) + (T)1/4*(3*rho_p_x_3 - 4*rho_p_x_4 + rho_p_x_5)*(3*rho_p_x_3 - 4*rho_p_x_4 + rho_p_x_5);

				// y-components
				TT s_m_y_1 = (T)13/12*(rho_m_y_1 - 2*rho_m_y_2 + rho_m_y_3)*(rho_m_y_1 - 2*rho_m_y_2 + rho_m_y_3) + (T)1/4*(rho_m_y_1 - 4*rho_m_y_2 + 3*rho_m_y_3)*(rho_m_y_1 - 4*rho_m_y_2 + 3*rho_m_y_3);
				TT s_m_y_2 = (T)13/12*(rho_m_y_2 - 2*rho_m_y_3 + rho_m_y_4)*(rho_m_y_2 - 2*rho_m_y_3 + rho_m_y_4) + (T)1/4*(rho_m_y_2 - rho_m_y_4)*(rho_m_y_2 - rho_m_y_4);
				TT s_m_y_3 = (T)13/12*(rho_m_y_3 - 2*rho_m_y_4 + rho_m_y_5)*(rho_m_y_3 - 2*rho_m_y_4 + rho_m_y_5) + (T)1/4*(3*rho_m_y_3 - 4*rho_m_y_4 + rho_m_y_5)*(3*rho_m_y_3 - 4*rho_m_y_4 + rho_m_y_5);

				TT s_p_y_1 = (T)13/12*(rho_p_y_1 - 2*rho_p_y_2 + rho_p_y_3)*(rho_p_y_1 - 2*rho_p_y_2 + rho_p_y_3) + (T)1/4*(rho_p_y_1 - 4*rho_p_y_2 + 3*rho_p_y_3)*(rho_p_y_1 - 4*rho_p_y_2 + 3*rho_p_y_3);
				TT s_p_y_2 = (T)13/12*(rho_p_y_2 - 2*rho_p_y_3 + rho_p_y_4)*(rho_p_y_2 - 2*rho_p_y_3 + rho_p_y_4) + (T)1/4*(rho_p_y_2 - rho_p_y_4)*(rho_p_y_2 - rho_p_y_4);
				TT s_p_y_3 = (T)13/12*(rho_p_y_3 - 2*rho_p_y_4 + rho_p_y_5)*(rho_p_y_3 - 2*rho_p_y_4 + rho_p_y_5) + (T)1/4*(3*rho_p_y_3 - 4*rho_p_y_4 + rho_p_y_5)*(3*rho_p_y_3 - 4*rho_p_y_4 + rho_p_y_5);

				// Weights
				// x-componets
				TT a_m_x_1 = (T)1/10*(T)1/((epsilon + s_m_x_1)*(epsilon + s_m_x_1));
				TT a_m_x_2 = (T)6/10*(T)1/((epsilon + s_m_x_2)*(epsilon + s_m_x_2));
				TT a_m_x_3 = (T)3/10*(T)1/((epsilon + s_m_x_3)*(epsilon + s_m_x_3));

				TT sum_m_x = a_m_x_1 + a_m_x_2 + a_m_x_3;

				TT w_m_x_1 = a_m_x_1/sum_m_x;
				TT w_m_x_2 = a_m_x_2/sum_m_x;
				TT w_m_x_3 = a_m_x_3/sum_m_x;

				TT a_p_x_1 = (T)1/10*(T)1/((epsilon + s_p_x_1)*(epsilon + s_p_x_1));
				TT a_p_x_2 = (T)6/10*(T)1/((epsilon + s_p_x_2)*(epsilon + s_p_x_2));
				TT a_p_x_3 = (T)3/10*(T)1/((epsilon + s_p_x_3)*(epsilon + s_p_x_3));

				TT sum_p_x = a_p_x_1 + a_p_x_2 + a_p_x_3;

				TT w_p_x_1 = a_p_x_1/sum_p_x;
				TT w_p_x_2 = a_p_x_2/sum_p_x;
				TT w_p_x_3 = a_p_x_3/sum_p_x;

				// y-component
				TT a_m_y_1 = (T)1/10*(T)1/((epsilon + s_m_y_1)*(epsilon + s_m_y_1));
				TT a_m_y_2 = (T)6/10*(T)1/((epsilon + s_m_y_2)*(epsilon + s_m_y_2));
				TT a_m_y_3 = (T)3/10*(T)1/((epsilon + s_m_y_3)*(epsilon + s_m_y_3));

				TT sum_m_y = a_m_y_1 + a_m_y_2 + a_m_y_3;

				TT w_m_y_1 = a_m_y_1/sum_m_y;
				TT w_m_y_2 = a_m_y_2/sum_m_y;
				TT w_m_y_3 = a_m_y_3/sum_m_y;

				TT a_p_y_1 = (T)1/10*(T)1/((epsilon + s_p_y_1)*(epsilon + s_p_y_1));
				TT a_p_y_2 = (T)6/10*(T)1/((epsilon + s_p_y_2)*(epsilon + s_p_y_2));
				TT a_p_y_3 = (T)3/10*(T)1/((epsilon + s_p_y_3)*(epsilon + s_p_y_3));

				TT sum_p_y = a_p_y_1 + a_p_y_2 + a_p_y_3;

				TT w_p_y_1 = a_p_y_1/sum_p_y;
				TT w_p_y_2 = a_p_y_2/sum_p_y;
				TT w_p_y_3 = a_p_y_3/sum_p_y;

				// Approximation of derivatives
				// rho_x
				TT rho_m_x = w_m_x_1*(rho_m_x_1*(T)1/3 - rho_m_x_2*(T)7/6 + rho_m_x_3*(T)11/6) + w_m_x_2*(rho_m_x_2*(-(T)1/6) + rho_m_x_3*(T)5/6 + rho_m_x_4*(T)1/3) + w_m_x_3*(rho_m_x_3*(T)1/3 + rho_m_x_4*(T)5/6 - rho_m_x_5*(T)1/6);
				TT rho_p_x = w_p_x_1*(rho_p_x_1*(T)1/3 - rho_p_x_2*(T)7/6 + rho_p_x_3*(T)11/6) + w_p_x_2*(rho_p_x_2*(-(T)1/6) + rho_p_x_3*(T)5/6 + rho_p_x_4*(T)1/3) + w_p_x_3*(rho_p_x_3*(T)1/3 + rho_p_x_4*(T)5/6 - rho_p_x_5*(T)1/6);

				// rho_y
				TT rho_m_y = w_m_y_1*(rho_m_y_1*(T)1/3 - rho_m_y_2*(T)7/6 + rho_m_y_3*(T)11/6) + w_m_y_2*(rho_m_y_2*(-(T)1/6) + rho_m_y_3*(T)5/6 + rho_m_y_4*(T)1/3) + w_m_y_3*(rho_m_y_3*(T)1/3 + rho_m_y_4*(T)5/6 - rho_m_y_5*(T)1/6);
				TT rho_p_y = w_p_y_1*(rho_p_y_1*(T)1/3 - rho_p_y_2*(T)7/6 + rho_p_y_3*(T)11/6) + w_p_y_2*(rho_p_y_2*(-(T)1/6) + rho_p_y_3*(T)5/6 + rho_p_y_4*(T)1/3) + w_p_y_3*(rho_p_y_3*(T)1/3 + rho_p_y_4*(T)5/6 - rho_p_y_5*(T)1/6);

				TT u_vel, v_vel;

				u_vel = velocity(i, j).x;
				v_vel = velocity(i, j).y;

				if (u_vel > (TT)0)
				{
					if (v_vel > (TT)0)
					{
						rho_array(i, j) = rho_array(i, j) - dt*(u_vel*rho_m_x + v_vel*rho_m_y);
					}
					else
					{
						rho_array(i, j) = rho_array(i, j) - dt*(u_vel*rho_m_x + v_vel*rho_p_y);
					}
				}
				else
				{
					if (v_vel > (TT)0)
					{
						rho_array(i, j) = rho_array(i, j) - dt*(u_vel*rho_p_x + v_vel*rho_m_y);
					}
					else
					{
						rho_array(i, j) = rho_array(i, j) - dt*(u_vel*rho_p_x + v_vel*rho_p_y);
					}
				}
			}
		}
	}

	// WENO 5th for MAC grid
	static void WENO5th(FIELD_STRUCTURE_2D<TT>& rho, FIELD_STRUCTURE_2D<TT>& rho_ghost, const FIELD_STRUCTURE_2D<TT>& velocity_x, const FIELD_STRUCTURE_2D<TT>& velocity_y, const T& dt, const T& epsilon)
	{
		const ARRAY_2D<TT>& velocity_array_x(velocity_x.array_for_this);
		const ARRAY_2D<TT>& velocity_array_y(velocity_y.array_for_this);
		ARRAY_2D<TT>& rho_array(rho.array_for_this);

		T one_over_dx(rho.one_over_dx), one_over_dy(rho.one_over_dy);

		// x-components
		for (int j = rho.grid.j_start; j <= rho.grid.j_end; j++)
		{
			for (int i = rho.grid.i_start; i <= rho.grid.i_end; i++)
			{
				TT rho_m_x_1 = one_over_dx*(rho_ghost(i - 2, j) - rho_ghost(i - 3, j));
				TT rho_m_x_2 = one_over_dx*(rho_ghost(i - 1, j) - rho_ghost(i - 2, j));
				TT rho_m_x_3 = one_over_dx*(rho_ghost(i, j) - rho_ghost(i - 1, j));
				TT rho_m_x_4 = one_over_dx*(rho_ghost(i + 1, j) - rho_ghost(i, j));
				TT rho_m_x_5 = one_over_dx*(rho_ghost(i + 2, j) - rho_ghost(i + 1, j));

				TT rho_p_x_1 = one_over_dx*(rho_ghost(i + 3, j) - rho_ghost(i + 2, j));
				TT rho_p_x_2 = one_over_dx*(rho_ghost(i + 2, j) - rho_ghost(i + 1, j));
				TT rho_p_x_3 = one_over_dx*(rho_ghost(i + 1, j) - rho_ghost(i , j));
				TT rho_p_x_4 = one_over_dx*(rho_ghost(i , j) - rho_ghost(i - 1, j));
				TT rho_p_x_5 = one_over_dx*(rho_ghost(i - 1, j) - rho_ghost(i - 2, j));

				// y-components
				TT rho_m_y_1 = one_over_dy*(rho_ghost(i, j - 2) - rho_ghost(i, j - 3));
				TT rho_m_y_2 = one_over_dy*(rho_ghost(i, j - 1) - rho_ghost(i, j - 2));
				TT rho_m_y_3 = one_over_dy*(rho_ghost(i, j) - rho_ghost(i, j - 1));
				TT rho_m_y_4 = one_over_dy*(rho_ghost(i, j + 1) - rho_ghost(i, j));
				TT rho_m_y_5 = one_over_dy*(rho_ghost(i, j + 2) - rho_ghost(i, j + 1));

				TT rho_p_y_1 = one_over_dy*(rho_ghost(i, j + 3) - rho_ghost(i, j + 2));
				TT rho_p_y_2 = one_over_dy*(rho_ghost(i, j + 2) - rho_ghost(i, j + 1));
				TT rho_p_y_3 = one_over_dy*(rho_ghost(i, j + 1) - rho_ghost(i, j));
				TT rho_p_y_4 = one_over_dy*(rho_ghost(i, j) - rho_ghost(i, j - 1));
				TT rho_p_y_5 = one_over_dy*(rho_ghost(i, j - 1) - rho_ghost(i, j - 2));
				
				// Smoothness
				// x-components
				TT s_m_x_1 = (T)13/12*(rho_m_x_1 - 2*rho_m_x_2 + rho_m_x_3)*(rho_m_x_1 - 2*rho_m_x_2 + rho_m_x_3) + (T)1/4*(rho_m_x_1 - 4*rho_m_x_2 + 3*rho_m_x_3)*(rho_m_x_1 - 4*rho_m_x_2 + 3*rho_m_x_3);
				TT s_m_x_2 = (T)13/12*(rho_m_x_2 - 2*rho_m_x_3 + rho_m_x_4)*(rho_m_x_2 - 2*rho_m_x_3 + rho_m_x_4) + (T)1/4*(rho_m_x_2 - rho_m_x_4)*(rho_m_x_2 - rho_m_x_4);
				TT s_m_x_3 = (T)13/12*(rho_m_x_3 - 2*rho_m_x_4 + rho_m_x_5)*(rho_m_x_3 - 2*rho_m_x_4 + rho_m_x_5) + (T)1/4*(3*rho_m_x_3 - 4*rho_m_x_4 + rho_m_x_5)*(3*rho_m_x_3 - 4*rho_m_x_4 + rho_m_x_5);

				TT s_p_x_1 = (T)13/12*(rho_p_x_1 - 2*rho_p_x_2 + rho_p_x_3)*(rho_p_x_1 - 2*rho_p_x_2 + rho_p_x_3) + (T)1/4*(rho_p_x_1 - 4*rho_p_x_2 + 3*rho_p_x_3)*(rho_p_x_1 - 4*rho_p_x_2 + 3*rho_p_x_3);
				TT s_p_x_2 = (T)13/12*(rho_p_x_2 - 2*rho_p_x_3 + rho_p_x_4)*(rho_p_x_2 - 2*rho_p_x_3 + rho_p_x_4) + (T)1/4*(rho_p_x_2 - rho_p_x_4)*(rho_p_x_2 - rho_p_x_4);
				TT s_p_x_3 = (T)13/12*(rho_p_x_3 - 2*rho_p_x_4 + rho_p_x_5)*(rho_p_x_3 - 2*rho_p_x_4 + rho_p_x_5) + (T)1/4*(3*rho_p_x_3 - 4*rho_p_x_4 + rho_p_x_5)*(3*rho_p_x_3 - 4*rho_p_x_4 + rho_p_x_5);

				// y-components
				TT s_m_y_1 = (T)13/12*(rho_m_y_1 - 2*rho_m_y_2 + rho_m_y_3)*(rho_m_y_1 - 2*rho_m_y_2 + rho_m_y_3) + (T)1/4*(rho_m_y_1 - 4*rho_m_y_2 + 3*rho_m_y_3)*(rho_m_y_1 - 4*rho_m_y_2 + 3*rho_m_y_3);
				TT s_m_y_2 = (T)13/12*(rho_m_y_2 - 2*rho_m_y_3 + rho_m_y_4)*(rho_m_y_2 - 2*rho_m_y_3 + rho_m_y_4) + (T)1/4*(rho_m_y_2 - rho_m_y_4)*(rho_m_y_2 - rho_m_y_4);
				TT s_m_y_3 = (T)13/12*(rho_m_y_3 - 2*rho_m_y_4 + rho_m_y_5)*(rho_m_y_3 - 2*rho_m_y_4 + rho_m_y_5) + (T)1/4*(3*rho_m_y_3 - 4*rho_m_y_4 + rho_m_y_5)*(3*rho_m_y_3 - 4*rho_m_y_4 + rho_m_y_5);

				TT s_p_y_1 = (T)13/12*(rho_p_y_1 - 2*rho_p_y_2 + rho_p_y_3)*(rho_p_y_1 - 2*rho_p_y_2 + rho_p_y_3) + (T)1/4*(rho_p_y_1 - 4*rho_p_y_2 + 3*rho_p_y_3)*(rho_p_y_1 - 4*rho_p_y_2 + 3*rho_p_y_3);
				TT s_p_y_2 = (T)13/12*(rho_p_y_2 - 2*rho_p_y_3 + rho_p_y_4)*(rho_p_y_2 - 2*rho_p_y_3 + rho_p_y_4) + (T)1/4*(rho_p_y_2 - rho_p_y_4)*(rho_p_y_2 - rho_p_y_4);
				TT s_p_y_3 = (T)13/12*(rho_p_y_3 - 2*rho_p_y_4 + rho_p_y_5)*(rho_p_y_3 - 2*rho_p_y_4 + rho_p_y_5) + (T)1/4*(3*rho_p_y_3 - 4*rho_p_y_4 + rho_p_y_5)*(3*rho_p_y_3 - 4*rho_p_y_4 + rho_p_y_5);

				// Weights
				// x-componets
				TT a_m_x_1 = (T)1/10*(T)1/((epsilon + s_m_x_1)*(epsilon + s_m_x_1));
				TT a_m_x_2 = (T)6/10*(T)1/((epsilon + s_m_x_2)*(epsilon + s_m_x_2));
				TT a_m_x_3 = (T)3/10*(T)1/((epsilon + s_m_x_3)*(epsilon + s_m_x_3));

				TT sum_m_x = a_m_x_1 + a_m_x_2 + a_m_x_3;

				TT w_m_x_1 = a_m_x_1/sum_m_x;
				TT w_m_x_2 = a_m_x_2/sum_m_x;
				TT w_m_x_3 = a_m_x_3/sum_m_x;

				TT a_p_x_1 = (T)1/10*(T)1/((epsilon + s_p_x_1)*(epsilon + s_p_x_1));
				TT a_p_x_2 = (T)6/10*(T)1/((epsilon + s_p_x_2)*(epsilon + s_p_x_2));
				TT a_p_x_3 = (T)3/10*(T)1/((epsilon + s_p_x_3)*(epsilon + s_p_x_3));

				TT sum_p_x = a_p_x_1 + a_p_x_2 + a_p_x_3;

				TT w_p_x_1 = a_p_x_1/sum_p_x;
				TT w_p_x_2 = a_p_x_2/sum_p_x;
				TT w_p_x_3 = a_p_x_3/sum_p_x;

				// y-component
				TT a_m_y_1 = (T)1/10*(T)1/((epsilon + s_m_y_1)*(epsilon + s_m_y_1));
				TT a_m_y_2 = (T)6/10*(T)1/((epsilon + s_m_y_2)*(epsilon + s_m_y_2));
				TT a_m_y_3 = (T)3/10*(T)1/((epsilon + s_m_y_3)*(epsilon + s_m_y_3));

				TT sum_m_y = a_m_y_1 + a_m_y_2 + a_m_y_3;

				TT w_m_y_1 = a_m_y_1/sum_m_y;
				TT w_m_y_2 = a_m_y_2/sum_m_y;
				TT w_m_y_3 = a_m_y_3/sum_m_y;

				TT a_p_y_1 = (T)1/10*(T)1/((epsilon + s_p_y_1)*(epsilon + s_p_y_1));
				TT a_p_y_2 = (T)6/10*(T)1/((epsilon + s_p_y_2)*(epsilon + s_p_y_2));
				TT a_p_y_3 = (T)3/10*(T)1/((epsilon + s_p_y_3)*(epsilon + s_p_y_3));

				TT sum_p_y = a_p_y_1 + a_p_y_2 + a_p_y_3;

				TT w_p_y_1 = a_p_y_1/sum_p_y;
				TT w_p_y_2 = a_p_y_2/sum_p_y;
				TT w_p_y_3 = a_p_y_3/sum_p_y;

				// Approximation of derivatives
				// rho_x
				TT rho_m_x = w_m_x_1*(rho_m_x_1*(T)1/3 - rho_m_x_2*(T)7/6 + rho_m_x_3*(T)11/6) + w_m_x_2*(rho_m_x_2*(-(T)1/6) + rho_m_x_3*(T)5/6 + rho_m_x_4*(T)1/3) + w_m_x_3*(rho_m_x_3*(T)1/3 + rho_m_x_4*(T)5/6 - rho_m_x_5*(T)1/6);
				TT rho_p_x = w_p_x_1*(rho_p_x_1*(T)1/3 - rho_p_x_2*(T)7/6 + rho_p_x_3*(T)11/6) + w_p_x_2*(rho_p_x_2*(-(T)1/6) + rho_p_x_3*(T)5/6 + rho_p_x_4*(T)1/3) + w_p_x_3*(rho_p_x_3*(T)1/3 + rho_p_x_4*(T)5/6 - rho_p_x_5*(T)1/6);

				// rho_y
				TT rho_m_y = w_m_y_1*(rho_m_y_1*(T)1/3 - rho_m_y_2*(T)7/6 + rho_m_y_3*(T)11/6) + w_m_y_2*(rho_m_y_2*(-(T)1/6) + rho_m_y_3*(T)5/6 + rho_m_y_4*(T)1/3) + w_m_y_3*(rho_m_y_3*(T)1/3 + rho_m_y_4*(T)5/6 - rho_m_y_5*(T)1/6);
				TT rho_p_y = w_p_y_1*(rho_p_y_1*(T)1/3 - rho_p_y_2*(T)7/6 + rho_p_y_3*(T)11/6) + w_p_y_2*(rho_p_y_2*(-(T)1/6) + rho_p_y_3*(T)5/6 + rho_p_y_4*(T)1/3) + w_p_y_3*(rho_p_y_3*(T)1/3 + rho_p_y_4*(T)5/6 - rho_p_y_5*(T)1/6);

				TT u_vel, v_vel;

				// Velocity interpolation - Using MAC grid
				if (rho.is_scalar == true)
				{
					u_vel = (velocity_array_x(i + 1, j) + velocity_array_x(i, j))*(T)0.5;
					v_vel = (velocity_array_y(i, j + 1) + velocity_array_y(i, j))*(T)0.5;
				}

				if (rho.is_x_component == true)
				{
					u_vel = velocity_array_x(i, j);
					v_vel = (velocity_array_y(i, j) + velocity_array_y(i, j + 1) + velocity_array_y(i - 1, j) + velocity_array_y(i - 1, j + 1))*(T)0.25;
				}
				else if (rho.is_y_component == true)
				{
					u_vel = (velocity_array_x(i, j - 1) + velocity_array_x(i + 1, j - 1) + velocity_array_x(i, j) + velocity_array_x(i + 1, j))*(T)0.25;
					v_vel = velocity_array_y(i, j);
				}
				
				if (u_vel > 0)
				{
					if (v_vel > 0)
					{
						rho_array(i, j) = rho_array(i, j) - dt*(u_vel*rho_m_x + v_vel*rho_m_y);
					}
					else
					{
						rho_array(i, j) = rho_array(i, j) - dt*(u_vel*rho_m_x + v_vel*rho_p_y);
					}
				}
				else
				{
					if (v_vel > 0)
					{
						rho_array(i, j) = rho_array(i, j) - dt*(u_vel*rho_p_x + v_vel*rho_m_y);
					}
					else
					{
						rho_array(i, j) = rho_array(i, j) - dt*(u_vel*rho_p_x + v_vel*rho_p_y);
					}
				}
			}
		}
	}

	static void WENO5th(FIELD_STRUCTURE_2D<TT>& rho, FIELD_STRUCTURE_2D<TT>& rho_ghost, const FIELD_STRUCTURE_2D<TT>& velocity_x, const FIELD_STRUCTURE_2D<TT>& velocity_y, const T& dt, const T& epsilon, MULTITHREADING& multithreading, const int& thread_id)
	{
		const ARRAY_2D<TT>& velocity_array_x(velocity_x.array_for_this);
		const ARRAY_2D<TT>& velocity_array_y(velocity_y.array_for_this);
		ARRAY_2D<TT>& rho_array(rho.array_for_this);

		rho_ghost.FillGhostCellsFrom(rho_array, true, thread_id);

		multithreading.Sync(thread_id);

		T one_over_dx(rho.one_over_dx), one_over_dy(rho.one_over_dy);
				
		GRID_ITERATION_2D(rho.partial_grids[thread_id])
		{
			// x-components
			TT rho_m_x_1 = one_over_dx*(rho_ghost(i - 2, j) - rho_ghost(i - 3, j));
			TT rho_m_x_2 = one_over_dx*(rho_ghost(i - 1, j) - rho_ghost(i - 2, j));
			TT rho_m_x_3 = one_over_dx*(rho_ghost(i, j) - rho_ghost(i - 1, j));
			TT rho_m_x_4 = one_over_dx*(rho_ghost(i + 1, j) - rho_ghost(i, j));
			TT rho_m_x_5 = one_over_dx*(rho_ghost(i + 2, j) - rho_ghost(i + 1, j));

			TT rho_p_x_1 = one_over_dx*(rho_ghost(i + 3, j) - rho_ghost(i + 2, j));
			TT rho_p_x_2 = one_over_dx*(rho_ghost(i + 2, j) - rho_ghost(i + 1, j));
			TT rho_p_x_3 = one_over_dx*(rho_ghost(i + 1, j) - rho_ghost(i , j));
			TT rho_p_x_4 = one_over_dx*(rho_ghost(i , j) - rho_ghost(i - 1, j));
			TT rho_p_x_5 = one_over_dx*(rho_ghost(i - 1, j) - rho_ghost(i - 2, j));

			// y-components
			TT rho_m_y_1 = one_over_dy*(rho_ghost(i, j - 2) - rho_ghost(i, j - 3));
			TT rho_m_y_2 = one_over_dy*(rho_ghost(i, j - 1) - rho_ghost(i, j - 2));
			TT rho_m_y_3 = one_over_dy*(rho_ghost(i, j) - rho_ghost(i, j - 1));
			TT rho_m_y_4 = one_over_dy*(rho_ghost(i, j + 1) - rho_ghost(i, j));
			TT rho_m_y_5 = one_over_dy*(rho_ghost(i, j + 2) - rho_ghost(i, j + 1));

			TT rho_p_y_1 = one_over_dy*(rho_ghost(i, j + 3) - rho_ghost(i, j + 2));
			TT rho_p_y_2 = one_over_dy*(rho_ghost(i, j + 2) - rho_ghost(i, j + 1));
			TT rho_p_y_3 = one_over_dy*(rho_ghost(i, j + 1) - rho_ghost(i, j));
			TT rho_p_y_4 = one_over_dy*(rho_ghost(i, j) - rho_ghost(i, j - 1));
			TT rho_p_y_5 = one_over_dy*(rho_ghost(i, j - 1) - rho_ghost(i, j - 2));
				
			// Smoothness
			// x-components
			TT s_m_x_1 = (T)13/12*(rho_m_x_1 - 2*rho_m_x_2 + rho_m_x_3)*(rho_m_x_1 - 2*rho_m_x_2 + rho_m_x_3) + (T)1/4*(rho_m_x_1 - 4*rho_m_x_2 + 3*rho_m_x_3)*(rho_m_x_1 - 4*rho_m_x_2 + 3*rho_m_x_3);
			TT s_m_x_2 = (T)13/12*(rho_m_x_2 - 2*rho_m_x_3 + rho_m_x_4)*(rho_m_x_2 - 2*rho_m_x_3 + rho_m_x_4) + (T)1/4*(rho_m_x_2 - rho_m_x_4)*(rho_m_x_2 - rho_m_x_4);
			TT s_m_x_3 = (T)13/12*(rho_m_x_3 - 2*rho_m_x_4 + rho_m_x_5)*(rho_m_x_3 - 2*rho_m_x_4 + rho_m_x_5) + (T)1/4*(3*rho_m_x_3 - 4*rho_m_x_4 + rho_m_x_5)*(3*rho_m_x_3 - 4*rho_m_x_4 + rho_m_x_5);

			TT s_p_x_1 = (T)13/12*(rho_p_x_1 - 2*rho_p_x_2 + rho_p_x_3)*(rho_p_x_1 - 2*rho_p_x_2 + rho_p_x_3) + (T)1/4*(rho_p_x_1 - 4*rho_p_x_2 + 3*rho_p_x_3)*(rho_p_x_1 - 4*rho_p_x_2 + 3*rho_p_x_3);
			TT s_p_x_2 = (T)13/12*(rho_p_x_2 - 2*rho_p_x_3 + rho_p_x_4)*(rho_p_x_2 - 2*rho_p_x_3 + rho_p_x_4) + (T)1/4*(rho_p_x_2 - rho_p_x_4)*(rho_p_x_2 - rho_p_x_4);
			TT s_p_x_3 = (T)13/12*(rho_p_x_3 - 2*rho_p_x_4 + rho_p_x_5)*(rho_p_x_3 - 2*rho_p_x_4 + rho_p_x_5) + (T)1/4*(3*rho_p_x_3 - 4*rho_p_x_4 + rho_p_x_5)*(3*rho_p_x_3 - 4*rho_p_x_4 + rho_p_x_5);

			// y-components
			TT s_m_y_1 = (T)13/12*(rho_m_y_1 - 2*rho_m_y_2 + rho_m_y_3)*(rho_m_y_1 - 2*rho_m_y_2 + rho_m_y_3) + (T)1/4*(rho_m_y_1 - 4*rho_m_y_2 + 3*rho_m_y_3)*(rho_m_y_1 - 4*rho_m_y_2 + 3*rho_m_y_3);
			TT s_m_y_2 = (T)13/12*(rho_m_y_2 - 2*rho_m_y_3 + rho_m_y_4)*(rho_m_y_2 - 2*rho_m_y_3 + rho_m_y_4) + (T)1/4*(rho_m_y_2 - rho_m_y_4)*(rho_m_y_2 - rho_m_y_4);
			TT s_m_y_3 = (T)13/12*(rho_m_y_3 - 2*rho_m_y_4 + rho_m_y_5)*(rho_m_y_3 - 2*rho_m_y_4 + rho_m_y_5) + (T)1/4*(3*rho_m_y_3 - 4*rho_m_y_4 + rho_m_y_5)*(3*rho_m_y_3 - 4*rho_m_y_4 + rho_m_y_5);

			TT s_p_y_1 = (T)13/12*(rho_p_y_1 - 2*rho_p_y_2 + rho_p_y_3)*(rho_p_y_1 - 2*rho_p_y_2 + rho_p_y_3) + (T)1/4*(rho_p_y_1 - 4*rho_p_y_2 + 3*rho_p_y_3)*(rho_p_y_1 - 4*rho_p_y_2 + 3*rho_p_y_3);
			TT s_p_y_2 = (T)13/12*(rho_p_y_2 - 2*rho_p_y_3 + rho_p_y_4)*(rho_p_y_2 - 2*rho_p_y_3 + rho_p_y_4) + (T)1/4*(rho_p_y_2 - rho_p_y_4)*(rho_p_y_2 - rho_p_y_4);
			TT s_p_y_3 = (T)13/12*(rho_p_y_3 - 2*rho_p_y_4 + rho_p_y_5)*(rho_p_y_3 - 2*rho_p_y_4 + rho_p_y_5) + (T)1/4*(3*rho_p_y_3 - 4*rho_p_y_4 + rho_p_y_5)*(3*rho_p_y_3 - 4*rho_p_y_4 + rho_p_y_5);

			// Weights
			// x-componets
			TT a_m_x_1 = (T)1/10*(T)1/((epsilon + s_m_x_1)*(epsilon + s_m_x_1));
			TT a_m_x_2 = (T)6/10*(T)1/((epsilon + s_m_x_2)*(epsilon + s_m_x_2));
			TT a_m_x_3 = (T)3/10*(T)1/((epsilon + s_m_x_3)*(epsilon + s_m_x_3));

			TT sum_m_x = a_m_x_1 + a_m_x_2 + a_m_x_3;

			TT w_m_x_1 = a_m_x_1/sum_m_x;
			TT w_m_x_2 = a_m_x_2/sum_m_x;
			TT w_m_x_3 = a_m_x_3/sum_m_x;

			TT a_p_x_1 = (T)1/10*(T)1/((epsilon + s_p_x_1)*(epsilon + s_p_x_1));
			TT a_p_x_2 = (T)6/10*(T)1/((epsilon + s_p_x_2)*(epsilon + s_p_x_2));
			TT a_p_x_3 = (T)3/10*(T)1/((epsilon + s_p_x_3)*(epsilon + s_p_x_3));

			TT sum_p_x = a_p_x_1 + a_p_x_2 + a_p_x_3;

			TT w_p_x_1 = a_p_x_1/sum_p_x;
			TT w_p_x_2 = a_p_x_2/sum_p_x;
			TT w_p_x_3 = a_p_x_3/sum_p_x;

			// y-component
			TT a_m_y_1 = (T)1/10*(T)1/((epsilon + s_m_y_1)*(epsilon + s_m_y_1));
			TT a_m_y_2 = (T)6/10*(T)1/((epsilon + s_m_y_2)*(epsilon + s_m_y_2));
			TT a_m_y_3 = (T)3/10*(T)1/((epsilon + s_m_y_3)*(epsilon + s_m_y_3));

			TT sum_m_y = a_m_y_1 + a_m_y_2 + a_m_y_3;

			TT w_m_y_1 = a_m_y_1/sum_m_y;
			TT w_m_y_2 = a_m_y_2/sum_m_y;
			TT w_m_y_3 = a_m_y_3/sum_m_y;

			TT a_p_y_1 = (T)1/10*(T)1/((epsilon + s_p_y_1)*(epsilon + s_p_y_1));
			TT a_p_y_2 = (T)6/10*(T)1/((epsilon + s_p_y_2)*(epsilon + s_p_y_2));
			TT a_p_y_3 = (T)3/10*(T)1/((epsilon + s_p_y_3)*(epsilon + s_p_y_3));

			TT sum_p_y = a_p_y_1 + a_p_y_2 + a_p_y_3;

			TT w_p_y_1 = a_p_y_1/sum_p_y;
			TT w_p_y_2 = a_p_y_2/sum_p_y;
			TT w_p_y_3 = a_p_y_3/sum_p_y;

			// Approximation of derivatives
			// rho_x
			TT rho_m_x = w_m_x_1*(rho_m_x_1*(T)1/3 - rho_m_x_2*(T)7/6 + rho_m_x_3*(T)11/6) + w_m_x_2*(rho_m_x_2*(-(T)1/6) + rho_m_x_3*(T)5/6 + rho_m_x_4*(T)1/3) + w_m_x_3*(rho_m_x_3*(T)1/3 + rho_m_x_4*(T)5/6 - rho_m_x_5*(T)1/6);
			TT rho_p_x = w_p_x_1*(rho_p_x_1*(T)1/3 - rho_p_x_2*(T)7/6 + rho_p_x_3*(T)11/6) + w_p_x_2*(rho_p_x_2*(-(T)1/6) + rho_p_x_3*(T)5/6 + rho_p_x_4*(T)1/3) + w_p_x_3*(rho_p_x_3*(T)1/3 + rho_p_x_4*(T)5/6 - rho_p_x_5*(T)1/6);

			// rho_y
			TT rho_m_y = w_m_y_1*(rho_m_y_1*(T)1/3 - rho_m_y_2*(T)7/6 + rho_m_y_3*(T)11/6) + w_m_y_2*(rho_m_y_2*(-(T)1/6) + rho_m_y_3*(T)5/6 + rho_m_y_4*(T)1/3) + w_m_y_3*(rho_m_y_3*(T)1/3 + rho_m_y_4*(T)5/6 - rho_m_y_5*(T)1/6);
			TT rho_p_y = w_p_y_1*(rho_p_y_1*(T)1/3 - rho_p_y_2*(T)7/6 + rho_p_y_3*(T)11/6) + w_p_y_2*(rho_p_y_2*(-(T)1/6) + rho_p_y_3*(T)5/6 + rho_p_y_4*(T)1/3) + w_p_y_3*(rho_p_y_3*(T)1/3 + rho_p_y_4*(T)5/6 - rho_p_y_5*(T)1/6);

			TT u_vel, v_vel;

			// Velocity interpolation - Using MAC grid
			if (rho.is_scalar == true)
			{
				u_vel = (velocity_array_x(i + 1, j) + velocity_array_x(i, j))*(T)0.5;
				v_vel = (velocity_array_y(i, j + 1) + velocity_array_y(i, j))*(T)0.5;
			}

			if (rho.is_x_component == true)
			{
				u_vel = velocity_array_x(i, j);
				v_vel = (velocity_array_y(i, j) + velocity_array_y(i, j + 1) + velocity_array_y(i - 1, j) + velocity_array_y(i - 1, j + 1))*(T)0.25;
			}
			else if (rho.is_y_component == true)
			{
				u_vel = (velocity_array_x(i, j - 1) + velocity_array_x(i + 1, j - 1) + velocity_array_x(i, j) + velocity_array_x(i + 1, j))*(T)0.25;
				v_vel = velocity_array_y(i, j);
			}
				
			if (u_vel > 0)
			{
				if (v_vel > 0)
				{
					rho_array(i, j) = rho_array(i, j) - dt*(u_vel*rho_m_x + v_vel*rho_m_y);
				}
				else
				{
					rho_array(i, j) = rho_array(i, j) - dt*(u_vel*rho_m_x + v_vel*rho_p_y);
				}
			}
			else
			{
				if (v_vel > 0)
				{
					rho_array(i, j) = rho_array(i, j) - dt*(u_vel*rho_p_x + v_vel*rho_m_y);
				}
				else
				{
					rho_array(i, j) = rho_array(i, j) - dt*(u_vel*rho_p_x + v_vel*rho_p_y);
				}
			}
		}
		multithreading.Sync(thread_id);
	}

	static void WENO5thReinitialization(FIELD_STRUCTURE_2D<TT>& rho, FIELD_STRUCTURE_2D<TT>& rho_ghost, const T& dt, const T& epsilon, const FIELD_STRUCTURE_2D<T>& sign_function)
	{
		ARRAY_2D<TT>& rho_array(rho.array_for_this);

		rho_ghost.FillGhostCellsFrom(rho_array, true);

		T one_over_dx(rho.one_over_dx), one_over_dy(rho.one_over_dy);

		// x-components
		for (int j = rho.grid.j_start; j <= rho.grid.j_end; j++)
		{
			for (int i = rho.grid.i_start; i <= rho.grid.i_end; i++)
			{
				TT rho_m_x_1 = one_over_dx*(rho_ghost(i - 2, j) - rho_ghost(i - 3, j));
				TT rho_m_x_2 = one_over_dx*(rho_ghost(i - 1, j) - rho_ghost(i - 2, j));
				TT rho_m_x_3 = one_over_dx*(rho_ghost(i, j) - rho_ghost(i - 1, j));
				TT rho_m_x_4 = one_over_dx*(rho_ghost(i + 1, j) - rho_ghost(i, j));
				TT rho_m_x_5 = one_over_dx*(rho_ghost(i + 2, j) - rho_ghost(i + 1, j));

				TT rho_p_x_1 = one_over_dx*(rho_ghost(i + 3, j) - rho_ghost(i + 2, j));
				TT rho_p_x_2 = one_over_dx*(rho_ghost(i + 2, j) - rho_ghost(i + 1, j));
				TT rho_p_x_3 = one_over_dx*(rho_ghost(i + 1, j) - rho_ghost(i , j));
				TT rho_p_x_4 = one_over_dx*(rho_ghost(i , j) - rho_ghost(i - 1, j));
				TT rho_p_x_5 = one_over_dx*(rho_ghost(i - 1, j) - rho_ghost(i - 2, j));

				// y-components
				TT rho_m_y_1 = one_over_dy*(rho_ghost(i, j - 2) - rho_ghost(i, j - 3));
				TT rho_m_y_2 = one_over_dy*(rho_ghost(i, j - 1) - rho_ghost(i, j - 2));
				TT rho_m_y_3 = one_over_dy*(rho_ghost(i, j) - rho_ghost(i, j - 1));
				TT rho_m_y_4 = one_over_dy*(rho_ghost(i, j + 1) - rho_ghost(i, j));
				TT rho_m_y_5 = one_over_dy*(rho_ghost(i, j + 2) - rho_ghost(i, j + 1));

				TT rho_p_y_1 = one_over_dy*(rho_ghost(i, j + 3) - rho_ghost(i, j + 2));
				TT rho_p_y_2 = one_over_dy*(rho_ghost(i, j + 2) - rho_ghost(i, j + 1));
				TT rho_p_y_3 = one_over_dy*(rho_ghost(i, j + 1) - rho_ghost(i, j));
				TT rho_p_y_4 = one_over_dy*(rho_ghost(i, j) - rho_ghost(i, j - 1));
				TT rho_p_y_5 = one_over_dy*(rho_ghost(i, j - 1) - rho_ghost(i, j - 2));
				
				// Smoothness
				// x-components
				TT s_m_x_1 = (T)13/12*(rho_m_x_1 - 2*rho_m_x_2 + rho_m_x_3)*(rho_m_x_1 - 2*rho_m_x_2 + rho_m_x_3) + (T)1/4*(rho_m_x_1 - 4*rho_m_x_2 + 3*rho_m_x_3)*(rho_m_x_1 - 4*rho_m_x_2 + 3*rho_m_x_3);
				TT s_m_x_2 = (T)13/12*(rho_m_x_2 - 2*rho_m_x_3 + rho_m_x_4)*(rho_m_x_2 - 2*rho_m_x_3 + rho_m_x_4) + (T)1/4*(rho_m_x_2 - rho_m_x_4)*(rho_m_x_2 - rho_m_x_4);
				TT s_m_x_3 = (T)13/12*(rho_m_x_3 - 2*rho_m_x_4 + rho_m_x_5)*(rho_m_x_3 - 2*rho_m_x_4 + rho_m_x_5) + (T)1/4*(3*rho_m_x_3 - 4*rho_m_x_4 + rho_m_x_5)*(3*rho_m_x_3 - 4*rho_m_x_4 + rho_m_x_5);

				TT s_p_x_1 = (T)13/12*(rho_p_x_1 - 2*rho_p_x_2 + rho_p_x_3)*(rho_p_x_1 - 2*rho_p_x_2 + rho_p_x_3) + (T)1/4*(rho_p_x_1 - 4*rho_p_x_2 + 3*rho_p_x_3)*(rho_p_x_1 - 4*rho_p_x_2 + 3*rho_p_x_3);
				TT s_p_x_2 = (T)13/12*(rho_p_x_2 - 2*rho_p_x_3 + rho_p_x_4)*(rho_p_x_2 - 2*rho_p_x_3 + rho_p_x_4) + (T)1/4*(rho_p_x_2 - rho_p_x_4)*(rho_p_x_2 - rho_p_x_4);
				TT s_p_x_3 = (T)13/12*(rho_p_x_3 - 2*rho_p_x_4 + rho_p_x_5)*(rho_p_x_3 - 2*rho_p_x_4 + rho_p_x_5) + (T)1/4*(3*rho_p_x_3 - 4*rho_p_x_4 + rho_p_x_5)*(3*rho_p_x_3 - 4*rho_p_x_4 + rho_p_x_5);

				// y-components
				TT s_m_y_1 = (T)13/12*(rho_m_y_1 - 2*rho_m_y_2 + rho_m_y_3)*(rho_m_y_1 - 2*rho_m_y_2 + rho_m_y_3) + (T)1/4*(rho_m_y_1 - 4*rho_m_y_2 + 3*rho_m_y_3)*(rho_m_y_1 - 4*rho_m_y_2 + 3*rho_m_y_3);
				TT s_m_y_2 = (T)13/12*(rho_m_y_2 - 2*rho_m_y_3 + rho_m_y_4)*(rho_m_y_2 - 2*rho_m_y_3 + rho_m_y_4) + (T)1/4*(rho_m_y_2 - rho_m_y_4)*(rho_m_y_2 - rho_m_y_4);
				TT s_m_y_3 = (T)13/12*(rho_m_y_3 - 2*rho_m_y_4 + rho_m_y_5)*(rho_m_y_3 - 2*rho_m_y_4 + rho_m_y_5) + (T)1/4*(3*rho_m_y_3 - 4*rho_m_y_4 + rho_m_y_5)*(3*rho_m_y_3 - 4*rho_m_y_4 + rho_m_y_5);

				TT s_p_y_1 = (T)13/12*(rho_p_y_1 - 2*rho_p_y_2 + rho_p_y_3)*(rho_p_y_1 - 2*rho_p_y_2 + rho_p_y_3) + (T)1/4*(rho_p_y_1 - 4*rho_p_y_2 + 3*rho_p_y_3)*(rho_p_y_1 - 4*rho_p_y_2 + 3*rho_p_y_3);
				TT s_p_y_2 = (T)13/12*(rho_p_y_2 - 2*rho_p_y_3 + rho_p_y_4)*(rho_p_y_2 - 2*rho_p_y_3 + rho_p_y_4) + (T)1/4*(rho_p_y_2 - rho_p_y_4)*(rho_p_y_2 - rho_p_y_4);
				TT s_p_y_3 = (T)13/12*(rho_p_y_3 - 2*rho_p_y_4 + rho_p_y_5)*(rho_p_y_3 - 2*rho_p_y_4 + rho_p_y_5) + (T)1/4*(3*rho_p_y_3 - 4*rho_p_y_4 + rho_p_y_5)*(3*rho_p_y_3 - 4*rho_p_y_4 + rho_p_y_5);

				// Weights
				// x-componets
				TT a_m_x_1 = (T)1/10*(T)1/((epsilon + s_m_x_1)*(epsilon + s_m_x_1));
				TT a_m_x_2 = (T)6/10*(T)1/((epsilon + s_m_x_2)*(epsilon + s_m_x_2));
				TT a_m_x_3 = (T)3/10*(T)1/((epsilon + s_m_x_3)*(epsilon + s_m_x_3));

				TT sum_m_x = a_m_x_1 + a_m_x_2 + a_m_x_3;

				TT w_m_x_1 = a_m_x_1/sum_m_x;
				TT w_m_x_2 = a_m_x_2/sum_m_x;
				TT w_m_x_3 = a_m_x_3/sum_m_x;
	
				TT a_p_x_1 = (T)1/10*(T)1/((epsilon + s_p_x_1)*(epsilon + s_p_x_1));
				TT a_p_x_2 = (T)6/10*(T)1/((epsilon + s_p_x_2)*(epsilon + s_p_x_2));
				TT a_p_x_3 = (T)3/10*(T)1/((epsilon + s_p_x_3)*(epsilon + s_p_x_3));

				TT sum_p_x = a_p_x_1 + a_p_x_2 + a_p_x_3;
	
				TT w_p_x_1 = a_p_x_1/sum_p_x;
				TT w_p_x_2 = a_p_x_2/sum_p_x;
				TT w_p_x_3 = a_p_x_3/sum_p_x;

				// y-component
				TT a_m_y_1 = (T)1/10*(T)1/((epsilon + s_m_y_1)*(epsilon + s_m_y_1));
				TT a_m_y_2 = (T)6/10*(T)1/((epsilon + s_m_y_2)*(epsilon + s_m_y_2));
				TT a_m_y_3 = (T)3/10*(T)1/((epsilon + s_m_y_3)*(epsilon + s_m_y_3));

				TT sum_m_y = a_m_y_1 + a_m_y_2 + a_m_y_3;
	
				TT w_m_y_1 = a_m_y_1/sum_m_y;
				TT w_m_y_2 = a_m_y_2/sum_m_y;
				TT w_m_y_3 = a_m_y_3/sum_m_y;
	
				TT a_p_y_1 = (T)1/10*(T)1/((epsilon + s_p_y_1)*(epsilon + s_p_y_1));
				TT a_p_y_2 = (T)6/10*(T)1/((epsilon + s_p_y_2)*(epsilon + s_p_y_2));
				TT a_p_y_3 = (T)3/10*(T)1/((epsilon + s_p_y_3)*(epsilon + s_p_y_3));

				TT sum_p_y = a_p_y_1 + a_p_y_2 + a_p_y_3;

				TT w_p_y_1 = a_p_y_1/sum_p_y;
				TT w_p_y_2 = a_p_y_2/sum_p_y;
				TT w_p_y_3 = a_p_y_3/sum_p_y;

				// Approximation of derivatives
				// rho_x
				T rho_m_x = w_m_x_1*(rho_m_x_1*(T)1/3 - rho_m_x_2*(T)7/6 + rho_m_x_3*(T)11/6) + w_m_x_2*(rho_m_x_2*(-(T)1/6) + rho_m_x_3*(T)5/6 + rho_m_x_4*(T)1/3) + w_m_x_3*(rho_m_x_3*(T)1/3 + rho_m_x_4*(T)5/6 - rho_m_x_5*(T)1/6);
				T rho_p_x = w_p_x_1*(rho_p_x_1*(T)1/3 - rho_p_x_2*(T)7/6 + rho_p_x_3*(T)11/6) + w_p_x_2*(rho_p_x_2*(-(T)1/6) + rho_p_x_3*(T)5/6 + rho_p_x_4*(T)1/3) + w_p_x_3*(rho_p_x_3*(T)1/3 + rho_p_x_4*(T)5/6 - rho_p_x_5*(T)1/6);
	
				// rho_y
				T rho_m_y = w_m_y_1*(rho_m_y_1*(T)1/3 - rho_m_y_2*(T)7/6 + rho_m_y_3*(T)11/6) + w_m_y_2*(rho_m_y_2*(-(T)1/6) + rho_m_y_3*(T)5/6 + rho_m_y_4*(T)1/3) + w_m_y_3*(rho_m_y_3*(T)1/3 + rho_m_y_4*(T)5/6 - rho_m_y_5*(T)1/6);
				T rho_p_y = w_p_y_1*(rho_p_y_1*(T)1/3 - rho_p_y_2*(T)7/6 + rho_p_y_3*(T)11/6) + w_p_y_2*(rho_p_y_2*(-(T)1/6) + rho_p_y_3*(T)5/6 + rho_p_y_4*(T)1/3) + w_p_y_3*(rho_p_y_3*(T)1/3 + rho_p_y_4*(T)5/6 - rho_p_y_5*(T)1/6);
			
				// Define the sign function
				//T s = rho(i, j)/sqrt(POW2(rho(i, j)) + POW2(rho.dx));

				T s = sign_function(i, j);

				// Define Sub functions
				T smx = s*rho_m_x, spx = s*rho_p_x, smy = s*rho_m_y, spy = s*rho_p_y;

				if ((spx >= 0) && (smx >= 0))
				{
					if ((spy >= 0) && (smy >= 0))
					{
						rho_array(i, j) = rho_array(i, j) - dt*s*(sqrt(POW2(rho_m_x) + POW2(rho_m_y)) - 1);
					}
					else if ((spy <= 0) && (smy <= 0))
					{
						rho_array(i, j) = rho_array(i, j) - dt*s*(sqrt(POW2(rho_m_x) + POW2(rho_p_y)) - 1);
					}
					else if ((spy > 0) && (smy < 0))
					{
						rho_array(i, j) = rho_array(i, j) - dt*s*(sqrt(POW2(rho_m_x)) - 1);
					}
					else if ((spy < 0) && (smy > 0))
					{
						T ss = s*(abs(rho_p_y) - abs(rho_m_y))/(rho_p_y - rho_m_y);
						if (ss > 0)
						{
							rho_array(i, j) = rho_array(i, j) - dt*s*(sqrt(POW2(rho_m_x) + POW2(rho_m_y)) - 1);
						}
						else
						{
							rho_array(i, j) = rho_array(i, j) - dt*s*(sqrt(POW2(rho_m_x) + POW2(rho_p_y)) - 1);
						}
					}
				}
				else if ((spx <= 0) && (smx <= 0))
				{
					if ((spy >= 0) && (smy >= 0))
					{
						rho_array(i, j) = rho_array(i, j) - dt*s*(sqrt(POW2(rho_p_x) + POW2(rho_m_y)) - 1);
					}
					else if ((spy <= 0) && (smy <= 0))
					{
						rho_array(i, j) = rho_array(i, j) - dt*s*(sqrt(POW2(rho_p_x) + POW2(rho_p_y)) - 1);
					}
					else if ((spy > 0) && (smy < 0))
					{
						rho_array(i, j) = rho_array(i, j) - dt*s*(sqrt(POW2(rho_p_x)) - 1);
					}
					else if ((spy < 0) && (smy > 0))
					{
						T ss = s*(abs(rho_p_y) - abs(rho_m_y))/(rho_p_y - rho_m_y);
						if (ss > 0)
						{
							rho_array(i, j) = rho_array(i, j) - dt*s*(sqrt(POW2(rho_p_x) + POW2(rho_m_y)) - 1);
						}
						else
						{
							rho_array(i, j) = rho_array(i, j) - dt*s*(sqrt(POW2(rho_p_x) + POW2(rho_p_y)) - 1);
						}
					}
				}
				else if ((spx > 0) && (smx < 0))
				{
					if ((spy >= 0) && (smy >= 0))
					{
						rho_array(i, j) = rho_array(i, j) - dt*s*(sqrt(POW2(rho_m_y)) - 1);
					}
					else if ((spy <= 0) && (smy <= 0))
					{
						rho_array(i, j) = rho_array(i, j) - dt*s*(sqrt(POW2(rho_p_y)) - 1);
					}
					else if ((spy > 0) && (smy < 0))
					{
						rho_array(i, j) = rho_array(i, j) - dt*s*(- 1);
					}
					else if ((spy < 0) && (smy > 0))
					{
						T ss = s*(abs(rho_p_y) - abs(rho_m_y))/(rho_p_y - rho_m_y);
						if (ss > 0)
						{
							rho_array(i, j) = rho_array(i, j) - dt*s*(sqrt(POW2(rho_m_y)) - 1);
						}
						else
						{
							rho_array(i, j) = rho_array(i, j) - dt*s*(sqrt(POW2(rho_p_y)) - 1);
						}
					}
				}
				else if ((spx < 0) && (smx > 0))
				{
					T sss = s*(abs(rho_p_x) - abs(rho_m_x))/(rho_p_x - rho_m_x);
						
					if (sss > 0)
					{
						if ((spy >= 0) && (smy >= 0))
						{
							rho_array(i, j) = rho_array(i, j) - dt*s*(sqrt(POW2(rho_m_x) + POW2(rho_m_y)) - 1);
						}
						else if ((spy <= 0) && (smy <= 0))
						{
							rho_array(i, j) = rho_array(i, j) - dt*s*(sqrt(POW2(rho_m_x) + POW2(rho_p_y)) - 1);
						}
						else if ((spy > 0) && (smy < 0))
						{
							rho_array(i, j) = rho_array(i, j) - dt*s*(sqrt(POW2(rho_m_x)) - 1);
						}
						else if ((spy < 0) && (smy > 0))
						{
							T ss = s*(abs(rho_p_y) - abs(rho_m_y))/(rho_p_y - rho_m_y);
							if (ss > 0)
							{
								rho_array(i, j) = rho_array(i, j) - dt*s*(sqrt(POW2(rho_m_x) + POW2(rho_m_y)) - 1);
							}
							else
							{
								rho_array(i, j) = rho_array(i, j) - dt*s*(sqrt(POW2(rho_m_x) + POW2(rho_p_y)) - 1);
							}
						}
					}
					else
					{
						if ((spy >= 0) && (smy >= 0))
						{
							rho_array(i, j) = rho_array(i, j) - dt*s*(sqrt(POW2(rho_p_x) + POW2(rho_m_y)) - 1);
						}
						else if ((spy <= 0) && (smy <= 0))
						{
							rho_array(i, j) = rho_array(i, j) - dt*s*(sqrt(POW2(rho_p_x) + POW2(rho_p_y)) - 1);
						}
						else if ((spy > 0) && (smy < 0))
						{
							rho_array(i, j) = rho_array(i, j) - dt*s*(sqrt(POW2(rho_p_x)) - 1);
						}
						else if ((spy < 0) && (smy > 0))
						{
							T ss = s*(abs(rho_p_y) - abs(rho_m_y))/(rho_p_y - rho_m_y);
							if (ss > 0)
							{
								rho_array(i, j) = rho_array(i, j) - dt*s*(sqrt(POW2(rho_p_x) + POW2(rho_m_y)) - 1);
							}
							else
							{
								rho_array(i, j) = rho_array(i, j) - dt*s*(sqrt(POW2(rho_p_x) + POW2(rho_p_y)) - 1);
							}
						}
					}
				}
				if (s*rho_array(i, j) < (T)0)
				{
					if (s < 0)
					{
						rho_array(i, j) = -epsilon;
					}
					if (s > 0)
					{
						rho_array(i, j) = epsilon;
					}
					//rho_array(i, j) = rho_ghost(i, j);
				}
			}	
		}
		
		
		//ARRAY_2D<TT>& rho_array(rho.array_for_this);

		////rho_ghost.FillGhostCellsFrom(rho_array, true);
		//rho_ghost.FillGhostCellsContinuousDerivativesFrom(rho_array, true);

		//T one_over_dx(rho.one_over_dx), one_over_dy(rho.one_over_dy);

		//
		//for (int j = rho.grid.j_start; j <= rho.grid.j_end; j++)
		//{
		//	for (int i = rho.grid.i_start; i <= rho.grid.i_end; i++)
		//	{
		//		// x-components
		//		TT rho_m_x_1 = one_over_dx*(rho_ghost(i - 2, j) - rho_ghost(i - 3, j));
		//		TT rho_m_x_2 = one_over_dx*(rho_ghost(i - 1, j) - rho_ghost(i - 2, j));
		//		TT rho_m_x_3 = one_over_dx*(rho_ghost(i, j) - rho_ghost(i - 1, j));
		//		TT rho_m_x_4 = one_over_dx*(rho_ghost(i + 1, j) - rho_ghost(i, j));
		//		TT rho_m_x_5 = one_over_dx*(rho_ghost(i + 2, j) - rho_ghost(i + 1, j));

		//		TT rho_p_x_1 = one_over_dx*(rho_ghost(i + 3, j) - rho_ghost(i + 2, j));
		//		TT rho_p_x_2 = one_over_dx*(rho_ghost(i + 2, j) - rho_ghost(i + 1, j));
		//		TT rho_p_x_3 = one_over_dx*(rho_ghost(i + 1, j) - rho_ghost(i , j));
		//		TT rho_p_x_4 = one_over_dx*(rho_ghost(i , j) - rho_ghost(i - 1, j));
		//		TT rho_p_x_5 = one_over_dx*(rho_ghost(i - 1, j) - rho_ghost(i - 2, j));

		//		// y-components
		//		TT rho_m_y_1 = one_over_dy*(rho_ghost(i, j - 2) - rho_ghost(i, j - 3));
		//		TT rho_m_y_2 = one_over_dy*(rho_ghost(i, j - 1) - rho_ghost(i, j - 2));
		//		TT rho_m_y_3 = one_over_dy*(rho_ghost(i, j) - rho_ghost(i, j - 1));
		//		TT rho_m_y_4 = one_over_dy*(rho_ghost(i, j + 1) - rho_ghost(i, j));
		//		TT rho_m_y_5 = one_over_dy*(rho_ghost(i, j + 2) - rho_ghost(i, j + 1));

		//		TT rho_p_y_1 = one_over_dy*(rho_ghost(i, j + 3) - rho_ghost(i, j + 2));
		//		TT rho_p_y_2 = one_over_dy*(rho_ghost(i, j + 2) - rho_ghost(i, j + 1));
		//		TT rho_p_y_3 = one_over_dy*(rho_ghost(i, j + 1) - rho_ghost(i, j));
		//		TT rho_p_y_4 = one_over_dy*(rho_ghost(i, j) - rho_ghost(i, j - 1));
		//		TT rho_p_y_5 = one_over_dy*(rho_ghost(i, j - 1) - rho_ghost(i, j - 2));
		//		
		//		// Smoothness
		//		// x-components
		//		TT s_m_x_1 = (T)13/12*(rho_m_x_1 - 2*rho_m_x_2 + rho_m_x_3)*(rho_m_x_1 - 2*rho_m_x_2 + rho_m_x_3) + (T)1/4*(rho_m_x_1 - 4*rho_m_x_2 + 3*rho_m_x_3)*(rho_m_x_1 - 4*rho_m_x_2 + 3*rho_m_x_3);
		//		TT s_m_x_2 = (T)13/12*(rho_m_x_2 - 2*rho_m_x_3 + rho_m_x_4)*(rho_m_x_2 - 2*rho_m_x_3 + rho_m_x_4) + (T)1/4*(rho_m_x_2 - rho_m_x_4)*(rho_m_x_2 - rho_m_x_4);
		//		TT s_m_x_3 = (T)13/12*(rho_m_x_3 - 2*rho_m_x_4 + rho_m_x_5)*(rho_m_x_3 - 2*rho_m_x_4 + rho_m_x_5) + (T)1/4*(3*rho_m_x_3 - 4*rho_m_x_4 + rho_m_x_5)*(3*rho_m_x_3 - 4*rho_m_x_4 + rho_m_x_5);

		//		TT s_p_x_1 = (T)13/12*(rho_p_x_1 - 2*rho_p_x_2 + rho_p_x_3)*(rho_p_x_1 - 2*rho_p_x_2 + rho_p_x_3) + (T)1/4*(rho_p_x_1 - 4*rho_p_x_2 + 3*rho_p_x_3)*(rho_p_x_1 - 4*rho_p_x_2 + 3*rho_p_x_3);
		//		TT s_p_x_2 = (T)13/12*(rho_p_x_2 - 2*rho_p_x_3 + rho_p_x_4)*(rho_p_x_2 - 2*rho_p_x_3 + rho_p_x_4) + (T)1/4*(rho_p_x_2 - rho_p_x_4)*(rho_p_x_2 - rho_p_x_4);
		//		TT s_p_x_3 = (T)13/12*(rho_p_x_3 - 2*rho_p_x_4 + rho_p_x_5)*(rho_p_x_3 - 2*rho_p_x_4 + rho_p_x_5) + (T)1/4*(3*rho_p_x_3 - 4*rho_p_x_4 + rho_p_x_5)*(3*rho_p_x_3 - 4*rho_p_x_4 + rho_p_x_5);

		//		// y-components
		//		TT s_m_y_1 = (T)13/12*(rho_m_y_1 - 2*rho_m_y_2 + rho_m_y_3)*(rho_m_y_1 - 2*rho_m_y_2 + rho_m_y_3) + (T)1/4*(rho_m_y_1 - 4*rho_m_y_2 + 3*rho_m_y_3)*(rho_m_y_1 - 4*rho_m_y_2 + 3*rho_m_y_3);
		//		TT s_m_y_2 = (T)13/12*(rho_m_y_2 - 2*rho_m_y_3 + rho_m_y_4)*(rho_m_y_2 - 2*rho_m_y_3 + rho_m_y_4) + (T)1/4*(rho_m_y_2 - rho_m_y_4)*(rho_m_y_2 - rho_m_y_4);
		//		TT s_m_y_3 = (T)13/12*(rho_m_y_3 - 2*rho_m_y_4 + rho_m_y_5)*(rho_m_y_3 - 2*rho_m_y_4 + rho_m_y_5) + (T)1/4*(3*rho_m_y_3 - 4*rho_m_y_4 + rho_m_y_5)*(3*rho_m_y_3 - 4*rho_m_y_4 + rho_m_y_5);

		//		TT s_p_y_1 = (T)13/12*(rho_p_y_1 - 2*rho_p_y_2 + rho_p_y_3)*(rho_p_y_1 - 2*rho_p_y_2 + rho_p_y_3) + (T)1/4*(rho_p_y_1 - 4*rho_p_y_2 + 3*rho_p_y_3)*(rho_p_y_1 - 4*rho_p_y_2 + 3*rho_p_y_3);
		//		TT s_p_y_2 = (T)13/12*(rho_p_y_2 - 2*rho_p_y_3 + rho_p_y_4)*(rho_p_y_2 - 2*rho_p_y_3 + rho_p_y_4) + (T)1/4*(rho_p_y_2 - rho_p_y_4)*(rho_p_y_2 - rho_p_y_4);
		//		TT s_p_y_3 = (T)13/12*(rho_p_y_3 - 2*rho_p_y_4 + rho_p_y_5)*(rho_p_y_3 - 2*rho_p_y_4 + rho_p_y_5) + (T)1/4*(3*rho_p_y_3 - 4*rho_p_y_4 + rho_p_y_5)*(3*rho_p_y_3 - 4*rho_p_y_4 + rho_p_y_5);

		//		// Weights
		//		// x-componets
		//		TT a_m_x_1 = (T)1/10*(T)1/((epsilon + s_m_x_1)*(epsilon + s_m_x_1));
		//		TT a_m_x_2 = (T)6/10*(T)1/((epsilon + s_m_x_2)*(epsilon + s_m_x_2));
		//		TT a_m_x_3 = (T)3/10*(T)1/((epsilon + s_m_x_3)*(epsilon + s_m_x_3));

		//		TT sum_m_x = a_m_x_1 + a_m_x_2 + a_m_x_3;

		//		TT w_m_x_1 = a_m_x_1/sum_m_x;
		//		TT w_m_x_2 = a_m_x_2/sum_m_x;
		//		TT w_m_x_3 = a_m_x_3/sum_m_x;
	
		//		TT a_p_x_1 = (T)1/10*(T)1/((epsilon + s_p_x_1)*(epsilon + s_p_x_1));
		//		TT a_p_x_2 = (T)6/10*(T)1/((epsilon + s_p_x_2)*(epsilon + s_p_x_2));
		//		TT a_p_x_3 = (T)3/10*(T)1/((epsilon + s_p_x_3)*(epsilon + s_p_x_3));

		//		TT sum_p_x = a_p_x_1 + a_p_x_2 + a_p_x_3;
	
		//		TT w_p_x_1 = a_p_x_1/sum_p_x;
		//		TT w_p_x_2 = a_p_x_2/sum_p_x;
		//		TT w_p_x_3 = a_p_x_3/sum_p_x;

		//		// y-component
		//		TT a_m_y_1 = (T)1/10*(T)1/((epsilon + s_m_y_1)*(epsilon + s_m_y_1));
		//		TT a_m_y_2 = (T)6/10*(T)1/((epsilon + s_m_y_2)*(epsilon + s_m_y_2));
		//		TT a_m_y_3 = (T)3/10*(T)1/((epsilon + s_m_y_3)*(epsilon + s_m_y_3));

		//		TT sum_m_y = a_m_y_1 + a_m_y_2 + a_m_y_3;
	
		//		TT w_m_y_1 = a_m_y_1/sum_m_y;
		//		TT w_m_y_2 = a_m_y_2/sum_m_y;
		//		TT w_m_y_3 = a_m_y_3/sum_m_y;
	
		//		TT a_p_y_1 = (T)1/10*(T)1/((epsilon + s_p_y_1)*(epsilon + s_p_y_1));
		//		TT a_p_y_2 = (T)6/10*(T)1/((epsilon + s_p_y_2)*(epsilon + s_p_y_2));
		//		TT a_p_y_3 = (T)3/10*(T)1/((epsilon + s_p_y_3)*(epsilon + s_p_y_3));

		//		TT sum_p_y = a_p_y_1 + a_p_y_2 + a_p_y_3;

		//		TT w_p_y_1 = a_p_y_1/sum_p_y;
		//		TT w_p_y_2 = a_p_y_2/sum_p_y;
		//		TT w_p_y_3 = a_p_y_3/sum_p_y;

		//		// Approximation of derivatives
		//		// rho_x
		//		T rho_m_x = w_m_x_1*(rho_m_x_1*(T)1/3 - rho_m_x_2*(T)7/6 + rho_m_x_3*(T)11/6) + w_m_x_2*(rho_m_x_2*(-(T)1/6) + rho_m_x_3*(T)5/6 + rho_m_x_4*(T)1/3) + w_m_x_3*(rho_m_x_3*(T)1/3 + rho_m_x_4*(T)5/6 - rho_m_x_5*(T)1/6);
		//		T rho_p_x = w_p_x_1*(rho_p_x_1*(T)1/3 - rho_p_x_2*(T)7/6 + rho_p_x_3*(T)11/6) + w_p_x_2*(rho_p_x_2*(-(T)1/6) + rho_p_x_3*(T)5/6 + rho_p_x_4*(T)1/3) + w_p_x_3*(rho_p_x_3*(T)1/3 + rho_p_x_4*(T)5/6 - rho_p_x_5*(T)1/6);
	
		//		// rho_y
		//		T rho_m_y = w_m_y_1*(rho_m_y_1*(T)1/3 - rho_m_y_2*(T)7/6 + rho_m_y_3*(T)11/6) + w_m_y_2*(rho_m_y_2*(-(T)1/6) + rho_m_y_3*(T)5/6 + rho_m_y_4*(T)1/3) + w_m_y_3*(rho_m_y_3*(T)1/3 + rho_m_y_4*(T)5/6 - rho_m_y_5*(T)1/6);
		//		T rho_p_y = w_p_y_1*(rho_p_y_1*(T)1/3 - rho_p_y_2*(T)7/6 + rho_p_y_3*(T)11/6) + w_p_y_2*(rho_p_y_2*(-(T)1/6) + rho_p_y_3*(T)5/6 + rho_p_y_4*(T)1/3) + w_p_y_3*(rho_p_y_3*(T)1/3 + rho_p_y_4*(T)5/6 - rho_p_y_5*(T)1/6);
		//	
		//		// Define the sign function
		//		TT s = sign_function(i, j);

		//		// Define Sub functions
		//		TT smx, spx, smy, spy;
		//		
		//		smx = s*rho_m_x;
		//		spx = s*rho_p_x;
		//		smy = s*rho_m_y;
		//		spy = s*rho_p_y;
		//		
		//		TT rho_x, rho_y;

		//		if ((spx >= 0) && (smx >= 0))
		//		{
		//			rho_x = rho_m_x;
		//		}
		//		if ((spx <= 0) && (smx <= 0))
		//		{
		//			rho_x = rho_p_x;
		//		}
		//		if ((spx > 0) && (smx < 0))
		//		{
		//			rho_x = 0;
		//		}
		//		if ((spx < 0) && (smx > 0))
		//		{
		//			TT ss = s*(abs(rho_p_x) - abs(rho_m_x))/(rho_p_x - rho_m_x);
		//			if (ss > 0)
		//			{
		//				rho_x = rho_m_x;
		//			}
		//			else
		//			{
		//				rho_x = rho_p_x;
		//			}
		//		}
		//		
		//		if ((spy >= 0) && (smy >= 0))
		//		{
		//			rho_y = rho_m_y;
		//		}
		//		if ((spy <= 0) && (smy <= 0))
		//		{
		//			rho_y = rho_p_y;
		//		}
		//		if ((spy > 0) && (smy < 0))
		//		{
		//			rho_y = 0;
		//		}
		//		if ((spy < 0) && (smy > 0))
		//		{
		//			TT ss = s*(abs(rho_p_y) - abs(rho_m_y))/(rho_p_y - rho_m_y);
		//			if (ss > 0)
		//			{
		//				rho_y = rho_m_y;
		//			}
		//			else
		//			{
		//				rho_y = rho_p_y;
		//			}
		//		}

		//		rho_array(i, j) = rho_array(i, j) - dt*s*(sqrt(POW2(rho_x) + POW2(rho_y)) - 1);
		//	}	
		//}
	}

	static void GodunovReinitialization(FIELD_STRUCTURE_2D<TT>& rho, FIELD_STRUCTURE_2D<TT>& rho_ghost, const T& dt, const T& epsilon, const FIELD_STRUCTURE_2D<T>& sign_function, const FIELD_STRUCTURE_2D<T>& phi_0)
	{
		ARRAY_2D<TT>& rho_array(rho.array_for_this);

		rho_ghost.FillGhostCellsContinuousDerivativesFrom(rho_array, true);

		T one_over_dx(rho.one_over_dx), one_over_dy(rho.one_over_dy);

		GRID_ITERATION_2D(rho.grid)
		{
			TT rho_m_x = one_over_dx*(rho_ghost(i, j) - rho_ghost(i - 1, j));
			TT rho_p_x = one_over_dx*(rho_ghost(i + 1, j) - rho_ghost(i, j));
			TT rho_m_y = one_over_dy*(rho_ghost(i, j) - rho_ghost(i, j - 1));
			TT rho_p_y = one_over_dy*(rho_ghost(i, j + 1) - rho_ghost(i, j));

			TT a_p, b_m, c_p, d_m, a_m, b_p, c_m, d_p;

			a_m = min(rho_m_x, 0);
			a_p = max(rho_m_x, 0);
			b_m = min(rho_p_x, 0);
			b_p = max(rho_p_x, 0);
			c_m = min(rho_m_y, 0);
			c_p = max(rho_m_y, 0);
			d_m = min(rho_p_y, 0);
			d_p = max(rho_p_y, 0);

			// Define the sign function
			TT s = sign_function(i, j);

			if (s > 0)
			{
				rho_array(i, j) = rho_ghost(i, j) - dt*s*(sqrt(max(POW2(a_p), POW2(b_m)) + max(POW2(c_p), POW2(d_m))) - (T)1);
				
				if (rho_array(i, j) < 0)
				{
					rho_array(i, j) = rho_ghost(i, j);
				}
			}
			else if (s < 0)
			{
				rho_array(i, j) = rho_ghost(i, j) - dt*s*(sqrt(max(POW2(a_m), POW2(b_p)) + max(POW2(c_m), POW2(d_p))) - (T)1);
				
				if (rho_array(i, j) > 0)
				{
					rho_array(i, j) = rho_ghost(i, j);
				}
			}
			else
			{
				rho_array(i, j) = rho_ghost(i, j);
			}
			
		}
	}

	static void SubcellFixedReinitialization(FIELD_STRUCTURE_2D<TT>& rho, FIELD_STRUCTURE_2D<TT>& rho_ghost, const T& dt, const T& epsilon, const FIELD_STRUCTURE_2D<T>& sign_function, const FIELD_STRUCTURE_2D<T>& phi_0)
	{
		ARRAY_2D<TT>& rho_array(rho.array_for_this);

		rho_ghost.FillGhostCellsContinuousDerivativesFrom(rho_array, true);

		T one_over_dx(rho.one_over_dx), one_over_dy(rho.one_over_dy);

		FIELD_STRUCTURE_2D<TT> D;
		D.Initialize(rho.grid, 2);

		GRID_ITERATION_2D(rho.grid)
		{
			T central = sqrt(POW2(phi_0(i + 1, j) - phi_0(i - 1, j)) + POW2(phi_0(i, j + 1) - phi_0(i, j - 1)))/(T)2;
			T one_side_p = sqrt(POW2(phi_0(i + 1, j) - phi_0(i, j)) + POW2(phi_0(i, j + 1) - phi_0(i, j)));
			T one_side_m = sqrt(POW2(phi_0(i, j) - phi_0(i - 1, j)) + POW2(phi_0(i, j) - phi_0(i, j - 1)));
			T max_value = max(max(central, one_side_p), max(one_side_m, (T)1e-6));
			D(i, j) = rho.grid.dx*phi_0(i, j)/max_value;
		}
		
		GRID_ITERATION_2D(rho.grid)
		{
			if (phi_0(i, j)*phi_0(i - 1, j) < 0 || phi_0(i, j)*phi_0(i + 1, j) < 0 || phi_0(i, j)*phi_0(i, j - 1) < 0 || phi_0(i, j)*phi_0(i, j + 1) < 0)
			{
				if (phi_0(i, j) < (T)0)
				{
					rho_array(i, j) = rho_ghost(i, j) - dt*one_over_dx*(-abs(rho_ghost(i, j)) - D(i, j));
					
					if (rho_array(i, j) > (T)0)
					{
						rho_array(i, j) = rho_ghost(i, j);
					}
				}
				else
				{
					rho_array(i, j) = rho_ghost(i, j) - dt*one_over_dx*(abs(rho_ghost(i, j)) - D(i, j));
					
					if (rho_array(i, j) < (T)0)
					{
						rho_array(i, j) = rho_ghost(i, j);
					}
				}
			}
			else
			{
				TT rho_m_x = one_over_dx*(rho_ghost(i, j) - rho_ghost(i - 1, j));
				TT rho_p_x = one_over_dx*(rho_ghost(i + 1, j) - rho_ghost(i, j));
				TT rho_m_y = one_over_dy*(rho_ghost(i, j) - rho_ghost(i, j - 1));
				TT rho_p_y = one_over_dy*(rho_ghost(i, j + 1) - rho_ghost(i, j));

				TT a_p, b_m, c_p, d_m, a_m, b_p, c_m, d_p;

				a_p = max(rho_m_x, 0);
				b_m = min(rho_p_x, 0);
				c_p = max(rho_m_y, 0);
				d_m = min(rho_p_y, 0);
				a_m = min(rho_m_x, 0);
				b_p = max(rho_p_x, 0);
				c_m = min(rho_m_y, 0);
				d_p = max(rho_p_y, 0);

				// Define the sign function
				
				T s = phi_0(i, j);
				
				if (s > 0)
				{
					rho_array(i, j) = rho_ghost(i, j) - dt*(sqrt(max(POW2(a_p), POW2(b_m)) + max(POW2(c_p), POW2(d_m))) - (T)1);
					
					if (rho_array(i, j) < 0)
					{
						rho_array(i, j) = rho_ghost(i, j);
					}
				}
				else if (s < 0)
				{
					rho_array(i, j) = rho_ghost(i, j) + dt*(sqrt(max(POW2(a_m), POW2(b_p)) + max(POW2(c_m), POW2(d_p))) - (T)1);

					if (rho_array(i, j) > 0)
					{
						rho_array(i, j) = rho_ghost(i, j);
					}
				}
				else
				{
					rho_array(i, j) = rho_ghost(i, j);
				}
			}
		}
	}

	static void WENO5thReinitialization(FIELD_STRUCTURE_2D<TT>& rho, FIELD_STRUCTURE_2D<TT>& rho_ghost, const T& dt, const T& epsilon, const FIELD_STRUCTURE_2D<T>& sign_function, MULTITHREADING& multithreading, const int& thread_id)
	{
		ARRAY_2D<TT>& rho_array(rho.array_for_this);

		rho_ghost.FillGhostCellsFrom(rho_array, true, thread_id);

		multithreading.Sync(thread_id);

		T one_over_dx(rho.one_over_dx), one_over_dy(rho.one_over_dy);

		GRID_ITERATION_2D(rho.partial_grids[thread_id])
		{
			// x-components
			TT rho_m_x_1 = one_over_dx*(rho_ghost(i - 2, j) - rho_ghost(i - 3, j));
			TT rho_m_x_2 = one_over_dx*(rho_ghost(i - 1, j) - rho_ghost(i - 2, j));
			TT rho_m_x_3 = one_over_dx*(rho_ghost(i, j) - rho_ghost(i - 1, j));
			TT rho_m_x_4 = one_over_dx*(rho_ghost(i + 1, j) - rho_ghost(i, j));
			TT rho_m_x_5 = one_over_dx*(rho_ghost(i + 2, j) - rho_ghost(i + 1, j));

			TT rho_p_x_1 = one_over_dx*(rho_ghost(i + 3, j) - rho_ghost(i + 2, j));
			TT rho_p_x_2 = one_over_dx*(rho_ghost(i + 2, j) - rho_ghost(i + 1, j));
			TT rho_p_x_3 = one_over_dx*(rho_ghost(i + 1, j) - rho_ghost(i , j));
			TT rho_p_x_4 = one_over_dx*(rho_ghost(i , j) - rho_ghost(i - 1, j));
			TT rho_p_x_5 = one_over_dx*(rho_ghost(i - 1, j) - rho_ghost(i - 2, j));

			// y-components
			TT rho_m_y_1 = one_over_dy*(rho_ghost(i, j - 2) - rho_ghost(i, j - 3));
			TT rho_m_y_2 = one_over_dy*(rho_ghost(i, j - 1) - rho_ghost(i, j - 2));
			TT rho_m_y_3 = one_over_dy*(rho_ghost(i, j) - rho_ghost(i, j - 1));
			TT rho_m_y_4 = one_over_dy*(rho_ghost(i, j + 1) - rho_ghost(i, j));
			TT rho_m_y_5 = one_over_dy*(rho_ghost(i, j + 2) - rho_ghost(i, j + 1));

			TT rho_p_y_1 = one_over_dy*(rho_ghost(i, j + 3) - rho_ghost(i, j + 2));
			TT rho_p_y_2 = one_over_dy*(rho_ghost(i, j + 2) - rho_ghost(i, j + 1));
			TT rho_p_y_3 = one_over_dy*(rho_ghost(i, j + 1) - rho_ghost(i, j));
			TT rho_p_y_4 = one_over_dy*(rho_ghost(i, j) - rho_ghost(i, j - 1));
			TT rho_p_y_5 = one_over_dy*(rho_ghost(i, j - 1) - rho_ghost(i, j - 2));
				
			// Smoothness
			// x-components
			TT s_m_x_1 = (T)13/12*(rho_m_x_1 - 2*rho_m_x_2 + rho_m_x_3)*(rho_m_x_1 - 2*rho_m_x_2 + rho_m_x_3) + (T)1/4*(rho_m_x_1 - 4*rho_m_x_2 + 3*rho_m_x_3)*(rho_m_x_1 - 4*rho_m_x_2 + 3*rho_m_x_3);
			TT s_m_x_2 = (T)13/12*(rho_m_x_2 - 2*rho_m_x_3 + rho_m_x_4)*(rho_m_x_2 - 2*rho_m_x_3 + rho_m_x_4) + (T)1/4*(rho_m_x_2 - rho_m_x_4)*(rho_m_x_2 - rho_m_x_4);
			TT s_m_x_3 = (T)13/12*(rho_m_x_3 - 2*rho_m_x_4 + rho_m_x_5)*(rho_m_x_3 - 2*rho_m_x_4 + rho_m_x_5) + (T)1/4*(3*rho_m_x_3 - 4*rho_m_x_4 + rho_m_x_5)*(3*rho_m_x_3 - 4*rho_m_x_4 + rho_m_x_5);

			TT s_p_x_1 = (T)13/12*(rho_p_x_1 - 2*rho_p_x_2 + rho_p_x_3)*(rho_p_x_1 - 2*rho_p_x_2 + rho_p_x_3) + (T)1/4*(rho_p_x_1 - 4*rho_p_x_2 + 3*rho_p_x_3)*(rho_p_x_1 - 4*rho_p_x_2 + 3*rho_p_x_3);
			TT s_p_x_2 = (T)13/12*(rho_p_x_2 - 2*rho_p_x_3 + rho_p_x_4)*(rho_p_x_2 - 2*rho_p_x_3 + rho_p_x_4) + (T)1/4*(rho_p_x_2 - rho_p_x_4)*(rho_p_x_2 - rho_p_x_4);
			TT s_p_x_3 = (T)13/12*(rho_p_x_3 - 2*rho_p_x_4 + rho_p_x_5)*(rho_p_x_3 - 2*rho_p_x_4 + rho_p_x_5) + (T)1/4*(3*rho_p_x_3 - 4*rho_p_x_4 + rho_p_x_5)*(3*rho_p_x_3 - 4*rho_p_x_4 + rho_p_x_5);

			// y-components
			TT s_m_y_1 = (T)13/12*(rho_m_y_1 - 2*rho_m_y_2 + rho_m_y_3)*(rho_m_y_1 - 2*rho_m_y_2 + rho_m_y_3) + (T)1/4*(rho_m_y_1 - 4*rho_m_y_2 + 3*rho_m_y_3)*(rho_m_y_1 - 4*rho_m_y_2 + 3*rho_m_y_3);
			TT s_m_y_2 = (T)13/12*(rho_m_y_2 - 2*rho_m_y_3 + rho_m_y_4)*(rho_m_y_2 - 2*rho_m_y_3 + rho_m_y_4) + (T)1/4*(rho_m_y_2 - rho_m_y_4)*(rho_m_y_2 - rho_m_y_4);
			TT s_m_y_3 = (T)13/12*(rho_m_y_3 - 2*rho_m_y_4 + rho_m_y_5)*(rho_m_y_3 - 2*rho_m_y_4 + rho_m_y_5) + (T)1/4*(3*rho_m_y_3 - 4*rho_m_y_4 + rho_m_y_5)*(3*rho_m_y_3 - 4*rho_m_y_4 + rho_m_y_5);

			TT s_p_y_1 = (T)13/12*(rho_p_y_1 - 2*rho_p_y_2 + rho_p_y_3)*(rho_p_y_1 - 2*rho_p_y_2 + rho_p_y_3) + (T)1/4*(rho_p_y_1 - 4*rho_p_y_2 + 3*rho_p_y_3)*(rho_p_y_1 - 4*rho_p_y_2 + 3*rho_p_y_3);
			TT s_p_y_2 = (T)13/12*(rho_p_y_2 - 2*rho_p_y_3 + rho_p_y_4)*(rho_p_y_2 - 2*rho_p_y_3 + rho_p_y_4) + (T)1/4*(rho_p_y_2 - rho_p_y_4)*(rho_p_y_2 - rho_p_y_4);
			TT s_p_y_3 = (T)13/12*(rho_p_y_3 - 2*rho_p_y_4 + rho_p_y_5)*(rho_p_y_3 - 2*rho_p_y_4 + rho_p_y_5) + (T)1/4*(3*rho_p_y_3 - 4*rho_p_y_4 + rho_p_y_5)*(3*rho_p_y_3 - 4*rho_p_y_4 + rho_p_y_5);

			// Weights
			// x-componets
			TT a_m_x_1 = (T)1/10*(T)1/((epsilon + s_m_x_1)*(epsilon + s_m_x_1));
			TT a_m_x_2 = (T)6/10*(T)1/((epsilon + s_m_x_2)*(epsilon + s_m_x_2));
			TT a_m_x_3 = (T)3/10*(T)1/((epsilon + s_m_x_3)*(epsilon + s_m_x_3));

			TT sum_m_x = a_m_x_1 + a_m_x_2 + a_m_x_3;

			TT w_m_x_1 = a_m_x_1/sum_m_x;
			TT w_m_x_2 = a_m_x_2/sum_m_x;
			TT w_m_x_3 = a_m_x_3/sum_m_x;
	
			TT a_p_x_1 = (T)1/10*(T)1/((epsilon + s_p_x_1)*(epsilon + s_p_x_1));
			TT a_p_x_2 = (T)6/10*(T)1/((epsilon + s_p_x_2)*(epsilon + s_p_x_2));
			TT a_p_x_3 = (T)3/10*(T)1/((epsilon + s_p_x_3)*(epsilon + s_p_x_3));

			TT sum_p_x = a_p_x_1 + a_p_x_2 + a_p_x_3;
	
			TT w_p_x_1 = a_p_x_1/sum_p_x;
			TT w_p_x_2 = a_p_x_2/sum_p_x;
			TT w_p_x_3 = a_p_x_3/sum_p_x;

			// y-component
			TT a_m_y_1 = (T)1/10*(T)1/((epsilon + s_m_y_1)*(epsilon + s_m_y_1));
			TT a_m_y_2 = (T)6/10*(T)1/((epsilon + s_m_y_2)*(epsilon + s_m_y_2));
			TT a_m_y_3 = (T)3/10*(T)1/((epsilon + s_m_y_3)*(epsilon + s_m_y_3));

			TT sum_m_y = a_m_y_1 + a_m_y_2 + a_m_y_3;
	
			TT w_m_y_1 = a_m_y_1/sum_m_y;
			TT w_m_y_2 = a_m_y_2/sum_m_y;
			TT w_m_y_3 = a_m_y_3/sum_m_y;
	
			TT a_p_y_1 = (T)1/10*(T)1/((epsilon + s_p_y_1)*(epsilon + s_p_y_1));
			TT a_p_y_2 = (T)6/10*(T)1/((epsilon + s_p_y_2)*(epsilon + s_p_y_2));
			TT a_p_y_3 = (T)3/10*(T)1/((epsilon + s_p_y_3)*(epsilon + s_p_y_3));

			TT sum_p_y = a_p_y_1 + a_p_y_2 + a_p_y_3;

			TT w_p_y_1 = a_p_y_1/sum_p_y;
			TT w_p_y_2 = a_p_y_2/sum_p_y;
			TT w_p_y_3 = a_p_y_3/sum_p_y;

			// Approximation of derivatives
			// rho_x
			T rho_m_x = w_m_x_1*(rho_m_x_1*(T)1/3 - rho_m_x_2*(T)7/6 + rho_m_x_3*(T)11/6) + w_m_x_2*(rho_m_x_2*(-(T)1/6) + rho_m_x_3*(T)5/6 + rho_m_x_4*(T)1/3) + w_m_x_3*(rho_m_x_3*(T)1/3 + rho_m_x_4*(T)5/6 - rho_m_x_5*(T)1/6);
			T rho_p_x = w_p_x_1*(rho_p_x_1*(T)1/3 - rho_p_x_2*(T)7/6 + rho_p_x_3*(T)11/6) + w_p_x_2*(rho_p_x_2*(-(T)1/6) + rho_p_x_3*(T)5/6 + rho_p_x_4*(T)1/3) + w_p_x_3*(rho_p_x_3*(T)1/3 + rho_p_x_4*(T)5/6 - rho_p_x_5*(T)1/6);
	
			// rho_y
			T rho_m_y = w_m_y_1*(rho_m_y_1*(T)1/3 - rho_m_y_2*(T)7/6 + rho_m_y_3*(T)11/6) + w_m_y_2*(rho_m_y_2*(-(T)1/6) + rho_m_y_3*(T)5/6 + rho_m_y_4*(T)1/3) + w_m_y_3*(rho_m_y_3*(T)1/3 + rho_m_y_4*(T)5/6 - rho_m_y_5*(T)1/6);
			T rho_p_y = w_p_y_1*(rho_p_y_1*(T)1/3 - rho_p_y_2*(T)7/6 + rho_p_y_3*(T)11/6) + w_p_y_2*(rho_p_y_2*(-(T)1/6) + rho_p_y_3*(T)5/6 + rho_p_y_4*(T)1/3) + w_p_y_3*(rho_p_y_3*(T)1/3 + rho_p_y_4*(T)5/6 - rho_p_y_5*(T)1/6);
			
			// Define the sign function
			T s = sign_function(i, j);

			// Define Sub functions
			T smx = s*rho_m_x, spx = s*rho_p_x, smy = s*rho_m_y, spy = s*rho_p_y;

			if ((spx >= 0) && (smx >= 0))
			{
				if ((spy >= 0) && (smy >= 0))
				{
					rho_array(i, j) = rho_array(i, j) - dt*s*(sqrt(POW2(rho_m_x) + POW2(rho_m_y)) - 1);
				}
				else if ((spy <= 0) && (smy <= 0))
				{
					rho_array(i, j) = rho_array(i, j) - dt*s*(sqrt(POW2(rho_m_x) + POW2(rho_p_y)) - 1);
				}
				else if ((spy > 0) && (smy < 0))
				{
					rho_array(i, j) = rho_array(i, j) - dt*s*(sqrt(POW2(rho_m_x)) - 1);
				}
				else if ((spy < 0) && (smy > 0))
				{
					T ss = s*(abs(rho_p_y) - abs(rho_m_y))/(rho_p_y - rho_m_y);
					if (ss > 0)
					{
						rho_array(i, j) = rho_array(i, j) - dt*s*(sqrt(POW2(rho_m_x) + POW2(rho_m_y)) - 1);
					}
					else
					{
						rho_array(i, j) = rho_array(i, j) - dt*s*(sqrt(POW2(rho_m_x) + POW2(rho_p_y)) - 1);
					}
				}
			}
			else if ((spx <= 0) && (smx <= 0))
			{
				if ((spy >= 0) && (smy >= 0))
				{
					rho_array(i, j) = rho_array(i, j) - dt*s*(sqrt(POW2(rho_p_x) + POW2(rho_m_y)) - 1);
				}
				else if ((spy <= 0) && (smy <= 0))
				{
					rho_array(i, j) = rho_array(i, j) - dt*s*(sqrt(POW2(rho_p_x) + POW2(rho_p_y)) - 1);
				}
				else if ((spy > 0) && (smy < 0))
				{
					rho_array(i, j) = rho_array(i, j) - dt*s*(sqrt(POW2(rho_p_x)) - 1);
				}
				else if ((spy < 0) && (smy > 0))
				{
					T ss = s*(abs(rho_p_y) - abs(rho_m_y))/(rho_p_y - rho_m_y);
					if (ss > 0)
					{
						rho_array(i, j) = rho_array(i, j) - dt*s*(sqrt(POW2(rho_p_x) + POW2(rho_m_y)) - 1);
					}
					else
					{
						rho_array(i, j) = rho_array(i, j) - dt*s*(sqrt(POW2(rho_p_x) + POW2(rho_p_y)) - 1);
					}
				}
			}
			else if ((spx > 0) && (smx < 0))
			{
				if ((spy >= 0) && (smy >= 0))
				{
					rho_array(i, j) = rho_array(i, j) - dt*s*(sqrt(POW2(rho_m_y)) - 1);
				}
				else if ((spy <= 0) && (smy <= 0))
				{
					rho_array(i, j) = rho_array(i, j) - dt*s*(sqrt(POW2(rho_p_y)) - 1);
				}
				else if ((spy > 0) && (smy < 0))
				{
					rho_array(i, j) = rho_array(i, j) - dt*s*(- 1);
				}
				else if ((spy < 0) && (smy > 0))
				{
					T ss = s*(abs(rho_p_y) - abs(rho_m_y))/(rho_p_y - rho_m_y);
					if (ss > 0)
					{
						rho_array(i, j) = rho_array(i, j) - dt*s*(sqrt(POW2(rho_m_y)) - 1);
					}
					else
					{
						rho_array(i, j) = rho_array(i, j) - dt*s*(sqrt(POW2(rho_p_y)) - 1);
					}
				}
			}
			else if ((spx < 0) && (smx > 0))
			{
				T sss = s*(abs(rho_p_x) - abs(rho_m_x))/(rho_p_x - rho_m_x);
					
				if (sss > 0)
				{
					if ((spy >= 0) && (smy >= 0))
					{
						rho_array(i, j) = rho_array(i, j) - dt*s*(sqrt(POW2(rho_m_x) + POW2(rho_m_y)) - 1);
					}
					else if ((spy <= 0) && (smy <= 0))
					{
						rho_array(i, j) = rho_array(i, j) - dt*s*(sqrt(POW2(rho_m_x) + POW2(rho_p_y)) - 1);
					}
					else if ((spy > 0) && (smy < 0))
					{
						rho_array(i, j) = rho_array(i, j) - dt*s*(sqrt(POW2(rho_m_x)) - 1);
					}
					else if ((spy < 0) && (smy > 0))
					{
						T ss = s*(abs(rho_p_y) - abs(rho_m_y))/(rho_p_y - rho_m_y);
						if (ss > 0)
						{
							rho_array(i, j) = rho_array(i, j) - dt*s*(sqrt(POW2(rho_m_x) + POW2(rho_m_y)) - 1);
						}
						else
						{
							rho_array(i, j) = rho_array(i, j) - dt*s*(sqrt(POW2(rho_m_x) + POW2(rho_p_y)) - 1);
						}
					}
				}
				else if (sss < 0)
				{
					if ((spy >= 0) && (smy >= 0))
					{
						rho_array(i, j) = rho_array(i, j) - dt*s*(sqrt(POW2(rho_p_x) + POW2(rho_m_y)) - 1);
					}
					else if ((spy <= 0) && (smy <= 0))
					{
						rho_array(i, j) = rho_array(i, j) - dt*s*(sqrt(POW2(rho_p_x) + POW2(rho_p_y)) - 1);
					}
					else if ((spy > 0) && (smy < 0))
					{
						rho_array(i, j) = rho_array(i, j) - dt*s*(sqrt(POW2(rho_p_x)) - 1);
					}
					else if ((spy < 0) && (smy > 0))
					{
						T ss = s*(abs(rho_p_y) - abs(rho_m_y))/(rho_p_y - rho_m_y);
						if (ss > 0)
						{
							rho_array(i, j) = rho_array(i, j) - dt*s*(sqrt(POW2(rho_p_x) + POW2(rho_m_y)) - 1);
						}
						else
						{
							rho_array(i, j) = rho_array(i, j) - dt*s*(sqrt(POW2(rho_p_x) + POW2(rho_p_y)) - 1);
						}
					}
				}
			}
		}
		multithreading.Sync(thread_id);
	}

	// ENO3rd Method
	static void ENO3rd(FIELD_STRUCTURE_2D<TT>& rho, FIELD_STRUCTURE_2D<TT>& rho_ghost, const FIELD_STRUCTURE_2D<VT>& velocity, const T& dt)
	{
		const ARRAY_2D<VT>& velocity_array(velocity.array_for_this);
		ARRAY_2D<TT>& rho_array(rho.array_for_this);

		rho_ghost.FillGhostCellsFrom(rho_array, true);

		T dx(rho.dx), dy(rho.dy);
		T one_over_dx(rho.one_over_dx), one_over_dy(rho.one_over_dy), one_over_2dx(rho.one_over_2dx), one_over_2dy(rho.one_over_2dy);
		T one_over_3dx(one_over_dx*(T)1/3), one_over_3dy(one_over_dy*(T)1/3);

		for (int j = rho.grid.j_start; j <= rho.grid.j_end; j++)
		{
			for (int i = rho.grid.i_start; i <= rho.grid.i_end; i++)
			{
				// x-component
				TT diff_1_x_n_3 = one_over_dx*(rho_ghost(i - 2, j) - rho_ghost(i - 3, j));
				TT diff_1_x_n_2 = one_over_dx*(rho_ghost(i - 1, j) - rho_ghost(i - 2, j));
				TT diff_1_x_n_1 = one_over_dx*(rho_ghost(i,     j) - rho_ghost(i - 1, j));
				TT diff_1_x_0   = one_over_dx*(rho_ghost(i + 1, j) - rho_ghost(i    , j));
				TT diff_1_x_p_1 = one_over_dx*(rho_ghost(i + 2, j) - rho_ghost(i + 1, j));
				TT diff_1_x_p_2 = one_over_dx*(rho_ghost(i + 3, j) - rho_ghost(i + 2, j));

				TT diff_2_x_n_2 = one_over_2dx*(diff_1_x_n_2 - diff_1_x_n_3);
				TT diff_2_x_n_1 = one_over_2dx*(diff_1_x_n_1 - diff_1_x_n_2);
				TT diff_2_x_0   = one_over_2dx*(diff_1_x_0   - diff_1_x_n_1);
				TT diff_2_x_p_1 = one_over_2dx*(diff_1_x_p_1 - diff_1_x_0  );
				TT diff_2_x_p_2 = one_over_2dx*(diff_1_x_p_2 - diff_1_x_p_1);

				TT diff_3_x_n_2 = one_over_3dx*(diff_2_x_n_1 - diff_2_x_n_2);
				TT diff_3_x_n_1 = one_over_3dx*(diff_2_x_0   - diff_2_x_n_1);
				TT diff_3_x_0   = one_over_3dx*(diff_2_x_p_1 - diff_2_x_0  );
				TT diff_3_x_p_1 = one_over_3dx*(diff_2_x_p_2 - diff_2_x_p_1);

				TT rho_m_x;
				if (abs(diff_2_x_n_1) <= abs(diff_2_x_0))
				{
					if (abs(diff_3_x_n_2) <= abs(diff_3_x_n_1))
					{
						rho_m_x = diff_1_x_n_1 + diff_2_x_n_1*dx + 2*diff_3_x_n_2*dx*dx;
					}
					else
					{
						rho_m_x = diff_1_x_n_1 + diff_2_x_n_1*dx + 2*diff_3_x_n_1*dx*dx;
					}
				}
				else
				{
					if (abs(diff_3_x_n_1) <= abs(diff_3_x_0))
					{
						rho_m_x = diff_1_x_n_1 + diff_2_x_0*dx - diff_3_x_n_1*dx*dx;
					}
					else
					{
						rho_m_x = diff_1_x_n_1 + diff_2_x_0*dx - diff_3_x_0*dx*dx;
					}
				}

				TT rho_p_x;
				if (abs(diff_2_x_0) <= abs(diff_2_x_p_1))
				{
					if (abs(diff_3_x_n_1) <= abs(diff_3_x_0))
					{
						rho_p_x = diff_1_x_0 - diff_2_x_0*dx - diff_3_x_n_1*dx*dx;
					}
					else
					{
						rho_p_x = diff_1_x_0 - diff_2_x_0*dx - diff_3_x_0*dx*dx;
					}
				}
				else
				{
					if (abs(diff_3_x_0) <= abs(diff_3_x_p_1))
					{
						rho_p_x = diff_1_x_0 - diff_2_x_p_1*dx + 2*diff_3_x_0*dx*dx;
					}
					else
					{
						rho_p_x = diff_1_x_0 - diff_2_x_p_1*dx + 2*diff_3_x_p_1*dx*dx;
					}
				}

				// y-component
				TT diff_1_y_n_3 = one_over_dy*(rho_ghost(i, j - 2) - rho_ghost(i, j - 3));
				TT diff_1_y_n_2 = one_over_dy*(rho_ghost(i, j - 1) - rho_ghost(i, j - 2));
				TT diff_1_y_n_1 = one_over_dy*(rho_ghost(i,     j) - rho_ghost(i, j - 1));
				TT diff_1_y_0   = one_over_dy*(rho_ghost(i, j + 1) - rho_ghost(i    , j));
				TT diff_1_y_p_1 = one_over_dy*(rho_ghost(i, j + 2) - rho_ghost(i, j + 1));
				TT diff_1_y_p_2 = one_over_dy*(rho_ghost(i, j + 3) - rho_ghost(i, j + 2));

				TT diff_2_y_n_2 = one_over_2dy*(diff_1_y_n_2 - diff_1_y_n_3);
				TT diff_2_y_n_1 = one_over_2dy*(diff_1_y_n_1 - diff_1_y_n_2);
				TT diff_2_y_0   = one_over_2dy*(diff_1_y_0   - diff_1_y_n_1);
				TT diff_2_y_p_1 = one_over_2dy*(diff_1_y_p_1 - diff_1_y_0  );
				TT diff_2_y_p_2 = one_over_2dy*(diff_1_y_p_2 - diff_1_y_p_1);

				TT diff_3_y_n_2 = one_over_3dy*(diff_2_y_n_1 - diff_2_y_n_2);
				TT diff_3_y_n_1 = one_over_3dy*(diff_2_y_0   - diff_2_y_n_1);
				TT diff_3_y_0   = one_over_3dy*(diff_2_y_p_1 - diff_2_y_0  );
				TT diff_3_y_p_1 = one_over_3dy*(diff_2_y_p_2 - diff_2_y_p_1);

				TT rho_m_y;
				if (abs(diff_2_y_n_1) <= abs(diff_2_y_0))
				{
					if (abs(diff_3_y_n_2) <= abs(diff_3_y_n_1))
					{
						rho_m_y = diff_1_y_n_1 + diff_2_y_n_1*dy + 2*diff_3_y_n_2*dy*dy;
					}
					else
					{
						rho_m_y = diff_1_y_n_1 + diff_2_y_n_1*dy + 2*diff_3_y_n_1*dy*dy;
					}
				}
				else
				{
					if (abs(diff_3_y_n_1) <= abs(diff_3_y_0))
					{
						rho_m_y = diff_1_y_n_1 + diff_2_y_0*dy - diff_3_y_n_1*dy*dy;
					}
					else
					{
						rho_m_y = diff_1_y_n_1 + diff_2_y_0*dy - diff_3_y_0*dy*dy;
					}
				}

				TT rho_p_y;
				if (abs(diff_2_y_0) <= abs(diff_2_y_p_1))
				{
					if (abs(diff_3_y_n_1) <= abs(diff_3_y_0))
					{
						rho_p_y = diff_1_y_0 - diff_2_y_0*dy - diff_3_y_n_1*dy*dy;
					}
					else
					{
						rho_p_y = diff_1_y_0 - diff_2_y_0*dy - diff_3_y_0*dy*dy;
					}
				}
				else
				{
					if (abs(diff_3_y_0) <= abs(diff_3_y_p_1))
					{
						rho_p_y = diff_1_y_0 - diff_2_y_p_1*dy + 2*diff_3_y_0*dy*dy;
					}
					else
					{
						rho_p_y = diff_1_y_0 - diff_2_y_p_1*dy + 2*diff_3_y_p_1*dy*dy;
					}
				}

				// Velocity Interpolation - Using MAC grid
				T u_vel = (velocity_array(i - 1, j).x + velocity_array(i + 1, j).x)*(T)0.5;
				T v_vel = (velocity_array(i, j - 1).y + velocity_array(i, j + 1).y)*(T)0.5;

				if (u_vel > 0)
				{
					if (v_vel > 0)
					{
						rho_array(i, j) = rho_array(i, j) - dt*(u_vel*rho_m_x + v_vel*rho_m_y);
					}
					else
					{
						rho_array(i, j) = rho_array(i, j) - dt*(u_vel*rho_m_x + v_vel*rho_p_y);
					}
				}
				else
				{
					if (v_vel > 0)
					{
						rho_array(i, j) = rho_array(i, j) - dt*(u_vel*rho_p_x + v_vel*rho_m_y);
					}
					else
					{
						rho_array(i, j) = rho_array(i, j) - dt*(u_vel*rho_p_x + v_vel*rho_p_y);
					}
				}
			}
		}
	}

	// For the velocity update
	static void ENO3rd(FIELD_STRUCTURE_2D<TT>& rho, FIELD_STRUCTURE_2D<TT>& rho_ghost, const FIELD_STRUCTURE_2D<TT>& velocity_x, const FIELD_STRUCTURE_2D<TT>& velocity_y, const T& dt)
	{
		const ARRAY_2D<TT>& velocity_array_x(velocity_x.array_for_this);
		const ARRAY_2D<TT>& velocity_array_y(velocity_y.array_for_this);

		ARRAY_2D<TT>& rho_array(rho.array_for_this);

		rho_ghost.FillGhostCellsFrom(rho_array, true);

		T dx(rho.dx), dy(rho.dy);
		T one_over_dx(rho.one_over_dx), one_over_dy(rho.one_over_dy), one_over_2dx(rho.one_over_2dx), one_over_2dy(rho.one_over_2dy);
		T one_over_3dx(one_over_dx*(T)1/3), one_over_3dy(one_over_dy*(T)1/3);

		for (int j = rho.grid.j_start; j <= rho.grid.j_end; j++)
		{
			for (int i = rho.grid.i_start; i <= rho.grid.i_end; i++)
			{
				// x-component
				TT diff_1_x_n_3 = one_over_dx*(rho_ghost(i - 2, j) - rho_ghost(i - 3, j));
				TT diff_1_x_n_2 = one_over_dx*(rho_ghost(i - 1, j) - rho_ghost(i - 2, j));
				TT diff_1_x_n_1 = one_over_dx*(rho_ghost(i,     j) - rho_ghost(i - 1, j));
				TT diff_1_x_0   = one_over_dx*(rho_ghost(i + 1, j) - rho_ghost(i    , j));
				TT diff_1_x_p_1 = one_over_dx*(rho_ghost(i + 2, j) - rho_ghost(i + 1, j));
				TT diff_1_x_p_2 = one_over_dx*(rho_ghost(i + 3, j) - rho_ghost(i + 2, j));

				TT diff_2_x_n_2 = one_over_2dx*(diff_1_x_n_2 - diff_1_x_n_3);
				TT diff_2_x_n_1 = one_over_2dx*(diff_1_x_n_1 - diff_1_x_n_2);
				TT diff_2_x_0   = one_over_2dx*(diff_1_x_0   - diff_1_x_n_1);
				TT diff_2_x_p_1 = one_over_2dx*(diff_1_x_p_1 - diff_1_x_0  );
				TT diff_2_x_p_2 = one_over_2dx*(diff_1_x_p_2 - diff_1_x_p_1);

				TT diff_3_x_n_2 = one_over_3dx*(diff_2_x_n_1 - diff_2_x_n_2);
				TT diff_3_x_n_1 = one_over_3dx*(diff_2_x_0   - diff_2_x_n_1);
				TT diff_3_x_0   = one_over_3dx*(diff_2_x_p_1 - diff_2_x_0  );
				TT diff_3_x_p_1 = one_over_3dx*(diff_2_x_p_2 - diff_2_x_p_1);

				TT rho_m_x;
				if (abs(diff_2_x_n_1) <= abs(diff_2_x_0))
				{
					if (abs(diff_3_x_n_2) <= abs(diff_3_x_n_1))
					{
						rho_m_x = diff_1_x_n_1 + diff_2_x_n_1*dx + 2*diff_3_x_n_2*dx*dx;
					}
					else
					{
						rho_m_x = diff_1_x_n_1 + diff_2_x_n_1*dx + 2*diff_3_x_n_1*dx*dx;
					}
				}
				else
				{
					if (abs(diff_3_x_n_1) <= abs(diff_3_x_0))
					{
						rho_m_x = diff_1_x_n_1 + diff_2_x_0*dx - diff_3_x_n_1*dx*dx;
					}
					else
					{
						rho_m_x = diff_1_x_n_1 + diff_2_x_0*dx - diff_3_x_0*dx*dx;
					}
				}

				TT rho_p_x;
				if (abs(diff_2_x_0) <= abs(diff_2_x_p_1))
				{
					if (abs(diff_3_x_n_1) <= abs(diff_3_x_0))
					{
						rho_p_x = diff_1_x_0 - diff_2_x_0*dx - diff_3_x_n_1*dx*dx;
					}
					else
					{
						rho_p_x = diff_1_x_0 - diff_2_x_0*dx - diff_3_x_0*dx*dx;
					}
				}
				else
				{
					if (abs(diff_3_x_0) <= abs(diff_3_x_p_1))
					{
						rho_p_x = diff_1_x_0 - diff_2_x_p_1*dx + 2*diff_3_x_0*dx*dx;
					}
					else
					{
						rho_p_x = diff_1_x_0 - diff_2_x_p_1*dx + 2*diff_3_x_p_1*dx*dx;
					}
				}

				// y-component
				TT diff_1_y_n_3 = one_over_dy*(rho_ghost(i, j - 2) - rho_ghost(i, j - 3));
				TT diff_1_y_n_2 = one_over_dy*(rho_ghost(i, j - 1) - rho_ghost(i, j - 2));
				TT diff_1_y_n_1 = one_over_dy*(rho_ghost(i,     j) - rho_ghost(i, j - 1));
				TT diff_1_y_0   = one_over_dy*(rho_ghost(i, j + 1) - rho_ghost(i    , j));
				TT diff_1_y_p_1 = one_over_dy*(rho_ghost(i, j + 2) - rho_ghost(i, j + 1));
				TT diff_1_y_p_2 = one_over_dy*(rho_ghost(i, j + 3) - rho_ghost(i, j + 2));

				TT diff_2_y_n_2 = one_over_2dy*(diff_1_y_n_2 - diff_1_y_n_3);
				TT diff_2_y_n_1 = one_over_2dy*(diff_1_y_n_1 - diff_1_y_n_2);
				TT diff_2_y_0   = one_over_2dy*(diff_1_y_0   - diff_1_y_n_1);
				TT diff_2_y_p_1 = one_over_2dy*(diff_1_y_p_1 - diff_1_y_0  );
				TT diff_2_y_p_2 = one_over_2dy*(diff_1_y_p_2 - diff_1_y_p_1);

				TT diff_3_y_n_2 = one_over_3dy*(diff_2_y_n_1 - diff_2_y_n_2);
				TT diff_3_y_n_1 = one_over_3dy*(diff_2_y_0   - diff_2_y_n_1);
				TT diff_3_y_0   = one_over_3dy*(diff_2_y_p_1 - diff_2_y_0  );
				TT diff_3_y_p_1 = one_over_3dy*(diff_2_y_p_2 - diff_2_y_p_1);

				TT rho_m_y;
				if (abs(diff_2_y_n_1) <= abs(diff_2_y_0))
				{
					if (abs(diff_3_y_n_2) <= abs(diff_3_y_n_1))
					{
						rho_m_y = diff_1_y_n_1 + diff_2_y_n_1*dy + 2*diff_3_y_n_2*dy*dy;
					}
					else
					{
						rho_m_y = diff_1_y_n_1 + diff_2_y_n_1*dy + 2*diff_3_y_n_1*dy*dy;
					}
				}
				else
				{
					if (abs(diff_3_y_n_1) <= abs(diff_3_y_0))
					{
						rho_m_y = diff_1_y_n_1 + diff_2_y_0*dy - diff_3_y_n_1*dy*dy;
					}
					else
					{
						rho_m_y = diff_1_y_n_1 + diff_2_y_0*dy - diff_3_y_0*dy*dy;
					}
				}

				TT rho_p_y;
				if (abs(diff_2_y_0) <= abs(diff_2_y_p_1))
				{
					if (abs(diff_3_y_n_1) <= abs(diff_3_y_0))
					{
						rho_p_y = diff_1_y_0 - diff_2_y_0*dy - diff_3_y_n_1*dy*dy;
					}
					else
					{
						rho_p_y = diff_1_y_0 - diff_2_y_0*dy - diff_3_y_0*dy*dy;
					}
				}
				else
				{
					if (abs(diff_3_y_0) <= abs(diff_3_y_p_1))
					{
						rho_p_y = diff_1_y_0 - diff_2_y_p_1*dy + 2*diff_3_y_0*dy*dy;
					}
					else
					{
						rho_p_y = diff_1_y_0 - diff_2_y_p_1*dy + 2*diff_3_y_p_1*dy*dy;
					}
				}

				T u_vel, v_vel;

				// Velocity Interpolation - Using MAC grid
				if (rho.is_scalar == true)
				{
					u_vel = (velocity_array_x(i + 1, j) + velocity_array_x(i, j))*(T)0.5;
					v_vel = (velocity_array_y(i, j + 1) + velocity_array_y(i, j))*(T)0.5;
				}

				if (rho.is_x_component == true)
				{
					u_vel = velocity_array_x(i, j);
					v_vel = (velocity_array_y(i, j) + velocity_array_y(i, j + 1) + velocity_array_y(i - 1, j) + velocity_array_y(i - 1, j + 1))*(T)0.25;
				}
				else if (rho.is_y_component == true)
				{
					u_vel = (velocity_array_x(i, j - 1) + velocity_array_x(i + 1, j - 1) + velocity_array_x(i, j) + velocity_array_x(i + 1, j))*(T)0.25;
					v_vel = velocity_array_y(i, j);
				}

				if (u_vel > 0)
				{
					if (v_vel > 0)
					{
						rho_array(i, j) = rho_array(i, j) - dt*(u_vel*rho_m_x + v_vel*rho_m_y);
					}
					else
					{
						rho_array(i, j) = rho_array(i, j) - dt*(u_vel*rho_m_x + v_vel*rho_p_y);
					}
				}
				else
				{
					if (v_vel > 0)
					{
						rho_array(i, j) = rho_array(i, j) - dt*(u_vel*rho_p_x + v_vel*rho_m_y);
					}
					else
					{
						rho_array(i, j) = rho_array(i, j) - dt*(u_vel*rho_p_x + v_vel*rho_p_y);
					}
				}
			}
		}
	}

	static void ENO3rd(FIELD_STRUCTURE_2D<TT>& rho, FIELD_STRUCTURE_2D<TT>& rho_ghost, const FIELD_STRUCTURE_2D<TT>& velocity_x, const FIELD_STRUCTURE_2D<TT>& velocity_y, const T& dt, MULTITHREADING& multithreading, const int& thread_id)
	{
		const ARRAY_2D<TT>& velocity_array_x(velocity_x.array_for_this);
		const ARRAY_2D<TT>& velocity_array_y(velocity_y.array_for_this);

		ARRAY_2D<TT>& rho_array(rho.array_for_this);

		rho_ghost.FillGhostCellsFrom(rho_array, true, thread_id);

		multithreading.Sync(thread_id);

		T dx(rho.dx), dy(rho.dy);
		T one_over_dx(rho.one_over_dx), one_over_dy(rho.one_over_dy), one_over_2dx(rho.one_over_2dx), one_over_2dy(rho.one_over_2dy);
		T one_over_3dx(one_over_dx*(T)1/3), one_over_3dy(one_over_dy*(T)1/3);

		GRID_ITERATION_2D(rho.partial_grids[thread_id])
		{
			// x-component
			TT diff_1_x_n_3 = one_over_dx*(rho_ghost(i - 2, j) - rho_ghost(i - 3, j));
			TT diff_1_x_n_2 = one_over_dx*(rho_ghost(i - 1, j) - rho_ghost(i - 2, j));
			TT diff_1_x_n_1 = one_over_dx*(rho_ghost(i,     j) - rho_ghost(i - 1, j));
			TT diff_1_x_0   = one_over_dx*(rho_ghost(i + 1, j) - rho_ghost(i    , j));
			TT diff_1_x_p_1 = one_over_dx*(rho_ghost(i + 2, j) - rho_ghost(i + 1, j));
			TT diff_1_x_p_2 = one_over_dx*(rho_ghost(i + 3, j) - rho_ghost(i + 2, j));

			TT diff_2_x_n_2 = one_over_2dx*(diff_1_x_n_2 - diff_1_x_n_3);
			TT diff_2_x_n_1 = one_over_2dx*(diff_1_x_n_1 - diff_1_x_n_2);
			TT diff_2_x_0   = one_over_2dx*(diff_1_x_0   - diff_1_x_n_1);
			TT diff_2_x_p_1 = one_over_2dx*(diff_1_x_p_1 - diff_1_x_0  );
			TT diff_2_x_p_2 = one_over_2dx*(diff_1_x_p_2 - diff_1_x_p_1);

			TT diff_3_x_n_2 = one_over_3dx*(diff_2_x_n_1 - diff_2_x_n_2);
			TT diff_3_x_n_1 = one_over_3dx*(diff_2_x_0   - diff_2_x_n_1);
			TT diff_3_x_0   = one_over_3dx*(diff_2_x_p_1 - diff_2_x_0  );
			TT diff_3_x_p_1 = one_over_3dx*(diff_2_x_p_2 - diff_2_x_p_1);

			TT rho_m_x;
			if (abs(diff_2_x_n_1) <= abs(diff_2_x_0))
			{
				if (abs(diff_3_x_n_2) <= abs(diff_3_x_n_1))
				{
					rho_m_x = diff_1_x_n_1 + diff_2_x_n_1*dx + 2*diff_3_x_n_2*dx*dx;
				}
				else
				{
					rho_m_x = diff_1_x_n_1 + diff_2_x_n_1*dx + 2*diff_3_x_n_1*dx*dx;
				}
			}
			else
			{
				if (abs(diff_3_x_n_1) <= abs(diff_3_x_0))
				{
					rho_m_x = diff_1_x_n_1 + diff_2_x_0*dx - diff_3_x_n_1*dx*dx;
				}
				else
				{
					rho_m_x = diff_1_x_n_1 + diff_2_x_0*dx - diff_3_x_0*dx*dx;
				}
			}

			TT rho_p_x;
			if (abs(diff_2_x_0) <= abs(diff_2_x_p_1))
			{
				if (abs(diff_3_x_n_1) <= abs(diff_3_x_0))
				{
					rho_p_x = diff_1_x_0 - diff_2_x_0*dx - diff_3_x_n_1*dx*dx;
				}
				else
				{
					rho_p_x = diff_1_x_0 - diff_2_x_0*dx - diff_3_x_0*dx*dx;
				}
			}
			else
			{
				if (abs(diff_3_x_0) <= abs(diff_3_x_p_1))
				{
					rho_p_x = diff_1_x_0 - diff_2_x_p_1*dx + 2*diff_3_x_0*dx*dx;
				}
				else
				{
					rho_p_x = diff_1_x_0 - diff_2_x_p_1*dx + 2*diff_3_x_p_1*dx*dx;
				}
			}
			
			// y-component
			TT diff_1_y_n_3 = one_over_dy*(rho_ghost(i, j - 2) - rho_ghost(i, j - 3));
			TT diff_1_y_n_2 = one_over_dy*(rho_ghost(i, j - 1) - rho_ghost(i, j - 2));
			TT diff_1_y_n_1 = one_over_dy*(rho_ghost(i,     j) - rho_ghost(i, j - 1));
			TT diff_1_y_0   = one_over_dy*(rho_ghost(i, j + 1) - rho_ghost(i    , j));
			TT diff_1_y_p_1 = one_over_dy*(rho_ghost(i, j + 2) - rho_ghost(i, j + 1));
			TT diff_1_y_p_2 = one_over_dy*(rho_ghost(i, j + 3) - rho_ghost(i, j + 2));
			
			TT diff_2_y_n_2 = one_over_2dy*(diff_1_y_n_2 - diff_1_y_n_3);
			TT diff_2_y_n_1 = one_over_2dy*(diff_1_y_n_1 - diff_1_y_n_2);
			TT diff_2_y_0   = one_over_2dy*(diff_1_y_0   - diff_1_y_n_1);
			TT diff_2_y_p_1 = one_over_2dy*(diff_1_y_p_1 - diff_1_y_0  );
			TT diff_2_y_p_2 = one_over_2dy*(diff_1_y_p_2 - diff_1_y_p_1);

			TT diff_3_y_n_2 = one_over_3dy*(diff_2_y_n_1 - diff_2_y_n_2);
			TT diff_3_y_n_1 = one_over_3dy*(diff_2_y_0   - diff_2_y_n_1);
			TT diff_3_y_0   = one_over_3dy*(diff_2_y_p_1 - diff_2_y_0  );
			TT diff_3_y_p_1 = one_over_3dy*(diff_2_y_p_2 - diff_2_y_p_1);

			TT rho_m_y;
			if (abs(diff_2_y_n_1) <= abs(diff_2_y_0))
			{
				if (abs(diff_3_y_n_2) <= abs(diff_3_y_n_1))
				{
					rho_m_y = diff_1_y_n_1 + diff_2_y_n_1*dy + 2*diff_3_y_n_2*dy*dy;
				}
				else
				{
					rho_m_y = diff_1_y_n_1 + diff_2_y_n_1*dy + 2*diff_3_y_n_1*dy*dy;
				}
			}
			else
			{
				if (abs(diff_3_y_n_1) <= abs(diff_3_y_0))
				{
					rho_m_y = diff_1_y_n_1 + diff_2_y_0*dy - diff_3_y_n_1*dy*dy;
				}
				else
				{
					rho_m_y = diff_1_y_n_1 + diff_2_y_0*dy - diff_3_y_0*dy*dy;
				}
			}

			TT rho_p_y;
			if (abs(diff_2_y_0) <= abs(diff_2_y_p_1))
			{
				if (abs(diff_3_y_n_1) <= abs(diff_3_y_0))
				{
					rho_p_y = diff_1_y_0 - diff_2_y_0*dy - diff_3_y_n_1*dy*dy;
				}
				else
				{
					rho_p_y = diff_1_y_0 - diff_2_y_0*dy - diff_3_y_0*dy*dy;
				}
			}
			else
			{
				if (abs(diff_3_y_0) <= abs(diff_3_y_p_1))
				{
					rho_p_y = diff_1_y_0 - diff_2_y_p_1*dy + 2*diff_3_y_0*dy*dy;
				}
				else
				{
					rho_p_y = diff_1_y_0 - diff_2_y_p_1*dy + 2*diff_3_y_p_1*dy*dy;
				}
			}

			T u_vel, v_vel;

			// Velocity Interpolation - Using MAC grid
			if (rho.is_scalar == true)
			{
				u_vel = (velocity_array_x(i + 1, j) + velocity_array_x(i, j))*(T)0.5;
				v_vel = (velocity_array_y(i, j + 1) + velocity_array_y(i, j))*(T)0.5;
			}
			if (rho.is_x_component == true)
			{
				u_vel = velocity_array_x(i, j);
				v_vel = (velocity_array_y(i, j) + velocity_array_y(i, j + 1) + velocity_array_y(i - 1, j) + velocity_array_y(i - 1, j + 1))*(T)0.25;
			}
			else if (rho.is_y_component == true)
			{
				u_vel = (velocity_array_x(i, j - 1) + velocity_array_x(i + 1, j - 1) + velocity_array_x(i, j) + velocity_array_x(i + 1, j))*(T)0.25;
				v_vel = velocity_array_y(i, j);
			}
			if (u_vel > 0)
			{
				if (v_vel > 0)
				{
					rho_array(i, j) = rho_array(i, j) - dt*(u_vel*rho_m_x + v_vel*rho_m_y);
				}
				else
				{
					rho_array(i, j) = rho_array(i, j) - dt*(u_vel*rho_m_x + v_vel*rho_p_y);
				}
			}
			else
			{
				if (v_vel > 0)
				{
					rho_array(i, j) = rho_array(i, j) - dt*(u_vel*rho_p_x + v_vel*rho_m_y);
				}
				else
				{
					rho_array(i, j) = rho_array(i, j) - dt*(u_vel*rho_p_x + v_vel*rho_p_y);
				}
			}
		}
		multithreading.Sync(thread_id);
	}
};
