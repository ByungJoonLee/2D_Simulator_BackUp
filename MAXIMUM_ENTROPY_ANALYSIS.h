#pragma once

#include "ARRAY.h"
#include "WORLD_DISCRETIZATION_2D.h"
#include "FIELD_STRUCTURE_1D.h"
#include "SCRIPT_READER.h"

class MAXIMUM_ENTROPY_ANALYSIS
{
public: // Essential Data
    WORLD_DISCRETIZATION_2D*	world_discretization;
	
	ARRAY<T>					x;
    ARRAY<T>					A;
	ARRAY<T>					tau, best;

	int							Nsamps;
	int							Ncoef;
    
    int							mfl;
	T							bptu;
	T							Power;
    int							ihr, nihr;

	T							min_f, max_f;

	// For drawing
	GRID_STRUCTURE_1D			base_grid;
	FIELD_STRUCTURE_1D<T>		power_spectral_density; 


public: // Constructor and Destructor
    MAXIMUM_ENTROPY_ANALYSIS(void)
	{}

    ~MAXIMUM_ENTROPY_ANALYSIS(void)
	{}

public: // Initialization Function
    void InitializeFromScriptBlock(const SCRIPT_BLOCK& signal_processing_block)
	{
	    bptu = signal_processing_block.GetFloat("bins_per_time_unit", (T)10);
        mfl = signal_processing_block.GetInteger("minimum_filter_length", (int)30);
        Ncoef = signal_processing_block.GetInteger("number_of_coefs", (int)10);
		
		x.Initialize(Nsamps);

		ifstream infile("MESA_Data.txt");

		string line;
		while (getline(infile, line))
		{
			istringstream iss(line);

			int i;
			T a;
			
			x[i] = a;
			
			if (!(iss >> i >> a))
			{
				break;
			}
		}
		
		A.Initialize(Ncoef);
		best.Initialize(Ncoef);
        tau.Initialize(nihr);

		//base_grid.Initialize(nihr, 0, min_f, max_f);
		power_spectral_density.Initialize(world_discretization->world_grid_1d);

        //Detren(Nsamps);
        //Lopass(Nsamps);
	}

public: // Member Function
    void Estimate()
	{
	    Power = 0.0;
			
		Coef(x, Power, A);

		for (int i = 0; i < best.length; i++)
		{
			best[i] = -A[i];
		}

        MESA(min_f, max_f, nihr, Ncoef, Power, bptu);
	}

	void Coef(const ARRAY<T>& data, T& xms, ARRAY<T>& coef)
	{
		int n = data.length, m = coef.length;

		// Initialize Power
		T p = 0.0;

		ARRAY<T> wk1(n), wk2(n), wkm(m);
		
		for (int i = 0; i < n; i++)
		{
			p += POW2(data[i]);
		}

		xms = p/n;

		wk1[0] = data[0];
		wk2[n -2] = data[n - 1];

		for (int i = 1; i < n - 1; i++)
		{
			wk1[i] = data[i];
			wk2[i - 1] = data[i];
		}
		
		for (int k = 0; k < m; k++)
		{
			T num = 0.0, denom = 0.0;
			
			for (int j = 0; j < (n - k - 1); j++)
			{
				num += (wk1[j]*wk2[j]);
				denom += (POW2(wk1[j]) + POW2(wk2[j]));
			}
			
			coef[k] = 2.0*num/denom;
			
			xms *= (1.0 - POW2(coef[k]));
			
			for (int i = 0; i < k; i++)
			{
				coef[i] = wkm[i] - coef[k]*wkm[k-1-i];
			}
			
			if (k == m-1)
			{
				return;
			}
			
			for (int i = 0; i <= k; i++)
			{
				wkm[i] = coef[i];
			}

			for (int i = 0; i < (n - k - 2); i++)
			{
				wk1[i] -= (wkm[k]*wk2[i]);
				wk2[i] = wk2[i + 1] - wkm[k]*wk1[i + 1];
			}
		}
		throw("never get here in coef");
	}


    void MESA(const T& min_f, const T & max_f, const int& nihr, const int& num_of_coefs, const T& p0, const T& bptu)
	{
	    const T pi2 = (T)2*PI;
        
		for (int i = 0; i < nihr; i++)
		{
            T freq = min_f + (i - 1)*(max_f - min_f)/(T)(nihr - 1);
            tau[i] = (T)1/freq;
            
			T ream = 1;
            T imag = 0;

			for (int j = 0; j < num_of_coefs; j++)
			{
                T omega = -pi2*freq*j;
                ream = ream + best[j]*cos(omega);
                imag = imag + best[j]*sin(omega);
			}

			power_spectral_density[i] = p0/(ream*ream + imag*imag);
			            
			tau[i] = tau[i]/bptu;
			
			if (tau[i] < (T)0.1)
			{
                break;
			}
		}
	}

    void Detren(const int& num_of_samps)
	{
	    ARRAY<T> ti;
        ti.Initialize(num_of_samps);

        ARRAY<T> res;
        res.Initialize(num_of_samps);

		T sum_x = 0;
        T sum_ti = 0;
        T sum_xt = 0;
        T sum_tsq = 0;
        T en = (T)num_of_samps;
		
		for (int i = 0; i < num_of_samps; i++)
		{
            ti[i] = float(i);
            sum_x = sum_x + x[i];
            sum_ti = sum_ti + ti[i];
            sum_tsq = sum_tsq + (ti[i]*ti[i]);
            sum_xt = sum_xt + (x[i]*x[i]);
		}

        T x_bar = sum_x/en;
        T t_bar = sum_ti/en;
        T q = (sum_xt - (sum_x*sum_ti)/en)/(sum_tsq - (sum_ti*sum_ti)/en);
        // x_int = x_bar + (q*(0-t_bar));
		for (int i = 0; i < num_of_samps; i++)
		{
            T x_hat = x_bar + q*(ti[i] - t_bar);
            res[i] = x[i] - x_hat;
            x[i] = res[i];
		}
	}

    void Lopass(const int& num_of_samps)
	{
	    ARRAY<T> y;
        y.Initialize(num_of_samps);

        // Two pole low pass Betterworth Filter
        y[0] = x[0];
        y[1] = x[1];

        T P = 9.656851;
        T Q = -3.414213;
        T C = 10.24264;

        // End Filter bank
		for (int i = 2; i < num_of_samps; i++)
		{
            y[i] = (x[i] + 2*x[i-1] + x[i-2] + P*y[i-1] + Q*y[i-2])/C;
		}

		for (int i = 0; i < num_of_samps; i++)
		{
            x[i] = y[i];
		}
	}
};