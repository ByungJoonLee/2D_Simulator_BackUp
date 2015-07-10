#pragma once

#include "ARRAY.h"
#include "SCRIPT_READER.h"

class MAXIMUM_ENTROPY_ANALYSIS
{
public: // Essential Data
    ARRAY<T>            x;
    ARRAY<T>            B;
	ARRAY<T>            R;
    ARRAY<T>            A;
    ARRAY<T>            A0;
    ARRAY<T>            tau, best, F;

	int                 Nsamps;
	int                 Ncoef;
    
    T                   mfl;
	T                   p0;
    T                   bptu;
	T                   Power;
    T                   Mfpe;
    T                   Fpe;
    T                   Nom, Denom;
    T                   Lfilt;

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
        Nsamps = signal_processing_block.GetInteger("number_of_samples", (int)100);
        Ncoef = Nsamps/2;

		B.Initialize(Nsamps);
        x.Initialize(Nsamps);
        tau.Initialize(Nsamps);
        best.Initialize(Nsamps);
        F.Initialize(Nsamps);

        Detren(Nsamps);
        Lopass(Nsamps);
	}

public: // Member Function
    void Estimate()
	{
	    Power = (T)0;
		Mfpe = (T)0;

		for (int i = 0; i < Nsamps; i++)
		{
            B[i] = x[i];
            Power = Power + x[i]*x[i];
		}

        Power = Power/(T)Nsamps;
        Fpe = (T)(Nsamps + 1)/(T)(Nsamps - 1)*Power;
        T Ftemp = Fpe;
        Fpe = 0.0;

        A[0] = -1;
        A0[0] = -1;

        int NN = Ncoef - 1;

		for (int k = 0; k < NN; k++)
		{
            int NMM = Nsamps - k;
            
			// Update B, B'
			for (int i = 0; i < NMM; i++)
			{
                x[i] = x[i] - A[k]*B[i];
				B[i] = B[i+1] - A[k]*x[i+1];
			}

            // Compute A(M,M)
            Nom = (T)0;
            Denom = (T)0;

			for (int i = 0; i < NMM; i++)
			{
                Nom = Nom + x[i]*B[i];
                Denom = Denom + x[i]*x[i] + B[i]*B[i];
			}

            A[k+1] = 2*Nom/Denom;

			Power = Power*(1-A[k+1]*A[k+1]);

            // Compute Akaike FPE
            int N0 = Nsamps/2;
			
			if (k == 0 ||  k >= N0)
			{
                // Doing nothing!
			}
			else
			{
                Fpe = (T)(Nsamps + k)/(T)(Nsamps - k)*Power;
                Fpe = Fpe/Ftemp;
                Fpe = log10(Fpe);
			}
			
			if (k == 0)
			{
                // Doint nothing!
			}
			else
			{
				for (int i = 1; i < k; i++)
				{
					A0[i] = A[i] - A[k+1]*A[k+2-i];
				}

				for (int i = 1; i < k; i++)
				{
					A[i] = A0[i];
				}

				A0[k+1] = A[k+1]; 
			}
            
			if (Fpe >= Mfpe)
			{
                continue;
			}

			if (k < mfl)
			{
                continue;
			}

			if (k == mfl)
			{
				Lfilt = k + 1;
				Mfpe = Fpe;

				p0 = Power;

				for (int i = 0; i < Lfilt; i++)
				{
					best[i] = A0[i];
				} 
			}
		}

        // Fill Array with Coeff. For Proper Length Filter
        
        Ncoef = Lfilt;
        A[0] = -1;
        
        // Fix the signs
		for (int i = 0; i < Ncoef; i++)
		{
            best[i] = -best[i];
		}

        T min_f = 2.0/(T)Nsamps;
        T max_f = 0.5;

        MESA(min_f, max_f, Nsamps, Ncoef, p0, bptu);
	}

    void MESA(const T& min_f, const T & max_f, const int& num_of_samps, const int& num_of_coefs, const T& p0, const T& bptu)
	{
	    const T pi2 = (T)2*PI;
        const int ihr = 32;
        const int nihr = num_of_samps*ihr;

		for (int i = 0; i < nihr; i++)
		{
            T freq = min_f + (i - 1)*(max_f - min_f)/(T)(nihr - 1);
            tau[i] = (T)1/freq;
            
			T ream = 0;
            T imag = 0;

			for (int j = 0; j < num_of_coefs; j++)
			{
                T omega = -pi2*freq*j;
                ream = ream + best[j]*cos(omega);
                imag = imag + best[j]*sin(omega);
			}

            F[i] = p0/(ream*ream + imag*imag);
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