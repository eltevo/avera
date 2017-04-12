#include <stdio.h>
#include <math.h>
#include "global_variables.h"

#ifdef USE_SINGLE_PRECISION
typedef float REAL;
#else
typedef double REAL;
#endif

//This file contains the functions that required for the simple boxing backreaction, and the local friedmann equation solver.

double friedman_solver_step(double a0, double h, double Omega_lambda, double Omega_r, double Omega_m, double Omega_k, double H0);

void density_field(REAL **x, double* RHO, int DENSITY_CELLS)
{
	int i, DENSITY_CELLS3;
	int h,k,l;
	int index;
	DENSITY_CELLS3 = pow(DENSITY_CELLS, 3);
	for(i=0;i<DENSITY_CELLS3;i++)
	{
		RHO[i] = 0;
	}
	for(i=0;i<N;i++)
	{
		h=(int)(floor(((double)x[i][0]/(double)L)*(double)DENSITY_CELLS));
		k=(int)(floor(((double)x[i][1]/(double)L)*(double)DENSITY_CELLS));
		l=(int)(floor(((double)x[i][2]/(double)L)*(double)DENSITY_CELLS));
		index = h+(DENSITY_CELLS*k)+(DENSITY_CELLS*DENSITY_CELLS*l);
		RHO[index]=RHO[index]+1.0;
	}
	int rho_min, rho_max;
	rho_min = 0;
	rho_max = 0;
	for(i=0;i<(DENSITY_CELLS3);i++)
	{
		RHO[i] = (double)M*RHO[i]/(pow(((double)L/((double)DENSITY_CELLS)), 3.0));
		RHO[i] = RHO[i]/rho_crit*pow(a_start/a_max,-3.0);
		if(RHO[rho_min] > RHO[i])
			rho_min = i;

		if(RHO[rho_max] < RHO[i])
			rho_max = i;
	}
	printf("The minimal and the maximal density (in Omega_m):\n Rho_min = %g\t Rho_max = %g\t\n", RHO[rho_min],RHO[rho_max]);
	return;
}

double nonis_friedmann(double* RHO, int DENSITY_CELLS)
{
	int i, DENSITY_CELLS3;
	double delta_a_local, delta_V_box, delta_a_effective, delta_a_homogeneous;
	DENSITY_CELLS3 = pow(DENSITY_CELLS, 3);
	delta_V_box = 0;
	for(i=0;i<(DENSITY_CELLS3);i++)
	{
		if(LOC_CURV == 1)
		{
			delta_a_local = friedman_solver_step(a, h, Omega_lambda, Omega_r, (Omega_b+RHO[i]), 1.0-Omega_lambda-Omega_r-(Omega_b+RHO[i]), H0)-a_prev;
		}
		else
		{
			delta_a_local = friedman_solver_step(a, h, Omega_lambda, Omega_r, (Omega_b+RHO[i]), Omega_k, H0)-a_prev;
		}
		//The full volumetric change of the box:
		delta_V_box += pow(L/((double) DENSITY_CELLS), 3.0)*(pow((delta_a_local+a_prev),3.0)-pow(a_prev,3.0));

	}
	delta_a_effective = pow((delta_V_box/pow(L, 3.0) + pow(a_prev,3.0)), 1.0/3.0)-a_prev; //The change of the effective scale factor:
	delta_a_homogeneous = friedman_solver_step(a, h, Omega_lambda, Omega_r, Omega_m, Omega_k, H0)-a_prev;
	printf("delta_a_effective = %lg\t delta_a_homogeneous = %lg\n", delta_a_effective, delta_a_homogeneous);
	return delta_a_effective;
}
