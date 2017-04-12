#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include "global_variables.h"

#include "voro++.hh"
using namespace voro;

#ifdef USE_SINGLE_PRECISION
typedef float REAL;
#else
typedef double REAL;
#endif

extern int N;
extern REAL L;
extern REAL M;
extern double *RHO;
extern REAL **x;

double friedman_solver_step(double a0, double h, double Omega_lambda, double Omega_r, double Omega_m, double Omega_k, double H0);

void get_voronoi()
{
	int i;
	// Set up the number of blocks that the container is divided into
	printf("Calculating Voronoi cells...\n");
	//Timing
	REAL voro_start_time = (REAL) clock () / (REAL) CLOCKS_PER_SEC;
	int n_x=6,n_y=6,n_z=6;
	// Create a container with the geometry given above, and make it
        // periodic in each of the three coordinates. Allocate space for
        // eight particles within each computational block
	container con(0,(double)L,0,(double)L,0,(double)L,n_x,n_y,n_z, true,true,true,64);
	for(i = 0; i < N; i++)
	{
		con.put(i,(double)x[i][0],(double)x[i][1],(double)x[i][2]);
	}
	//Calculating volumes of all cells
	i = 0;
	c_loop_all cla(con);
	voronoicell c;
	printf("...done.\n Calculating cell volumes...\n");
	double vvol = 0;
	int vvol_min, vvol_max;
	vvol_min = 0;
	vvol_max = 0;
	if(cla.start()) do if (con.compute_cell(c,cla))
	{
			vol_cell[i] = c.volume();
			vvol += vol_cell[i];
			RHO[i] = (double)M/(vol_cell[i]*pow(a_start/a_max, 3.0)*rho_crit);
			if(vol_cell[vvol_min] > vol_cell[i])
				vvol_min = i;

			if(vol_cell[vvol_max] < vol_cell[i])
                        	vvol_max = i;
			i++;
	} while (cla.inc());
	printf("...done.\n");
	printf("Full box volume:\t\t %g\n"
               "Sum volume of the Voronoi cells:\t %g\n"
               "Difference:\t\t\t %g\n\n",(double)(L*L*L),vvol,vvol-double(L*L*L));
	printf(" Volume of the smallest cell:\t %f\n Volume of the largest cell:\t %f\n",vol_cell[vvol_min],vol_cell[vvol_max]);
	printf("The minimal and the maximal density(Omega_m):\n Rho_max = %g\t Rho_min = %g\t\n", RHO[vvol_min],RHO[vvol_max]);
	//Timing
	REAL voro_end_time = (REAL) clock () / (REAL) CLOCKS_PER_SEC;
	printf("Voro++ CPU time = %lfs\n", voro_end_time-voro_start_time);
	return;
}

double nonis_friedmann_voronoi(double* RHO, double a_prev)
{
        int i;
        double new_a_local, delta_V_box, a_effective, delta_a_effective, delta_a_homogeneous;
        delta_V_box = 0;
        for(i=0;i<N;i++)
        {
                if(LOC_CURV == 1)
                {
                        	new_a_local = friedman_solver_step(a, h, Omega_lambda, Omega_r, (Omega_b+RHO[i]), 1.0-Omega_lambda-Omega_r-(Omega_b+RHO[i]), H0);
                }
                else
                {
                        new_a_local = friedman_solver_step(a, h, Omega_lambda, Omega_r, (Omega_b+RHO[i]), Omega_k, H0);
                }

                delta_V_box += vol_cell[i]*(pow(new_a_local,3.0) - pow(a_prev, 3.0));

        }
        a_effective = pow((delta_V_box/pow((double)L, 3.0) + pow(a_prev,3.0)), 1.0/3.0); //Calculating the effecive scale factor
	delta_a_effective = a_effective - a_prev;
        delta_a_homogeneous = friedman_solver_step(a, h, Omega_lambda, Omega_r, Omega_m, Omega_k, H0)-a_prev;
        printf("delta_a_effective = %g\t delta_a_homogeneous = %g\n", delta_a_effective, delta_a_homogeneous);
        return a_effective;
}
