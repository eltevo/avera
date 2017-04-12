#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>
#include <omp.h>
#include <time.h>
#include<DTFE.h>
#include "global_variables.h"

#ifdef USE_SINGLE_PRECISION
typedef float REAL;
#else
typedef double REAL;
#endif

//density field estimation with DTFE method

void DTFE_density(REAL** x)
{
printf("DTFE starting...\n");
//időmérés
REAL start_time = (REAL) clock () / (REAL) CLOCKS_PER_SEC;
REAL omp_start_time = omp_get_wtime();
//időmérés

char NO_grid[14];
sprintf(NO_grid, "%i", DENSITY_CELLS);
char L_str[19];
sprintf(L_str, "%.13f", L);
char box_min_str[4];
sprintf(box_min_str, "0.0");
char grid_str[7];
sprintf(grid_str, "--grid");
char periodic_str[11];
sprintf(periodic_str, "--periodic");
char box_str[6];
sprintf(box_str, "--box");
char verbose_str[7];
sprintf(verbose_str, "--verbose");
char DTFE_str[5];
sprintf(DTFE_str, "DTFE");
char one_str[2];
sprintf(one_str, "1");
//read the program options
int argc;
char *argv_DTFE[] = {DTFE_str, grid_str, NO_grid, periodic_str, box_str, box_min_str, L_str, box_min_str, L_str, box_min_str, L_str, verbose_str, one_str};
char sph_str[6];
sprintf(sph_str, "--SPH");
char neighbors_str[3];
sprintf(neighbors_str, "40");
char *argv_SPH[] = {DTFE_str, grid_str, NO_grid, periodic_str, box_str, box_min_str, L_str, box_min_str, L_str, box_min_str, L_str, verbose_str, one_str, sph_str, neighbors_str};

User_options userOptions;
if(NONISOTROPIC_EXPANSION == 3)
{
argc = 13;
userOptions.readOptions(argc, argv_DTFE, false, true);
}
else
{
//if NONISOTROPIC_EXPANSION == 4, the code will use SPH field estimator
printf("SPH density estimation...\n");
sprintf(one_str, "1");
argc = 15;
userOptions.readOptions(argc, argv_SPH, false, true);
}

// define variables to store the particle data and the output results
std::vector<Particle_data> particles; // vector for particle data
std::vector<Sample_point> samplingCoordinates;// sampling coordinates for non-regular grid
Quantities uquantities;           // class to keep the results of the interpolation
Quantities aquantities;		// structure that will store the volume averaged output quantities - i.e. the fields volume averaged over the sampling cells (NOTE: "aQuantities" stands for "averaged quantities")

/* read the input particle data into the DTFE data structure */
/* here: pos - vector to particle positions, vel-particle velocities, etc ... */
unsigned N_un = (unsigned) N;
particles.reserve(N); // reserve memory for ’noParticles’ particles
double M_new = (double) M/rho_crit*pow(a_max/a_start, 3.0);
for (size_t i=0; i<N_un; ++i)
{
	Particle_data temp;
	for (int j=0; j<3; ++j)		// read particle i-th position
		temp.position(j) = (double) x[i][j];
	temp.weight() = M_new;			// read particle i-th weight(mass)
	particles.push_back(temp);
}
// compute the DTFE interpolation
DTFE(&particles, samplingCoordinates, userOptions, &uquantities, &aquantities);
// the above function deletes the data contained in the ’particles’ vector
int DENSITY_CELLS3 = pow(DENSITY_CELLS,3);
//printf("DENSITY_CELLS3 = %i\n", DENSITY_CELLS3);
int rho_min, rho_max;
rho_min = 0;
rho_max = 0;
for(int i=0; i<DENSITY_CELLS3; i++)
{
	RHO[i] = ((double) aquantities.density[i]);//Cell density in Omega_m
	if(RHO[rho_min] > RHO[i])
		rho_min = i;

	if(RHO[rho_max] < RHO[i])
		rho_max = i;

}
//timing
REAL end_time = (REAL) clock () / (REAL) CLOCKS_PER_SEC;
REAL omp_end_time = omp_get_wtime();
//timing

printf("The minimal and the maximal density (in Omega_m):\n Rho_min = %g\t Rho_max = %g\t\n", RHO[rho_min],RHO[rho_max]);
printf("...DTFE done.\n");
if(NONISOTROPIC_EXPANSION == 3)
{
	printf("DTFE CPU time = %lfs\n", end_time-start_time);
	printf("DTFE RUN time = %lfs\n", omp_end_time-omp_start_time);
}
else
{
	printf("SPH CPU time = %lfs\n", end_time-start_time);
	printf("SPH RUN time = %lfs\n", omp_end_time-omp_start_time);
}
return;
}
