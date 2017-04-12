#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <omp.h>
#include <time.h>
#include "global_variables.h"

#ifdef USE_SINGLE_PRECISION
typedef float REAL;
#else
typedef double REAL;
#endif

int t,N,el,hl;
int e[2202][4];
int H[2202][4];
REAL SOFT_CONST[8];
REAL w[3];
double a_max;
REAL** x;
REAL** F;
double h, h_min, h_max, rcut, Ekin, Epot, T, rho, t_next, t_bigbang;
REAL mean_err;
double FIRST_T_OUT, H_OUT; //First output time, output frequency in Gy
double rho_crit; //Critical density

REAL c21,c31,c32,c41,c42,c43,c51,c52,c53,c54,c61,c62,c63,c64,c65,a1,a3,a4,a5,b1,b3,b4,b5,b6;
REAL x4, err, errmax;
REAL beta, ParticleRadi;

int IS_PERIODIC, COSMOLOGY;
REAL L;
char IC_FILE[1024];
char OUT_DIR[1024];
extern char __BUILD_DATE;
int IC_FORMAT; // 0: ascii, 1:GADGET

double Omega_b,Omega_lambda,Omega_dm,Omega_r,Omega_k,Omega_m,H0,H0_start,Hubble_param, Decel_param, delta_Hubble_param, Hubble_tmp; //Cosmologycal parameters

double epsilon=1;
double sigma=1;
REAL G;//Newtonian gravitational constant
REAL M;//Particle mass
double a, a_start,a_prev,a_tmp;//Scalefactor, scalefactor at the starting time, previous scalefactor
double Omega_m_eff; //Effective Omega_m
double delta_a, a_prev1, a_prev2, h_prev;
int NONISOTROPIC_EXPANSION; //inhomogeneous expansion (Backreaction) (0=standard cosmology, 1=naiv backreaction with boxes, 2=backreaction with voronoi-cells, 3=DTFE backreaction, 4=naive SPH backreaction)

int RESTART; //Restarted simulation(0=no, 1=yes)
double T_RESTART; //Time of restart
double A_RESTART; //Scalefactor at the time of restart
double H_RESTART; //Hubble-parameter at the time of restart


//Variables for the inhomogeneous "Friedmann-equation" integrator
int DENSITY_CELLS;
double* RHO;
double* vol_cell;
float DTFE_MEMORY;
int LOC_CURV; // Using local (=1) or global (=0) curvature in backreaction simulations

//Functions for reading GADGET2 format IC
int gadget_format_conversion(void);
int load_snapshot(char *fname, int files);
int allocate_memory(void);
int reordering(void);

//Functions for the naive box backreaction
void density_field(REAL **x, double* RHO, int DENSITY_CELLS);
double nonis_friedmann(double* RHO, int DENSITY_CELLS);
//Function for the Voronoi-cell backreaction
void get_voronoi();
double nonis_friedmann_voronoi(double* RHO, double a_prev);
//Functions for the DTFE backreaction
void DTFE_density(REAL** x);


void read_ic(FILE *ic_file, int N);
void read_param(FILE *param_file);
void step(REAL** x, REAL** F);
void kiiras(REAL** x);
void Log_write(REAL** x);
void forces_old(REAL** x, REAL** F);
void forces_old_periodic(REAL**x, REAL**F);
void forces_EWALD(REAL** x, REAL** F);
double friedmann_solver_start(double a0, double t0, double h, double Omega_lambda, double Omega_r, double Omega_m, double H0, double a_start);
double friedman_solver_step(double a0, double h, double Omega_lambda, double Omega_r, double Omega_m, double Omega_k, double H0);
int ewald_space(REAL R, int ewald_index[2102][4]);
//Functions for rescaling
void rescaling();
double CALCULATE_decel_param(double a, double a_prev1, double a_prev2, double h, double h_prev);


void read_ic(FILE *ic_file, int N)
{
int i,j;

x = (REAL**)malloc(N*sizeof(REAL*)); //Allocating memory
for(i = 0; i < N; i++)
	{
		x[i] = (REAL*)malloc(6*sizeof(REAL));
	}

F = (REAL**)malloc(N*sizeof(REAL*)); 
for(i = 0; i < N; i++)
{
	F[i] = (REAL*)malloc(3*sizeof(REAL));
}



printf("\nReading IC from the %s file...\n", IC_FILE);
for(i=0; i<N; i++) //reading
{
	for(j=0; j<6; j++)
	{
		#ifdef USE_SINGLE_PRECISION
		fscanf(ic_file, "%f", & x[i][j]);
		#else
		fscanf(ic_file, "%lf", & x[i][j]);
		#endif

	}

}
printf("...done.\n\n");
fclose(ic_file);
return;
}


void kiiras(REAL** x)
{
	int i,k;
	char A[20];
	if(COSMOLOGY == 1)
	{
		sprintf(A, "%d", (int)(round(100*t_next*47.1482347621227)));
	}
	else
	{
		sprintf(A, "%d", (int)(round(100*t_next)));
	}
	char filename[0x100];
	snprintf(filename, sizeof(filename), "%st%s.dat", OUT_DIR, A);
	if(COSMOLOGY == 0)
	{
		printf("Saving: t= %f, file: \"%st%s.dat\" \n", t_next, OUT_DIR, A);
	}
	else
	{
		printf("Saving: t= %f, file: \"%st%s.dat\" \n", t_next*47.1482347621227, OUT_DIR, A);
	}
	FILE *coordinate_file;
	if(t < 1)
	{
		coordinate_file = fopen(filename, "w");
	}
	else
	{
		coordinate_file = fopen(filename, "a");
	}

	for(i=0; i<N; i++)
	{
		for(k=0; k<6; k++)
		{
			fprintf(coordinate_file, "%.16f\t",x[i][k]);
		}
		fprintf(coordinate_file, "\n");
	}

	fclose(coordinate_file);
}

void Log_write(REAL** x) //Writing logfile
{
	FILE *LOGFILE;
	char A[] = "Logfile.dat";
	char filename[0x100];
	snprintf(filename, sizeof(filename), "%s%s", OUT_DIR, A);
	LOGFILE = fopen(filename, "a");
	fprintf(LOGFILE, "%.15f\t%e\t%e\t%.15f\t%.15f\t%.15f\t%.15f\t%.10f\n", T*47.1482347621227, errmax, h*47.1482347621227, a, a_max/a-1, Hubble_param*20.7386814448645, Decel_param, Omega_m_eff);
	fclose(LOGFILE);
}


int main(int argc, char *argv[])
{
	printf("-------------------------------------------------------------------\nCCLEA v0.7.0.0\n (Cosmological Code with Local Expansion and Averaging)\n\n Gabor Racz, 2016\n Department of Physics of Complex Systems, Eötvös Loránd University\n\n");
	printf("Build date: %zu\n-------------------------------------------------------------------\n\n", (unsigned long) &__BUILD_DATE);
	int i;
	RESTART = 0;
	T_RESTART = 0;
	if( argc != 2)
	{
		fprintf(stderr, "Missing parameter file!\n");
		fprintf(stderr, "Call with: ./CCLEA  <parameter file>\n");
		return (-1);
	}
	FILE *param_file = fopen(argv[1], "r");
	read_param(param_file);
	if(IS_PERIODIC>1)
	{
		el = ewald_space(3.6,e);
		if(IS_PERIODIC>2)
		{
			hl = ewald_space(8.0,H);
		}
	}
	if(IC_FORMAT != 0 && IC_FORMAT != 1)
        {
                fprintf(stderr, "Error: bad IC format!\nExiting.\n");
                return (-1);
        }
	if(IC_FORMAT == 0)
	{
		FILE *ic_file = fopen(IC_FILE, "r");
		read_ic(ic_file, N);
	}
	if(IC_FORMAT == 1)
	{
		int files;
		printf("The IC file is in Gadget format.\nThe IC determines the box size.\n");
		files = 1;      /* number of files per snapshot */
		x = (REAL**)malloc(N*sizeof(REAL*)); //Allocating memory
		for(i = 0; i < N; i++)
		{
			x[i] = (REAL*)malloc(6*sizeof(REAL));
		}
		F = (REAL**)malloc(N*sizeof(REAL*));
		for(i = 0; i < N; i++)
		{
			F[i] = (REAL*)malloc(3*sizeof(REAL));
		}
		load_snapshot(IC_FILE, files);
		reordering();
		gadget_format_conversion();
	}

	//Allocating memory for the density field
	int DENSITY_CELLS3;
	if(NONISOTROPIC_EXPANSION == 2)
	{
        	DENSITY_CELLS3 = N;
		vol_cell = (double *)malloc(N*sizeof(double));
	}
        else
        {
		DENSITY_CELLS3 = pow(DENSITY_CELLS,3);
	}
	if(NONISOTROPIC_EXPANSION == 3)
	{
		DTFE_MEMORY = 500.0*(((float) N) / 1000000.0) + 4.0*5.0*(((float) N) / 1000000.0) + 4.0*5.0*(((float) DENSITY_CELLS3) / 1000000.0);
		printf("Using DTFE.\nMemory needed for DTFE:%.2fMB\n", DTFE_MEMORY);
	}
	RHO = (double*)malloc(DENSITY_CELLS3*sizeof(double));
	//Critical density, Gravitational constant and particle masses
	if(COSMOLOGY ==1)
	{
		Omega_m = Omega_b+Omega_dm;
                Omega_k = 1.-Omega_m-Omega_lambda-Omega_r;
		G = 1;
		rho_crit = 3*H0*H0/(8*pi*G);
                //rescaling();
		M = Omega_dm*rho_crit*pow(L, 3.0)/((REAL) N);
		rescaling();
	}
	else
	{
		printf("Particle mass:\tM = %.10f\n\n", M);
		G = 1;
	}
	T=0.0;
	beta = ParticleRadi;
	
	recalculate_softening();
	t_next = 0.;
	T = 0;
	REAL Delta_T_out;
	a=1;//scalefactor
	//Calculating initial time
	if(COSMOLOGY == 1)
	{
		a = a_start;
		a_tmp = a;
		printf("a_start/a_ma=%.8f\tz=%.8f\n", a, a_max/a-1);
		
		if(RESTART == 0)
		{
		T = friedmann_solver_start(1,0,h_min*0.00021,Omega_lambda,Omega_r,Omega_m,H0, a_start);
		}
		else
		{
			T = T_RESTART/47.1482347621227; //if the simulation is restarted
			Hubble_param = H_RESTART;
			a = A_RESTART;
			recalculate_softening();
		}
		Delta_T_out = H_OUT/47.1482347621227; //Output frequency
		if(FIRST_T_OUT >= T) //Calculating first output time
		{
			t_next = FIRST_T_OUT/47.1482347621227;
		}
		else
		{
			t_next = T+Delta_T_out;
		}
		printf("Initial time:\tt_start = %.10fGy\nInitial scalefactor:\ta_start = %.8f\nMaximal scalefactor:\ta_max=%.8f\n\n", T*47.1482347621227, a, a_max);
	}
	else
	{
		a = a_max;
		Hubble_param = 0;
		T = 0.0; //If we do not running cosmological simulations, the initial time will be 0.
		printf("t_start = %f\tt_max = %f\n", T, a_max);
		a_tmp = 0;
		Delta_T_out = H_OUT;
		t_next = T+Delta_T_out;
	}
	//Timing
	REAL SIM_start_time = (REAL) clock () / (REAL) CLOCKS_PER_SEC;
	REAL SIM_omp_start_time = omp_get_wtime();
	//Timing

	//Initial force calculation
	if(IS_PERIODIC < 2)
	{
		forces_old(x, F);
	}
	if(IS_PERIODIC == 2)
	{
		forces_old_periodic(x, F);
	}
	
	//The simulation is starting...
	if(COSMOLOGY == 1)
	{
	a_prev1 = friedman_solver_step(a, -1*h, Omega_lambda, Omega_r, Omega_m, Omega_k, H0);
	a_prev2 = friedman_solver_step(a_prev1, -1*h, Omega_lambda, Omega_r, Omega_m, Omega_k, H0);
	}
	else
	{
	a_prev1 = a;
	a_prev2 = a;
	}
	h_prev = h;
	Hubble_tmp = H0_start;
	Hubble_param = H0_start;
	if(COSMOLOGY == 0)
	{
		Hubble_param = 0;
	}
	printf("The simulation is starting...\n");
	REAL T_prev,Hubble_param_prev;
	T_prev = T;
	Hubble_param_prev = Hubble_param;
	for(t=0; a_tmp<a_max; t++)
	{
		printf("\n\n----------------------------------------------------------------------------------------------\n");
		if(COSMOLOGY == 1)
                {
                        printf("Timestep %i, t=%.8fGy, h=%fGy, a=%.8f, H=%.8fkm/s/Mpc, z=%.8f:\n", t, T*47.1482347621227, h*47.1482347621227, a, Hubble_param*20.7386814448645, a_max/a-1.0);
                }
                else
                {
                        printf("Timestep %i, t=%f, h=%f, a=%f:\n", t, T, h, a);
                }
		Hubble_param_prev = Hubble_param;
		T_prev = T;
		T = T+h;
		step(x, F);
		Log_write(x);	//Writing logfile

		if(T > t_next)
		{
			kiiras(x);
			t_next=t_next+Delta_T_out;
			if(COSMOLOGY == 1)
			{
				printf("t = %f Gy\n\th=%f Gy\n", T*47.1482347621227, h*47.1482347621227);
			}
			else
			{
				printf("t = %f\n\terr_max = %e\th=%f\n", T, errmax, h);
			}
		}

		//Changing timestep length
		h_prev = h;
		if(h<h_max || h>h_min || mean_err/errmax>1)
		{
			h = (double) pow(2*mean_err*beta/errmax, 0.5);
		}

		if(h<h_min)
		{
			h=h_min;
		}
		if(h>h_max)
		{
			h=h_max;
		}
	}
	kiiras(x); //writing output
	printf("\n\n----------------------------------------------------------------------------------------------\n");
	printf("The simulation ended. The final state:\n");
	if(COSMOLOGY == 1)
        {
	printf("Timestep %i, t=%.8fGy, h=%f, a=%.8f, H=%.8f, z=%.8f\n", t, T*47.1482347621227, h*47.1482347621227, a, Hubble_param*20.7386814448645, a_max/a-1.0);

	double a_end, b_end;
	a_end = (Hubble_param - Hubble_param_prev)/(a-a_prev);
	b_end = Hubble_param_prev-a_end*a_prev;
	double H_end = a_max*a_end+b_end;
	a_end = (T - T_prev)/(a-a_prev);
        b_end = T_prev-a_end*a_prev;
	double T_end = a_max*a_end+b_end;
	printf("\nAt a = %lf state, with linear interpolation:\n",a_max);
	printf("t=%.8fGy, a=%.8f, H=%.8fkm/s/Mpc\n\n", T_end*47.1482347621227, a_max, H_end*20.7386814448645);
	}
	else
	{
		printf("Timestep %i, t=%f, h=%f, a=%f:\n", t, T, h, a);
	}
	printf("Running time of the simulation:\n");
	//Timing
	REAL SIM_end_time = (REAL) clock () / (REAL) CLOCKS_PER_SEC;
	REAL SIM_omp_end_time = omp_get_wtime();
	//Timing
	printf("CPU time = %lfs\n", SIM_end_time-SIM_start_time);
	printf("RUN time = %lfs\n", SIM_omp_end_time-SIM_omp_start_time);
	return 0;
}
