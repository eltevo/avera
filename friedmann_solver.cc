#include <stdio.h>
#include <math.h>
#include "global_variables.h"

#ifdef USE_SINGLE_PRECISION
typedef float REAL;
#else
typedef double REAL;
#endif

//Friedman-equation integrator. We use 4th order Runge-Kutta integrator

double friedman_solver_step(double a0, double h, double Omega_lambda, double Omega_r, double Omega_m, double Omega_k, double H0);

double friedmann_solver_start(double a0, double t0, double h, double Omega_lambda, double Omega_r, double Omega_m, double H0, double a_start)
{
	//In this function we calculate the time of Big Bang, and the initial time for the simulation.
	printf("Calculating time for the initial scale factor...\n");
	double b, b_tmp, t_cosm_tmp;
	double t_cosm = t0;
	double t_start;
	double t_start_err=1e-20;
	double h_var;
	b = a0;	//t=1-ben a skálafaktor
	//double t_bigbang; //???
	printf("h=%eGy\n", h*47.1482347621227);
	double Omega_k = 1.-Omega_m-Omega_lambda-Omega_r;
	//Solving the "da/dt = a*H0*sqrt(Omega_m*pow(a, -3)+Omega_r*pow(a, -4)+Omega_lambda)" differential equation
	b_tmp = b;
	while(0<b)
	{
		b_tmp = b;
		b = friedman_solver_step(b, -h, Omega_lambda, Omega_r, Omega_m, Omega_k, H0);
		t_cosm -= h;
	}
	t_bigbang=t_cosm+h; //rough estimation for t_bibgang.
	b = b_tmp;
	printf("First guess: %.12f Gy\n\n", -t_bigbang*47.1482347621227);
	//Searching for t_start.
	h_var = -0.5*h;
	while(fabs(h_var)>t_start_err)
	{
		b_tmp = b;
		b = friedman_solver_step(b, h_var, Omega_lambda, Omega_r, Omega_m, Omega_k, H0);
		t_cosm_tmp = t_cosm;
		t_cosm=t_cosm+h_var;
		if(b>0)
		{
			//printf("After BB: t=%.15e\th_var=%e\n", t_cosm*47.1482347621227, h_var);
			//h_var=0.5*h_var;
		}
		else
		{
			//printf("Before BB: t=%.15e\th_var=%e\n", t_cosm*47.1482347621227, h_var);
			b = b_tmp;
			t_cosm = t_cosm_tmp;
			h_var=0.5*h_var;
		}
	}
	t_bigbang = t_cosm;
	//Setting t=0 to Big Bang:
	t_start =  -1.0 * t_bigbang;
	return t_start;
}

double friedman_solver_step(double a0, double h, double Omega_lambda, double Omega_r, double Omega_m, double Omega_k, double H0)
{
	int collapse;
	double b,k1,k2,k3,k4,K;
	double j,l,m,n;

	a_prev = a0;
	if(H0>0)
	{
		collapse=0;
	}
	else
	{
		collapse=1;
	}
	b = a0;
	if(fabs(Omega_k) < 1e-8)
	{
		k1 = b*H0*sqrt(Omega_m*pow(b, -3.0)+Omega_r*pow(b, -4.0)+Omega_lambda);
		k2 = (b+h*k1/2.0)*H0*sqrt(Omega_m*pow((b+h*k1/2.0), -3.0)+Omega_r*pow((b+h*k1/2), -4.0)+Omega_lambda);
		k3 = (b+h*k2/2.0)*H0*sqrt(Omega_m*pow((b+h*k2/2.0), -3.0)+Omega_r*pow((b+h*k2/2), -4.0)+Omega_lambda);
		k4 = (b+h*k3)*H0*sqrt(Omega_m*pow((b+h*k3), -3.0)+Omega_r*pow((b+h*k3), -4.0)+Omega_lambda);
		K = h*(k1+k2*2.0+k3*2.0+k4)/6.0;
		b += K;
	}
	else
	{
		j= Omega_m*pow(b, -3.0)+Omega_r*pow(b, -4.0)+Omega_lambda+Omega_k*pow(b, -2.0);
		k1 = b*H0*sqrt(fabs(j));
		l = Omega_m*pow((b+h*k1/2.0), -3.0)+Omega_r*pow((b+h*k1/2.0), -4.0)+Omega_lambda+Omega_k*pow((b+h*k1/2.0), -2.0);
		k2 = (b+h*k1/2.0)*H0*sqrt(fabs(l));
		m = Omega_m*pow((b+h*k2/2.0), -3.0)+Omega_r*pow((b+h*k2/2.0), -4.0)+Omega_lambda+Omega_k*pow((b+h*k2/2.0), -2.0);
		k3 = (b+h*k2/2.0)*H0*sqrt(fabs(m));
		n=Omega_m*pow((b+h*k3), -3.0)+Omega_r*pow((b+h*k3), -4.0)+Omega_lambda+Omega_k*pow((b+h*k3), -2.0);
		k4 = (b+h*k3)*H0*sqrt(fabs(n));
		if(j<0&&l<0&&m<0&&n<0)
		{
			collapse = 1;
		}

		if(collapse==0)
		{
			K = h*(k1+k2*2.0+k3*2.0+k4)/6.0;
		}
		else
		{
			K = -h*(k1+k2*2.0+k3*2.0+k4)/6.0;
		}
		b = b + K;
	}
	return b;
}

void rescaling()
{
	int i;
	printf("Cosmological parameters before rescaling:\n");
	printf("a_start=\t%lf\na_max=\t\t%lf\nH0=\t\t%lf\nOmega_b=\t%lf\nOmega_dm=\t%lf\nOmega_m=\t%lf\nOmega_lambda=\t%lf\nOmega_k=\t%lf\nM=\t\t%lf\n", a_start, a_max, H0*20.7386814448645, Omega_b, Omega_dm, Omega_m, Omega_lambda, Omega_k, M);
	printf("\nRescaling...\n");
	//Saving old cosmological parameters
	double rho_crit_old = rho_crit;
	double Omega_b_old = Omega_b;
	double Omega_dm_old = Omega_dm;
	double Omega_m_old = Omega_m;
	double Omega_r_old = Omega_r;
	double Omega_lambda_old = Omega_lambda;
	double Omega_k_old = Omega_k;
	double a_start_old = a_start;
	double a_max_old = a_max;
	double H0_old = H0;
	//Here we rescale the cosmological parameters, so that a = 1 will be a_start
	//We must start this with rho_crit:
	rho_crit = rho_crit_old*(Omega_m_old*pow((a_max/a_start), 3.0) + Omega_r_old*pow((a_max/a_start), 4.0) + Omega_lambda_old + Omega_k_old*pow((a_max/a_start), 2.0));
	//The new Omega parameters:
	Omega_m = rho_crit_old*Omega_m_old*pow((a_max/a_start), 3.0)/rho_crit;
	Omega_b = rho_crit_old*Omega_b_old*pow((a_max/a_start), 3.0)/rho_crit;
	Omega_dm = rho_crit_old*Omega_dm_old*pow((a_max/a_start), 3.0)/rho_crit;
	Omega_r = rho_crit_old*Omega_r_old*pow((a_max/a_start), 4.0)/rho_crit;
	Omega_lambda = rho_crit_old*Omega_lambda_old/rho_crit;
	Omega_k = 1.0-Omega_m-Omega_lambda-Omega_r;

	//H0 and a_start:
	a_max = a_max_old/a_start_old;
	a_start = 1.0;
	a = 1.0;
	if(NONISOTROPIC_EXPANSION == 0)
	{
		H0 = H0_old*sqrt(Omega_m_old*pow(a_max/a_start,3.0) + Omega_r_old*pow(a_max/a_start,4.0) + Omega_lambda_old + Omega_k_old*pow(a_max/a_start,2.0));
		printf("The rescaled Hubble-constant: H(a=%f) =  %f km/s/Mpc\n", a_start, H0*20.7386814448645);
		M = Omega_dm*rho_crit*pow(L*(a_start/a_max), 3.0)/((double) N);
	}
	else
	{
		H0 = H0_start; //it is in the parameter file.
		M = Omega_dm*rho_crit*pow(L*(a_start/a_max), 3.0)/((double) N);
		printf("The particle masses with the rescaled critical density:\nM=\t\t%lf\n", M);
		//Because we set a different H0, than we must change the critcal density and the particle masses
		rho_crit = 3*pow(H0,2)/(8*pi*G);
		M = Omega_dm*rho_crit*pow(L*(a_start/a_max), 3.0)/((double) N);
		printf("With H0_start:\n The critical density:\trho_crit=%lf 10e11M_sol/Mpc^3\n particles mass:\tM=%lf\n", rho_crit, M);
		
	}
	//Rescaling speeds. If one uses Gadget format: http://wwwmpa.mpa-garching.mpg.de/gadget/gadget-list/0113.html
	if(RESTART == 0)
	{
	for(i=0;i<N;i++)
	{
		x[i][3] = x[i][3]*sqrt(a_max/a_start);
		x[i][4] = x[i][4]*sqrt(a_max/a_start);
		x[i][5] = x[i][5]*sqrt(a_max/a_start);
	}
	}
	printf("Cosmological parameters after the rescaling:\n");
        printf("a_start=\t%lf\na_max=\t\t%lf\nH0=\t\t%lfkm/s/Mpc\nOmega_b=\t%lf\nOmega_dm=\t%lf\nOmega_m=\t%lf\nOmega_lambda=\t%lf\nOmega_k=\t%lf\nM=\t\t%lf*10e+11M_sol\n", a_start, a_max, H0*20.7386814448645, Omega_b, Omega_dm, Omega_m, Omega_lambda, Omega_k, M);
	return;	
}

double CALCULATE_decel_param(double a, double a_prev1, double a_prev2, double h, double h_prev)
{
	double Decel_param_out;
	if(NONISOTROPIC_EXPANSION == 0)
	{
		double SUM_OMEGA_TMP = Omega_m*pow((a/a_start), -3.0) + Omega_r*pow((a/a_start), -4.0) + Omega_lambda + Omega_k*pow((a/a_start), -2.0);
		double Omega_m_tmp = Omega_m*pow((a/a_start), -3.0)/SUM_OMEGA_TMP;
		double Omega_r_tmp = Omega_r*pow((a/a_start), -4.0)/SUM_OMEGA_TMP;
		double Omega_lambda_tmp = Omega_lambda/SUM_OMEGA_TMP;
		double Omega_sum_tmp = Omega_m_tmp + Omega_r_tmp + Omega_lambda_tmp;
		Decel_param_out = Omega_sum_tmp/2 + Omega_r_tmp - 1.5*Omega_lambda_tmp;
	}
	else
	{
		//Using the first and the second numeric derivative, we calculate the deceleration parameter:
		double diff2_a = 2*(h*a_prev2-(h+h_prev)*a_prev1+h_prev*a)/(h*h*h_prev+h*h_prev*h_prev);
		double diff1_a = -1*(-1*h*h*a_prev2+pow((h+h_prev),2)*a_prev1 - (2*h*h_prev+pow(h_prev,2))*a)/(h*h*h_prev+h*h_prev*h_prev);
		Decel_param_out = -1*(diff2_a*a)/pow(diff1_a, 2);
	}
	return Decel_param_out;
}
