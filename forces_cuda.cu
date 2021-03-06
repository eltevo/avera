#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <omp.h>
#include <time.h>
#include "global_variables.h"
#include "cuda_runtime.h"
#include "device_launch_parameters.h"

#define BLOCKSIZE 256

#ifdef USE_SINGLE_PRECISION
	typedef float REAL;
#else
	typedef double REAL;
#endif

extern int H[2202][4];
extern int e[2202][4];
extern REAL w[3];
extern REAL G, M;
extern int N, hl, el;

int ewald_space(REAL R, int ewald_index[2102][4]);


cudaError_t forces_old_cuda(REAL**x, REAL**F);
cudaError_t forces_old_periodic_cuda(REAL**x, REAL**F);

void forces_old(REAL**x, REAL**F)
{
	forces_old_cuda(x, F);
	return;
}

void forces_old_periodic(REAL**x, REAL**F)
{
	forces_old_periodic_cuda(x, F);
	return;
}


void recalculate_softening();

__global__ void ForceKernel_old(int n, int N, const REAL *xx, const REAL *xy, const REAL *xz, REAL *F, int IS_PERIODIC, REAL M, REAL L, REAL *SOFT_CONST, REAL beta)
{
	REAL Fx_tmp, Fy_tmp, Fz_tmp;
	REAL r, dx, dy, dz, wij;
	int i, j, id;

	id = blockIdx.x * blockDim.x + threadIdx.x;
	Fx_tmp = Fy_tmp = Fz_tmp = 0.0;
	for (i = id; i<N; i+=n)
		{
			for (j = 0; j<N; j++)
			{
				//calculating particle distances
				dx = (xx[j] - xx[i]);
				dy = (xy[j] - xy[i]);
				dz = (xz[j] - xz[i]);
				//in this function we use only the nearest image
				if (IS_PERIODIC == 1)
				{
					if (fabs(dx)>0.5*L)
						dx = dx - L*dx / fabs(dx);
					if (fabs(dy)>0.5*L)
						dy = dy - L*dy / fabs(dy);
					if (fabs(dz)>0.5*L)
						dz = dz - L*dz / fabs(dz);
				}

				r = sqrt(pow(dx, 2) + pow(dy, 2) + pow(dz, 2));
				wij = 0.0;
				if (r <= beta)
				{
					wij = M*(SOFT_CONST[0] * pow(r, 3) + SOFT_CONST[1] * pow(r, 2) + SOFT_CONST[2]);
				}
				if (r > beta && r <= 2 * beta)
				{
					wij = M*(SOFT_CONST[3] * pow(r, 3) + SOFT_CONST[4] * pow(r, 2) + SOFT_CONST[5] * r + SOFT_CONST[6] + SOFT_CONST[7] / pow(r, 3));
				}
				if (r > 2 * beta)
				{
					wij = M / (pow(r, 3));
				}
				Fx_tmp += wij*(dx);
				Fy_tmp += wij*(dy);
				Fz_tmp += wij*(dz);

			}
			F[3*i] += Fx_tmp;
			F[3*i+1] += Fy_tmp;
			F[3*i+2] += Fz_tmp;
			Fx_tmp = Fy_tmp = Fz_tmp = 0.0;
			
		}
}

__global__ void ForceKernel_old_periodic(int n, int N, const REAL *xx, const REAL *xy, const REAL *xz, REAL *F, int IS_PERIODIC, REAL M, REAL L, REAL *SOFT_CONST, REAL beta, int *e, int el)
{
	REAL Fx_tmp, Fy_tmp, Fz_tmp;
	REAL r, dx, dy, dz, wij;
	int i, j, m, id;
	id = blockIdx.x * blockDim.x + threadIdx.x;
	Fx_tmp = Fy_tmp = Fz_tmp = 0;
	for (i = id; i<N; i=i+n)
	{
		for (j = 0; j<N; j++)
		{
			//calculating particle distances
			dx = (xx[j] - xx[i]);
			dy = (xy[j] - xy[i]);
			dz = (xz[j] - xz[i]);
			//in this function we use multiple images
			for (m = 0; m < 3*el; m = m+3)
			{
				r = sqrt(pow((dx - ((REAL)e[m])*L), 2) + pow((dy - ((REAL)e[m+1])*L), 2) + pow((dz-((REAL)e[m+2])*L), 2));
				wij = 0.0;
				if (r <= beta)
				{
					wij = M*(SOFT_CONST[0] * pow(r, 3) + SOFT_CONST[1] * pow(r, 2) + SOFT_CONST[2]);
				}
				if (r > beta && r <= 2 * beta)
				{
					wij = M*(SOFT_CONST[3] * pow(r, 3) + SOFT_CONST[4] * pow(r, 2) + SOFT_CONST[5] * r + SOFT_CONST[6] + SOFT_CONST[7] / pow(r, 3));
				}
				if (r > 2 * beta && r < 2.6*L)
				{
					wij = M / (pow(r, 3));
				}
				if (wij != 0)
				{
					Fx_tmp += wij*(dx - ((REAL)e[m])*L);
					Fy_tmp += wij*(dy - ((REAL)e[m + 1])*L);
					Fz_tmp += wij*(dz - ((REAL)e[m + 2])*L);
				}
			}
			
		}
		F[3 * i] += Fx_tmp;
		F[3 * i + 1] += Fy_tmp;
		F[3 * i + 2] += Fz_tmp;
		Fx_tmp = Fy_tmp = Fz_tmp = 0;
	}

}


void recalculate_softening()
{
	printf("Calculating smoothing length...\n");
	beta = ParticleRadi;
	if(COSMOLOGY ==1)
	{
	beta = ParticleRadi*(REAL)(a_max/a);
	}
	SOFT_CONST[0] = 32.0/pow(2.0*beta, 6);
        SOFT_CONST[1] = -38.4/pow(2.0*beta, 5);
        SOFT_CONST[2] = 32.0/(3.0*pow(2.0*beta, 3));

        SOFT_CONST[3] = -32.0/(3.0*pow(2*beta, 6));
        SOFT_CONST[4] = 38.4/pow(2.0*beta, 5);
        SOFT_CONST[5] = -48.0/pow(2.0*beta, 4);
        SOFT_CONST[6] = 64.0/(3.0*pow(2.0*beta, 3));
        SOFT_CONST[7] = -1.0/15.0;
	printf("The new smoothing length:\tbeta = %lf\n", beta);
}

cudaError_t forces_old_cuda(REAL**x, REAL**F) //Force calculation on GPU
{
	int i, j;
	int mprocessors;
	cudaDeviceGetAttribute(&mprocessors, cudaDevAttrMultiProcessorCount, 0);
	printf("GPU force calculation.\n Number of threads: %i\n", 32*mprocessors*BLOCKSIZE);
	REAL start_time, end_time;
	REAL *xx_tmp, *xy_tmp, *xz_tmp, *F_tmp;
	REAL *dev_xx= 0;
	REAL *dev_xy= 0;
	REAL *dev_xz= 0;
	REAL *dev_F = 0;
	REAL *dev_SOFT_CONST;

	//Checking for the GPU
	cudaError_t cudaStatus;
	cudaStatus = cudaSetDevice(0);
	if (cudaStatus != cudaSuccess) {
		fprintf(stderr, "cudaSetDevice failed!  Do you have a CUDA-capable GPU installed?");
		goto Error;
	}
	//timing
	start_time = (REAL) clock () / (REAL) CLOCKS_PER_SEC;
	//timing
	// Allocate GPU buffers for coordinate vectors
	cudaStatus = cudaMalloc((void**)&dev_xx, N * sizeof(REAL));
	if (cudaStatus != cudaSuccess) {
		fprintf(stderr, "cudaMalloc failed!");
		goto Error;
	}
	cudaStatus = cudaMalloc((void**)&dev_xy, N * sizeof(REAL));
        if (cudaStatus != cudaSuccess) {
                fprintf(stderr, "cudaMalloc failed!");
                goto Error;
        }
	cudaStatus = cudaMalloc((void**)&dev_xz, N * sizeof(REAL));
        if (cudaStatus != cudaSuccess) {
                fprintf(stderr, "cudaMalloc failed!");
                goto Error;
        }
	// Allocate GPU buffers for force vectors
	cudaStatus = cudaMalloc((void**)&dev_F, 3 * N * sizeof(REAL));
	if (cudaStatus != cudaSuccess) {
		fprintf(stderr, "cudaMalloc failed!");
		goto Error;
	}
	// Allocate GPU buffers for SOFT_CONST vector
	cudaStatus = cudaMalloc((void**)&dev_SOFT_CONST, 8 * sizeof(REAL));
	if (cudaStatus != cudaSuccess) {
		fprintf(stderr, "cudaMalloc failed!");
		goto Error;
	}
	//Converting the Nx3 matrix to 3Nx1 vector.
	xx_tmp = (REAL*)malloc(N*sizeof(REAL));
	xy_tmp = (REAL*)malloc(N*sizeof(REAL));
	xz_tmp = (REAL*)malloc(N*sizeof(REAL));
	F_tmp = (REAL*)malloc(3 * N*sizeof(REAL));
	for(i = 0; i < N; i++)
	{
		xx_tmp[i] = x[i][0];
		xy_tmp[i] = x[i][1];
		xz_tmp[i] = x[i][2];
		for(j=0; j<3; j++)
			F_tmp[3*i + j] = 0.0f;
	}
	// Copy input vectors from host memory to GPU buffers.
	cudaStatus = cudaMemcpy(dev_xx, xx_tmp, N * sizeof(REAL), cudaMemcpyHostToDevice);
	if (cudaStatus != cudaSuccess) {
		fprintf(stderr, "cudaMemcpy in failed!");
		goto Error;
	}
	cudaStatus = cudaMemcpy(dev_xy, xy_tmp, N * sizeof(REAL), cudaMemcpyHostToDevice);
        if (cudaStatus != cudaSuccess) {
                fprintf(stderr, "cudaMemcpy in failed!");
                goto Error;
        }
	cudaStatus = cudaMemcpy(dev_xz, xz_tmp, N * sizeof(REAL), cudaMemcpyHostToDevice);
        if (cudaStatus != cudaSuccess) {
                fprintf(stderr, "cudaMemcpy in failed!");
                goto Error;
        }
	cudaStatus = cudaMemcpy(dev_F, F_tmp, 3 * N * sizeof(REAL), cudaMemcpyHostToDevice);
	if (cudaStatus != cudaSuccess) {
		fprintf(stderr, "cudaMemcpy in failed!");
		goto Error;
	}
	cudaStatus = cudaMemcpy(dev_SOFT_CONST, SOFT_CONST, 8 * sizeof(REAL), cudaMemcpyHostToDevice);
	if (cudaStatus != cudaSuccess) {
		fprintf(stderr, "cudaMemcpy in failed!");
		goto Error;
	}
	// Launch a kernel on the GPU
	ForceKernel_old<<<32*mprocessors, BLOCKSIZE>>>(32 * mprocessors * BLOCKSIZE, N, dev_xx, dev_xy, dev_xz, dev_F, IS_PERIODIC, M, L, dev_SOFT_CONST, beta);
	
	// Check for any errors launching the kernel
	cudaStatus = cudaGetLastError();
	if (cudaStatus != cudaSuccess) {
		fprintf(stderr, "ForceKernel launch failed: %s\n", cudaGetErrorString(cudaStatus));
		goto Error;
	}
	// cudaDeviceSynchronize waits for the kernel to finish, and returns
	// any errors encountered during the launch.
	cudaStatus = cudaDeviceSynchronize();
	if (cudaStatus != cudaSuccess) {
		fprintf(stderr, "cudaDeviceSynchronize returned error code %d after launching ForceKernel!\n", cudaStatus);
		goto Error;
	}
	// Copy output vector from GPU buffer to host memory.
	cudaStatus = cudaMemcpy(F_tmp, dev_F, 3 * N * sizeof(REAL), cudaMemcpyDeviceToHost);
	if (cudaStatus != cudaSuccess) {
		fprintf(stderr, "cudaMemcpy out failed!");
		goto Error;
	}
	//converting back the 3Nx1 vector ot Nx3 matrix.
	for (i = 0; i < N; i++)
	{
		for (j = 0; j < 3; j++)
		{
			F[i][j] = F_tmp[3 * i + j];
		}
	}
	free(F_tmp);
	free(xx_tmp);
	free(xy_tmp);
	free(xz_tmp);
	//timing
	end_time = (REAL) clock () / (REAL) CLOCKS_PER_SEC;
	//timing
	printf("Force calculation finished.\n");
	printf("Force calculation GPU time = %lfs\n", end_time-start_time);
	cudaFree(dev_xx);
	cudaFree(dev_xy);
	cudaFree(dev_xz);
	cudaFree(dev_F);
	cudaFree(dev_SOFT_CONST);
	cudaThreadExit();
Error:
	cudaFree(dev_xx);
	cudaFree(dev_xy);
	cudaFree(dev_xz);
	cudaFree(dev_F);
	cudaFree(dev_SOFT_CONST);
	cudaThreadExit();

	return cudaStatus;
}


cudaError_t forces_old_periodic_cuda(REAL**x, REAL**F) //Force calculation with multiple images on GPU
{
	int i, j;
	int mprocessors;
	cudaDeviceGetAttribute(&mprocessors, cudaDevAttrMultiProcessorCount, 0);
	printf("GPU force calculation.\n Number of threads: %i\n", 32*mprocessors*BLOCKSIZE);
	REAL start_time, end_time;
	REAL *xx_tmp, *xy_tmp, *xz_tmp, *F_tmp;
	REAL *dev_xx= 0;
	REAL *dev_xy= 0;
	REAL *dev_xz= 0;
	REAL *dev_F = 0;
	REAL *dev_SOFT_CONST;
	int *dev_e;
	int e_tmp[6606];

	//Checking for the GPU
	cudaError_t cudaStatus;
	cudaStatus = cudaSetDevice(0);
	if (cudaStatus != cudaSuccess) {
		fprintf(stderr, "cudaSetDevice failed!  Do you have a CUDA-capable GPU installed?");
		goto Error;
	}
	//timing
	start_time = (REAL)clock() / (REAL)CLOCKS_PER_SEC;
	//timing
	// Allocate GPU buffers for coordinate vectors
	cudaStatus = cudaMalloc((void**)&dev_xx, N * sizeof(REAL));
	if (cudaStatus != cudaSuccess) {
		fprintf(stderr, "cudaMalloc failed!");
		goto Error;
	}
	cudaStatus = cudaMalloc((void**)&dev_xy, N * sizeof(REAL));
        if (cudaStatus != cudaSuccess) {
                fprintf(stderr, "cudaMalloc failed!");
                goto Error;
        }
	cudaStatus = cudaMalloc((void**)&dev_xz, N * sizeof(REAL));
        if (cudaStatus != cudaSuccess) {
                fprintf(stderr, "cudaMalloc failed!");
                goto Error;
        }
	// Allocate GPU buffers for force vectors
	cudaStatus = cudaMalloc((void**)&dev_F, 3 * N * sizeof(REAL));
	if (cudaStatus != cudaSuccess) {
		fprintf(stderr, "cudaMalloc failed!");
		goto Error;
	}
	// Allocate GPU buffers for SOFT_CONST vector
	cudaStatus = cudaMalloc((void**)&dev_SOFT_CONST, 8 * sizeof(REAL));
	if (cudaStatus != cudaSuccess) {
		fprintf(stderr, "cudaMalloc failed!");
		goto Error;
	}
	// Allocate GPU buffers for e matrix
	cudaStatus = cudaMalloc((void**)&dev_e, 6606 * sizeof(int));
	if (cudaStatus != cudaSuccess) {
		fprintf(stderr, "cudaMalloc failed!");
		goto Error;
	}
	//Converting e matrix into a vector
	for (i = 0; i < 2202; i++)
	{
		for (j = 0; j < 3; j++)
		{
			e_tmp[3 * i + j] = e[i][j];
		}
	}
	//Converting the Nx3 matrix to 3Nx1 vector.
        xx_tmp = (REAL*)malloc(N*sizeof(REAL));
        xy_tmp = (REAL*)malloc(N*sizeof(REAL));
        xz_tmp = (REAL*)malloc(N*sizeof(REAL));
        F_tmp = (REAL*)malloc(3 * N*sizeof(REAL));
        for(i = 0; i < N; i++)
        {
                xx_tmp[i] = x[i][0];
                xy_tmp[i] = x[i][1];
                xz_tmp[i] = x[i][2];
                for(j=0; j<3; j++)
                        F_tmp[3*i + j] = 0;
        }
        // Copy input vectors from host memory to GPU buffers.
        cudaStatus = cudaMemcpy(dev_xx, xx_tmp, N * sizeof(REAL), cudaMemcpyHostToDevice);
        if (cudaStatus != cudaSuccess) {
                fprintf(stderr, "cudaMemcpy in failed!");
                goto Error;
        }
        cudaStatus = cudaMemcpy(dev_xy, xy_tmp, N * sizeof(REAL), cudaMemcpyHostToDevice);
        if (cudaStatus != cudaSuccess) {
                fprintf(stderr, "cudaMemcpy in failed!");
                goto Error;
        }
        cudaStatus = cudaMemcpy(dev_xz, xz_tmp, N * sizeof(REAL), cudaMemcpyHostToDevice);
        if (cudaStatus != cudaSuccess) {
                fprintf(stderr, "cudaMemcpy in failed!");
                goto Error;
        }
	cudaStatus = cudaMemcpy(dev_F, F_tmp, 3 * N * sizeof(REAL), cudaMemcpyHostToDevice);
	if (cudaStatus != cudaSuccess) {
		fprintf(stderr, "cudaMemcpy F failed!");
		goto Error;
	}
	cudaStatus = cudaMemcpy(dev_SOFT_CONST, SOFT_CONST, 8 * sizeof(REAL), cudaMemcpyHostToDevice);
	if (cudaStatus != cudaSuccess) {
		fprintf(stderr, "cudaMemcpy SOFT_CONST failed!");
		goto Error;
	}
	cudaStatus = cudaMemcpy(dev_e, e_tmp, 6606 * sizeof(int), cudaMemcpyHostToDevice);
	if (cudaStatus != cudaSuccess) {
		fprintf(stderr, "cudaMemcpy e failed!");
		goto Error;
	}
	// Launch a kernel on the GPU with one thread for each element.
	ForceKernel_old_periodic << <32*mprocessors, BLOCKSIZE>> >(32*mprocessors * BLOCKSIZE, N, dev_xx, dev_xy, dev_xz, dev_F, IS_PERIODIC, M, L, dev_SOFT_CONST, beta, dev_e, el);

	// Check for any errors launching the kernel
	cudaStatus = cudaGetLastError();
	if (cudaStatus != cudaSuccess) {
		fprintf(stderr, "ForceKernel_periodic launch failed: %s\n", cudaGetErrorString(cudaStatus));
		goto Error;
	}
	// cudaDeviceSynchronize waits for the kernel to finish, and returns
	// any errors encountered during the launch.
	cudaStatus = cudaDeviceSynchronize();
	if (cudaStatus != cudaSuccess) {
		fprintf(stderr, "cudaDeviceSynchronize returned error code %d after launching ForceKernel_periodic!\n", cudaStatus);
		goto Error;
	}
	// Copy output vector from GPU buffer to host memory.
	cudaStatus = cudaMemcpy(F_tmp, dev_F, 3 * N * sizeof(REAL), cudaMemcpyDeviceToHost);
	if (cudaStatus != cudaSuccess) {
		fprintf(stderr, "cudaMemcpy out failed!");
		goto Error;
	}
	//converting back the 3Nx1 vector ot Nx3 matrix.
	for (i = 0; i < N; i++)
	{
		for (j = 0; j < 3; j++)
		{
			F[i][j] = F_tmp[3 * i + j];
		}
	}
	free(F_tmp);
	free(xx_tmp);
	free(xy_tmp);
	free(xz_tmp);
	//timing
	end_time = (REAL)clock() / (REAL)CLOCKS_PER_SEC;
	//timing
	printf("Force calculation finished.\n");
	printf("Force calculation GPU time = %lfs\n", end_time - start_time);
	cudaFree(dev_xx);
	cudaFree(dev_xy);
	cudaFree(dev_xz);
	cudaFree(dev_F);
	cudaFree(dev_e);
	cudaFree(dev_SOFT_CONST);
	cudaThreadExit();
Error:
	cudaFree(dev_xx);
	cudaFree(dev_xy);
	cudaFree(dev_xz);
	cudaFree(dev_F);
	cudaFree(dev_e);
	cudaFree(dev_SOFT_CONST);
	cudaThreadExit();

	return cudaStatus;
}
