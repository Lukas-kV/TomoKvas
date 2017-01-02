#ifndef KERNEL_H
#define KERNEL_H

#include <time.h>
#include "cuda_runtime.h"
#include "device_launch_parameters.h"
//#include "math_functions.h"

#include <stdio.h>
#include <iostream>
#include <algorithm>
#include <vector>
#include <sstream>

#include "cufft.h"


using namespace std;

//#define PI (2 * asin(1.))
#define PI 3.14159265358979323846

typedef struct _rec_par
{
	float *dev_reconstruction;
	size_t rec_p;
	size_t rec_X;
	float dX;
	size_t rec_Y;
	float dY;
	size_t rec_Z;
	float dZ;
	int ZBlockSize;
	float rot; /// otoceni rekonstrukce v rad
	unsigned int dirX, dirY, dirZ; /// smer os
	float scale;

} rec_par, *prec_par;

const int POINTERS = 10;

typedef struct _sin_par
{
	float *dev_sinogram[POINTERS];
	size_t sin_in_block;
	float *dev_dsin;
	float *dev_dcos;
	size_t sin_p;
	size_t sin_X;
	float cX;
	size_t sin_Y;
//	float cY;
	size_t sin_Z;
	float cZ;
	float delta_f;
	float scaling; // FFT normalisation = numer of samples in FFT
	float distance;
	float sampling;
} sin_par, *psin_par;


string GetDeltaTime( clock_t _start_time, clock_t _now);

//***** device functions ******//

__device__ float fQ3D(float _p, float _f, size_t _beta, size_t _sin_p, size_t _sin_X, size_t _sin_Y, size_t _sin_Z, float _cZ, float _cX, float **_dev_sinogram, float _delta_f , int sin_per_block); ///< interpolace a vypocet hodnoty ze sinogramu

//***** kernels *****//

cudaError_t RunAbsKernel(size_t _data_X, cufftComplex *_dev_data, float *_dev_filter, float _sampling);
cudaError_t RunFiltrationKernel(size_t _data_X, size_t _data_Y, cufftComplex *_dev_data, float *_dev_filter);
cudaError_t RunFBPConeKernel(rec_par reconstruction, sin_par sinogram);
cudaError_t RunRB_comma3DKernel(sin_par sinogram);
cudaError_t RunCollectionKernel(rec_par reconstruction, float *_dev_collection, size_t _collection_p);


#endif //KERNEL_H