#pragma once

#include <time.h>
#include "cuda_runtime.h"
#include "device_launch_parameters.h"
//#include "math_functions.h"

#include <stdio.h>
#include <iostream>
#include <algorithm>
#include <vector>
#include <sstream>

#include "mucat_file.h"

#include "cufft.h"
//#include "sm_13_float_functions.h"


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

#define POINTERS 10

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

__global__ void AbsKernel(size_t _data_X, cufftComplex *_dev_data, float *_dev_filter, float _sampling);

__global__ void FiltrationKernel(size_t _data_X, size_t _data_Y, cufftComplex *_dev_data, float *_dev_filter);

__global__ void FBPConeKernel(rec_par reconstruction, sin_par sinogram);

__global__ void RB_comma3DKernel(size_t _sinogram_X, size_t _sinogram_Y, size_t _sinogram_Z, int _sZ, float *_dev_sinogram, float _distance);

__global__ void CollectionKernel(rec_par reconstruction, float *_dev_collection, size_t _collection_p);

class FilteredBackProjection 
{
public:
	FilteredBackProjection();
	~FilteredBackProjection();


	void SetSinogramProperties( float *_sinogram, const int &_sinogram_X, const float &_centre_X, const int &_sinogram_Y, const int &_sinogram_Z, const float &_centre_Z, const float &_delta_fi, float _distance, float _sampling);
	//void SetReconstructionProperties(float *_reconstruction, const int &_reconstruction_X, const float &_centre_X, const int &_reconstruction_Y, const float &_centre_Y, const int &_reconstruction_Z, const float &_centre_Z, const float &_rot);

	void SetReconstructionProperties(bin_siz *_output_file, const int &_reconstruction_X, const float &_centre_X, const int &_reconstruction_Y, const float &_centre_Y, const int &_reconstruction_Z, const float &_centre_Z, const float &_rot);

	void SetAxisDirection(const bool &_X = false, const bool &_Y = false, const bool &_Z = false);

	void SetGPUid(const int &_ID = 0);
	void SetBlockSize(const int &_sz = 16);
	void SetAllocationLimit(const long long int &_memory_limit = (long long int)1024 * 1024 * 1024);
	void SetCollecting(const bool &_collecting = false);
	void SetShowInfo(const bool &_show = false);
	void SetScaleXYZ(const float &_scale);

	cudaError_t Initialize();
	cudaError_t CalcSinCos();
	cudaError_t CalculateReconstruction();
	void fShowInfo();

	static size_t GetFreeMemInfo(int _GPUId);

private:
	int devID; ///< GPU ID
	cudaDeviceProp deviceProps; ///< GPU properties
	cudaError_t cudaStatus;

	rec_par recn;
	sin_par sing;

	float *dev_collecting;
	size_t collecting_p; 
	bool collecting; // zapne kolektovani dat
	bool show; // zapne zobrazeni dat

	float *dev_block;
	float *sinogram;
	float dZ;
	float *reconstruction;

	int size_of_Z_block;
	int nr_of_Z_blocks;
	float dcorZ;

	unsigned long long memory_limit; ///< velikost limitni pameti
	unsigned long long sinograms_per_limit; ///< pocet sinogramu na limitni pamet
	unsigned long long memory_sinograms; ///< velikost vsech sinogramu v B
	unsigned long long memory_sin; ///< velikost jedineho sinogramu v B
	unsigned long long memory_rest; ///< zbytek sinogramu pres celociselny pocet limitnich oblasti
	int  sinograms_per_rest; ///< pocet sinogramu ve zbyvajici pameti

	int nr_of_sin_blocks;

	cudaError_t fRBcomma();
	cudaError_t fConvolution();
	void FillFilter(float *_filter, int size);
	float RamLak(int _value, float _sm);

	clock_t start, global_start, pr;

	size_t free, total;

	bin_siz *output_file;
};

