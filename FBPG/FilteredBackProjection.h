#ifndef FILTEREDBACKPROJECTION_H
#define FILTEREDBACKPROJECTION_H

#include "cuda_runtime.h"
#include "mucat_file.h"
#include "kernel.cuh"

class FilteredBackProjection
{
public:
	FilteredBackProjection();
	~FilteredBackProjection();


	void SetSinogramProperties(float *_sinogram, const int &_sinogram_X, const float &_centre_X, const int &_sinogram_Y, const int &_sinogram_Z, const float &_centre_Z, const float &_delta_fi, float _distance, float _sampling);

	void SetReconstructionProperties(bin_siz *_output_file, const int &_reconstruction_X, const float &_centre_X, const int &_reconstruction_Y, const float &_centre_Y, const int &_reconstruction_Z, const float &_centre_Z, const float &_rot);

	void SetAxisDirection(const bool &_X = false, const bool &_Y = false, const bool &_Z = false);

	void SetGPUid(const int &_ID = 0) { devID = _ID; };
	void SetBlockSize(const int &_sz = 16) { size_of_Z_block = _sz; };
	void SetAllocationLimit(const long long int &_memory_limit = (long long int)1024 * 1024 * 1024) { memory_limit = _memory_limit; };
	void SetCollecting(const bool &_collecting = false) { collecting = _collecting; };
	void SetShowInfo(const bool &_show = false) { show = _show; };
	void SetScaleXYZ(const float &_scale) { recn.scale = _scale; };

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

#endif //FILTEREDBACKPROJECTION_H
