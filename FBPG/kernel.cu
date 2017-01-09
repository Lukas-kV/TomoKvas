#include "kernel.cuh"

std::string GetDeltaTime( clock_t _start_time, clock_t _now)
{
	double delta = _now - _start_time;	
	int t = 0;
	std::ostringstream os;	

	if(delta > 1000)
		t = 1;		

	if(delta > 1000 * 60)
		t = 2;
		

	switch (t)
	{
		case 0:
			{
				os << delta << " ms";
				break;
			}
		case 1:
			{
				os << delta / 1000. << " s";
				break;
			}
		case 2:
			{
				os << delta / 1000. / 60. << " min";
				break;	
			}
	}
	std::string result(os.str());
	return result;
}

///* convolution *///

__global__ void FiltrationKernel(size_t _data_X, size_t _data_Y, cufftComplex *_dev_data, float *_dev_filter)
{
	unsigned int irow = blockIdx.x * blockDim.x + threadIdx.x; 

	if (irow < _data_Y)
	{
		cufftComplex *row = (cufftComplex* )((char*)_dev_data + (unsigned long long)irow * _data_X * sizeof(_dev_data[0]));
		for(unsigned int x = 0; x < _data_X; x++)
		{
			row[x].x *= _dev_filter[x];
			row[x].y *= _dev_filter[x];
		}
	}
}

__global__ void AbsKernel(size_t _data_X, cufftComplex *_dev_data, float *_dev_filter, float _sampling)
{
	unsigned int px = blockIdx.x * blockDim.x + threadIdx.x; 
	if (px < _data_X)
	{
		float X = _dev_data[px].x;
		float Y = _dev_data[px].y;
		_dev_filter[px] = _sampling * sqrt(X * X + Y * Y);
	}
}


/*** 3D reconstaructions *///

__global__ void FBPConeKernel(rec_par reconstruction, sin_par sinogram)
{
	unsigned int x = blockIdx.x * blockDim.x + threadIdx.x; 
	unsigned int y = blockIdx.y * blockDim.y + threadIdx.y;
	unsigned int z = blockIdx.z * blockDim.z + threadIdx.z;


	if (x < reconstruction.rec_X && y < reconstruction.rec_Y && z < reconstruction.ZBlockSize)
	{

		//float aX = x + reconstruction.dX; /// poloha vzhledem k pocatku ve stredu rekonstrukce
		//float aY = y + reconstruction.dY; 
		//float aZ = z + reconstruction.dZ; 

		float aX = reconstruction.scale * (x + reconstruction.dX); /// poloha vzhledem k pocatku ve stredu rekonstrukce
		float aY = reconstruction.scale * (y + reconstruction.dY); 
		float aZ = reconstruction.scale * (z + reconstruction.dZ); 

		(*(unsigned int *)&aX) ^= reconstruction.dirX;
		(*(unsigned int *)&aY) ^= reconstruction.dirY;
		(*(unsigned int *)&aZ) ^= reconstruction.dirZ;

		int beta; ///< celocislena beta = cislo radku v sinogramu

		float result = 0; ///< promena pro sumaci vysledku

		float U;
		float S;
		float T;
		float dsin;
		float dcos;

		if(sqrt(aX*aX + aY*aY) < sinogram.sin_X / 2. )
		{
			for( beta = 0; beta < sinogram.sin_Y; beta++)
			{

				dsin = sinogram.dev_dsin[beta];
				dcos = sinogram.dev_dcos[beta];
				S = aX * dcos + aY * dsin;
				T = aY * dcos - aX * dsin;
				U = sinogram.distance / (sinogram.distance - S);

				result +=  U * U * fQ3D( T * U, aZ * U, beta, sinogram.sin_p, sinogram.sin_X, sinogram.sin_Y, sinogram.sin_Z, sinogram.cZ, sinogram.cX, sinogram.dev_sinogram, sinogram.delta_f, sinogram.sin_in_block);
			}	
			result *= sinogram.delta_f * sinogram.scaling;
		}
	
		float* row = (float*)((char*)reconstruction.dev_reconstruction + (unsigned long long)reconstruction.rec_p * ( y + reconstruction.rec_Y * z)); /// vypocteme  ukazatel na radek rekonstrukce
		row[x] = result;
	}
}

__device__ float fQ3D(float _p, float _f, size_t _beta, size_t _sin_p, size_t _sin_X, size_t _sin_Y, size_t _sin_Z, float _cZ, float _cX, float **_dev_sinogram, float _delta_f , int sin_per_block) ///< interpolace a vypocet hodnoty ze sinogramu
{

	float scP = _p + _cX;
	float scF = _f + _cZ;

	if((scP < 0) || (scP > (_sin_X - 1)))  /// jestli je poloha v jenom radku mimo p rozmer sinogramu result += 0;
		return 0;		

	if((scF < 0) || ((scF) > (_sin_Z - 1)))  /// jestli je poloha v jenom radku mimo f rozmer sinogramu result += 0;
		return 0;		

	int af = ceil(scF); // nejblizsi vyzsi cele cislo v ose f pro vypocet radku
	int bf = floor(scF); // nejblizsi nizsi cele cislo v ose f 
	/// zde vybrat podle af a bf kterou pouzit promenou alokace

	int ablock = af / sin_per_block;
	int air = af - ablock * sin_per_block;

	int bblock = bf / sin_per_block;
	int bir = bf - bblock * sin_per_block;

	int a = ceil(scP); // nejblizsi vyzsi cele cislo v radku
	int b = floor(scP); // nejblizsi nizsi cele cislo v radku

	if(af != bf)
	{
		float va,vb;

		//* interpolace mezi radky
		float* rowa = (float*)((char*)(_dev_sinogram[ablock]) + (unsigned long long)_sin_p * (_beta  + _sin_Y * air)); /// vypocteme  ukazatel na radek sinogramu s uhlem beta
		float* rowb = (float*)((char*)(_dev_sinogram[bblock]) + (unsigned long long)_sin_p * (_beta + _sin_Y * bir)); 

		if(a != b)
		{
			va = (rowa[a] + (rowa[b] - rowa[a]) * (-scP + a));	
			vb = (rowb[a] + (rowb[b] - rowb[a]) * (-scP + a));	/// interpolace b - a = -1
//			vb = (rowb[a] + (rowb[b] - rowb[a]) * (_p + cP - a)/(b - a));	
		}
		else
		{
			va = rowa[a];
			vb = rowb[b];
		}

		return (va + (vb - va) * (-scF + af));		
	}
	else
	{
		float* row = (float*)((char*)(_dev_sinogram[ablock]) + (unsigned long long)_sin_p * (_beta  + _sin_Y * air)); 

		if(a != b)
			return (row[a] + (row[b] - row[a]) * (-scP + a));	

		return row[a];
	}
}

__global__ void RB_comma3DKernel(sin_par sinogram)
{
	unsigned int irow = blockIdx.x * blockDim.x + threadIdx.x; 
	unsigned int isin = blockIdx.y * blockDim.y + threadIdx.y; 

	if (irow < sinogram.sin_Y && isin < sinogram.sin_Z)
	{
		//int ce = (sinogram.sin_X - 1) / 2;
		float cs = isin - sinogram.cZ;
		float ce;

		int block = isin / sinogram.sin_in_block;
		int ir = isin - block * sinogram.sin_in_block;
	
		float *row = (float* )((char*)sinogram.dev_sinogram[block] +  (unsigned long long)sinogram.sin_p * (irow + sinogram.sin_Y * ir));

		for(int x = 0; x < sinogram.sin_X; x++)
		{
			ce = x - sinogram.cX;
			row[x] *= sinogram.distance / sqrt(ce*ce + sinogram.distance * sinogram.distance + cs*cs);	
		}
	}
}

__global__ void CollectionKernel(rec_par reconstruction, float *_dev_collection, size_t _collection_p)
{
	unsigned int x = blockIdx.x * blockDim.x + threadIdx.x; 
	unsigned int y = blockIdx.y * blockDim.y + threadIdx.y;
	unsigned int z = blockIdx.z * blockDim.z + threadIdx.z;

	if (x < reconstruction.rec_X && y < reconstruction.rec_Y && z < reconstruction.ZBlockSize)
	{
			float* R_row = (float*)((char*)reconstruction.dev_reconstruction + (unsigned long long)reconstruction.rec_p * (y + reconstruction.rec_Y * z)); /// vypocteme  ukazatel na radek rekonstrukce
			float* C_row = (float*)((char*)_dev_collection + (unsigned long long)_collection_p * (y + reconstruction.rec_Y * z)); /// vypocteme  ukazatel na radek rekonstrukce
			R_row[x] += C_row[x];
	}
}




cudaError_t RunAbsKernel(size_t _data_X, cufftComplex *_dev_data, float *_dev_filter, float _sampling)
{
	dim3 threadsPerBlock(16);
	dim3 numBlocks(ceil((float)_data_X / threadsPerBlock.x));
	AbsKernel <<<numBlocks, threadsPerBlock >>>(_data_X, _dev_data, _dev_filter, _sampling);

	cudaError_t cudaStatus = cudaDeviceSynchronize();
	if (cudaStatus != cudaSuccess)
	{		
		std::cerr << "convolution: filter calculation cudaDeviceSynchronize returned " << cudaGetErrorString(cudaStatus) << "after launching AbsKernel! " << std::endl;
	}
	return cudaStatus;
}

cudaError_t RunFiltrationKernel(size_t _data_X, size_t _data_Y, cufftComplex *_dev_data, float *_dev_filter)
{
	dim3 threadsPerBlock(16);
	dim3 numBlocks(ceil((float)_data_Y / threadsPerBlock.x));
	FiltrationKernel <<<numBlocks, threadsPerBlock >>>(_data_X, _data_Y, _dev_data, _dev_filter);

	cudaError_t cudaStatus = cudaDeviceSynchronize();
	if (cudaStatus != cudaSuccess)
	{
		std::cerr << "convolution: filtration cudaDeviceSynchronize returned " << cudaGetErrorString(cudaStatus) << "after launching FiltrationKernel! " << std::endl;
	}

	return cudaStatus;
}

cudaError_t RunFBPConeKernel(rec_par recn, sin_par sing)
{
	dim3 threadsPerBlock(8, 8, 8); ///< zvolim velikotak jst bloku 8 x 8 x 8 = 512 vlaken
	dim3 numBlocks(ceil((float)recn.rec_X / threadsPerBlock.x), ceil((float)recn.rec_Y / threadsPerBlock.y), ceil((float)recn.ZBlockSize / threadsPerBlock.z));  ///< spocitam potrebny pocet bloku na vypocet cele velikosti rekonstrukce
	FBPConeKernel <<<numBlocks, threadsPerBlock >>>(recn, sing);

	// cudaDeviceSynchronize waits for the kernel to finish, and returns
	// any errors encountered during the launch.
	cudaError_t cudaStatus = cudaDeviceSynchronize();
	if (cudaStatus != cudaSuccess)
	{
		std::cerr << "reconstruction: cudaDeviceSynchronize returned " << cudaGetErrorString(cudaStatus) << " after launching FBPKernel!" << std::endl;
	}
	return cudaStatus;
}

cudaError_t RunRB_comma3DKernel(sin_par sing)
{
	dim3 threadsPerBlock(16, 16);
	dim3 numBlocks(ceil((float)sing.sin_Y / threadsPerBlock.x), ceil((float)sing.sin_Z / threadsPerBlock.y));
	RB_comma3DKernel <<<numBlocks, threadsPerBlock >>>(sing);

	cudaError_t cudaStatus = cudaDeviceSynchronize();
	if (cudaStatus != cudaSuccess)
	{
		std::cerr << "RB': cudaDeviceSynchronize returned " << cudaGetErrorString(cudaStatus) << " after launching RB_comma3DKernel!" << std::endl;
	}

	return cudaStatus;
}

cudaError_t RunCollectionKernel(rec_par recn, float *_dev_collection, size_t _collection_p)
{
	dim3 threadsPerBlock(8, 8, 8); ///< zvolim velikost bloku 8 x 8 x 8 = 512 vlaken
	dim3 numBlocks(ceil((float)recn.rec_X / threadsPerBlock.x), ceil((float)recn.rec_Y / threadsPerBlock.y), ceil((float)recn.ZBlockSize / threadsPerBlock.z));  ///< spocitam potrebny pocet bloku na vypocet cele velikosti rekonstrukce
	CollectionKernel <<<numBlocks, threadsPerBlock >>>(recn, _dev_collection, _collection_p);

	cudaError_t cudaStatus = cudaDeviceSynchronize();
	if (cudaStatus != cudaSuccess)
	{
		std::cerr << "reconstruction: cudaDeviceSynchronize returned " << cudaGetErrorString(cudaStatus) << " after launching CollectingKernel!" << std::endl;
	}

	return cudaStatus;
}
