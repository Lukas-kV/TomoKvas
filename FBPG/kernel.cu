#include "kernel.cuh"

string GetDeltaTime( clock_t _start_time, clock_t _now)
{
	float delta = _now - _start_time;	
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


 //**  class implementation **//
FilteredBackProjection::FilteredBackProjection()
{
	for(int i = 0; i < POINTERS; i++)
		sing.dev_sinogram[i] = NULL;
	sing.dev_dcos = NULL;
	sing.dev_dsin = NULL;
	dev_block = NULL;
	recn.dev_reconstruction = NULL;
	dev_collecting = NULL;
	collecting = false;
	show = false;

	sing.scaling = 1;

	recn.dirX = 0;
	recn.dirY = 0;
	recn.dirZ = 0;

	recn.scale = 1;

	cudaStatus = cudaSuccess;

	devID = 0;

	size_of_Z_block = 16;
	memory_limit = (long long int)1024 * 1024 * 1024; ///< 1 GiB

	memory_sinograms = 0;
	nr_of_sin_blocks = 1;
	sinograms_per_limit = 1;

	pr = 0;

	output_file = NULL;
}

FilteredBackProjection::~FilteredBackProjection()
{

	for(int i = 0; i < POINTERS; i++)
		if(sing.dev_sinogram[i]) 
			cudaFree(sing.dev_sinogram[i]);

	if(sing.dev_dsin) 
		cudaFree(sing.dev_dsin);

	if(sing.dev_dcos) 
		cudaFree(sing.dev_dcos);

	if(recn.dev_reconstruction)
		cudaFree(recn.dev_reconstruction);

	if(dev_collecting)
		cudaFree(dev_collecting);

    cudaStatus = cudaDeviceReset();
    if (cudaStatus != cudaSuccess) 
	{
        cerr << "cudaDeviceReset failed!" << endl;
    }

	cout <<  "reconstruction complete: " << GetDeltaTime(global_start, clock()) << endl;
}

void FilteredBackProjection::SetGPUid(const int &_ID)
{
	devID = _ID;
}

void FilteredBackProjection::SetBlockSize(const int &_sz)
{
	size_of_Z_block = _sz;
}

void FilteredBackProjection::SetCollecting(const bool &_collecting)
{
	collecting = _collecting;
}

void FilteredBackProjection::SetShowInfo(const bool &_show)
{
	show = _show;
}

void FilteredBackProjection::SetAllocationLimit(const long long int &_memory_limit)
{
	memory_limit = _memory_limit;
}

void FilteredBackProjection::SetScaleXYZ(const float &_scale)
{
	recn.scale = _scale;
}

size_t FilteredBackProjection::GetFreeMemInfo(int _GPUId)
{
	cudaDeviceProp deviceProps; 
	size_t free, total;

	cudaGetDeviceProperties(&deviceProps, _GPUId);
	cudaMemGetInfo(&free, &total); /// velikost pameti volna / celkova	

	return free;
}

cudaError_t FilteredBackProjection::Initialize()
{

	global_start = start = clock(); 

	cudaGetDeviceProperties(&deviceProps, devID);
	cudaMemGetInfo(&free, &total); /// velikost pameti volna / celkova

	if(!collecting && !show) 
		fShowInfo();

	cudaStatus = cudaSetDevice(devID);
    if (cudaStatus != cudaSuccess) 
	{
        cerr << "cudaSetDevice failed!  Do you have a CUDA-capable GPU installed?" << endl;
		return cudaStatus;
    }

	sing.sin_in_block = sinograms_per_limit = memory_limit / memory_sin; /// celociselne vypoctu pocet sinogramu do memory limitu
	memory_limit = sinograms_per_limit * memory_sin; /// zpetne vypoctu limit tak aby se do nej vesel celocisleny pocet sinogramu

	nr_of_sin_blocks = ceil((float)memory_sinograms / (float)memory_limit); /// vypocte nejmensi vyssi pocet potrebnych bloku
	memory_rest = memory_limit - (memory_limit * nr_of_sin_blocks - memory_sinograms); /// vypocte se zbyvajici pamet
	sinograms_per_rest = memory_rest / memory_sin; /// vypocte se pocet sinogramu ve zbyvajici pameti

	unsigned long long int cnt = nr_of_sin_blocks;
	if(memory_rest)
	{
		cnt--;
		cudaStatus = cudaMallocPitch( (void**)&sing.dev_sinogram[cnt], &sing.sin_p, sing.sin_X * sizeof(sing.dev_sinogram[0][0]), sing.sin_Y * sinograms_per_rest );
		if (cudaStatus != cudaSuccess) 
		{	
			cerr << "requested memsize: " << memory_sin * sinograms_per_rest / 1024. / 1024. << " MiB" << endl;
			cerr << "sinogram cudaMallocPitch failed!" << endl;
			return cudaStatus;
		}
	}

	for(int i = 0; i < cnt; i++)
	{
		cudaStatus = cudaMallocPitch( (void**)&sing.dev_sinogram[i], &sing.sin_p, sing.sin_X * sizeof(sing.dev_sinogram[0][0]), sing.sin_Y * sinograms_per_limit);
		if (cudaStatus != cudaSuccess) 
		{	
			cerr << "requested memsize: " << memory_sin * sinograms_per_limit / 1024. / 1024. << " MiB" << endl;
			cerr << "sinogram cudaMallocPitch failed!" << endl;
			return cudaStatus;
		}
	}

	cout << "copying " << memory_sinograms / 1024. / 1024. << " MiB of sinograms to GPU memory in " << nr_of_sin_blocks << " memory block(s)" << endl;
	for(unsigned long long int i = 0; i < cnt; i++)
	{
		float *data_set = sinogram + i * sinograms_per_limit * sing.sin_X * sing.sin_Y;
		cudaStatus = cudaMemcpy2D(sing.dev_sinogram[i], sing.sin_p, data_set, sing.sin_X * sizeof(sinogram[0]), sing.sin_X * sizeof(sinogram[0]), sing.sin_Y * sinograms_per_limit, cudaMemcpyHostToDevice);
		if (cudaStatus != cudaSuccess) 
		{
			cerr << "sinogram cudaMemcpy2D failed!" << endl;
			cudaFree(sing.dev_sinogram[i]);
			return cudaStatus;
		}
	}

	if(memory_rest)
	{
		float *data_set = sinogram + cnt * sinograms_per_limit  * sing.sin_X * sing.sin_Y;
		cudaStatus = cudaMemcpy2D(sing.dev_sinogram[cnt], sing.sin_p, data_set, sing.sin_X * sizeof(sinogram[0]), sing.sin_X * sizeof(sinogram[0]), sing.sin_Y * sinograms_per_rest, cudaMemcpyHostToDevice);
		if (cudaStatus != cudaSuccess) 
		{
			cerr << "sinogram cudaMemcpy2D failed!" << endl;
			cudaFree(sing.dev_sinogram[cnt]);
			return cudaStatus;
		}	
	}
	
	float nr_Z = (float)(recn.rec_Z - 1)/ (float)size_of_Z_block;
	nr_of_Z_blocks = ceil(nr_Z);
	dcorZ = (nr_of_Z_blocks - nr_Z) * size_of_Z_block / 2;

	cudaStatus = cudaMallocPitch( (void**)&recn.dev_reconstruction, &recn.rec_p, recn.rec_X * sizeof(recn.dev_reconstruction[0]), recn.rec_Y * size_of_Z_block);
    if (cudaStatus != cudaSuccess) 
	{
        cerr << "reconstruction cudaMallocPitch failed!" << endl;
		return cudaStatus; 
    }

	if(collecting)
	{
		cudaStatus = cudaMallocPitch( (void**)&dev_collecting, &collecting_p, recn.rec_X * sizeof(recn.dev_reconstruction[0]), recn.rec_Y * size_of_Z_block);
		if (cudaStatus != cudaSuccess) 
		{
	        cerr << "collecting cudaMallocPitch failed!" << endl;
			return cudaStatus; 
		}
	}

	cout <<  "initialization: " << GetDeltaTime(start, clock()) << endl;

	CalcSinCos();

	cudaMemGetInfo(&free, &total);
	cout << "GPU reconstruction block: " << (unsigned long long)size_of_Z_block * recn.rec_X * recn.rec_Y * sizeof(float) / 1024 / 1024 << "MiB" << endl;
	cout <<"Free GPU RAM: " << free / 1024 / 1024 << "MiB" << endl;

	return cudaStatus;
}

cudaError_t FilteredBackProjection::CalcSinCos()
{
	float *dsin = new float[sing.sin_Y];
	float *dcos = new float[sing.sin_Y];
	double beta_rad;

	for( int beta = 0; beta < sing.sin_Y; beta++)
	{
		beta_rad = beta * sing.delta_f + recn.rot;
		dsin[beta] = sin(beta_rad);
		dcos[beta] = cos(beta_rad);
	}

	cudaStatus = cudaMalloc( (void**)&sing.dev_dsin, sing.sin_Y * sizeof(sing.dev_dsin[0]));
    if (cudaStatus != cudaSuccess) 
	{
        cerr << "sin array cudaMalloc failed!" << endl;
		return cudaStatus;
    }

	cudaStatus = cudaMemcpy(sing.dev_dsin, dsin, sing.sin_Y * sizeof(dsin[0]), cudaMemcpyHostToDevice);
    if (cudaStatus != cudaSuccess) 
	{
        cerr << "sin array cudaMemcpy failed!" << endl;
		cudaFree(sing.dev_dsin);
		return cudaStatus;
    }

	cudaStatus = cudaMalloc( (void**)&sing.dev_dcos, sing.sin_Y * sizeof(sing.dev_dcos[0]));
    if (cudaStatus != cudaSuccess) 
	{
        cerr << "cos array cudaMalloc failed!" << endl;
		return cudaStatus;
    }

	cudaStatus = cudaMemcpy(sing.dev_dcos, dcos, sing.sin_Y * sizeof(dcos[0]), cudaMemcpyHostToDevice);
    if (cudaStatus != cudaSuccess) 
	{
        cerr << "cos array cudaMemcpy failed!" << endl;
		cudaFree(sing.dev_dcos);
		return cudaStatus;
    }

	delete[] dsin;
	delete[] dcos;

	return cudaStatus;
}

void FilteredBackProjection::SetSinogramProperties( float *_sinogram, const int &_sinogram_X, const float &_centre_X, const int &_sinogram_Y, const int &_sinogram_Z, const float &_centre_Z, const float &_delta_fi, float _distance, float _sampling)
{
	sinogram = _sinogram;

	sing.sin_X = _sinogram_X;
	sing.sin_Y = _sinogram_Y;
	sing.sin_Z = _sinogram_Z;
	sing.delta_f = _delta_fi;
	sing.distance = _distance;	
	sing.cX = _centre_X;
	sing.cZ = _centre_Z;
	sing.sampling = _sampling;

	memory_sin = (unsigned long long)sing.sin_X * sizeof(_sinogram[0]) * sing.sin_Y;
	memory_sinograms =  memory_sin * sing.sin_Z;



}

void FilteredBackProjection::SetAxisDirection(const bool &_X, const bool &_Y, const bool &_Z)
{
	if(_X)
		recn.dirX = (unsigned int)1 << 31;
	if(_Y)
		recn.dirY = (unsigned int)1 << 31; /// < maska pro zmenu znamenka osy pomoci XOR osa je v float
	if(_Z)
		recn.dirZ = (unsigned int)1 << 31;
}

/*void FilteredBackProjection::SetReconstructionProperties(float *_reconstruction, const int &_reconstruction_X, const float &_centre_X, const int &_reconstruction_Y, const float &_centre_Y, const int &_reconstruction_Z, const float &_centre_Z, const float &_rot)
{
	reconstruction = _reconstruction;

	recn.rec_X = _reconstruction_X;
	recn.rec_Y = _reconstruction_Y;
	recn.rec_Z = _reconstruction_Z;

	recn.dX = -_centre_X;
	recn.dY = -_centre_Y;
	recn.dZ = dZ = -_centre_Z;
	recn.rot = _rot;
}
*/

void FilteredBackProjection::SetReconstructionProperties(bin_siz *_output_file, const int &_reconstruction_X, const float &_centre_X, const int &_reconstruction_Y, const float &_centre_Y, const int &_reconstruction_Z, const float &_centre_Z, const float &_rot)
{

	output_file = _output_file;

	recn.rec_X = _reconstruction_X;
	recn.rec_Y = _reconstruction_Y;
	recn.rec_Z = _reconstruction_Z;

	recn.dX = -_centre_X;
	recn.dY = -_centre_Y;
	recn.dZ = dZ = -_centre_Z;
	recn.rot = _rot;
}

cudaError_t FilteredBackProjection::fRBcomma()
{
	start = clock();

	dim3 threadsPerBlock(16,16);
	dim3 numBlocks( ceil((float)sing.sin_Y / threadsPerBlock.x), ceil((float)sing.sin_Z / threadsPerBlock.y));
	RB_comma3DKernel<<<numBlocks, threadsPerBlock>>>(sing);

	cudaStatus = cudaDeviceSynchronize();
	if (cudaStatus != cudaSuccess) 
	{
		cerr << "RB': cudaDeviceSynchronize returned error code %d after launching addKernel!" << endl;
	}
	
	cout <<  "calculation of RB': " << GetDeltaTime(start, clock()) << endl;
	return cudaStatus;
}

cudaError_t FilteredBackProjection::fConvolution()
{

	start = clock();

	cufftResult cufftStatus;

	cufftHandle plan = NULL;
	cufftComplex *dev_data = NULL;
	cufftComplex *dev_filter = NULL;

	float *dev_abs_filter = NULL;

	int p = 2;
	while(p < sing.sin_X)
		p *= 2;	
	p *= 2;
	
	sing.scaling /= p;

	cout << "sinogram size x: " << sing.sin_X << endl << "zero padding to: " << p << endl;

	size_t comp_size = p / 2 + 1;

	float *filter = new float[comp_size * 2];
	for(int x = 0; x < comp_size * 2; x++)
		filter[x] = 0;

	FillFilter(filter, sing.sin_X); // zero padded filter

  /*filter[0] = filter[1] = 0;
	for(int i = 1; i < comp_size; i++)
	{
		float t = (float)i / (float)comp_size;
		filter[i * 2] = filter[i * 2 + 1] = t;
	}

	FILE *f = fopen("test.raw", "wb");
	fwrite(filter, sizeof(filter[0]), comp_size * 2, f);
	fclose(f);*/

	cudaStatus = cudaMalloc((void**)&dev_filter, sizeof(dev_filter[0]) * comp_size);
    if (cudaStatus != cudaSuccess) 
	{
        cerr << "convolution: filter malloc failed!" << endl;
		return cudaStatus; 
    }
	
	cudaStatus = cudaMalloc((void**)&dev_abs_filter, sizeof(dev_abs_filter[0]) * comp_size);
    if (cudaStatus != cudaSuccess) 
	{
        cerr << "convolution: filter abs malloc failed!" << endl;
		return cudaStatus; 
    }

	cudaStatus = cudaMemcpy(dev_filter, filter, p * sizeof(filter[0]), cudaMemcpyHostToDevice);
	if (cudaStatus != cudaSuccess) 
	{
		cerr << "convolution: filter cudaMemcpy2D failed! " << endl;
		cudaFree(dev_filter);
		return cudaStatus;
	}

	cufftStatus = cufftPlan1d(&plan, p, CUFFT_R2C, 1);
	if (cufftStatus != CUFFT_SUCCESS) 
	{
		cerr << "convolution: Filter FFT plan failed! " << endl;
		return (cudaError_t)cufftStatus;
	}
	
	cufftStatus = cufftSetCompatibilityMode(plan, CUFFT_COMPATIBILITY_NATIVE);
	if (cufftStatus != CUFFT_SUCCESS) 
	{
		cerr << "convolution: FFT set native failed! " << endl;
		return (cudaError_t)cufftStatus;
	}

	cufftStatus = cufftExecR2C(plan, (cufftReal*)dev_filter, dev_filter);
	if (cufftStatus != CUFFT_SUCCESS) 
	{
        cerr << "convolution: Filter FFT execute R2C failed! " << endl;
		cufftDestroy(plan);
		return (cudaError_t)cufftStatus;
	}

	cufftDestroy(plan);

	dim3 threadsPerBlock(16);
	dim3 numBlocks( ceil((float)comp_size / threadsPerBlock.x));
	AbsKernel<<<numBlocks, threadsPerBlock>>>( comp_size, dev_filter, dev_abs_filter, sing.sampling);

	cudaStatus = cudaDeviceSynchronize();
	if (cudaStatus != cudaSuccess) 
	{
		cerr << "convolution: filter calculation cudaDeviceSynchronize returned error code " << cudaStatus << "after launching addKernel! "  << endl;
		return cudaStatus; 
	}

	cudaFree(dev_filter);

	//* FFT *//
	cudaStatus = cudaMalloc((void**)&dev_data, sizeof(dev_data[0]) * comp_size * sing.sin_Y);
    if (cudaStatus != cudaSuccess) 
	{
        cerr << "convolution: FFT malloc failed!" << endl;
		return cudaStatus; 
    }

	for(int n = 0; n < sing.sin_Z; n++)
	{

		int block = n / sinograms_per_limit;
		int ir = n - block * sinograms_per_limit;


			/* Create a 1D FFT plan. */
		cufftStatus = cufftPlan1d(&plan, p, CUFFT_R2C, sing.sin_Y);
		if (cufftStatus != CUFFT_SUCCESS) 
		{
			cerr << "convolution: FFT plan failed! " << n << endl;
			cudaFree(dev_data);
			return (cudaError_t)cufftStatus;
		}
	
		cufftStatus = cufftSetCompatibilityMode(plan, CUFFT_COMPATIBILITY_NATIVE);
		if (cufftStatus != CUFFT_SUCCESS) 
		{
			cerr << "convolution: FFT set native failed! " << n << endl;
			cudaFree(dev_data);
			return (cudaError_t)cufftStatus;
		}

		dev_block = (float*)((char*)sing.dev_sinogram[block] + ir * ( sing.sin_p * sing.sin_Y ));

		cudaStatus = cudaMemcpy2D(dev_data, comp_size * sizeof(dev_data[0]), dev_block, sing.sin_p, sing.sin_X * sizeof(dev_block[0]), sing.sin_Y, cudaMemcpyDeviceToDevice);
		if (cudaStatus != cudaSuccess) 
		{
			cerr << "convolution: sinogram to block cudaMemcpy2D failed! "  << n << endl;
			cufftDestroy(plan);
			cudaFree(dev_data);
			return cudaStatus;
		}

		cudaStatus = cudaMemset2D((float *)dev_data + sing.sin_X, comp_size * sizeof(dev_data[0]), 0, sizeof(float) * (2 * comp_size - sing.sin_X), sing.sin_Y);
		if (cudaStatus != cudaSuccess) 
		{
			cerr << "convolution: zeroing pitch cudaMemset2D failed! "  << n << endl;
			cufftDestroy(plan);
			cudaFree(dev_data);
			return cudaStatus;
		}
		
		/* Use the CUFFT plan to transform the signal in place. */
		cufftStatus = cufftExecR2C(plan, (cufftReal*)dev_data, dev_data);
		if (cufftStatus != CUFFT_SUCCESS) 
		{
	        cerr << "convolution: FFT execute R2C failed! " << n << endl;
			cufftDestroy(plan);
			cudaFree(dev_data);
			return (cudaError_t)cufftStatus;
		}

		cufftDestroy(plan);
	
		//* filtration *//
		dim3 threadsPerBlock(16);
		dim3 numBlocks( ceil((float)sing.sin_Y / threadsPerBlock.x));
		FiltrationKernel<<<numBlocks, threadsPerBlock>>>( comp_size, sing.sin_Y, dev_data, dev_abs_filter);

		cudaStatus = cudaDeviceSynchronize();
		if (cudaStatus != cudaSuccess) 
		{
			cerr << "convolution: filtration cudaDeviceSynchronize returned error code " << cudaStatus << "after launching addKernel! "  << n << endl;
			cudaFree(dev_data);
			cufftDestroy(plan);
			return cudaStatus; 
		}

		cufftStatus = cufftPlan1d(&plan, p, CUFFT_C2R, sing.sin_Y);
		if (cufftStatus != CUFFT_SUCCESS) 
		{
			cerr << "convolution: FFT plan failed! "  << n << endl;
			cudaFree(dev_data);
			return (cudaError_t)cufftStatus;
		} 

		cufftStatus = cufftSetCompatibilityMode(plan, CUFFT_COMPATIBILITY_NATIVE);
		if (cufftStatus != CUFFT_SUCCESS) 
		{
			cerr << "convolution: FFT set native failed! " << n << endl;
			cudaFree(dev_data);
			return (cudaError_t)cufftStatus;
		}

		//* Backward FFT *//

		cufftStatus = cufftExecC2R(plan, dev_data, (cufftReal*)dev_data);
		if (cufftStatus != CUFFT_SUCCESS) 
		{
			cerr << "FFT-1 execute C2R failed! "  << n << endl;
			cudaFree(dev_data);
			cufftDestroy(plan);
			return (cudaError_t)cufftStatus;
		}

		cudaStatus = cudaMemcpy2D(dev_block, sing.sin_p, dev_data,  comp_size * sizeof(dev_data[0]), sing.sin_X * sizeof(dev_block[0]), sing.sin_Y, cudaMemcpyDeviceToDevice);
		if (cudaStatus != cudaSuccess) 
		{
			cerr << "convolution: block to sinogram cudaMemcpy2D failed! "  << n << endl;
			cudaFree(dev_data);
			cufftDestroy(plan);
			return cudaStatus;
		}

		cufftDestroy(plan); /* Destroy the CUFFT plan. */
	}

	cudaFree(dev_data); /* uvolneni komplexnich dat */
	cudaFree(dev_abs_filter);

	delete[] filter;

	cout <<  "convolution: " << GetDeltaTime(start, clock()) << endl;
	return  cudaStatus;
}

void FilteredBackProjection::FillFilter(float *_filter, int size)
{
	// naplneni 1/2 h(nt), KaK[3.6 (176)]
	for(int i = 0; i < size; i++)
	{
		_filter[i] = RamLak(i - size / 2, sing.sampling) / 2;
	}
}

float FilteredBackProjection::RamLak(int _value, float _sm)
{
	if(_value == 0)
		return 1./(4. * _sm * _sm);

	if((_value % 2) == 0)
		return 0.;

	return -1./(_value * _value * PI * PI * _sm * _sm);

}

cudaError_t FilteredBackProjection::CalculateReconstruction()
{

	fRBcomma();
	fConvolution();

	start = clock();
	
	float p = nr_of_Z_blocks / 100.;
	printf("calculating reconstruction %.2f %%      \r",0);

	recn.ZBlockSize = size_of_Z_block;

	//output_file->SetSlicesPerBlock(size_of_Z_block);

	for(int n = 0; n < nr_of_Z_blocks; n++)
	{
		recn.dZ = dZ + ( n - (float)nr_of_Z_blocks / 2. ) * size_of_Z_block + dcorZ;
		clock_t kt = clock();


		dim3 threadsPerBlock(8, 8, 8); ///< zvolim velikotak jst bloku 8 x 8 x 8 = 512 vlaken
		dim3 numBlocks(ceil((float)recn.rec_X / threadsPerBlock.x), ceil((float)recn.rec_Y / threadsPerBlock.y), ceil((float)size_of_Z_block / threadsPerBlock.z));  ///< spocitam potrebny pocet bloku na vypocet cele velikosti rekonstrukce
		FBPConeKernel<<<numBlocks, threadsPerBlock>>>(recn, sing); 
		
		// cudaDeviceSynchronize waits for the kernel to finish, and returns
		// any errors encountered during the launch.
		cudaStatus = cudaDeviceSynchronize();
		if (cudaStatus != cudaSuccess) 
		{
			 cerr << "reconstruction: cudaDeviceSynchronize returned error code " << cudaStatus << " after launching FBPKernel! in block " << n << endl;
			return cudaStatus; 
		}

		int SizeOfZBlock = size_of_Z_block;
		if((n + 1) * size_of_Z_block > recn.rec_Z)
			SizeOfZBlock = recn.rec_Z - n * size_of_Z_block;

		//float *block = reconstruction + n * recn.rec_X * recn.rec_Y * size_of_Z_block; /// nahradit objektem ktery bude vracet memory block, nejspise stejne misto v pameti, ale s jinym obsahem podle n
		float *block = NULL;
 		block = output_file->GetBlock( n, collecting);

		pr += clock() - kt;

	
		if(collecting)
		{			
			cudaStatus = cudaMemcpy2D(dev_collecting, collecting_p, block, recn.rec_X * sizeof(reconstruction[0]), recn.rec_X * sizeof(recn.dev_reconstruction[0]), recn.rec_Y * SizeOfZBlock, cudaMemcpyHostToDevice);
			if (cudaStatus != cudaSuccess) 
			{
				cerr << "reconstruction: collecting cudaMemcpy failed!" << endl;
				return cudaStatus; 
			}


			/*FILE *f = fopen("d:\\tt.raw","ab");
			fwrite(block,  recn.rec_X * recn.rec_Y, sizeof(float), f); 
			fclose(f);
			*/
			
			dim3 threadsPerBlock(8, 8, 8); ///< zvolim velikost bloku 8 x 8 x 8 = 512 vlaken
			dim3 numBlocks(ceil((float)recn.rec_X / threadsPerBlock.x), ceil((float)recn.rec_Y / threadsPerBlock.y), ceil((float)size_of_Z_block / threadsPerBlock.z));  ///< spocitam potrebny pocet bloku na vypocet cele velikosti rekonstrukce
			CollectionKernel<<<numBlocks, threadsPerBlock>>>(recn, dev_collecting, collecting_p); 
		
			// cudaDeviceSynchronize waits for the kernel to finish, and returns
			// any errors encountered during the launch.
			cudaStatus = cudaDeviceSynchronize();
			if (cudaStatus != cudaSuccess) 
			{
				cerr << "reconstruction: cudaDeviceSynchronize returned error code " << cudaStatus << " after launching CollectingKernel! in block " << n << endl;
				return cudaStatus; 
			}
		}
		//clock_t copy = clock();
		// Copy output matrix from GPU buffer to host memory.	  
		cudaStatus = cudaMemcpy2D(block, recn.rec_X * sizeof(reconstruction[0]), recn.dev_reconstruction, recn.rec_p, recn.rec_X * sizeof(recn.dev_reconstruction[0]), recn.rec_Y * SizeOfZBlock, cudaMemcpyDeviceToHost);
		if (cudaStatus != cudaSuccess) 
		{
			 cerr << "reconstruction: result cudaMemcpy failed!" << endl;
			return cudaStatus; 
		}

		//cout << "copy time: " << GetDeltaTime(copy, clock()) << endl;
		printf("calculating reconstruction %.2f %%      \r", (n + 1) / p);
	}
	//output_file->LastBlock();

	cout << endl;
	cout << "BP kernel:" << GetDeltaTime(-pr, 0) << endl <<  "reconstruction: " << GetDeltaTime(start, clock()) << endl;

	pr = 0;

	return  cudaStatus;	
}

void FilteredBackProjection::fShowInfo()
{
	cout << endl;
	cout << "GPU:" << endl;
	cout << "  CUDA device " << deviceProps.name << " has " << deviceProps.multiProcessorCount << " Multi-Processors" << endl;
	cout << "  Each Multi-Processor has " << deviceProps.maxThreadsPerMultiProcessor <<" threads" << endl;
	cout << "  Max threads per block " << deviceProps.maxThreadsPerBlock << endl;
	cout << "  Max threads per dim x " << deviceProps.maxThreadsDim[0] << endl;
	cout << "                      y " << deviceProps.maxThreadsDim[1] << endl;
	cout << "                      z " << deviceProps.maxThreadsDim[2] << endl;
	cout << "  Total     RAM         " << deviceProps.totalGlobalMem / 1024 / 1024 << "MiB" << endl << endl;
	cout << "  Free      RAM         " << free / 1024 / 1024 << "MiB" << endl << endl;
}


///d:\xmt\jp_tooth5b.cra
///d:\CONETEST\CONETEST.CRA
//d:\large_data\dm_scroll_03.cra