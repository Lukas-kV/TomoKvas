#include "rsq_load.h"

//#define PI (2 * asin(1.))
#define PI 3.14159265358979323846

rsq::rsq(string _filename)
{
	fstream = NULL;
	idata = NULL;
	filename = _filename;

	header_offset = 3584;
	file_size = 0;
	sin_x = 1024;  
	sin_y = 2000; 
	sin_y_fc = sin_y + 2; 
	sin_x_cpl = 0;
	sin_y_cpl = 0;

	detector_height = 248;

	axis_move = 21;
	distance = 125810. / 24.; ///< default distance for scanco uCT
	delta_fi = 2 * PI / sin_y; ///< default value for fullscan scanco uCT
	flip = false;
	invert = false;
	complete = false;
	minimal = false;
	parallel = false;

	bHeader = false;

	num_of_sin = 0; 
	num_of_readed_sin = 0;

	zsize = 0;
}

rsq::~rsq()
{
	if(idata)
		delete[] idata;
}

void rsq::SetOffset(unsigned long int offset)
{
	header_offset = offset;
}

void rsq::GetSize(unsigned int &_sin_x, unsigned int &_sin_y)
{
	_sin_x = sin_x_cpl;
	_sin_y = sin_y_cpl;
}

unsigned int rsq::GetNum()
{
	return num_of_sin;
}

unsigned int rsq::GetNumLoaded()
{
	return num_of_readed_sin;
}

void rsq::CalculateFileSize()
{
	
	_fseeki64( fstream, 0, SEEK_END ); 
	file_size = _ftelli64(fstream); 
	rewind(fstream); 

	//file_size = iHeader.nr_of_bytes;

}

void rsq::CalculateNumOfSinograms()
{
	num_of_sin = (file_size - header_offset) / (sin_x * sin_y_fc * 2);
}

void rsq::CalculateSizeOfSinogram()
{
	if(complete)
	{
		if(parallel)
		{
			sin_x_cpl = 2 * sin_x - 1 - 2 * axis_move;
			sin_y_cpl = sin_y / 2;	
			return;
		}

		if(minimal)
		{
			sin_x_cpl = 2 * sin_x - 1 - 2 * axis_move;
			sin_y_cpl = ceil((PI + 2*atan(sin_x_cpl / 2. / distance)) / delta_fi);
			return;
		}

		sin_x_cpl = 2 * sin_x - 1 - 2 * axis_move;	
		sin_y_cpl = sin_y;
		return;
	}

	sin_x_cpl = sin_x - axis_move;
	sin_y_cpl = sin_y;
	return;
}

void rsq::SetAxisMove(unsigned int _axis)
{
	axis_move = _axis;
}

void rsq::SetFlip(bool _flip)
{
	flip = _flip;
}

void rsq::SetParallel(bool _para)
{
	parallel = _para;
}

void rsq::SetMinimal(bool _min)
{
	minimal = _min;
}

void rsq::SetInvert(bool _invert)
{
	invert = _invert;
}

void rsq::SetDistance(double _distance)
{
	distance = _distance;
}

void rsq::SetDeltaFi(double _delta_fi)
{
	delta_fi = _delta_fi;
}

void rsq::SetDetectorHeight(unsigned int _height)
{
	detector_height = _height;
}

int	 rsq::GetDetectorHeight()
{
	return detector_height;
}

void rsq::SetCompleting(bool _complete)
{
	complete = _complete;
}

int rsq::LoadImages(unsigned int _NumOfWholeDetector)
{
	unsigned int first = _NumOfWholeDetector * detector_height;
	unsigned int last = detector_height - 1 + _NumOfWholeDetector * detector_height;
	return LoadImages(first, last);
}

int rsq::LoadImagesUpper(unsigned int _NumOfWholeDetector)
{
	unsigned int first = (detector_height - 1) / 2  + _NumOfWholeDetector * detector_height;
	unsigned int last = detector_height - 1 + _NumOfWholeDetector * detector_height;
	return LoadImages(first, last);	
}

int rsq::LoadImagesLower(unsigned int _NumOfWholeDetector)
{
	unsigned int first = _NumOfWholeDetector * detector_height;
	unsigned int last = (detector_height - 1) / 2  + _NumOfWholeDetector * detector_height;
	return LoadImages(first, last);
}

int rsq::LoadImages(unsigned int first, unsigned int last)
{
	int error = 0;

	if( error = fopen_s(&fstream, filename.c_str(), "rb")) 
		return error;

	LoadHeader();
	detector_height = iHeader.nr_of_projections;
	zsize = iHeader.dimz_p;

	CalculateFileSize();
	CalculateNumOfSinograms();
	CalculateSizeOfSinogram();

	if(last > num_of_sin || last < first)
		return - 1;

	if(idata)
	{
		delete[] idata;
		idata = NULL;		
	}

	num_of_readed_sin = last + 1 - first;
	unsigned int size_of_sin = sin_x * sin_y_fc;
	unsigned int size_of_block =  size_of_sin * num_of_readed_sin;

	idata = new unsigned short int[size_of_block];

	unsigned long long move =  header_offset + (unsigned long long)size_of_sin * 2 * first;
	_fseeki64(fstream, move , 0);

	if( fread(idata, 2, size_of_block, fstream) !=  size_of_block)
	{
		fclose(fstream);
		fstream = NULL;
		num_of_readed_sin = 0;
		return -1;
	}

/*	int max = 0,min = 0;
	for(int a = 0; a < size_of_block; a++)
	{
		if(idata[a] > max) max = idata[a];
		if(idata[a] < min) min = idata[a];
	}
	cout << max << "  " << min << endl;
*/
	fclose(fstream);
	fstream = NULL;
	return error;
}

int  rsq::LoadImages()
{
	int error;

	if( error = fopen_s(&fstream, filename.c_str(), "rb")) 
		return error;

	CalculateFileSize();
	CalculateNumOfSinograms();

	if(idata)
	{
		delete[] idata;
		idata = NULL;
	}

	unsigned int size_of_block = sin_x * sin_y_fc * num_of_sin;

	idata = new unsigned short int[size_of_block];

	_fseeki64(fstream, header_offset, 0);

	if( fread(idata, 2, size_of_block, fstream) !=  size_of_block)
	{
		fclose(fstream);
		return -1;
	}

	fclose(fstream);
	return error;

}

int rsq::ProccessAndCopy(double *_post_processed_data)
{
	double *ffdata = new double[sin_x * sin_y];
	unsigned short int *iblock = NULL;
	double *dblock = NULL;

	double p = num_of_readed_sin;
	p /= 100;

	for(unsigned int n = 0; n < num_of_readed_sin; n++)
	{
		iblock = idata + n * sin_y_fc * sin_x;
		dblock = _post_processed_data + n * sin_y_cpl * sin_x_cpl;

		FlatFieldCorrection(iblock, ffdata);
		if(complete)
			if(parallel)
				Half2Para(ffdata, dblock);
			else			
				Half2FullFan(ffdata, dblock);
			
		else
			CutHalfSin(ffdata, dblock);

		printf("preparing sinograms %.2f %%      \r", (n + 1) / p);
	}

	delete[] ffdata;
	return 0;
}

void rsq::FlatFieldCorrection(unsigned short *image, double *result)
{
	
	for(int y = 2; y < sin_y_fc; y++)
	{
		double *drow = result + (y - 2) * sin_x; 
		unsigned short int *irow = image + y * sin_x;
		double value = 0;

		if(invert)
			for(int x = 0; x < sin_x; x++)
			{
				value = (irow[x] - image[x]) / (double)(image[sin_x + x] - image[x]);
				value = BeamH.A * log(value) + BeamH.B;
				drow[x] = -value;
			}
		else
			for(int x = 0; x < sin_x; x++)
			{
				value = (irow[x] - image[x]) / (double)(image[sin_x + x] - image[x]);
				value = BeamH.A * log(value) + BeamH.B;
				drow[x] = value;
			}
	}

	return;
}

void rsq::Half2FullFan(double *image, double *result)
{
	double *atan_vect = new double[sin_x - axis_move];
	for(int t = -sin_x + 1 + axis_move; t < 0; t++)
	{
		atan_vect[-t] = atan(t / distance);
	}

	int start = 0;
	int stop = sin_y_cpl;
	if (flip)
	{
		start = sin_y - sin_y_cpl;
		stop = sin_y;
	}


	for(int f = start; f < stop; f++)
	{
		double *row_p;
		if (flip)
			row_p = result + (sin_y_cpl - 1 - f + start) * sin_x_cpl;
		else
			row_p = result + f * sin_x_cpl;

		double fi_rad =  f * delta_fi;

		memcpy(row_p + sin_x - 1 - axis_move, image + f * sin_x + axis_move, (sin_x - axis_move) * sizeof(image[0]));

		for(int t = -sin_x + 1 + axis_move; t < 0; t++)
		{
			double at = atan_vect[-t];

			double new_f = (fi_rad - 2 * at - PI) / delta_fi;

			if(new_f > (sin_y - 1) || new_f < 0)
				new_f = (fi_rad - 2 * at + PI) / delta_fi;

			int a = ceil(new_f);
			int b = floor(new_f);

			double A;

			if(a > sin_y - 1)
			{					
				A = image[(a - sin_y)* sin_x + (-t) + axis_move];
			}
			else
				A = image[a * sin_x + (-t) + axis_move];

			double B = image[b * sin_x + (-t) + axis_move];


			if(a != b)
				row_p[t + sin_x - 1 - axis_move] = (A + (B - A) * (-new_f + a));	///interpolace protoze b - a = -1 stale
			else
				row_p[t + sin_x - 1 - axis_move] = A;

		}
	}

	delete[] atan_vect;
}

void rsq::CutHalfSin(double *image, double *result)
{
	if (flip)
		for(int y = 0; y < sin_y; y++)
		{
			double *irow = image + y * sin_x + axis_move;
			double *rrow = result + (sin_y - 1 - y)* sin_x_cpl;
			memcpy(rrow, irow, sin_x_cpl * sizeof(rrow[0]));
		}
	else
		for(int y = 0; y < sin_y; y++)
		{
			double *irow = image + y * sin_x + axis_move;
			double *rrow = result + y * sin_x_cpl;
			memcpy(rrow, irow, sin_x_cpl * sizeof(rrow[0]));
		}
}

void rsq::Half2Para(double *image, double *result)
{
	for(int f = 0; f < sin_y_cpl; f++)
	{
		double *rrow = result + f * sin_x_cpl;
		for(int u = 0; u < sin_x_cpl; u++)
		{	
			rrow[u] = Fan2Para(image, f, u);
		}
	}
		
}

double rsq::Fan2Para(double *image, unsigned int _fi_index, unsigned int _u_index)
{
	unsigned int center_index = sin_x - 1 - axis_move;
	double u_index = ((double)_u_index - (double)center_index);
	double beta_index = 0;
	
	if(u_index >= 0.)
		beta_index = _fi_index + asin(u_index / distance)/ delta_fi;
	else
		beta_index = _fi_index - (asin(u_index / distance) - PI)/ delta_fi;

	double s_index = u_index; // * distance * distance  /( R * sqrt(distance * distance - u_index * u_index));

	if(beta_index < .0)
		beta_index += sin_y;

	if(beta_index > sin_y - 1)
		beta_index -= sin_y;

	if(u_index < 0.)
		s_index *= -1;

	int ba = ceil(beta_index);
	int bb = floor(beta_index);

	int sa = ceil(s_index);
	int sb = floor(s_index);
	int msa = sa + axis_move;
	int msb = sa + axis_move;

	if(ba != bb)
	{
		double va,vb;

		//* interpolace mezi radky
		double* rowa = image + ba * sin_x ; /// vypocteme  ukazatel na radek sinogramu s uhlem beta
		double* rowb = NULL;

		if(bb >= 0)
			rowb = image + bb * sin_x ; 
		else
			rowb = image + (bb + sin_y) * sin_x ; 

		if(sa != sb)
		{
			va = (rowa[msa] + (rowa[msb] - rowa[msa]) * (- s_index + sa));	
			vb = (rowb[msa] + (rowb[msb] - rowb[msa]) * (- s_index + sa));	/// interpolace b - a = -1
		}
		else
		{
			va = rowa[msa];
			vb = rowb[msb];
		}

		return (va + (vb - va) * (- beta_index + ba));		
	}
	else
	{
		double* row = image + ba * sin_x;

		if(sa != sb)
			return (row[msa] + (row[msb] - row[msa]) * (s_index - sa));	

		return row[msa];
	}
}

void rsq::LoadHeader()
{
	if(bHeader) 
		return;

	if(!fstream)
		return;

	cout << "size of header struct: " << sizeof(iHeader) << endl;

	fread(&iHeader, 1, sizeof(iHeader), fstream);

	char tmp[17];
	memcpy(tmp,iHeader.check, 16 * sizeof(tmp[0]));
	tmp[16] = 0;

	char bh[152];
	fseek(fstream, 3432, 0);
	fread(bh, 1, 152, fstream);

	cout << "string: " << tmp << endl;
	cout << "data offset: " << 512 + iHeader.data_offset * 512 << " B" << endl;
	cout << "filesize: " << iHeader.nr_of_bytes << " B" << endl;
	cout << "filesize blocks 512B: " << iHeader.nr_of_blocks << endl;
	cout << "patient index: " << iHeader.patient_index << endl;
	cout << "scaner ID: " << iHeader.scanner_id << endl;
	cout << "scaner type: " << iHeader.scanner_type << endl;
	
	cout << "creation date: " << (time_t)iHeader.creation_date << endl;

	cout << "dim x: " << iHeader.dimx_p << " px" << endl;
	cout << "dim y: " << iHeader.dimy_p << " px" << endl;
	cout << "dim z: " << iHeader.dimz_p  << " px" << endl;

	cout << "dim x: " << iHeader.dimx_um << " um" << endl;
	cout << "dim y: " << iHeader.dimy_um << " um" << endl;
	cout << "dim z: " << iHeader.dimz_um << " um" << endl;

	cout << "number of samples: " << iHeader.nr_of_samples << endl;
	cout << "number of projections: " << iHeader.nr_of_projections << endl;
	cout << "scan distance: " << iHeader.scandist_um << " um"<< endl;

	cout << "reference line: " << iHeader.reference_line_um << " um"<< endl;

	//cout << "name: " << iHeader.name << endl;

	cout << "beam hardening correction: " << bh << endl;

	sscanf(bh, "Corr at %d %d %d:  %lf %lf %lf", &BeamH.X[0], &BeamH.X[1], &BeamH.X[2], &BeamH.Y[0], &BeamH.Y[1], &BeamH.Y[2]);

	cout << "beam hardening correction points (X,Y):" << endl;
	for(int i = 0; i < 3; i++)
		cout << "	" << BeamH.X[i] << "	" << BeamH.Y[i] << endl;

	BeamH.A = (BeamH.Y[2] - BeamH.Y[0]) / (BeamH.X[2] - BeamH.X[0]);
	BeamH.B = BeamH.Y[0] - BeamH.A * BeamH.X[0];

	cout << endl << "beam hardening correction curve:" << endl;
	cout << "	A = " << BeamH.A << endl;
	cout << "	B = " << BeamH.B << endl;

	bHeader = true;
}

int	 rsq::GetNoOfWholeBlocks()
{
	return zsize / detector_height;
}