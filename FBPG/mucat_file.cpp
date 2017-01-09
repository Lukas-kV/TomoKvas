#include "mucat_file.h"

volatile bool sync_flag;


cra_con::cra_con(string _path, string _filename)
{
	filename = _filename;
	filename.append(".CRA");
	path = _path;
	fstream = nullptr;

	prs[0].assign("FILE NAME:              ");
	prs[1].assign("SOURCE NAME:            ");
	prs[2].assign("SAMPLE SPACING:         ");
	prs[3].assign("DISTANCE:               ");
	prs[4].assign("DETECTOR X PIXELS:      ");
	prs[5].assign("DETECTOR Z PIXELS:      ");
	prs[6].assign("NUMBER OF PROJECTIONS:  ");
	prs[7].assign("X DIMENSION:            ");
	prs[8].assign("Z DIMENSION:            ");
	prs[9].assign("ROTATE ANGLE:           ");
	prs[10].assign("CENTRE:                 ");
	prs[11].assign("NUMBER OF BLOCKS:       ");
	prs[12].assign("CLOCKWISE:");
	prs[13].assign("SCALE:                  ");
	
	header.sample_spacing = 0;
	header.distance = 0;
	header.detector_x = 0;
	header.detector_z = 0;
	header.num_of_projections = 0;
	header.x_dimension = 0;
	header.z_dimension = 0;
	header.rotate_angle = 0;
	header.centre_of_rotation[0] = 0;
	header.num_of_blocks = 1;
	header.clockwise = false;
	header.distance_unit = 0x63;
	header.spacing_unit = 0x75;
	header.scale = 1.;

	start = clock();

	ActiveBlock = 0;

	NoOfSinBlocks = 1;
	SinActiveBlock = 0;

	distance = 0;
	delta_fi = 0;

}

cra_con::~cra_con()
{
	//cout << "preparing sinograms: " << GetDeltaTime(start, clock()) << endl;
}

int cra_con::GetMetaData()
{
	int error = 0;

	string fn;
	fn.assign(path);
	fn.append("\\");
	fn.append(filename);

	ifstream file(fn);
	if(!file.is_open())
	{
		error = -1;
		return error;
	}
	string line;


	int centre_cnt = 0;

	while( getline( file, line ) )
	{
		for(int i = 0; i < STRINGS; i++)	
			if(!line.find(prs[i]))
			{
				line.erase(0,prs[i].length());
				istringstream stm(line);
				switch (i)
				{
					case 0:
						header.raw_file = line;
						break;
					case 1:
						header.source_name = line;
						break;
					case 2:
						stm >> header.sample_spacing;
						break;
					case 3:
						stm >> header.distance;
						break;
					case 4:
						stm >> header.detector_x;
						break;
					case 5:
						stm >> header.detector_z;
						break;
					case 6:
						stm >> header.num_of_projections;
						break;
					case 7:
						stm >> header.x_dimension;
						break;
					case 8:
						stm >> header.z_dimension;
						break;
					case 9:
						stm >> header.rotate_angle;
						break;
					case 10:
						stm >> header.centre_of_rotation[centre_cnt];
						centre_cnt++;
						break;
					case 11:
						stm >> header.num_of_blocks;
						break;			
					case 12:
						header.clockwise = true;
						break;
					case 13:
						stm >> header.scale;
						break;			

				}
			}

	}
	
	file.close();

	distance = (header.distance * 10000. / header.sample_spacing );
	delta_fi = 2 * PI / header.num_of_projections;

	unsigned int tmp;

	float phi = atan(header.detector_x / 2 / distance);
	tmp = 2 * ((int)(distance * sin(phi) + 1));
	
	if(header.x_dimension > tmp || header.x_dimension < 1)
		header.x_dimension = tmp;

	
	tmp = (distance - header.x_dimension / 2) / distance * header.detector_z;

	tmp &= ~1;

	if(header.z_dimension > tmp)
		header.z_dimension = tmp;

	if(header.z_dimension < 1)
		header.z_dimension = 1;


	return error;
}

void cra_con::PrintMetaData()
{
	cout << endl << "MetaData: " << endl;
	cout << "   file name:             " << header.raw_file << endl;
	cout << "   source name:           " << header.source_name << endl;
	cout << "   sample spacing:        " << header.sample_spacing << " " << header.spacing_unit << endl;
	cout << "   distance:              " << header.distance << " " << header.distance_unit << endl;
	cout << "   detector x pixels:     " << header.detector_x << endl;
	cout << "   detector y pixels:     " << header.detector_z << endl;
	cout << "   number of projections: " << header.num_of_projections << endl;
	cout << "   number of blocks:      " << header.num_of_blocks << endl;
	cout << "   x dimension:           " << header.x_dimension << endl;
	cout << "   z dimension:           " << header.z_dimension << endl;
	cout << "   rotate angle:          " << header.rotate_angle << endl;
	cout << "   clockwise:             " << header.clockwise << endl;
	for(int i = 0; i < header.num_of_blocks; i++)
	cout << "   block[" << i << "] centre:       " << header.centre_of_rotation[i] << endl;
	cout << "   reconstruction scale:  " << header.scale << endl;



}

int cra_con::LoadData(float *_data)
{
	int error = 0;
	
	string fn;
	fn.assign(path);
	fn.append("\\");
	fn.append(header.raw_file);
	
	block = (unsigned long long int)header.detector_x * header.detector_z * header.num_of_projections;

	if( error = fopen_s(&fstream, fn.c_str(), "rb")) 
		return error;

	int proj_size = header.detector_x * header.detector_z; /// velikost projekce

	float *data = new float[proj_size];

	unsigned int num_of_projections = GetNoProjection(); // pocet projekci k nacteni

	unsigned long long int sin_size = header.detector_x * num_of_projections;  /// velikost sinogramu

	unsigned long long int offset = SinActiveBlock * SinNoOfProjections * proj_size + block * ActiveBlock;
	offset *= sizeof(data[0]);

	_fseeki64(fstream, offset, 0);

	float p = num_of_projections / 100.;
	printf("\nloading data %.2f %%      \r",0.);

	if(header.clockwise)
	{
		/*{
			// nacteni 0 uhlu, protoze ten je bran jako startovni uhel
			if( fread(data, sizeof(data[0]), proj_size, fstream) !=  proj_size)
			{
				return -1;	
			}

			for(int s = 0; s < header.detector_z; s++) /// pres vsechny sinogramy
			{
				float *row = data + s * header.detector_x;
				float *new_row = _data + s * sin_size + 0 * header.detector_x;
				memcpy(new_row, row, header.detector_x * sizeof(row[0]));
			}
		}*/

		for(long long a = num_of_projections - 1; a >= 0 ; a--) /// pres vsechny uhly
		{
			
			if( fread(data, sizeof(data[0]), proj_size, fstream) !=  proj_size)
			{
				return -1;	
			}

			for(int s = 0; s < header.detector_z; s++) /// pres vsechny sinogramy
			{
				float *row = data + s * header.detector_x;
				float *new_row = _data + sin_size * s + a * header.detector_x;
				memcpy(new_row, row, header.detector_x * sizeof(row[0]));
			}
			printf("loading data %.2f %%      \r",(num_of_projections - a) / p);
		}
	}
	else	
	{
		for(long long a = 0; a < num_of_projections; a++) /// pres vsechny uhly
		{
			if( fread(data, sizeof(data[0]), proj_size, fstream) !=  proj_size)
			{
				return -1;	
			}

			for(int s = 0; s < header.detector_z; s++) /// pres vsechny sinogramy
			{
				float *row = data + s * header.detector_x;
				float *new_row = _data + sin_size * s + a * header.detector_x;
				memcpy(new_row, row, header.detector_x * sizeof(row[0]));
			}
			printf("loading data %.2f %%      \r",(a + 1) / p);
		}
	}

	cout << endl;
	delete[] data;
	fclose(fstream);
	fstream = nullptr;
	return error;
}

void cra_con::GetSize(unsigned int &_sin_x, unsigned int &_sin_y, unsigned int &_sin_z)
{
	_sin_x = header.detector_x;
	_sin_z = header.detector_z;
	_sin_y = header.num_of_projections;
}

void cra_con::GetRecSize(unsigned int &_rec_x, unsigned int &_rec_y, unsigned int &_rec_z)
{
	if(header.x_dimension == 0 )
		_rec_x = _rec_y = header.detector_x;
	else
		_rec_x = _rec_y = header.x_dimension;

	_rec_z = header.z_dimension;
}

double cra_con::GetCentreOfRotation()
{
	return header.centre_of_rotation[ActiveBlock];
}

double cra_con::GetDistance()
{
	return distance;
}

double cra_con::GetDFI()
{
	return delta_fi;
}

double cra_con::GetRotRad()
{	
	double result = 0;

	if(header.clockwise) // korekce zacatku sinogramu v pripade clockwise rotace
	{
		result += delta_fi;

		if(NoOfSinBlocks > 1) // pokud je vice nez jeden blok musi se zacit pootoceno o spravny pocet rad 
		{
			int x = 0;
			if(SinRestOfProjection > 0)
				x = NoOfSinBlocks - 2;
			else
				x = NoOfSinBlocks - 1;

			result += delta_fi * SinRestOfProjection + delta_fi * x * SinNoOfProjections;
			
			if(SinActiveBlock == NoOfSinBlocks - 1 && SinRestOfProjection > 0)
				result -= (SinActiveBlock - 1) * SinNoOfProjections * delta_fi + SinRestOfProjection * delta_fi;   /// < minus pro clockwise + pro anti-clockwise
			else
				result -= SinActiveBlock * SinNoOfProjections * delta_fi;   
		}
	}
	else
		result += SinActiveBlock * SinNoOfProjections * delta_fi;			

	//return (header.rotate_angle * PI / 180.) + SinActiveBlock * 2. * PI / header.num_of_projections;
	return result + (header.rotate_angle * PI / 180.);
}

double cra_con::GetSampling()
{
	return header.sample_spacing;
}

double cra_con::GetScale()
{
	return header.scale;
}

int cra_con::GetNumberOfBlocks()
{
	return header.num_of_blocks;
}

int cra_con::SetActiveBlock(int _block)
{
	if(_block < header.num_of_blocks)
		ActiveBlock = _block;
	else
		ActiveBlock = 0;

	return ActiveBlock;
}

long long cra_con::SetSinBlockMaxSize(long long int _size)
{
	SinBlockSize = _size;
	return SinBlockSize;
}

void cra_con::CalculateSinBlocks()
{
	SinogramSize = header.detector_x;
	SinogramSize *= header.detector_z;
	SinogramSize *= header.num_of_projections;
	SinogramSize *= sizeof(float);

	NoOfSinBlocks = SinogramSize / SinBlockSize;

	long long int tmp = header.detector_x;
	tmp *= header.detector_z;
	tmp *= sizeof(float);

	SinNoOfProjections =  SinBlockSize / tmp;
	SinRestOfProjection = 0;

	if(SinogramSize % SinBlockSize > 0)
	{
		SinRestOfProjection = header.num_of_projections - NoOfSinBlocks * SinNoOfProjections;
		NoOfSinBlocks += 1;
	}
}

int cra_con::GetSinNoBlocks()
{
	return NoOfSinBlocks;
}

void cra_con::GetSinBlockSize(unsigned int &_sin_x, unsigned int &_sin_y, unsigned int &_sin_z)
{
	_sin_x = header.detector_x;
	_sin_z = header.detector_z;
	_sin_y = GetNoProjection();
}

int cra_con::SetSinActiveBlock(int _block)
{
	if(_block < NoOfSinBlocks)
		SinActiveBlock = _block;
	else
		SinActiveBlock = 0;

	return SinActiveBlock;
}

int cra_con::GetNoProjection()
{
	int result = 0;
	if( SinActiveBlock < NoOfSinBlocks)
		result = SinNoOfProjections;

	if( (SinActiveBlock == (NoOfSinBlocks - 1)) && (SinRestOfProjection > 0))
		result = SinRestOfProjection;

	return result;
}

///******************************* results *****************************/////

bin_siz::bin_siz(string _path, string _filename)
{
	im_x = 0; 
	im_y = 0;
	im_z = 0;

	filename = _filename;
	path = _path;
	fstream = nullptr;
	im_size = 0;

	block_size = 0;
	slices_per_block = 16;
	slices_rest = 0;
	scan_block_offset =  0;
	scan_block_offset_memory = 0;
	scan_blocks_count = 1;

	h_file.assign(path);
	h_file.append("\\");
	h_file.append(filename);
	h_file.append(".siz");

	b_file.assign(path);
	b_file.append("\\");
	b_file.append(filename);
	b_file.append(".bin");

	ThreadIO = nullptr;

	sync_flag = false;

	for(int i = 0; i < ALMEM; i++)
	{
		data[i].data = nullptr;
		data[i].collecting = false;
		data[i].block = 0;
		data[i].filepos = 0;
		data[i].path = b_file;
		data[i].data_size = 0;
		data[i].type = FREE;
	}

//	data[0].type = WRITE;
//	data[1].type = WORK;
//	data[2].type = READ;

	counter.SetNumber(ALMEM);

}

bin_siz::~bin_siz()
{
	for(int i = 0; i < ALMEM; i++)
		if(data[i].data)
			delete[] data[i].data;

	SaveHeader();
}

void bin_siz::SetImageSize(int _x, int _y, int _z)
{
	im_x = _x; 
	im_y = _y;
	im_z = _z;
	im_size = _x * _y;

	CreateBinFile();
}

void bin_siz::SetSlicesPerBlock(int _z)
{
	slices_per_block = _z;
	block_size = slices_per_block * im_size;

	blocks_per_reconstruction = ceil((float)im_z / (float)slices_per_block);
	slices_rest = im_z - (im_z / slices_per_block)  * slices_per_block;

	for(int i = 0; i < ALMEM; i++)
		data[i].data = new float[block_size];

}

void bin_siz::SetScanBlockOffset(int _offset)
{
	scan_block_offset = _offset;
	scan_block_offset_memory = im_size * scan_block_offset * im_z * sizeof(float); 
}

void bin_siz::SetScanBlocksCount(int _sb)
{
	scan_blocks_count = _sb;
}

float bin_siz::GetBufferSize()
{
	return (float)ALMEM * block_size * sizeof(float) / 1024. / 1024.; 
}

int bin_siz::SaveHeader()
{
		if( error = fopen_s(&fstream, h_file.c_str(), "w")) 
			return error;
		fprintf(fstream,"%d %d %d 0", im_x, im_y, im_z * scan_blocks_count);
		fclose(fstream);

		return error;
}

int bin_siz::CreateBinFile()
{
	
	LARGE_INTEGER new_size;
	new_size.QuadPart = (long long int)scan_blocks_count * im_size * im_z * sizeof(float);

	HANDLE hFile = CreateFileA( b_file.c_str(), GENERIC_WRITE, 0, nullptr, CREATE_ALWAYS, FILE_ATTRIBUTE_NORMAL, nullptr);
	if(hFile == INVALID_HANDLE_VALUE)
	{
		error = GetLastError();
		cerr << "file creation error " << error << endl;
		return error;
	}

	SetFilePointerEx(hFile, new_size, nullptr, FILE_BEGIN);
	SetEndOfFile(hFile);
	CloseHandle(hFile);

	
	cout <<"reconstruction file created: " << b_file << endl << endl;
	
	/*if( error = fopen_s(&fstream, b_file.c_str(), "wb")) 
			return error;

		_fseeki64(fstream, (unsigned long long)scan_blocks_count * im_size * im_z * sizeof(float)-1, SEEK_SET);
		fwrite("0",1,1, fstream);

		fclose(fstream);

		cout <<"reconstruction file created: " << b_file << endl;
		*/

	return error;
}

float *bin_siz::GetBlock(int _block_n, bool _collecting)
{
	float *pointer = nullptr;

	if(sync_flag)
	{
		if(WaitForSingleObject(ThreadIO, 60000) == WAIT_TIMEOUT)
			cerr << "IO thread timeout" << endl;
		CloseHandle(ThreadIO);
		ThreadIO = nullptr;
		sync_flag = false;
	}

//	if(data[counter + 0].block != _block_n)
	counter++;
	//pointer = data[counter++].data;

	IO *buffer;
	buffer = data + (counter - 1);
	buffer->type = WRITE;

	buffer = data + (counter + 0);
	buffer->type = WORK;
	buffer->collecting = _collecting;
	buffer->block = _block_n;
	buffer->filepos = scan_block_offset_memory + im_size * _block_n * slices_per_block * sizeof(float);
	pointer = buffer->data;

	if(blocks_per_reconstruction - 1 == _block_n && slices_rest) // last block
		buffer->data_size = im_size * slices_rest * sizeof(float);
	else
		buffer->data_size = block_size * sizeof(float);


	buffer = data + (counter + 1);
	buffer->type = READ;
	if(blocks_per_reconstruction - 1 == _block_n)
	{
		buffer->collecting = true;
		buffer->block = 0;
		buffer->filepos = scan_block_offset_memory; // prvni blok v souboru
		buffer->data_size = block_size * sizeof(float);
	}
	else
	{
		buffer->collecting = _collecting;
		buffer->filepos = scan_block_offset_memory + (long long)(_block_n + 1) * im_size * slices_per_block * sizeof(float);
		buffer->data_size = block_size * sizeof(float);
	}

	//cout << endl << endl;
	//cout << "block: " << _block_n << endl;
	//cout << data[0].data_size << " " << data[0].type << " " << data[0].filepos << " " << 0 << endl ;
	//cout << data[1].data_size << " " << data[1].type << " " << data[1].filepos << " " << 1 << endl ;
	//cout << data[2].data_size << " " << data[2].type << " " << data[2].filepos << " " << 2 << endl ;


	ThreadIO = CreateThread(nullptr, 0, FileIO, &data, 0, nullptr);
    if (!ThreadIO)
	{
		unsigned int exit = 0;
		ExitProcess(exit);
		cerr << "IO thread exited with the flag: " << exit << endl;
	}
 
	return pointer;
}

void bin_siz::FlushBuffer()
{
	GetBlock(0);
	int cnt = 0;
	while(!sync_flag && cnt++ != 60)
		Sleep(1000);

	if(sync_flag)
	{
		WaitForSingleObject(ThreadIO, 60000);
		CloseHandle(ThreadIO);
		ThreadIO = nullptr;
		sync_flag = false;
	}

	for(int i = 0; i < ALMEM; i++)
	{
		data[i].collecting = false;
		data[i].data_size = 0;
	}

}


float *bin_siz::GetBlock_t(int _actual_block, bool _collecting)
{
	float *pointer = nullptr;
	IO *buffer;

	if(sync_flag)
	{
		if(WaitForSingleObject(ThreadIO, 60000) == WAIT_TIMEOUT)
			cerr << "IO thread timeout" << endl;
		CloseHandle(ThreadIO);
		ThreadIO = nullptr;
		sync_flag = false;
	}

/*	buffer = data + (counter + 0);
	if(buffer->block == _actual_block)		
		return buffer->data;
		*/

	counter++;

	buffer = data + (counter - 1);
	buffer->type = WRITE;

	buffer = data + (counter + 0);
	buffer->type = WORK;
	buffer->collecting = _collecting;
	buffer->block = _actual_block;
	buffer->filepos = scan_block_offset_memory + im_size * _actual_block * slices_per_block * sizeof(float);
	pointer = buffer->data;

	if(blocks_per_reconstruction - 1 == _actual_block && slices_rest) // last block
		buffer->data_size = im_size * slices_rest * sizeof(float);
	else
		buffer->data_size = block_size * sizeof(float);


	buffer = data + (counter + 1);
	buffer->type = READ;

	if(blocks_per_reconstruction - 1 == _actual_block)
	{
		buffer->collecting = true;
		buffer->block = 0;
		buffer->filepos = scan_block_offset_memory; // prvni blok v souboru
		buffer->data_size = block_size * sizeof(float);
	}
	else
	{
		buffer->collecting = _collecting;
		buffer->filepos = scan_block_offset_memory + (long long)(_actual_block + 1) * im_size * slices_per_block * sizeof(float);
		buffer->data_size = block_size * sizeof(float);
	}

	ThreadIO = CreateThread(nullptr, 0, FileIO, &data, 0, nullptr);
	if (!ThreadIO)
	{
		unsigned int exit;
		ExitProcess(exit);
		cerr << "IO thread exited with the flag: " << exit << endl;
	}

	return pointer;
}

void bin_siz::LastBlock()
{

}


DWORD WINAPI bin_siz::FileIO( LPVOID lpParam )
{
	sync_flag = true;

	pIO data = (pIO) lpParam;

	FILE *fstream = nullptr;
	if( fopen_s(&fstream, data[0].path.c_str(), "r+b"))
			cerr << "file open error" << endl;
	int rd = 0;

	for(int i = 0; i < ALMEM; i++)
	{
		if(!data[i].data_size) /// pokud je velikost dat 0 tak pokracuju
		{
			//cout << data[i].data_size << " " << i << " null" << endl ;
			continue;
		}

		_fseeki64(fstream, data[i].filepos, SEEK_SET);	
		//cout << data[i].type << endl;
		switch ( data[i].type )
		{
			case WRITE:
				{
					if(rd = fwrite(data[i].data, 1, data[i].data_size, fstream) != data[i].data_size)
							cerr << endl << "output file write error " << rd << endl << endl;
					
					//cout << "write: " << data[i].data_size << " " << data[i].filepos << " " << endl ;

					data[i].type = FINISH;
					break;
				}

			case READ:
				{
					if(data[i].collecting) 
						if(rd = fread(data[i].data, 1, data[i].data_size, fstream) != data[i].data_size)
							cerr << endl << "output file read error "<< rd << endl << endl;
					break;
				}
			case WORK:
				{
					break;
				}			
			case FINISH:
				{
					break;
				}
		}
	}

	fclose(fstream);
	return 0; 
}


//******* cyclic counter ***************//

cyclic_counter::cyclic_counter(unsigned int _numbers)
{
	numbers = _numbers;
	state = 0;
}

cyclic_counter::cyclic_counter(const cyclic_counter &_value)
{
	this->numbers = _value.numbers;
	this->state = _value.state;
}

cyclic_counter::~cyclic_counter()
{

}

unsigned int  cyclic_counter::operator++ (int)
{
	state++;
	if(state == numbers)	
		state = 0;
	return state;
}

unsigned int cyclic_counter::GetState()
{
	return state;
}

cyclic_counter cyclic_counter::operator= (const cyclic_counter &_value)
{
	this->state = _value.state;
	this->numbers = _value.numbers;
	return _value;
}

cyclic_counter cyclic_counter::operator= (const unsigned int &_value)
{
	if(_value >= numbers)
		state = _value % numbers;
	else
		state = _value;
	
	return *this;
}

void cyclic_counter::SetNumber(const unsigned int &_numbers)
{
	numbers = _numbers;
}

unsigned int cyclic_counter::operator+ (const int &_value)
{
	int result = state + _value;
	if(result >= numbers)
		result %= numbers;

	return result;
}

unsigned int cyclic_counter::operator- (const int &_value)
{
	int result = state - _value;
	if(result >= numbers)
		result %= numbers;

	if(result < 0) 
	{
		result %= numbers;
		result += numbers;
	}

	return result;

}
