#pragma once

#include <stdio.h>
#include <iostream>
#include <algorithm>
#include <time.h>
#include <sstream>
#include <fstream>
#include "tiffio.h"
#include <time.h>

#include <windows.h> 

#define PI 3.14159265358979323846

using namespace std;


typedef struct _cra_header
{
	string raw_file;
	string source_name;
	double sample_spacing; 
	char spacing_unit; // u = um
	double distance; 
	char distance_unit; // c = cm
	int detector_x; // px
	int detector_z; // px
	int num_of_projections;
	int num_of_blocks;
	int x_dimension;
	int z_dimension;
	double rotate_angle; // degrees
	bool clockwise;
	double centre_of_rotation[64];
	double scale;

} cra_header, *pcra_header;


#define STRINGS 14

class cra_con
{
public:
	cra_con(string _path, string _filename);
	~cra_con();


	int LoadData(float *_data);
	int GetMetaData();
	void PrintMetaData();
	void GetSize(unsigned int &_sin_x, unsigned int &_sin_y, unsigned int &_sin_z);
	void GetRecSize(unsigned int &_rec_x, unsigned int &_rec_y, unsigned int &_rec_z);
	double GetCentreOfRotation();
	double GetDistance();
	double GetDFI();
	double GetRotRad();
	double GetSampling();
	double GetScale();

	int SetActiveBlock(int _block);
	int GetNumberOfBlocks();

	long long SetSinBlockMaxSize(long long int _size = (long long int)5 * 1024 * 1024 * 1024); /// default na 5 GiB
	int GetSinNoBlocks(); /// vrati pocet bloku na rozdeleni sinogramu vypocitanych z velikosti maximalni velikosti pameti
	void GetSinBlockSize(unsigned int &_sin_x, unsigned int &_sin_y, unsigned int &_sin_z); /// vraci rozmery casti sinogramu (bloku) pro aktivni block
	int SetSinActiveBlock(int _block); /// nastavi aktivni block v pripade deleni sinogramu na casti
	void CalculateSinBlocks(); /// vypocte rozdeleni sinogramu
	

private:
	string filename;
	string path;
	FILE *fstream;
	cra_header header;
	unsigned long long int block;
	float distance;
	float delta_fi;

	int ActiveBlock; /// cislo bloku skenovani

	unsigned long long int SinogramSize; // celkova velikost sinogramu
	unsigned long long int SinBlockSize; // max velikost sinogramu nactena do pameti
	int SinNoOfProjections; // pocet projekci na jeden block
	int SinRestOfProjection; // zbyvajici pocet projekci, pokud nejde sinogram rozdelit presne
	int NoOfSinBlocks; // pocet bloku do kterych bude sinogram rozdelen
	unsigned long long SinActiveBlock; /// block sinogramu ktery se bude nacitat do alokovane pameti

	int GetNoProjection(); /// vrati pocet projekci pro aktivni sin block

	string prs[STRINGS];
	clock_t start;

};


class cyclic_counter
{
public:
	cyclic_counter(unsigned int _numbers = 10);
	cyclic_counter(const cyclic_counter &_value);
	~cyclic_counter();

	unsigned int GetState();
	void SetNumber(const unsigned int &_numbers);

	unsigned int operator++ (int);
	unsigned int operator+ (const int &_value);
	unsigned int operator- (const int &_value);
	cyclic_counter operator= (const cyclic_counter &_value);
	cyclic_counter operator= (const unsigned int &_value);

private:
	int numbers;
	int state;
};


#define WRITE 0
#define READ 1
#define WORK 2
#define FINISH 3
#define FREE 4;
#define RW 5;

typedef struct _IO
{
	float *data;
	bool	collecting; 
	int		type;
	int		block;
	unsigned long long data_size;
	unsigned long long filepos;
	string path;

}	IO, *pIO;

#define ALMEM 3

class bin_siz
{
public:
	bin_siz(string _path, string _filename);
	~bin_siz();

	void SetImageSize(int _x, int _y, int _z);
	void SetScanBlocksCount(int _sb);
	
	void SetSlicesPerBlock(int _z = 16);
	void SetScanBlockOffset(int _offset);
	float GetBufferSize();
	void FlushBuffer();

	float *GetBlock(int _block_n, bool _collecting = false);

	float *GetBlock_t(int _actual_block, bool _collecting = false);
	void LastBlock();

	static DWORD WINAPI FileIO( LPVOID lpParam );

private:
	string filename;
	string h_file;
	string b_file;
	string path;
	FILE *fstream;
	IO data[ALMEM];

	cyclic_counter counter;

	unsigned int im_x, im_y, im_z;

	int scan_block_offset;
	int scan_blocks_count;

	unsigned long long int scan_block_offset_memory;
	unsigned long long int block_size;

	int slices_per_block;
	int slices_rest;

	int blocks_per_reconstruction;

	int error;
	unsigned long long int im_size;


	HANDLE ThreadIO;
	int SaveHeader();
	int CreateFile();

};

