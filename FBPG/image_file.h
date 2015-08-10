#pragma once

#include <stdio.h>
#include <iostream>
#include <algorithm>
#include <time.h>
#include <sstream>
#include "tiffio.h"


using namespace std;

#define BW16 1
#define BW8  2
#define RGB8 3



class image_file
{
public:
	image_file(std::string _tmp_filename, std::string _output_dir);
	~image_file();

	int AppendImages(float *_data, int _cnt);
	void SetImageSize(int _x, int _y);
	void SetTiffType(int _type);

	int ConvertToTiff();

private:
	std::string filename, output_dir;
	FILE *fstream;
	int im_size_x;
	int im_size_y;
	int block_size;
	int im_count;
	double max, min;
	bool bfirst;

	int tiffType;

	void *image;

	int TiffWrite (const std::string &_filename);

};