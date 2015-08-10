#include "image_file.h"



image_file::image_file(std::string _filename, std::string _output_dir)
{
	filename = _filename;
	output_dir = _output_dir;

	fstream = NULL;

	im_size_x = 0;
	im_size_y = 0;
	im_count = 0;
	block_size = 0;

	bfirst = true;

	tiffType = BW16;
}


image_file::~image_file()
{
	remove(filename.c_str());
}


int image_file::AppendImages(float *_data, int _cnt)
{
	if(bfirst) 
	{
		max = min = _data[0];
		bfirst = false;
	}

	unsigned long int size = block_size * _cnt;
	for(unsigned long int s = 0; s < size; s++)
	{
		if(_data[s] > max) max = _data[s];
		if(_data[s] < min) min = _data[s];		
	}
		
	int error;

//	cout << filename.c_str() << endl;
	if( error = fopen_s(&fstream, filename.c_str(), "ab")) 
		return error;
	
	double p = _cnt / 100.;
	for(int i = 0; i < _cnt; i++)
	{
		float *block = _data + i * block_size;
		fwrite(block, sizeof(_data[0]), block_size, fstream);
		printf("appending results %.2f %%      \r", (i + 1) / p);
	}

	im_count += _cnt;

	fclose(fstream);

	return error;
}

void image_file::SetImageSize(int _x, int _y)
{
	im_size_x = _x;
	im_size_y = _y;
	block_size = _x *_y;
}

void image_file::SetTiffType(int _type)
{
	tiffType = _type;
}


int image_file::ConvertToTiff()
{
	int error;

	if( error = fopen_s(&fstream, filename.c_str(), "rb")) 
		return error;

	switch(tiffType)
	{
		case BW16:
			image = new unsigned short[block_size];
			break;
		
		case BW8:
			image = new char[block_size];
			break;

		case RGB8:
			image = new char[block_size * 3];
			break;
	}

	float *data = new float[block_size];
	char rname[1024];


	double p = im_count / 100.;

	cout << endl << "results max: " << max << endl << "results min: " << min << endl;
	for(int c = 0; c < im_count; c++)
	{
		fread(data, sizeof(data[0]), block_size, fstream);

		unsigned short * dstPtr = static_cast< unsigned short * >(image);

		if(max != min)
			for(int x = 0; x < block_size; x++)
				dstPtr[x] = (unsigned short)(65535 * (data[x] - min)/(max - min));

		sprintf(rname, "%s\\%08d.tif",output_dir.c_str(), c);
		TiffWrite(rname);
	
		printf("converting results %.2f %%      \r", (c + 1) / p);
	}

	fclose(fstream);
	delete[] data;
	delete[] image;

	return error;
}

int image_file::TiffWrite (const std::string &_filename)
{
	unsigned long stripBytes, numStrips, imageWidth, imageHeight;
	unsigned short bitsPerSample;
	unsigned short samplesPerPixel;
	unsigned short photometric;
	struct	tiff * m_tiff;


	if ((m_tiff = TIFFOpen(_filename.c_str(), "w" )) != NULL )
	{
	
		switch (tiffType) 
		{
			default:
				
			case BW16:
				{
					bitsPerSample = 16;
					imageHeight = numStrips = im_size_y;
					stripBytes = im_size_x * 2;
					imageWidth = im_size_x;
					samplesPerPixel = 1;
					photometric = PHOTOMETRIC_MINISBLACK;
					break;
				}
			case BW8:
				{
					bitsPerSample = 8;
					imageHeight = numStrips = im_size_y;
					imageWidth = stripBytes = im_size_x;
					samplesPerPixel = 1;
					photometric = PHOTOMETRIC_MINISBLACK;
					break;
				}
			case RGB8:
				{
					bitsPerSample = 8;
					imageHeight = numStrips = im_size_y;
					stripBytes = im_size_x * 3;
					imageWidth = im_size_x;
					samplesPerPixel = 3;
					photometric = PHOTOMETRIC_RGB;
					break;
				}
	
		}
		// All required fields for RGB and grayscale, except Xres, Yres, and res-unit.
		TIFFSetField( m_tiff, TIFFTAG_ROWSPERSTRIP, 1 );
		TIFFSetField( m_tiff, TIFFTAG_BITSPERSAMPLE, bitsPerSample );
		TIFFSetField( m_tiff, TIFFTAG_SAMPLESPERPIXEL, samplesPerPixel );
		TIFFSetField( m_tiff, TIFFTAG_IMAGEWIDTH,  imageWidth );
		TIFFSetField( m_tiff, TIFFTAG_IMAGELENGTH, imageHeight );
		TIFFSetField( m_tiff, TIFFTAG_PHOTOMETRIC, photometric );
		TIFFSetField( m_tiff, TIFFTAG_PLANARCONFIG, 1 );	// RGBRGB
		TIFFSetField( m_tiff, TIFFTAG_COMPRESSION,  1 );	// None

		switch (tiffType) 
		{
			case BW16: 
				{
					 // Grayscale 16 bit
					unsigned short * srcPtr = static_cast< unsigned short * >(image);
					unsigned long	stripWords;			// unsigned short per strip
					stripWords = stripBytes >> 1;
					for ( unsigned long strip = 0; strip < numStrips; strip++ ) {
						TIFFWriteEncodedStrip( m_tiff, strip, srcPtr, stripBytes );
						srcPtr += im_size_x;
					 }
					break;
				}

			case BW8: 
				{
					unsigned char * srcPtr = static_cast< unsigned char * >(image);
					for ( unsigned long strip = 0; strip < numStrips; strip++ ) {
						TIFFWriteEncodedStrip( m_tiff, strip, srcPtr, stripBytes );
						srcPtr += im_size_x;
					}
					break;
				}

			case RGB8: 
				{
					unsigned char * srcPtr = static_cast< unsigned char * >(image);
					for ( unsigned long strip = 0; strip < numStrips; strip++ ) {
						TIFFWriteEncodedStrip( m_tiff, strip, srcPtr, stripBytes );
						srcPtr += stripBytes;//_resolution.GetX();
					}
					break;
				}

		}

		TIFFClose( m_tiff );
	}
	else
		return 1;

	return 0;
}
