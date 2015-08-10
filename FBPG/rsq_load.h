#pragma once

#include <stdio.h>
#include <iostream>
#include <algorithm>
#include <time.h>
#include <sstream>


using namespace std;


typedef struct s_ima_data_type {
/*---------------------------------------------*/
char check[16];
int data_type;
int nr_of_bytes; /* either one of them */
int nr_of_blocks; /* or both, but min. of 1 */
int patient_index; /* 1 block = 512 bytes */
int scanner_id;
int creation_date[2];
/*---------------------------------------------*/
int dimx_p;
int dimy_p;
int dimz_p;
int dimx_um;
int dimy_um;
int dimz_um;
int slice_thickness_um;
int slice_increment_um;
int slice_1_pos_um;
int min_data_value;
int max_data_value;
int mu_scaling; /* p(x,y,z)/mu_scaling = value [1/cm] */
int nr_of_samples;
int nr_of_projections;
int scandist_um;
int scanner_type;
int sampletime_us;
int index_measurement;
int site; /* Coded value */
int reference_line_um;
int recon_alg; /* Coded value */
char name[40];
int energy; /*V */
int intensity; /* uA */
int fill[83];
/*---------------------------------------------*/
int data_offset; /* in 512-byte-blocks */
} ima_data_type, *ima_data_typeP; 


typedef struct s_beam_hardening_correction 
{
	//char format[152];
	int X[3];
	double Y[3];

	double A; ///< y = A * log(x) + B
	double B;

} beam_hardening;

class rsq
{
public:
	rsq(string _filename);
	~rsq();

	void SetOffset(unsigned long int offset); ///< set reading offset
	int LoadImagesUpper(unsigned int _NumOfWholeDetector);
	int LoadImagesLower(unsigned int _NumOfWholeDetector);
	int LoadImages(unsigned int _NumOfWholeDetector); ///< load images for one cone beam
	int LoadImages(); ///< load all images to RAM
	int LoadImages(unsigned int first, unsigned int last); ///< load selected range of images to ram, first and last index of sinogram. 
	void GetSize(unsigned int &_sin_x, unsigned int &_sin_y); ///< gets size of one full sinogram
	unsigned int GetNum(); ///< gets No. of sinograms
	unsigned int GetNumLoaded(); ///< gets No. of loaded sinograms
	void SetAxisMove(unsigned int _axis);
	int ProccessAndCopy(double *_post_processed_data); ///< do calculating of flat-field correction and completing to full sinogram and copy to memory

	void SetParallel(bool _para = true); 
	void SetMinimal(bool _min = true); 
	void SetCompleting(bool _complete = true); ///< does calculating of whole sinogram
	void SetFlip(bool _flip = true); ///< flip sinogram
	void SetInvert(bool _invert = true); ///< invert sinogram
	void SetDistance(double _distance); ///< set distance in pixels
	void SetDeltaFi(double _delta_fi); ///< set delata fi
	void SetDetectorHeight(unsigned int _height); ///< set height of detector
	int	 GetDetectorHeight(); ///< get height of detector
	int	 GetNoOfWholeBlocks(); ///< get number of blocks of size of detector

private:
	string filename; ///< name and path of rsq file
	FILE *fstream; ///< file stream
	unsigned short int *idata; ///< 16 bit data from rsq file
	ima_data_type iHeader;
	beam_hardening BeamH;


	unsigned long long int header_offset; ///< offset for reading sinograms
	long long int file_size; ///< filesize
	int sin_x; ///< size of sinogram x
	int sin_x_cpl; ///< size of sinogram x after it was completed by Half2FullFan
	int sin_y; ///< size of sinogram y
	int sin_y_cpl; ///<
	int sin_y_fc; ///< size of sinogram y plus two rows for flat-field correction
	unsigned int num_of_sin; ///< calculated number of sinograms 
	unsigned int num_of_readed_sin; ///< readed number of sinograms 
	unsigned int detector_height; ///< height of detector 
	unsigned int zsize;
	bool bHeader;

	unsigned int axis_move; ///< move central axis for half sinogram completing
	double distance; ///< distance from source to detector for Half2FullFan
	double delta_fi; ///< delta of angle between rows of sinogram for Half2FullFan
	bool flip;  /// flip sinogram final sinogram
	bool invert;  /// invert sinogram final sinogram
	bool complete; ///< true for calculate full sinogram from half
	bool minimal; ///< true for calculate minimal sinogram
	bool parallel; ///< true for calculate paralel sinogram

	void CalculateFileSize(); ///< calculate file size
	void CalculateNumOfSinograms(); 
	void CalculateSizeOfSinogram(); ///< clacculate size of complete sinogram
	void Half2FullFan(double *image, double *result);
	void FlatFieldCorrection(unsigned short *image, double *result);
	void CutHalfSin(double *image, double *result);
	void Half2Para(double *image, double *result);
	double Fan2Para(double *image,unsigned  int _fi_index, unsigned int _u_index);
	void LoadHeader();


};

