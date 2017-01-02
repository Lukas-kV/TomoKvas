#ifndef PARATC_H
#define PARATC_H

#include "kernel.cuh"
#include "rsq_load.h"
#include "image_file.h"
#include "mucat_file.h"
#include "FilteredBackProjection.h"

#include "tiffio.h"
#include <time.h>
#include "tinyxml.h"

using namespace std;

typedef struct _setting
{
	int GPUmemLimit;
	int GPUId;
	int GPUAllocLimit;
	//int HostMemLimit;
	double VerticalCenter; ///< posunuti stredu snimace vertikalne 
	int NoOfSlicesInOneBlock; /// pocet slices v jednom bloku ktere budou rekonstruopvany naraz
	bool x_direction, y_direction, z_direction;
	bool WaitForEnter; 

} setting, *psetting;


bool LoadSetting(setting *set);
void SaveSetting(setting *set);
void BoolAttribute(TiXmlElement* _element, const char *name, bool *_value);

#endif // !PARATC_H
