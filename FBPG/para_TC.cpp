#include "para_TC.h"



int main(int argc, char* argv[])
{
	cout << "CUDA accelerated cone beam reconstruction r45" << endl;
	cout << "Lukas Kvasnica, kvasnica@ufi.fme.vutbr.cz, ikv@vranet.cz" << endl;

	if(argc != 2)
	{
		string fn;
		fn.assign(argv[0]);
		int sl = fn.rfind("\\");
		fn = fn.substr(sl + 1, fn.length());		
		cerr << "usage: " << fn << " file.cra" << endl;
		cout << "press enter" << endl;
		cin.get();
		return -1;
	}

	string fn;
	fn.assign(argv[1]);
	string path;
	string filename;
	
	int sl = fn.rfind("\\");
	int dt = fn.rfind(".");

	if (sl==string::npos)
		path = ".\\";
	else
		path = fn.substr( 0, sl);

	filename = fn.substr( sl + 1, dt - sl - 1);

	unsigned int sx,sy,sz;
	unsigned int szb; // z size of block of sinograms
	int  ScanBlock = 0; /// Scanovaci block 

	unsigned int rec_x, rec_y, rec_z;
	
	setting set; // parametry z xml souboru

	if(!LoadSetting(&set))
		SaveSetting(&set);

	size_t AvailableMemory = FilteredBackProjection::GetFreeMemInfo(set.GPUId);

	cra_con *datafile = new cra_con(path, filename);

	datafile->GetMetaData();
	datafile->PrintMetaData();
	datafile->GetSize(sx, sy, sz);
	datafile->GetRecSize(rec_x, rec_y, rec_z);
	double ReconstructionScale = datafile->GetScale();

	rec_x *= ReconstructionScale;
	rec_y *= ReconstructionScale;
	rec_z *= ReconstructionScale;

	double sinRam = (long long int)sz * sx * sy * sizeof(float) / 1024. / 1024.;
	double recRam = (long long int)rec_x * rec_y * rec_z * sizeof(float) / 1024. / 1024.;
	double recSlicesRam = (long long int)rec_x * rec_y * set.NoOfSlicesInOneBlock * sizeof(float) / 1024. / 1024.;

	set.GPUmemLimit = (AvailableMemory / 1024. / 1024. - recSlicesRam) * 0.9; // 90 procent volne pameti
	cout << endl << "GPU memory" << endl;
	cout << "  available : " <<  AvailableMemory / 1024. / 1024. << " MiB" << endl;
	cout << "  limit (90% of free memory): " << set.GPUmemLimit << " MiB" << endl;
	cout << "  allocation limit: " << set.GPUAllocLimit << "MiB" << endl;


	cout << endl;
	cout << "reconstruction:" << endl;
	cout << "  size x: " << rec_x << endl << "  size y: " << rec_y << endl << "  size z: " << rec_z << endl;
	cout << "  size: " << recRam << " MiB per scan block" << endl;
	cout << "  memory for " << set.NoOfSlicesInOneBlock << " slices: " << recSlicesRam << " MiB" << endl << endl;

	datafile->SetActiveBlock(ScanBlock);

	datafile->SetSinBlockMaxSize((long long int)set.GPUmemLimit * 1024 *1024); /// nastaveni maximalni velikosti sinogramu - podle velikosti gpu memory
	datafile->CalculateSinBlocks();

	float delta_fi = datafile->GetDFI();
	float distance = datafile->GetDistance();
	float sampling = datafile->GetSampling();
	int ScanBlocks = datafile->GetNumberOfBlocks();
	int sinBlocks = datafile->GetSinNoBlocks();

	cout << "scan blocks: " << ScanBlocks << endl << endl;
	cout << "sinograms:" << endl;
	cout << "  size x: " << sx << endl << "  size y: " << sy << endl << "  size z: " << sz << endl;
	cout << "  memory: " << sinRam << " MiB per scan block" << endl;
	cout << "  distance: " << distance << " px"<< endl;
	cout << "  vertical offset: " << set.VerticalCenter << " px"<< endl;
	cout << "  sample spacing: " << sampling << " um"<< endl;
	cout << "  delta fi: " << delta_fi << " rad"<< endl;
	cout << "  split to " << sinBlocks << " block(s)"<< endl << endl;

	bin_siz *output_file = new bin_siz(path, filename);
	output_file->SetScanBlocksCount(ScanBlocks);
	output_file->SetImageSize(rec_x, rec_y, rec_z);
	output_file->SetSlicesPerBlock(set.NoOfSlicesInOneBlock);

	cout << "allocated IO reconstruction buffer " << output_file->GetBufferSize() << " MiB for "<< set.NoOfSlicesInOneBlock << " slices" << endl << endl;

	clock_t Sblock_start = clock();
	for(ScanBlock = 0; ScanBlock < ScanBlocks; ScanBlock++)
	{
		datafile->SetActiveBlock(ScanBlock);
		output_file->SetScanBlockOffset(ScanBlock);
		float centre = datafile->GetCentreOfRotation();
		cout << endl << "*****************" << endl;
		cout << "calculating scan block: " << ScanBlock + 1 << " / " << ScanBlocks << endl;
		cout << "  centre of rotation: " << centre << " px" <<endl << endl;

		clock_t block_start = clock();

		for(int SinBlock = 0; SinBlock < sinBlocks; SinBlock++ )
		{
			datafile->SetSinActiveBlock(SinBlock);
		//	datafile->SetSinActiveBlock(0);
			datafile->GetSinBlockSize(sx, sy, szb);

			float *data = new float[(unsigned long long)sx*sy*szb];
			datafile->LoadData(data);

			float rot = datafile->GetRotRad();
			cout << "calculating sinogram block " << SinBlock + 1 << " / " << sinBlocks << endl;
			cout << "  rotation correction: " << rot << " rad" <<endl;

			FilteredBackProjection *FBP = new FilteredBackProjection();
			FBP->SetGPUid(set.GPUId);
			FBP->SetBlockSize(set.NoOfSlicesInOneBlock);
			FBP->SetAllocationLimit((long long int)set.GPUAllocLimit * 1024 * 1024);

			FBP->SetCollecting(SinBlock); /// pri nule je colecting nastaveno na false jinak true
			FBP->SetShowInfo(ScanBlock); /// informace o GPU se zobrazi jen pri nultem bloku
			FBP->SetAxisDirection(set.x_direction, set.y_direction, set.z_direction);

			FBP->SetSinogramProperties(data, sx, centre, sy, sz , ((double)sz - 1) / 2. + set.VerticalCenter, delta_fi, distance, sampling * 0.0001); // samplovani na cm
			//FBP->SetReconstructionProperties(reconstruction, rec_x, ((double)rec_x - 1) / 2., rec_y, ((double)rec_y - 1) / 2., rec_z, 0, rot);
			FBP->SetReconstructionProperties(output_file, rec_x, ((double)rec_x - 1) / 2., rec_y, ((double)rec_y - 1) / 2., rec_z, 0, rot);
			FBP->SetScaleXYZ(1./ReconstructionScale);

			if(FBP->Initialize() != cudaSuccess)
			{
				delete FBP;
				cin.get();
				return -1;
			} 
			delete[] data; // data jsou v GPU

			FBP->CalculateReconstruction();
			delete FBP;

		}

		output_file->FlushBuffer();

		cout << endl << "scan block " <<  ScanBlock + 1 << " complete time: " << GetDeltaTime(block_start, clock()) << endl;
	}

	delete datafile;
	delete output_file;

	cout << endl << "total time: " << GetDeltaTime(Sblock_start, clock()) << endl;

	//SaveSetting(&set);

	if(set.WaitForEnter)
	{
		cout << endl << "press enter" << endl;
		cin.get();
	}

    return 0;

}


void BoolAttribute(TiXmlElement* _element, const char *name, bool *_value)
{
	int i = 0;
	_element->Attribute(name,&i);
	if(i == 0) *_value = false;
	else *_value = true;
}

bool LoadSetting(setting *set)
{
	TiXmlDocument document = TiXmlDocument("setting.xml");
	document.LoadFile();
	TiXmlElement* root = document.FirstChildElement( "GPU_FBP" );

	TiXmlElement* element;

	if(root)
	{
		element = root->FirstChildElement("GPU");
		if(element)
		{
			element->Attribute("ID",&set->GPUId); 
			element->Attribute("MemLimit_MiB",&set->GPUmemLimit); 
			element->Attribute("AllocLimit_MiB",&set->GPUAllocLimit); 
		}
		else 
		{
			set->GPUId = 0;
			set->GPUmemLimit = 5 * 1024; // 5 GiB
			set->GPUAllocLimit = 1024; // 1 GiB
		}

/*		element = root->FirstChildElement("HOST");
		if(element)
		{
			element->Attribute("MemLimit_MiB",&set->HostMemLimit); 
		}
		else 
		{
			set->HostMemLimit = 15 * 1024; // 5 GiB
		}*/

		element = root->FirstChildElement("RECONSTRUCTION");
		if(element)
		{
			element->Attribute("SimultaneouslyReconstructedSlices",&set->NoOfSlicesInOneBlock); 
			BoolAttribute(element, "X_direction", &set->x_direction);
			BoolAttribute(element, "Y_direction", &set->y_direction);
			BoolAttribute(element, "Z_direction", &set->z_direction);
			element->Attribute("DetectorVerticalCenterOffset", &set->VerticalCenter);
			
		}
		else {
			set->NoOfSlicesInOneBlock = 16;
			set->x_direction = false;
			set->y_direction = false;
			set->z_direction = false;
			set->VerticalCenter = 0;
		}

		element = root->FirstChildElement("GENERAL");
		if(element)
		{
			BoolAttribute(element, "WaitForEnter", &set->WaitForEnter);
		}
		else
		{
			set->WaitForEnter = true;
		}
	
		return true;
	}
	else
	{
		set->GPUId = 0;
		set->GPUmemLimit = 5 * 1024; // 5 GiB
		set->GPUAllocLimit = 1024; // 1 GiB
//		set->HostMemLimit = 15 * 1024; // 5 GiB
		set->NoOfSlicesInOneBlock = 16;
		set->x_direction = false;
		set->y_direction = false;
		set->z_direction = false;	
		set->WaitForEnter = true;
		set->VerticalCenter = 0;

		return false;
	}

}

void SaveSetting(setting *set)
{

	TiXmlDocument document = TiXmlDocument("setting.xml");
	TiXmlElement* root;
	TiXmlElement* element;

	root = new TiXmlElement("GPU_FBP");
	document.LinkEndChild(root);

	element = new TiXmlElement("GPU");
	element->SetAttribute("ID",set->GPUId);
	element->SetAttribute("MemLimit_MiB",set->GPUmemLimit);
	element->SetAttribute("AllocLimit_MiB",set->GPUAllocLimit);
	root->LinkEndChild(element);

/*	element = new TiXmlElement("HOST");
	element->SetAttribute("MemLimit_MiB",set->HostMemLimit);
	root->LinkEndChild(element);*/

	element = new TiXmlElement("RECONSTRUCTION");
	element->SetAttribute("SimultaneouslyReconstructedSlices",set->NoOfSlicesInOneBlock);
	element->SetAttribute("X_direction",(int)set->x_direction);
	element->SetAttribute("Y_direction",(int)set->y_direction);
	element->SetAttribute("Z_direction",(int)set->z_direction);
	element->SetAttribute("DetectorVerticalCenterOffset", set->VerticalCenter);
	root->LinkEndChild(element);

	element = new TiXmlElement("GENERAL");
	element->SetAttribute("WaitForEnter",(int)set->WaitForEnter);
	root->LinkEndChild(element);

	document.SaveFile();
}
