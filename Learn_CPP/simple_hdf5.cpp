#include <iostream>
#include <string>
#include "H5Cpp.h"
using namespace H5;

const H5std_string	FILE_NAME("h5tutr_dset.h5");
const H5std_string	DATASET_NAME("dset");
const int	 NX = 4;                     // dataset dimensions
const int	 NY = 6;
const int	 RANK = 2;

int main (void)
{



	H5File file(FILE_NAME, H5F_ACC_TRUNC);


	hsize_t dims[2];               // dataset dimensions
	dims[0] = NX;
	dims[1] = NY;
	DataSpace dataspace(RANK, dims);

	// Create the dataset.
	DataSet dataset = file.createDataSet(DATASET_NAME, PredType::STD_I32BE, dataspace);


    return 0;  // successfully terminated
}
