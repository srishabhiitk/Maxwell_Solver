#include <iostream>
#include <string>
#include "H5Cpp.h"
using namespace H5;

const H5std_string	FILE_NAME("h5tutr_dset.h5");
const H5std_string	DATASET_NAME("dset");
const int 	DIM0 = 4;	               // dataset dimensions
const int 	DIM1 = 6;

int main (void)
{

    // Data initialization.

    int i, j;
    int data[DIM0][DIM1];	    // buffer for data to write

    for (j = 0; j < DIM0; j++)
	for (i = 0; i < DIM1; i++)
	 data[j][i] = i * 6 + j + 1;

    // Try block to detect exceptions raised by any of the calls inside it

	H5File file(FILE_NAME, H5F_ACC_RDWR);
	DataSet dataset = file.openDataSet(DATASET_NAME);


	dataset.write(data, PredType::NATIVE_INT);

    

    return 0;  // successfully terminated
}
