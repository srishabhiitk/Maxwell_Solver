#ifndef WRITER_H
#define WRITER_H

#include <blitz/array.h>
#include <iostream>
#include <string>
#include <sstream>
#include "field.h"
#include "grid.h"
#include "hdf5.h"
//using namespace H5;

class writer {
    private:
        /** A const reference to the global variables stored in the grid class to access grid parameters */
        const grid &mesh;

        /** A const reference to the scalar field data to be written into output file */
        const field *outField;


        hsize_t gloSize[3];
        hsize_t locSize[3];

        hsize_t gloOffset[3];
        hsize_t locOffset[3];


        hid_t plist_id;

        hid_t fileHandle;

        hid_t dataSet;

        hid_t sourceDSpace;
        hid_t targetDSpace;

        herr_t status;

        std::stringstream temp;
        std::string datasetName;

    public:
        writer(const char *fileName, bool read, const grid &mesh, field *iField);

        void writeHDF5(double t);
        void readHDF5(std::string datasetName);
        void closeWriter();
};

/**
 ********************************************************************************************************************************************
 *  \class writer writer.h "lib/writer.h"
 *  \brief Class for all the global variables and functions related to writing output data of the solver.
 *
 *  The computational data from the solver is written in ASCII format in a .dat file.
 *  The data is written to be compatible with Tecplot and hence writes the corresponding headers in the solution files.
 *  MPI-IO is used to write the data in parallel, and hence all the associated MPI derived data-types are stored in this class.
 *  The class allows for both collocated and staggered grid data to be written in separate output files.
 ********************************************************************************************************************************************
 */

#endif
