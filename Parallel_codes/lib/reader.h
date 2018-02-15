#ifndef READER_H
#define READER_H

#include <string>
#include <fstream>
#include <blitz/array.h>
#include <yaml-cpp/yaml.h>

class reader {
    public:
        int npX;
        int npY;
        double dx;
        double dy;
        double dz;
        double S; //courant factor
        int num_timesteps;
        int Nx;
        int Ny;
        int Nz;
        int n_threads;


        reader();
        void readYAML();
        void checkData();

};

/**
 ********************************************************************************************************************************************
 *  \class reader reader.h "lib/reader.h"
 *  \brief  Contains all the global variables set by the user through the yaml file
 *
 *  The class reads the paramters.yaml file and stores all the simulation paramters in publicly accessible constants.
 *  The class also has a function to check the consistency of the user set paramters and throw exceptions.
 *  The class is best initialized as a constant to prevent inadvertent tampering of the global variables it contains.
 ********************************************************************************************************************************************
 */

#endif
