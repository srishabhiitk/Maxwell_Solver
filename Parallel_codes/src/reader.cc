#include <iostream>
#include "reader.h"

reader::reader() {
    readYAML();
    checkData();
}

/**
 ********************************************************************************************************************************************
 * \brief   Function to open the yaml file and read the parameters
 *
 *          The function opens the parameters.yaml file and reads the simulation parameters into its member variables that are publicly
 *          accessible.
 ********************************************************************************************************************************************
 */
void reader::readYAML() {
    std::ifstream inFile;
    inFile.open("parameters.yaml", std::ifstream::in);

    YAML::Node yamlNode;
    YAML::Parser parser(inFile);

    parser.GetNextDocument(yamlNode);
    yamlNode["X Length"] >> Nx;
    yamlNode["Y Length"] >> Ny;
    yamlNode["Z Length"] >> Nz;
    yamlNode["X Resolution"] >> dx;
    yamlNode["Y Resolution"] >> dy;
    yamlNode["Z Resolution"] >> dz;
    yamlNode["Courant factor"] >> S;
    yamlNode["Number of timesteps"] >> num_timesteps;

    parser.GetNextDocument(yamlNode);
    yamlNode["X Number of Procs"] >> npX;
    yamlNode["Y Number of Procs"] >> npY;

    inFile.close();
}

/**
 ********************************************************************************************************************************************
 * \brief   Function to perform a check on the consistency of user-set parameters
 *
 *          In order to catch potential errors early on, a few basic checks are performed here to validate the paramters set
 *          by the user.
 *          Additional checks to be performed on the paramters can be added to this function if necessary.
 ********************************************************************************************************************************************
 */
 void reader::checkData() {

     if ((Nx-2)%npX!=0) {
         std::cout << "ERROR: Condition-- (Number of points - 2) mod Number of processors --is not satisfied in X-direction. Aborting" << std::endl;
         exit(0);
     }
     if ((Ny-2)%npY!=0) {
         std::cout << "ERROR: Condition-- (Number of points - 2) mod Number of processors --is not satisfied in Y-direction. Aborting" << std::endl;
         exit(0);
     }


 }
