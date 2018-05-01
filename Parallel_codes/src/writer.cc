#include "writer.h"

/**
 ********************************************************************************************************************************************
 * \brief   Constructor of the reader class
 *
 *          The constructor initializes two modes of parallel file writing through calls to their
 *          respective initialization functions.
 *
 * \param   mesh is a const reference to the global data contained in the grid class
 * \param   iField is a const reference to the scalar field whose data is to be written
 ********************************************************************************************************************************************
 */
writer::writer(const char *fileName, bool read, const grid &mesh, field *iField): mesh(mesh), outField(iField) {
    // Create a property list for collectively opening a file by all processors
    plist_id = H5Pcreate(H5P_FILE_ACCESS);
    H5Pset_fapl_mpio(plist_id, MPI_COMM_WORLD, MPI_INFO_NULL);

    // First create a file handle with the path to the file
    if (read){
        fileHandle = H5Fopen(fileName, H5F_ACC_RDONLY, plist_id);
    }
    else{
        fileHandle = H5Fcreate(fileName, H5F_ACC_TRUNC, H5P_DEFAULT, plist_id);
    }


    // Close the property list for later reuse
    H5Pclose(plist_id);




    //Preparing source Dataspace that will be written to bigger Dataspace

    locSize[0] = outField->local_Nx;
    locSize[1] = outField->local_Ny;
    locSize[2] = outField->local_Nz;

    sourceDSpace = H5Screate_simple(3, locSize, NULL);

    locOffset[0] = 0;
    locOffset[1] = 0;
    locOffset[2] = 0;

    if (outField->xStag){
        locSize[0] = locSize[0]-1; //if field is stagged
        if (mesh.rankData.xRank==(mesh.rankData.npX-1)){
            locSize[0] = locSize[0]+1;
        }
    }
    else{
        locSize[0] = locSize[0]-2; //if field is colloquated
        locOffset[0] = 1;
        if (mesh.rankData.xRank==(mesh.rankData.npX-1)){
            locSize[0] = locSize[0]+1;
        }
        if (mesh.rankData.xRank==0){
            locSize[0] = locSize[0]+1;
            locOffset[0] = 0;
        }
    }

    //in y direction
    if (outField->yStag){
        locSize[1] = locSize[1]-1;
        if (mesh.rankData.yRank==(mesh.rankData.npY-1)){
            locSize[1] = locSize[1]+1;
        }
    }
    else{
        locSize[1] = locSize[1]-2;
        locOffset[1] = 1;
        if (mesh.rankData.yRank==(mesh.rankData.npY-1)){
            locSize[1] = locSize[1]+1;
        }
        if (mesh.rankData.yRank==0){
            locSize[1] = locSize[1]+1;
            locOffset[1] = 0;
        }
    }

    status = H5Sselect_hyperslab(sourceDSpace, H5S_SELECT_SET, locOffset, NULL, locSize, NULL);
    if (status) {
        if (mesh.rankData.rank == 0) {
            std::cout << "Error in creating hyperslab while writing data. Aborting" << std::endl;
        }
        exit(0);
    }



    //Now preparing the bigger dataspace in which the data from source will be written to
    gloSize[0] = mesh.inputData.Nx;
    gloSize[1] = mesh.inputData.Ny;
    gloSize[2] = mesh.inputData.Nz;

    if (outField->xStag){
        gloSize[0] = gloSize[0]-1;
    }
    if (outField->yStag){
        gloSize[1] = gloSize[1]-1;
    }
    if (outField->zStag){
        gloSize[2] = gloSize[2]-1;
    }

    targetDSpace = H5Screate_simple(3, gloSize, NULL);

    gloOffset[0] = mesh.colloq_start_index_x + locOffset[0];
    gloOffset[1] = mesh.colloq_start_index_y + locOffset[1];
    gloOffset[2] = 0;

    status = H5Sselect_hyperslab(targetDSpace, H5S_SELECT_SET, gloOffset, NULL, locSize, NULL);
    if (status) {
        if (mesh.rankData.rank == 0) {
            std::cout << "Error in creating hyperslab while writing data. Aborting" << std::endl;
        }
        exit(0);
    }
    // std::cout<<"Processor:"<<mesh.rankData.rank<<", Global size: "<<gloSize[0]<<","<<gloSize[1]<<","<<gloSize[2]<<std::endl;
    // std::cout<<"Processor:"<<mesh.rankData.rank<<", Local size: "<<locSize[0]<<","<<locSize[1]<<","<<locSize[2]<<std::endl;
    // std::cout<<"Processor:"<<mesh.rankData.rank<<", Global offset: "<<gloOffset[0]<<","<<gloOffset[1]<<","<<gloOffset[2]<<std::endl;
    // std::cout<<"Processor:"<<mesh.rankData.rank<<", Local offset: "<<locOffset[0]<<","<<locOffset[1]<<","<<locOffset[2]<<std::endl;
}




void writer::writeHDF5(double t) {

    temp<<"t="<<t;
    datasetName = temp.str();

    // Create the dataset *for the file*, linking it to the file handle.
    // Correspondingly, it will use the *core* dataspace, as only the core has to be written excluding the pads
    dataSet = H5Dcreate2(fileHandle, datasetName.c_str(), H5T_NATIVE_DOUBLE, targetDSpace, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);

    // Create a property list to use collective data write
    plist_id = H5Pcreate(H5P_DATASET_XFER);
    H5Pset_dxpl_mpio(plist_id, H5FD_MPIO_COLLECTIVE);

    // Write the dataset. Most important thing to note is that the 3rd and 4th arguments represent the *source* and *destination* dataspaces.
    // The source here is the sourceDSpace pointing to the memory buffer. Note that its view has been adjusted using hyperslab.
    // The destination is the targetDSpace. Though the targetDSpace is smaller than the sourceDSpace,
    // only the appropriate hyperslab within the sourceDSpace is transferred to the destination.
    status = H5Dwrite(dataSet, H5T_NATIVE_DOUBLE, sourceDSpace, targetDSpace, plist_id, outField->F.dataFirst());
    if (status) {
        if (mesh.rankData.rank == 0) {
            std::cout << "Error in writing output to HDF file. Aborting" << std::endl;
        }
        exit(0);
    }

    temp.str("");
    H5Pclose(plist_id);
    H5Dclose(dataSet);
}


void writer::readHDF5(std::string datasetName) {

    // Create the dataset *for the file*, linking it to the file handle.
    // Correspondingly, it will use the *core* dataspace, as only the core has to be written excluding the pads
    dataSet = H5Dopen(fileHandle, datasetName.c_str(), H5P_DEFAULT);

    // Create a property list to use collective data write
    plist_id = H5Pcreate(H5P_DATASET_XFER);
    H5Pset_dxpl_mpio(plist_id, H5FD_MPIO_COLLECTIVE);

    // Write the dataset. Most important thing to note is that the 3rd and 4th arguments represent the *source* and *destination* dataspaces.
    // The source here is the sourceDSpace pointing to the memory buffer. Note that its view has been adjusted using hyperslab.
    // The destination is the targetDSpace. Though the targetDSpace is smaller than the sourceDSpace,
    // only the appropriate hyperslab within the sourceDSpace is transferred to the destination.
    status = H5Dread(dataSet, H5T_NATIVE_DOUBLE, sourceDSpace, targetDSpace, plist_id, (void*)outField->F.data());
    if (status) {
        if (mesh.rankData.rank == 0) {
            std::cout << "Error in writing output to HDF file. Aborting" << std::endl;
        }
        exit(0);
    }

    temp.str("");
    H5Pclose(plist_id);
    H5Dclose(dataSet);
}


void writer::closeWriter(){
    // CLOSE/RELEASE RESOURCES
    H5Sclose(sourceDSpace);
    H5Sclose(targetDSpace);
    H5Fclose(fileHandle);
}
