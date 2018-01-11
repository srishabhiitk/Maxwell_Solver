#include "mpidata.h"
#include <mpi.h>

/**
 ********************************************************************************************************************************************
 * \brief   Constructor of the mpidata class
 *
 *          The short constructor of mpidata class merely resizes the array of MPI_Status and MPI_Request datatypes.
 *          The former is used in non-blocking communication of MPI_Irecv, while the later is used in the MPI_Waitall
 *          function to complete the non-blocking communication call.
 *
 * \param   inputArray is the blitz array whose sub-arrays have to be created and synchronised across processors
 * \param   parallelData is a const reference to the global data contained in the parallel class
 ********************************************************************************************************************************************
 */
mpidata::mpidata(blitz::Array<double, 3> inputArray, bool xStag, bool yStag, const parallel &parallelData, const grid &gridData_): dataField(inputArray), rankData(parallelData), gridData(gridData_) {
    shareRanks.resize(4);
    shareRanks=rankData.nearRanks;
    if (xStag){
        shareRanks(0) = MPI_PROC_NULL;
        shareRanks(1) = MPI_PROC_NULL;
    }
    if (yStag){
        shareRanks(2) = MPI_PROC_NULL;
        shareRanks(3) = MPI_PROC_NULL;
    }

    recvStatus.resize(4);
    recvRequest.resize(4);
}

/**
 ********************************************************************************************************************************************
 * \brief   Function to create the subarray MPI_Datatypes
 *
 *          Must be called only after the grid class has been initialized.
 *          The subarray data-types cannot be created within the constructor of the parallel class as it needs the grid parameters for
 *          setting the limits of the subarrays.
 *          For this, the grid class will have to be included in the parallel class.
 *
 *          However, the grid object cannot be passed to the parallel class as the grid class already includes the parallel object
 *          within itself, and a reverse include will raise cyclic dependency error.
 *          As a result, the mpidata class offers an additional layer over the parallel class for grid specific data transfer functions.
 *
 * \param   globSize stores the global size of a sub-domain - including core and pads
 * \param   coreSize stores the size of the core of the sub-domain and is similar to the collocCoreSize variable in the grid class
 * \param   padWidth contains the widths of pads along the 3 directions, namely padWidths TinyVector from the grid class
 ********************************************************************************************************************************************
 */
void mpidata::createSubarrays(const bool xStag, const bool yStag, const bool zStag) {

    blitz::TinyVector<int, 3> subSize;
    blitz::TinyVector<int, 3> saStarts;
    blitz::TinyVector<int, 3> globSize;

    globSize(0) = gridData.local_colloq_x;
    globSize(1) = gridData.local_colloq_y;
    globSize(2) = gridData.local_colloq_z;

    if (xStag){
        globSize(0) = globSize(0) - 1;
    }
    if (yStag){
        globSize(1) = globSize(1) - 1;
    }
    if (zStag){
        globSize(2) = globSize(2) - 1;
    }
    saStarts = 0;
    subSize = globSize;



    // Subarray for left side
    subSize(0) = 1;
    saStarts(0) = 1;
    MPI_Type_create_subarray(3, globSize.data(), subSize.data(), saStarts.data(), MPI_ORDER_C, MPI_DOUBLE_PRECISION, &sendSubarrayX0);
    MPI_Type_commit(&sendSubarrayX0);

    saStarts(0) = 0;
    MPI_Type_create_subarray(3, globSize.data(), subSize.data(), saStarts.data(), MPI_ORDER_C, MPI_DOUBLE_PRECISION, &recvSubarrayX0);
    MPI_Type_commit(&recvSubarrayX0);

    // Subarray for right side
    subSize(0) = 1;
    saStarts(0) = globSize(0) - 2;
    MPI_Type_create_subarray(3, globSize.data(), subSize.data(), saStarts.data(), MPI_ORDER_C, MPI_DOUBLE_PRECISION, &sendSubarrayX1);
    MPI_Type_commit(&sendSubarrayX1);

    saStarts(0) = globSize(0) - 1;
    MPI_Type_create_subarray(3, globSize.data(), subSize.data(), saStarts.data(), MPI_ORDER_C, MPI_DOUBLE_PRECISION, &recvSubarrayX1);
    MPI_Type_commit(&recvSubarrayX1);



    subSize = globSize;
    saStarts = 0;


    // Subarray for front side
    subSize(1) = 1;
    saStarts(1) = 1;
    MPI_Type_create_subarray(3, globSize.data(), subSize.data(), saStarts.data(), MPI_ORDER_C, MPI_DOUBLE_PRECISION, &sendSubarrayY0);
    MPI_Type_commit(&sendSubarrayY0);

    saStarts(1) = 0;
    MPI_Type_create_subarray(3, globSize.data(), subSize.data(), saStarts.data(), MPI_ORDER_C, MPI_DOUBLE_PRECISION, &recvSubarrayY0);
    MPI_Type_commit(&recvSubarrayY0);

    // Subarray for back side
    subSize(1) = 1;
    saStarts(1) = globSize(1) - 2;
    MPI_Type_create_subarray(3, globSize.data(), subSize.data(), saStarts.data(), MPI_ORDER_C, MPI_DOUBLE_PRECISION, &sendSubarrayY1);
    MPI_Type_commit(&sendSubarrayY1);

    saStarts(1) = globSize(1) - 1;
    MPI_Type_create_subarray(3, globSize.data(), subSize.data(), saStarts.data(), MPI_ORDER_C, MPI_DOUBLE_PRECISION, &recvSubarrayY1);
    MPI_Type_commit(&recvSubarrayY1);

}

/**
 ********************************************************************************************************************************************
 * \brief   Function to send data across all sub-domain faces
 *
 *          This is the core function of the mpidata class.
 *          The end slices of each sub-domain recieves data from their corresponding neighbouring sub-domains,
 *          while the interior slices of each sub-domain sends data to their corresponding neighbouring sub-domains.
 *
 *          All the data slices are send as subarray MPI derived data-types created in the \ref createSubarrays function.
 *          As a result, \ref syncData must be called only after the subarrays have been created.
 *
 *          The data transfer is implemented here with a mixture of blocking and non-blocking communication calls.
 *          The recieves are non-blocking, while the sends are blocking. This combination prevents inter-processor deadlock.
 ********************************************************************************************************************************************
 */
void mpidata::syncData() {
    recvRequest = MPI_REQUEST_NULL;

    MPI_Irecv(dataField.dataFirst(), 1, recvSubarrayX0, shareRanks(0), 1, MPI_COMM_WORLD, &recvRequest(0));
    MPI_Irecv(dataField.dataFirst(), 1, recvSubarrayX1, shareRanks(1), 1, MPI_COMM_WORLD, &recvRequest(1));
    MPI_Irecv(dataField.dataFirst(), 1, recvSubarrayY0, shareRanks(2), 1, MPI_COMM_WORLD, &recvRequest(2));
    MPI_Irecv(dataField.dataFirst(), 1, recvSubarrayY1, shareRanks(3), 1, MPI_COMM_WORLD, &recvRequest(3));

    MPI_Send(dataField.dataFirst(), 1, sendSubarrayX0, shareRanks(0), 1, MPI_COMM_WORLD);
    MPI_Send(dataField.dataFirst(), 1, sendSubarrayX1, shareRanks(1), 1, MPI_COMM_WORLD);
    MPI_Send(dataField.dataFirst(), 1, sendSubarrayY0, shareRanks(2), 1, MPI_COMM_WORLD);
    MPI_Send(dataField.dataFirst(), 1, sendSubarrayY1, shareRanks(3), 1, MPI_COMM_WORLD);

    MPI_Waitall(4, recvRequest.dataFirst(), recvStatus.dataFirst());
}
