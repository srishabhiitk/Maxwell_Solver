#include <iostream>
#include <math.h>
#include <blitz/array.h>
#include "parallel.h"
#include "reader.h"
#include "grid.h"
#include "vfield.h"
#include "writer.h"
#include <sstream>
#include <sys/time.h>
#include <maxwell.h>





int main() {
    // INITIALIZE MPI
    MPI_Init(NULL, NULL);

    // ALL PROCESSES READ THE INPUT PARAMETERS
    reader inputData;

    // INITIALIZE PARALLELIZATION DATA
    parallel mpi(inputData);

    grid gridData(inputData, mpi);

    maxwell solver(inputData, mpi, gridData);
    solver.solve();


    MPI_Finalize();

    return 0;
}
