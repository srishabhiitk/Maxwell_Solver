#ifndef GRID_H
#define GRID_H

#include "reader.h"
#include "parallel.h"

class grid {
    public:
        /** A const reference to the global variables stored in the reader class to access user set parameters */
        const reader &inputData;

        /** A const reference to the global variables stored in the parallel class to access MPI related parameters */
        const parallel &rankData;

        int local_colloq_x;
        int local_colloq_y;
        int local_colloq_z;


        int colloq_start_index_x;
        int colloq_start_index_y;
        int colloq_start_index_z;

        int colloq_end_index_x;
        int colloq_end_index_y;
        int colloq_end_index_z;

        bool isPlanar;




        grid(const reader &solParam, parallel &parallelData);
};

/**
 ********************************************************************************************************************************************
 *  \class grid grid.h "lib/grid.h"
 *  \brief  Contains all the global variables related to the grid, its slices, limits, and grid derivatives used
 *          throughout the solver
 ********************************************************************************************************************************************
 */

#endif
