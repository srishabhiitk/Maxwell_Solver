#ifndef PARALLEL_H
#define PARALLEL_H

#include <blitz/array.h>
#include <mpi.h>
#include "reader.h"

class parallel {
    private:
        void assignRanks();
        void getNeighbours();

    public:
        // ALL THE INTEGERS USED BELOW ARE POSITIVE. STILL IT IS BETTER TO USE int INSTEAD OF unsigned int [1]
        /** The MPI rank of each sub-domain */
        int rank;

        /** The total number of cores available for computation */
        int nProc;

        /** npX and npY indicates the number of sub-domain divisions along the X and Y directions respectively */
        //@{
        const int npX, npY;
        //@}

        /** xRank and yRank indicates the rank in terms of sub-domain divisions along the X and Y directions respectively.
         *  Like the global rank variable, these values also start from 0 to npX -1 and npY - 1 respectively. */
        //@{
        int xRank, yRank;
        //@}

        /** Array of ranks of the 4 neighbouring sub-domains - Left, Right, Front, Back */
        blitz::Array<int, 1> nearRanks;

        parallel(const reader &iDat);

        inline int findRank(int xR, int yR);
};



#endif
