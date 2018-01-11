#include "grid.h"


grid::grid(const reader &solParam, parallel &parallelData): inputData(solParam),
                                                            rankData(parallelData) {


    int Nx = inputData.Nx;
    int Ny = inputData.Ny;
    int npX = inputData.npX;
    int npY = inputData.npY;

    int xRank = rankData.xRank;
    int yRank = rankData.yRank;

    colloq_start_index_x = xRank*(Nx-2)/npX;
    colloq_start_index_y = yRank*(Ny-2)/npY;
    colloq_end_index_x = colloq_start_index_x+(Nx-2)/npX+1;
    colloq_end_index_y = colloq_start_index_y+(Ny-2)/npY+1;

    local_colloq_x = (Nx-2)/npX+2;
    local_colloq_y = (Ny-2)/npY+2;
    local_colloq_z = inputData.Nz;


}
