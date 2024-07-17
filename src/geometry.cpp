#include "./headers/lbm.hpp"
#include <omp.h>

void cylinder_generator(LBM &lb, int D)
{
    int Nx = lb.get_Nx();
    int Ny = lb.get_Ny();
    int Nz = lb.get_Nz();
    double rad=D/2;         // radius of cylinder
    double cx=(Nx/5);    // x-center of cylinder
    double cy=(Ny/2);    // y-center of cylinder

    #pragma omp parallel for
    for(int i=0; i<Nx; i++){
        for(int j=0; j<Ny; j++){
            for(int k=0; k<Nz; ++k)
            {
                double rx=i-0.5;
                double ry=j-0.5;
                if(sqrt(pow(cx-rx,2)+pow(cy-ry,2))<=rad){
                    lb.mixture[i][j][k].type=TYPE_A;
                    // lb.mixture[i][j][k].temp=0.3;
                }
            }
        }
    }
};
