#include "lbm.hpp"

void LBM::fill_BC()
{
    for(int i = 0; i < Nx; ++i){
        for(int j = 0; j < Ny; ++j){
            for(int k = 0; k < Nz; ++k){
                if(mixture[i][j][k].type == TYPE_F){
                    int i_nb = i - 1;
                    int i_pb = i + 1;
                    int j_nb = j - 1;
                    int j_pb = j + 1;
                    int k_nb = k - 1;
                    int k_pb = k + 1;

                    // Adiabatic Boundary Condition
                    if(mixture[i_nb][j][k].type == TYPE_A){
                        mixture[i_nb][j][k].rho = mixture[i][j][k].rho;
                        mixture[i_nb][j][k].u   = mixture[i][j][k].u;
                        mixture[i_nb][j][k].v   = mixture[i][j][k].v;
                        mixture[i_nb][j][k].w   = mixture[i][j][k].w;
                        mixture[i_nb][j][k].temp= mixture[i][j][k].temp;

                        #ifdef MULTICOMP
                        for(size_t a = 0; a < nSpecies; ++a)
                            species[a][i_nb][j][k].rho = species[a][i][j][k].rho;
                        #endif
                    }
                    if(mixture[i_pb][j][k].type == TYPE_A){
                        mixture[i_pb][j][k].rho = mixture[i][j][k].rho;
                        mixture[i_pb][j][k].u   = mixture[i][j][k].u;
                        mixture[i_pb][j][k].v   = mixture[i][j][k].v;
                        mixture[i_pb][j][k].w   = mixture[i][j][k].w;
                        mixture[i_pb][j][k].temp= mixture[i][j][k].temp;
                        #ifdef MULTICOMP
                        for(size_t a = 0; a < nSpecies; ++a)
                            species[a][i_pb][j][k].rho = species[a][i][j][k].rho;
                        #endif
                    }
                    if(mixture[i][j_nb][k].type == TYPE_A){
                        mixture[i][j_nb][k].rho = mixture[i][j][k].rho;
                        mixture[i][j_nb][k].u   = mixture[i][j][k].u;
                        mixture[i][j_nb][k].v   = mixture[i][j][k].v;
                        mixture[i][j_nb][k].w   = mixture[i][j][k].w;
                        mixture[i][j_nb][k].temp= mixture[i][j][k].temp;
                        #ifdef MULTICOMP
                        for(size_t a = 0; a < nSpecies; ++a)
                            species[a][i][j_nb][k].rho = species[a][i][j][k].rho;
                        #endif
                    }
                    if(mixture[i][j_pb][k].type == TYPE_A){
                        mixture[i][j_pb][k].rho = mixture[i][j][k].rho;
                        mixture[i][j_pb][k].u   = mixture[i][j][k].u;
                        mixture[i][j_pb][k].v   = mixture[i][j][k].v;
                        mixture[i][j_pb][k].w   = mixture[i][j][k].w;
                        mixture[i][j_pb][k].temp= mixture[i][j][k].temp;
                        #ifdef MULTICOMP
                        for(size_t a = 0; a < nSpecies; ++a)
                            species[a][i][j_pb][k].rho = species[a][i][j][k].rho;
                        #endif
                    }
                    if(mixture[i][j][k_nb].type == TYPE_A){
                        mixture[i][j][k_nb].rho = mixture[i][j][k].rho;
                        mixture[i][j][k_nb].u   = mixture[i][j][k].u;
                        mixture[i][j][k_nb].v   = mixture[i][j][k].v;
                        mixture[i][j][k_nb].w   = mixture[i][j][k].w;
                        mixture[i][j][k_nb].temp= mixture[i][j][k].temp;
                        #ifdef MULTICOMP
                        for(size_t a = 0; a < nSpecies; ++a)
                            species[a][i][j][k_nb].rho = species[a][i][j][k].rho;
                        #endif
                    }
                    if(mixture[i][j][k_pb].type == TYPE_A){
                        mixture[i][j][k_pb].rho = mixture[i][j][k].rho;
                        mixture[i][j][k_pb].u   = mixture[i][j][k].u;
                        mixture[i][j][k_pb].v   = mixture[i][j][k].v;
                        mixture[i][j][k_pb].w   = mixture[i][j][k].w;
                        mixture[i][j][k_pb].temp= mixture[i][j][k].temp;
                        #ifdef MULTICOMP
                        for(size_t a = 0; a < nSpecies; ++a)
                            species[a][i][j][k_pb].rho = species[a][i][j][k].rho;
                        #endif
                    }
                    
                    // Periodic Boundary Condition 
                    if(mixture[i_nb][j][k].type == TYPE_P){
                        mixture[i_nb][j][k].rho = mixture[Nx-2][j][k].rho;
                        mixture[i_nb][j][k].u   = mixture[Nx-2][j][k].u;
                        mixture[i_nb][j][k].v   = mixture[Nx-2][j][k].v;
                        mixture[i_nb][j][k].w   = mixture[Nx-2][j][k].w;
                        mixture[i_nb][j][k].temp= mixture[Nx-2][j][k].temp;
                        #ifdef MULTICOMP
                        for(size_t a = 0; a < nSpecies; ++a)
                            species[a][i_nb][j][k].rho = species[a][Nx-2][j][k].rho;
                        #endif
                    }
                    if(mixture[i_pb][j][k].type == TYPE_P){
                        mixture[i_pb][j][k].rho = mixture[1][j][k].rho;
                        mixture[i_pb][j][k].u   = mixture[1][j][k].u;
                        mixture[i_pb][j][k].v   = mixture[1][j][k].v;
                        mixture[i_pb][j][k].w   = mixture[1][j][k].w;
                        mixture[i_pb][j][k].temp= mixture[1][j][k].temp;
                        #ifdef MULTICOMP
                        for(size_t a = 0; a < nSpecies; ++a)
                            species[a][i_pb][j][k].rho = species[a][1][j][k].rho;
                        #endif
                    }
                    if(mixture[i][j_nb][k].type == TYPE_P){
                        mixture[i][j_nb][k].rho = mixture[i][Ny-2][k].rho;
                        mixture[i][j_nb][k].u   = mixture[i][Ny-2][k].u;
                        mixture[i][j_nb][k].v   = mixture[i][Ny-2][k].v;
                        mixture[i][j_nb][k].w   = mixture[i][Ny-2][k].w;
                        mixture[i][j_nb][k].temp= mixture[i][Ny-2][k].temp;
                        #ifdef MULTICOMP
                        for(size_t a = 0; a < nSpecies; ++a)
                            species[a][i][j_nb][k].rho = species[a][i][Ny-2][k].rho;
                        #endif
                    }
                    if(mixture[i][j_pb][k].type == TYPE_P){
                        mixture[i][j_pb][k].rho = mixture[i][1][k].rho;
                        mixture[i][j_pb][k].u   = mixture[i][1][k].u;
                        mixture[i][j_pb][k].v   = mixture[i][1][k].v;
                        mixture[i][j_pb][k].w   = mixture[i][1][k].w;
                        mixture[i][j_pb][k].temp= mixture[i][1][k].temp;
                        #ifdef MULTICOMP
                        for(size_t a = 0; a < nSpecies; ++a)
                            species[a][i][j_pb][k].rho = species[a][i][1][k].rho;
                        #endif
                    }
                    if(mixture[i][j][k_nb].type == TYPE_P){
                        mixture[i][j][k_nb].rho = mixture[i][j][Nz-2].rho;
                        mixture[i][j][k_nb].u   = mixture[i][j][Nz-2].u;
                        mixture[i][j][k_nb].v   = mixture[i][j][Nz-2].v;
                        mixture[i][j][k_nb].w   = mixture[i][j][Nz-2].w;
                        mixture[i][j][k_nb].temp= mixture[i][j][Nz-2].temp;
                        #ifdef MULTICOMP
                        for(size_t a = 0; a < nSpecies; ++a)
                            species[a][i][j][k_nb].rho = species[a][i][j][Nz-2].rho;
                        #endif
                    }
                    if(mixture[i][j][k_pb].type == TYPE_P){
                        mixture[i][j][k_pb].rho = mixture[i][j][1].rho;
                        mixture[i][j][k_pb].u   = mixture[i][j][1].u;
                        mixture[i][j][k_pb].v   = mixture[i][j][1].v;
                        mixture[i][j][k_pb].w   = mixture[i][j][1].w;
                        mixture[i][j][k_pb].temp= mixture[i][j][1].temp;
                        #ifdef MULTICOMP
                        for(size_t a = 0; a < nSpecies; ++a)
                            species[a][i][j][k_pb].rho = species[a][i][j][1].rho;
                        #endif
                    }


                }
            }
        }
    }
}