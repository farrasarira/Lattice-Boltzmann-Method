#include "lbm.hpp"

void LBM::Streaming()
{
    #ifdef PARALLEL 
        #pragma omp parallel for schedule(static, 1) 
    #endif
    for(int i=0; i<Nx; ++i)
    {
        for(int j=0; j<Ny; ++j)
        {
            for(int k = 0; k<Nz; ++k)
            {
                if(mixture[i][j][k].type==TYPE_F)
                {
                    for (int l=0; l < npop; ++l)
                    {
                        int i_nb, j_nb, k_nb;
                        i_nb = i - cx[l];
                        j_nb = j - cy[l];
                        k_nb = k - cz[l];

                        //---- Solid Boundary Condition ----------------------
                        if(mixture[i_nb][j_nb][k_nb].type==TYPE_S)
                        {
                            #ifndef MULTICOMP
                            mixture[i][j][k].f[l] = mixture[i][j][k].fpc[opposite[l]];
                            #else
                                for(size_t a = 0; a < nSpecies; ++a)
                                    species[a][i][j][k].f[l] = species[a][i][j][k].fpc[opposite[l]];
                            #endif
                        }
                        //---- Adiabatic Wall --------------------------
                        else if (mixture[i_nb][j_nb][k_nb].type==TYPE_A)
                        {
                            #ifndef MULTICOMP
                            mixture[i][j][k].f[l] = mixture[i][j][k].fpc[opposite[l]];
                            #else
                                for(size_t a = 0; a < nSpecies; ++a)
                                    species[a][i][j][k].f[l] = species[a][i][j][k].fpc[opposite[l]];
                            #endif
                            mixture[i][j][k].g[l] = mixture[i][j][k].gpc[opposite[l]];
                        }
                        //---- Inlet/Outlet Boundary Condition ---------------
                        else if (mixture[i_nb][j_nb][k_nb].type==TYPE_E)
                        {
                            // mixture[i][j][k].f[l] = mixture[i_nb][j_nb][k_nb].fpc[l];
                            mixture[i][j][k].g[l] = mixture[i_nb][j_nb][k_nb].gpc[l];
                        }
                        else //---- Periodic Boundary Condition --------------------
                        {
                            /* Alternative Periodic Code
                            if (i_nb < 1) i_nb = Nx-2;
                            else if(i_nb > Nx-2) i_nb = 1;

                            if (j_nb < 1) j_nb = Ny-2;
                            else if(j_nb > Ny-2) j_nb = 1;

                            if (k_nb < 1) k_nb = Nz-2;
                            else if(k_nb > Nz-2) k_nb = 1;
                            */

                            
                            i_nb = ((i_nb - 1 + (Nx-2)) % (Nx-2)) + 1;
                            j_nb = ((j_nb - 1 + (Ny-2)) % (Ny-2)) + 1;
                            k_nb = ((k_nb - 1 + (Nz-2)) % (Nz-2)) + 1;
                            
                            #ifndef MULTICOMP
                            mixture[i][j][k].f[l] = mixture[i_nb][j_nb][k_nb].fpc[l];
                            #else
                                for(size_t a = 0; a < nSpecies; ++a)
                                    species[a][i][j][k].f[l] = species[a][i_nb][j_nb][k_nb].fpc[l];
                            #endif
                            mixture[i][j][k].g[l] = mixture[i_nb][j_nb][k_nb].gpc[l];                            
                        }
                    }
                }
            }
        }
    }
}