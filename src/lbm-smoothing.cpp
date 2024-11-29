// #include "lbm.hpp"
// #include "FD.hpp"


// #ifndef MULTICOMP




// #elif defined MULTICOMP


// void LBM::Smoothing()
// {
//     double alpha = 0.1;

//     #ifdef PARALLEL 
//     #pragma omp parallel for schedule(static, 1) 
//     #endif
//     for(int i = 0; i < Nx; ++i){
//         for(int j = 0; j < Ny; ++j){
//             for(int k = 0; k < Nz; ++k){
                
//                 if(mixture[i][j][k].type == TYPE_F || mixture[i][j][k].type == TYPE_I){
//                     for(size_t a = 0; a < nSpecies; ++a){
//                         double dX_threshold = 1E-5;

//                         double dXa_dx = 0.0;
//                         double dXa_dy = 0.0;
//                         double dXa_dz = 0.0;
                        
//                         if(i-1 >=0 && i+1 < Nx)
//                             dXa_dx = fd_central(species[a][i-1][j][k].X, species[a][i][j][k].X, species[a][i+1][j][k].X, dx, species[a][i][j][k].u, mixture[i-1][j][k].type, mixture[i+1][j][k].type) ;
//                         if(j-1 >=0 && j+1 < Ny)
//                             dXa_dy = fd_central(species[a][i][j-1][k].X, species[a][i][j][k].X, species[a][i][j+1][k].X, dx, species[a][i][j][k].v, mixture[i][j-1][k].type, mixture[i][j+1][k].type) ;
//                         if(k-1 >=0 && k+1 < Nz)
//                             dXa_dz = fd_central(species[a][i][j][k-1].X, species[a][i][j][k].X, species[a][i][j][k+1].X, dx, species[a][i][j][k].w, mixture[i][j][k-1].type, mixture[i][j][k+1].type) ;

//                         if (abs(dXa_dx) > dX_threshold || abs(dXa_dy) > dX_threshold || abs(dXa_dz) > dX_threshold){
//                         // if (species[a][i][j][k].X < dX_threshold){
//                             // std::cout <<  i << " " << j << " " << k << " " << a << " | " << dXa_dx << " " << dXa_dy << " " << dXa_dz << std::endl;
//                             for(size_t l = 0; l < npop; ++l){
//                                 double stc_c   = 0.0;
//                                 double stc_x_l = 0.0;
//                                 double stc_x_r = 0.0;
//                                 short type_x_l = TYPE_O;
//                                 short type_x_r = TYPE_O;
//                                 double stc_y_l = 0.0;
//                                 double stc_y_r = 0.0;
//                                 short type_y_l = TYPE_O;
//                                 short type_y_r = TYPE_O;
//                                 double stc_z_l = 0.0;
//                                 double stc_z_r = 0.0;
//                                 short type_z_l = TYPE_O;
//                                 short type_z_r = TYPE_O;

//                                 stc_c   = species[a][i][j][k].f[l];

//                                 if(i-1 >=0 && i+1 < Nx){
//                                     stc_x_l = species[a][i-1][j][k].f[l];
//                                     stc_x_r = species[a][i+1][j][k].f[l];
//                                     type_x_l = mixture[i-1][j][k].type;
//                                     type_x_r = mixture[i+1][j][k].type;
//                                 }
//                                 if(j-1 >=0 && j+1 < Ny){
//                                     stc_y_l = species[a][i][j-1][k].f[l];
//                                     stc_y_r = species[a][i][j+1][k].f[l];
//                                     type_y_l = mixture[i][j-1][k].type;
//                                     type_y_r = mixture[i][j+1][k].type;
//                                 }
//                                 if(k-1 >=0 && k+1 < Nz){
//                                     stc_z_l = species[a][i][j][k-1].f[l];
//                                     stc_z_r = species[a][i][j][k+1].f[l];
//                                     type_z_l = mixture[i][j][k-1].type ;
//                                     type_z_r = mixture[i][j][k+1].type;
//                                 }
  
//                                 species[a][i][j][k].f[l] = species[a][i][j][k].f[l] + alpha*(fd_laplace(stc_c, stc_x_l, stc_x_r, stc_y_l, stc_y_r, stc_z_l, stc_z_r, dx, dy, dz, type_x_l, type_x_r, type_y_l, type_y_r, type_z_l, type_z_r));
//                             }
//                         }

//                     }
//                 }

//             }
//         }
//     }

// }
// #endif

#include "lbm.hpp"
#include "FD.hpp"


void LBM::Smoothing()
{
    double alpha = 0.1;

    #ifdef PARALLEL 
    #pragma omp parallel for schedule(dynamic) 
    #endif
    for(int i = 0; i < Nx; ++i){
        for(int j = 0; j < Ny; ++j){
            for(int k = 0; k < Nz; ++k){
                
                if(mixture[i][j][k].type == TYPE_F || mixture[i][j][k].type == TYPE_I){
                    for(size_t a = 0; a < nSpecies; ++a){
                        double dX_threshold = 1E-2;

                        double dXa_dx = 0.0;
                        double dXa_dy = 0.0;
                        double dXa_dz = 0.0;
                        
                        if(i-1 >=0 && i+1 < Nx)
                            dXa_dx = fd_central(species[a][i-1][j][k].X, species[a][i][j][k].X, species[a][i+1][j][k].X, dx, species[a][i][j][k].u, mixture[i-1][j][k].type, mixture[i+1][j][k].type) ;
                        if(j-1 >=0 && j+1 < Ny)
                            dXa_dy = fd_central(species[a][i][j-1][k].X, species[a][i][j][k].X, species[a][i][j+1][k].X, dx, species[a][i][j][k].v, mixture[i][j-1][k].type, mixture[i][j+1][k].type) ;
                        if(k-1 >=0 && k+1 < Nz)   
                            dXa_dz = fd_central(species[a][i][j][k-1].X, species[a][i][j][k].X, species[a][i][j][k+1].X, dx, species[a][i][j][k].w, mixture[i][j][k-1].type, mixture[i][j][k+1].type) ;

                        if (abs(dXa_dx) > dX_threshold || abs(dXa_dy) > dX_threshold || abs(dXa_dz) > dX_threshold){
                        // if (species[a][i][j][k].X < dX_threshold){
                            // std::cout <<  i << " " << j << " " << k << " " << a << " | " << dXa_dx << " " << dXa_dy << " " << dXa_dz << std::endl;
                            for(size_t l = 0; l < npop; ++l){
                                double stc_c   = 0.0;
                                double stc_x_l = 0.0;
                                double stc_x_r = 0.0;
                                double stc_y_l = 0.0;
                                double stc_y_r = 0.0;
                                double stc_z_l = 0.0;
                                double stc_z_r = 0.0;
                                short type_x_l = 0.0;
                                short type_x_r = 0.0;
                                short type_y_l = 0.0;
                                short type_y_r = 0.0;
                                short type_z_l = 0.0;
                                short type_z_r = 0.0;

                                stc_c   = species[a][i][j][k].f[l];
                                if(i-1 >=0 && i+1 < Nx){
                                    stc_x_l = species[a][i-1][j][k].f[l];
                                    stc_x_r = species[a][i+1][j][k].f[l];
                                    type_x_l = mixture[i-1][j][k].type;
                                    type_x_r = mixture[i+1][j][k].type;
                                }
                                if(j-1 >=0 && j+1 < Ny){
                                    stc_y_l = species[a][i][j-1][k].f[l];
                                    stc_y_r = species[a][i][j+1][k].f[l];
                                    type_y_l = mixture[i][j-1][k].type;
                                    type_y_r = mixture[i][j+1][k].type;
                                }
                                if(k-1 >=0 && k+1 < Nz){
                                    stc_z_l = species[a][i][j][k-1].f[l];
                                    stc_z_r = species[a][i][j][k+1].f[l];                                
                                    type_z_l = mixture[i][j][k-1].type;
                                    type_z_r = mixture[i][j][k+1].type;
                                }

                                species[a][i][j][k].f[l] = species[a][i][j][k].f[l] + alpha*(fd_laplace(stc_c, stc_x_l, stc_x_r, stc_y_l, stc_y_r, stc_z_l, stc_z_r, dx, dy, dz, type_x_l, type_x_r, type_y_l, type_y_r, type_z_l, type_z_r));
                                
                            }
                        }

                    }
                }

            }
        }
    }

}