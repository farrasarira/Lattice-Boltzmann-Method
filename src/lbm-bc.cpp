#include "lbm.hpp"
#include "math_util.hpp"
#include <omp.h>
#include "units.hpp"
#include "FD.hpp"

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
                    
                    // Periodic Boundary Condition 
                    if(mixture[i_nb][j][k].type == TYPE_P){
                        mixture[i_nb][j][k].rho = mixture[Nx-2][j][k].rho;
                        mixture[i_nb][j][k].u   = mixture[Nx-2][j][k].u;
                        mixture[i_nb][j][k].v   = mixture[Nx-2][j][k].v;
                        mixture[i_nb][j][k].w   = mixture[Nx-2][j][k].w;
                        mixture[i_nb][j][k].temp= mixture[Nx-2][j][k].temp;
                        #ifdef MULTICOMP
                        for(size_t a = 0; a < nSpecies; ++a){
                            species[a][i_nb][j][k].rho = species[a][Nx-2][j][k].rho;
                            species[a][i_nb][j][k].u = species[a][Nx-2][j][k].u;
                            species[a][i_nb][j][k].v = species[a][Nx-2][j][k].v;
                            species[a][i_nb][j][k].w = species[a][Nx-2][j][k].w;
                        }
                        #endif
                    }
                    if(mixture[i_pb][j][k].type == TYPE_P){
                        mixture[i_pb][j][k].rho = mixture[1][j][k].rho;
                        mixture[i_pb][j][k].u   = mixture[1][j][k].u;
                        mixture[i_pb][j][k].v   = mixture[1][j][k].v;
                        mixture[i_pb][j][k].w   = mixture[1][j][k].w;
                        mixture[i_pb][j][k].temp= mixture[1][j][k].temp;
                        #ifdef MULTICOMP
                        for(size_t a = 0; a < nSpecies; ++a){
                            species[a][i_pb][j][k].rho = species[a][1][j][k].rho;
                            species[a][i_pb][j][k].u = species[a][1][j][k].u;
                            species[a][i_pb][j][k].v = species[a][1][j][k].v;
                            species[a][i_pb][j][k].w = species[a][1][j][k].w;
                        }
                        #endif
                    }
                    if(mixture[i][j_nb][k].type == TYPE_P){
                        mixture[i][j_nb][k].rho = mixture[i][Ny-2][k].rho;
                        mixture[i][j_nb][k].u   = mixture[i][Ny-2][k].u;
                        mixture[i][j_nb][k].v   = mixture[i][Ny-2][k].v;
                        mixture[i][j_nb][k].w   = mixture[i][Ny-2][k].w;
                        mixture[i][j_nb][k].temp= mixture[i][Ny-2][k].temp;
                        #ifdef MULTICOMP
                        for(size_t a = 0; a < nSpecies; ++a){
                            species[a][i][j_nb][k].rho = species[a][i][Ny-2][k].rho;
                            species[a][i][j_nb][k].u = species[a][i][Ny-2][k].u;
                            species[a][i][j_nb][k].v = species[a][i][Ny-2][k].v;
                            species[a][i][j_nb][k].w = species[a][i][Ny-2][k].w;

                        }
                        #endif
                    }
                    if(mixture[i][j_pb][k].type == TYPE_P){
                        mixture[i][j_pb][k].rho = mixture[i][1][k].rho;
                        mixture[i][j_pb][k].u   = mixture[i][1][k].u;
                        mixture[i][j_pb][k].v   = mixture[i][1][k].v;
                        mixture[i][j_pb][k].w   = mixture[i][1][k].w;
                        mixture[i][j_pb][k].temp= mixture[i][1][k].temp;
                        #ifdef MULTICOMP
                        for(size_t a = 0; a < nSpecies; ++a){
                            species[a][i][j_pb][k].rho = species[a][i][1][k].rho;
                            species[a][i][j_pb][k].u = species[a][i][1][k].u;
                            species[a][i][j_pb][k].v = species[a][i][1][k].v;
                            species[a][i][j_pb][k].w = species[a][i][1][k].w;
                        
                        }
                        #endif
                    }
                    if(mixture[i][j][k_nb].type == TYPE_P){
                        mixture[i][j][k_nb].rho = mixture[i][j][Nz-2].rho;
                        mixture[i][j][k_nb].u   = mixture[i][j][Nz-2].u;
                        mixture[i][j][k_nb].v   = mixture[i][j][Nz-2].v;
                        mixture[i][j][k_nb].w   = mixture[i][j][Nz-2].w;
                        mixture[i][j][k_nb].temp= mixture[i][j][Nz-2].temp;
                        #ifdef MULTICOMP
                        for(size_t a = 0; a < nSpecies; ++a){
                            species[a][i][j][k_nb].rho = species[a][i][j][Nz-2].rho;
                            species[a][i][j][k_nb].u = species[a][i][j][Nz-2].u;
                            species[a][i][j][k_nb].v = species[a][i][j][Nz-2].v;
                            species[a][i][j][k_nb].w = species[a][i][j][Nz-2].w;
                        }
                        #endif
                    }
                    if(mixture[i][j][k_pb].type == TYPE_P){
                        mixture[i][j][k_pb].rho = mixture[i][j][1].rho;
                        mixture[i][j][k_pb].u   = mixture[i][j][1].u;
                        mixture[i][j][k_pb].v   = mixture[i][j][1].v;
                        mixture[i][j][k_pb].w   = mixture[i][j][1].w;
                        mixture[i][j][k_pb].temp= mixture[i][j][1].temp;
                        #ifdef MULTICOMP
                        for(size_t a = 0; a < nSpecies; ++a){
                            species[a][i][j][k_pb].rho = species[a][i][j][1].rho;
                            species[a][i][j][k_pb].u = species[a][i][j][1].u;
                            species[a][i][j][k_pb].v = species[a][i][j][1].v;
                            species[a][i][j][k_pb].w = species[a][i][j][1].w;
                        }
                        #endif
                    }


                }
            }
        }
    }
}

void LBM::TMS_BC()
{
    
    // Outflow =================================================================================================================================================
    #if defined MULTICOMP
    #ifdef PARALLEL 
        #pragma omp parallel for schedule(dynamic) 
    #endif
    for(int i=0; i<Nx; ++i){
        for(int j=0; j<Ny; ++j){
            for(int k = 0; k<Nz; ++k){
                if(mixture[i][j][k].type==TYPE_F){
                    int i_nb, j_nb, k_nb;
                    double fa_out[nSpecies][npop];
                    double g_out[npop];
                    double fa_loc[nSpecies][npop];
                    double g_loc[npop];
                    int l_interface = 999;

                    // check interface node                         
                    for (int l=0; l < npop; ++l){
                        i_nb = i - cx[l];
                        j_nb = j - cy[l];
                        k_nb = k - cz[l];

                        if( (mixture[i_nb][j_nb][k_nb].type==TYPE_O || mixture[i_nb][j_nb][k_nb].type==TYPE_O_C) && mixture[int (i+cx[l])][int (j+cy[l])][int (k+cz[l])].type==TYPE_F){    
                            l_interface = l;
                            break;
                        }                            
                    }

                    // check if there is the node is the interface with boundary
                    if (l_interface == 999 || !onlyOne_1(cx[l_interface], cy[l_interface], cz[l_interface]))
                        continue;

                    double rho_out = 0.0;
                    double vel_out[3] = {0.0};
                    double T_out = 0.0;
                    double rhoa_out[nSpecies] = {0.0};
                    double vela_out[nSpecies][3] = {0.0};
                    double p_out = 0.0;

                    int i_bdr = i-cx[l_interface];
                    int j_bdr = j-cy[l_interface];
                    int k_bdr = k-cz[l_interface];

                    if ( mixture[i_bdr][j_bdr][k_bdr].type == TYPE_O){
                        calculate_moment_point((int)(i+cx[l_interface]), (int)(j+cy[l_interface]), (int)(k+cz[l_interface]));

                        vel_out[0] = mixture[(int)(i+cx[l_interface])][(int)(j+cy[l_interface])][(int)(k+cz[l_interface])].u;
                        vel_out[1] = mixture[(int)(i+cx[l_interface])][(int)(j+cy[l_interface])][(int)(k+cz[l_interface])].v;
                        vel_out[2] = mixture[(int)(i+cx[l_interface])][(int)(j+cy[l_interface])][(int)(k+cz[l_interface])].w;
                        T_out      = mixture[(int)(i+cx[l_interface])][(int)(j+cy[l_interface])][(int)(k+cz[l_interface])].temp;
                        rho_out    = mixture[(int)(i+cx[l_interface])][(int)(j+cy[l_interface])][(int)(k+cz[l_interface])].rho;
                        for(size_t a = 0; a < nSpecies; ++a){
                            rhoa_out[a]    = species[a][(int)(i+cx[l_interface])][(int)(j+cy[l_interface])][(int)(k+cz[l_interface])].rho;
                            vela_out[a][0] = species[a][(int)(i+cx[l_interface])][(int)(j+cy[l_interface])][(int)(k+cz[l_interface])].u;
                            vela_out[a][1] = species[a][(int)(i+cx[l_interface])][(int)(j+cy[l_interface])][(int)(k+cz[l_interface])].v;
                            vela_out[a][2] = species[a][(int)(i+cx[l_interface])][(int)(j+cy[l_interface])][(int)(k+cz[l_interface])].w;
                        }
                        
                        p_out = mixture[(int)(i-cx[l_interface])][(int)(j-cy[l_interface])][(int)(k-cz[l_interface])].p;
                        // p_out = mixture[(int)(i+cx[l_interface])][(int)(j+cy[l_interface])][(int)(k+cz[l_interface])].p;

                        int rank = omp_get_thread_num();
                        auto gas = sols[rank]->thermo();   
                        std::vector <double> Y (gas->nSpecies());
                        for(size_t a = 0; a < nSpecies; ++a) Y[gas->speciesIndex(speciesName[a])] = rhoa_out[a] / rho_out;
                        gas->setMassFractions(&Y[0]);   
                        gas->setState_TP(units.si_temp(T_out), units.si_p(p_out));

                        rho_out = units.rho(gas->density());
                        for(size_t a = 0; a < nSpecies; ++a) rhoa_out[a] = Y[gas->speciesIndex(speciesName[a])] * rho_out;
                    }
                    else if( mixture[i_bdr][j_bdr][k_bdr].type == TYPE_O_C){
                        double sigma = 5.0;

                        int rank = omp_get_thread_num();
                        auto gas = sols[rank]->thermo();   
                        std::vector <double> Y (gas->nSpecies());
                        for(size_t a = 0; a < nSpecies; ++a) Y[gas->speciesIndex(speciesName[a])] = species[a][i][j][k].rho / mixture[i][j][k].rho;
                        gas->setMassFractions(&Y[0]);
                        gas->setState_DP(units.si_rho(mixture[i][j][k].rho), units.si_p(mixture[i][j][k].p));

                        int i_1 = i + cx[l_interface];
                        int j_1 = j + cy[l_interface];
                        int k_1 = k + cz[l_interface];
                        int i_2 = i + 2*cx[l_interface];
                        int j_2 = j + 2*cy[l_interface];
                        int k_2 = k + 2*cz[l_interface];
                        double spd_sound = units.u(gas->soundSpeed());

                        if (cx[l_interface] != 0.0){
                            double drho_dx  = fd_uw(mixture[i][j][k].rho   , mixture[i_1][j_1][k_1].rho, mixture[i_2][j_2][k_2].rho, dx, cx[l_interface]) ;
                            double du_dx    = fd_uw(mixture[i][j][k].u     , mixture[i_1][j_1][k_1].u  , mixture[i_2][j_2][k_2].u  , dx, cx[l_interface]) ;
                            double dv_dx    = fd_uw(mixture[i][j][k].v     , mixture[i_1][j_1][k_1].v  , mixture[i_2][j_2][k_2].v  , dx, cx[l_interface]) ;
                            double dw_dx    = fd_uw(mixture[i][j][k].w     , mixture[i_1][j_1][k_1].w  , mixture[i_2][j_2][k_2].w  , dx, cx[l_interface]) ;
                            double dp_dx    = fd_uw(mixture[i][j][k].p     , mixture[i_1][j_1][k_1].p  , mixture[i_2][j_2][k_2].p  , dx, cx[l_interface]) ;
                            double dYa_dx[nSpecies] = {};
                            for(size_t a = 0; a < nSpecies; ++a){
                                dYa_dx[a] = fd_uw(species[a][i][j][k].rho/mixture[i][j][k].rho, species[a][i_1][j_1][k_1].rho/mixture[i_1][j_1][k_1].rho, species[a][i_2][j_2][k_2].rho/mixture[i_2][j_2][k_2].rho, dx, cx[l_interface]);                            
                            }

                            double lambda_1 = mixture[i][j][k].u - spd_sound;
                            double lambda_2 = mixture[i][j][k].u;
                            double lambda_3 = mixture[i][j][k].u;
                            double lambda_4 = mixture[i][j][k].u;
                            double lambda_5 = mixture[i][j][k].u + spd_sound;

                            double L1 = lambda_1*(dp_dx - mixture[i][j][k].rho*spd_sound*du_dx);
                            double L2 = lambda_2*(spd_sound*spd_sound*drho_dx - dp_dx);
                            double L3 = lambda_3*(dv_dx);
                            double L4 = lambda_4*(dw_dx);
                            double L5 = lambda_5*(dp_dx + mixture[i][j][k].rho*spd_sound*du_dx);

                            double relax_K = sigma*spd_sound*(1.0 - (mixture[i][j][k].u/spd_sound)*(mixture[i][j][k].u/spd_sound))/Nx;

                            if(cx[l_interface] > 0) // left boundary
                                L5 = relax_K * ( mixture[i][j][k].p -  mixture[i_bdr][j_bdr][k_bdr].p);
                            else // right boundary
                                L1 = relax_K * ( mixture[i][j][k].p -  mixture[i_bdr][j_bdr][k_bdr].p);

                            double d1 = 1.0/(spd_sound*spd_sound) * (L2+0.5*(L5+L1));
                            double d2 = 1.0/(2.0*mixture[i][j][k].rho*spd_sound) * (L5-L1);
                            double d3 = L3;
                            double d4 = L4;
                            double d5 = 0.5*(L5+L1);
                            // double d6 = mixture[i][j][k].temp/(mixture[i][j][k].rho*spd_sound*spd_sound) * (-L2+0.5*(gas->cp_mass()/gas->cv_mass()-1.0)*(L5+L1) );

                            // rho_out    = mixture[i_1][j_1][k_1].rho;
                            // vel_out[0] = mixture[i_1][j_1][k_1].u;
                            // vel_out[1] = mixture[i_1][j_1][k_1].v;
                            // vel_out[2] = mixture[i_1][j_1][k_1].w;
                            // for(size_t a = 0; a < nSpecies; ++a)
                            //     rhoa_out[a]    = species[a][i_1][j_1][k_1].rho;
                            // p_out = mixture[i_bdr][j_bdr][k_bdr].p;

                            rho_out    = mixture[i][j][k].rho - dt_sim*d1;
                            vel_out[0] = mixture[i][j][k].u - dt_sim*d2;
                            vel_out[1] = mixture[i][j][k].v - dt_sim*d3;
                            vel_out[2] = mixture[i][j][k].w - dt_sim*d4;
                            p_out      = mixture[i][j][k].p - dt_sim*d5;
                            // T_out      = mixture[i][j][k].temp - dt_sim*d6;
                            for(size_t a = 0; a < nSpecies; ++a){
                                rhoa_out[a] = (species[a][i][j][k].rho/mixture[i][j][k].rho - dt_sim*mixture[i][j][k].u*dYa_dx[a]) * rho_out;
                                // rho_out = rho_out + dt_sim*mixture[i][j][k].u*dYa_dx[a]*mixture[i][j][k].rho;
                            }

                        }
                        else if (cy[l_interface] != 0.0){
                            double drho_dx  = fd_uw(mixture[i][j][k].rho   , mixture[i_1][j_1][k_1].rho, mixture[i_2][j_2][k_2].rho, dy, cy[l_interface]) ;
                            double du_dx    = fd_uw(mixture[i][j][k].u     , mixture[i_1][j_1][k_1].u  , mixture[i_2][j_2][k_2].u  , dy, cy[l_interface]) ;
                            double dv_dx    = fd_uw(mixture[i][j][k].v     , mixture[i_1][j_1][k_1].v  , mixture[i_2][j_2][k_2].v  , dy, cy[l_interface]) ;
                            double dw_dx    = fd_uw(mixture[i][j][k].w     , mixture[i_1][j_1][k_1].w  , mixture[i_2][j_2][k_2].w  , dy, cy[l_interface]) ;
                            double dp_dx    = fd_uw(mixture[i][j][k].p     , mixture[i_1][j_1][k_1].p  , mixture[i_2][j_2][k_2].p  , dy, cy[l_interface]) ;
                            double dYa_dx[nSpecies] = {};
                            for(size_t a = 0; a < nSpecies; ++a){
                                dYa_dx[a] = fd_uw(species[a][i][j][k].rho/mixture[i][j][k].rho, species[a][i_1][j_1][k_1].rho/mixture[i_1][j_1][k_1].rho, species[a][i_2][j_2][k_2].rho/mixture[i_2][j_2][k_2].rho, dy, cy[l_interface]);
                            }

                            double lambda_1 = mixture[i][j][k].v - spd_sound;
                            double lambda_2 = mixture[i][j][k].v;
                            double lambda_3 = mixture[i][j][k].v;
                            double lambda_4 = mixture[i][j][k].v;
                            double lambda_5 = mixture[i][j][k].v + spd_sound;

                            double L1 = lambda_1*(dp_dx - mixture[i][j][k].rho*spd_sound*dv_dx);
                            double L2 = lambda_2*(spd_sound*spd_sound*drho_dx - dp_dx);
                            double L3 = lambda_3*(dw_dx);
                            double L4 = lambda_4*(du_dx);
                            double L5 = lambda_5*(dp_dx + mixture[i][j][k].rho*spd_sound*dv_dx);

                            double relax_K = sigma*spd_sound*(1.0 - (mixture[i][j][k].v/spd_sound)*(mixture[i][j][k].v/spd_sound))/Ny;

                            if(cy[l_interface] > 0) // left boundary
                                L5 = relax_K * ( mixture[i][j][k].p -  mixture[i_bdr][j_bdr][k_bdr].p);
                            else // right boundary
                                L1 = relax_K * ( mixture[i][j][k].p -  mixture[i_bdr][j_bdr][k_bdr].p);

                            double d1 = 1.0/(spd_sound*spd_sound) * (L2+0.5*(L5+L1));
                            double d2 = 1.0/(2.0*mixture[i][j][k].rho*spd_sound) * (L5-L1);
                            double d3 = L3;
                            double d4 = L4;
                            double d5 = 0.5*(L5+L1);

                            // rho_out    = mixture[i_1][j_1][k_1].rho;
                            // vel_out[0] = mixture[i_1][j_1][k_1].u;
                            // vel_out[1] = mixture[i_1][j_1][k_1].v;
                            // vel_out[2] = mixture[i_1][j_1][k_1].w;
                            // for(size_t a = 0; a < nSpecies; ++a)
                            //     rhoa_out[a]    = species[a][i_1][j_1][k_1].rho;
                            // p_out = mixture[i_bdr][j_bdr][k_bdr].p;

                            rho_out    = mixture[i][j][k].rho - dt_sim*d1;
                            vel_out[1] = mixture[i][j][k].v - dt_sim*d2;
                            vel_out[2] = mixture[i][j][k].w - dt_sim*d3;
                            vel_out[0] = mixture[i][j][k].u - dt_sim*d4;
                            p_out      = mixture[i][j][k].p - dt_sim*d5;
                            for(size_t a = 0; a < nSpecies; ++a)
                                rhoa_out[a] = (species[a][i][j][k].rho/mixture[i][j][k].rho - dt_sim*mixture[i][j][k].v*dYa_dx[a]) * rho_out;
                        }
                        else if (cz[l_interface] != 0.0){
                            double drho_dx  = fd_uw(mixture[i][j][k].rho   , mixture[i_1][j_1][k_1].rho, mixture[i_2][j_2][k_2].rho, dz, cz[l_interface]) ;
                            double du_dx    = fd_uw(mixture[i][j][k].u     , mixture[i_1][j_1][k_1].u  , mixture[i_2][j_2][k_2].u  , dz, cz[l_interface]) ;
                            double dv_dx    = fd_uw(mixture[i][j][k].v     , mixture[i_1][j_1][k_1].v  , mixture[i_2][j_2][k_2].v  , dz, cz[l_interface]) ;
                            double dw_dx    = fd_uw(mixture[i][j][k].w     , mixture[i_1][j_1][k_1].w  , mixture[i_2][j_2][k_2].w  , dz, cz[l_interface]) ;
                            double dp_dx    = fd_uw(mixture[i][j][k].p     , mixture[i_1][j_1][k_1].p  , mixture[i_2][j_2][k_2].p  , dz, cz[l_interface]) ;
                            double dYa_dx[nSpecies] = {};
                            for(size_t a = 0; a < nSpecies; ++a){
                                dYa_dx[a] = fd_uw(species[a][i][j][k].rho/mixture[i][j][k].rho, species[a][i_1][j_1][k_1].rho/mixture[i_1][j_1][k_1].rho, species[a][i_2][j_2][k_2].rho/mixture[i_2][j_2][k_2].rho, dz, cz[l_interface]);
                            }

                            double lambda_1 = mixture[i][j][k].w - spd_sound;
                            double lambda_2 = mixture[i][j][k].w;
                            double lambda_3 = mixture[i][j][k].w;
                            double lambda_4 = mixture[i][j][k].w;
                            double lambda_5 = mixture[i][j][k].w + spd_sound;

                            double L1 = lambda_1*(dp_dx - mixture[i][j][k].rho*spd_sound*dw_dx);
                            double L2 = lambda_2*(spd_sound*spd_sound*drho_dx - dp_dx);
                            double L3 = lambda_3*(du_dx);
                            double L4 = lambda_4*(dv_dx);
                            double L5 = lambda_5*(dp_dx + mixture[i][j][k].rho*spd_sound*dw_dx);

                            double relax_K = sigma*spd_sound*(1.0 - (mixture[i][j][k].w/spd_sound)*(mixture[i][j][k].w/spd_sound))/Nz;

                            if(cy[l_interface] > 0) // left boundary
                                L5 = relax_K * ( mixture[i][j][k].p -  mixture[i_bdr][j_bdr][k_bdr].p);
                            else // right boundary
                                L1 = relax_K * ( mixture[i][j][k].p -  mixture[i_bdr][j_bdr][k_bdr].p);

                            double d1 = 1.0/(spd_sound*spd_sound) * (L2+0.5*(L5+L1));
                            double d2 = 1.0/(2.0*mixture[i][j][k].rho*spd_sound) * (L5-L1);
                            double d3 = L3;
                            double d4 = L4;
                            double d5 = 0.5*(L5+L1);

                            // rho_out    = mixture[i_1][j_1][k_1].rho;
                            // vel_out[0] = mixture[i_1][j_1][k_1].u;
                            // vel_out[1] = mixture[i_1][j_1][k_1].v;
                            // vel_out[2] = mixture[i_1][j_1][k_1].w;
                            // for(size_t a = 0; a < nSpecies; ++a)
                            //     rhoa_out[a]    = species[a][i_1][j_1][k_1].rho;
                            // p_out = mixture[i_bdr][j_bdr][k_bdr].p;

                            rho_out    = mixture[i][j][k].rho - dt_sim*d1;
                            vel_out[2] = mixture[i][j][k].w - dt_sim*d2;
                            vel_out[0] = mixture[i][j][k].u - dt_sim*d3;
                            vel_out[1] = mixture[i][j][k].v - dt_sim*d4;
                            p_out      = mixture[i][j][k].p - dt_sim*d5;
                            for(size_t a = 0; a < nSpecies; ++a)
                                rhoa_out[a] = (species[a][i][j][k].rho/mixture[i][j][k].rho - dt_sim*mixture[i][j][k].w*dYa_dx[a]) * rho_out;
                        }

                        std::fill(Y.begin(), Y.end(), 0.0);
                        for(size_t a = 0; a < nSpecies; ++a) Y[gas->speciesIndex(speciesName[a])] = rhoa_out[a] / rho_out;
                        gas->setMassFractions(&Y[0]);   
                        gas->setState_DP(units.si_rho(rho_out), units.si_p(p_out));
                        T_out = units.temp(gas->temperature());                        
                    }

                    // calculate_feq_geq(fa_out, g_out, rho_out, rhoa_out, vel_out, vela_out, T_out);
                    calculate_feq_geq(fa_out, g_out, rho_out, rhoa_out, vel_out, T_out);

                    // step 2: calculate f_loc, g_loc, and fa_loc (local distribution function)
                    double rho_loc = 0.0;
                    double vel_loc[3] = {0.0};
                    double rhoe_loc = 0.0;
                    double rhoa_loc[nSpecies] = {0.0};
                    double vela_loc[nSpecies][3] = {0.0};
                    
                    for (int l=0; l < npop; ++l){
                        i_nb = i - cx[l];
                        j_nb = j - cy[l];
                        k_nb = k - cz[l];

                        if (mixture[i_nb][j_nb][k_nb].type==TYPE_O || mixture[i_nb][j_nb][k_nb].type==TYPE_O_C){
                            for (size_t a = 0; a < nSpecies; ++a){
                                rhoa_loc[a] += fa_out[a][l];
                                vela_loc[a][0] += fa_out[a][l]*cx[l];
                                vela_loc[a][1] += fa_out[a][l]*cy[l];
                                vela_loc[a][2] += fa_out[a][l]*cz[l];
                            }
                            #ifndef ISOTHERM
                            rhoe_loc += g_out[l];
                            #endif
                        }
                        else{
                            for (size_t a = 0; a < nSpecies; ++a){
                                rhoa_loc[a] += species[a][i][j][k].f[l];
                                vela_loc[a][0] += species[a][i][j][k].f[l]*cx[l];
                                vela_loc[a][1] += species[a][i][j][k].f[l]*cy[l];
                                vela_loc[a][2] += species[a][i][j][k].f[l]*cz[l];
                            }
                            #ifndef ISOTHERM
                            rhoe_loc += mixture[i][j][k].g[l];
                            #endif
                        }
                    }

                    for (size_t a = 0; a < nSpecies; ++a){
                        rho_loc += rhoa_loc[a];
                        vel_loc[0] += vela_loc[a][0];
                        vel_loc[1] += vela_loc[a][1];
                        vel_loc[2] += vela_loc[a][2];

                        if(rhoa_loc[a] > 0.0){
                            vela_loc[a][0] = vela_loc[a][0] / rhoa_loc[a];
                            vela_loc[a][1] = vela_loc[a][1] / rhoa_loc[a];
                            vela_loc[a][2] = vela_loc[a][2] / rhoa_loc[a];
                        }
                        else{
                            vela_loc[a][0] = 0.0;
                            vela_loc[a][1] = 0.0;
                            vela_loc[a][2] = 0.0;
                        }
                    }
                    if(rho_loc > 0){
                        vel_loc[0] = vel_loc[0] / rho_loc;
                        vel_loc[1] = vel_loc[1] / rho_loc;
                        vel_loc[2] = vel_loc[2] / rho_loc;
                    }
                    else{
                        vel_loc[0] = 0.0;
                        vel_loc[1] = 0.0;
                        vel_loc[2] = 0.0;
                    }

                    #ifndef ISOTHERM
                    double internalEnergy=rhoe_loc/rho_loc;// - 0.5*v_sqr(vel_loc[0], vel_loc[1], vel_loc[2]); 
                    for (size_t a = 0; a < nSpecies; ++a){
                        internalEnergy = internalEnergy - 0.5 * rhoa_loc[a]/rho_loc * v_sqr(vela_loc[a][0], vela_loc[a][1], vela_loc[a][2]);
                    }     
                    // double internalEnergy=rhoe_loc/rho_loc - 0.5*v_sqr(vel_loc[0], vel_loc[1], vel_loc[2]); 
                    double T_loc = calculate_temp(internalEnergy, rho_loc, rhoa_loc);
                    #else
                    double T_loc = T_out;
                    #endif
                    // calculate_feq_geq(fa_loc, g_loc, rho_loc, rhoa_loc, vel_loc, vela_loc, T_out);
                    calculate_feq_geq(fa_loc, g_loc, rho_loc, rhoa_loc, vel_loc, vela_loc, T_loc);

                    for (int l=0; l < npop; ++l){
                        i_nb = i - cx[l];
                        j_nb = j - cy[l];
                        k_nb = k - cz[l];

                        if (mixture[i_nb][j_nb][k_nb].type==TYPE_O || mixture[i_nb][j_nb][k_nb].type==TYPE_O_C){
                            for(size_t a = 0; a < nSpecies; ++a)
                                species[a][i][j][k].f[l] = 2*fa_out[a][l] - fa_loc[a][l];
                                // species[a][i][j][k].f[l] = fa_out[a][l];

                            #ifndef ISOTHERM
                            mixture[i][j][k].g[l] = 2*g_out[l] - g_loc[l];
                            // mixture[i][j][k].g[l] = g_out[l];
                            #endif
                        }
                        else{
                            for(size_t a = 0; a < nSpecies; ++a)
                                species[a][i][j][k].f[l] = fa_out[a][l] + species[a][i][j][k].f[l] - fa_loc[a][l];
                                // species[a][i][j][k].f[l] = fa_out[a][l];

                            #ifndef ISOTHERM
                            mixture[i][j][k].g[l] = g_out[l] + mixture[i][j][k].g[l] - g_loc[l];
                            // mixture[i][j][k].g[l] = g_out[l];
                            #endif
                        }
                    }     

                }
            }
        }
    }
    #else
    #ifdef PARALLEL 
    #pragma omp parallel for schedule(dynamic) 
    #endif
    for(int i=0; i<Nx; ++i){
        for(int j=0; j<Ny; ++j){
            for(int k = 0; k<Nz; ++k){
                if(mixture[i][j][k].type==TYPE_F){
                    int i_nb, j_nb, k_nb;
                    double f_out[npop];
                    double g_out[npop];
                    double f_loc[npop];
                    double g_loc[npop];
                    bool interface_nodes[npop] = {false}; 
                    int n_d = 0;

                    // check interface node                         
                    for (int l=0; l < npop; ++l){
                        i_nb = i - cx[l];
                        j_nb = j - cy[l];
                        k_nb = k - cz[l];

                        if(mixture[i_nb][j_nb][k_nb].type==TYPE_O && mixture[int (i+cx[l])][int (j+cy[l])][int (k+cz[l])].type==TYPE_F){    
                            interface_nodes[l] = true;
                            n_d = n_d + 1;
                            break;
                        }                            
                    }

                    // check if there is the node is the interface with boundary
                    if (n_d == 0)
                        continue;

                    double vel_out[3] = {0.0};
                    double T_out = 0.0;
                    double rho_out = 0.0;
                    double p_out = 0.0;
                    for (int l=0; l < npop; ++l){
                        if (interface_nodes[l] == true){
                            vel_out[0]  = mixture[(int)(i+cx[l])][(int)(j+cy[l])][(int)(k+cz[l])].u;
                            vel_out[1]  = mixture[(int)(i+cx[l])][(int)(j+cy[l])][(int)(k+cz[l])].v;
                            vel_out[2]  = mixture[(int)(i+cx[l])][(int)(j+cy[l])][(int)(k+cz[l])].w;
                            T_out       = mixture[(int)(i+cx[l])][(int)(j+cy[l])][(int)(k+cz[l])].temp;
                            rho_out     = mixture[(int)(i+cx[l])][(int)(j+cy[l])][(int)(k+cz[l])].rho;
                            p_out = mixture[(int)(i-cx[l])][(int)(j-cy[l])][(int)(k-cz[l])].p;
                            break;
                        }
                    }

                    rho_out = p_out / (gas_const * T_out);

                    // std::cout << i << " | " << vel_in[0] << " | " << T_in << std::endl;
                    calculate_feq_geq(f_out, g_out, rho_out, vel_out, T_out);
                    
                    // step 2: calculate f_loc, g_loc, and fa_loc (local distribution function)
                    double rho_loc = 0.0;
                    double rhou_loc = 0.0;
                    double rhov_loc = 0.0;
                    double rhow_loc = 0.0;
                    double rhoe_loc = 0.0;

                    for (int l=0; l < npop; ++l){
                        i_nb = i - cx[l];
                        j_nb = j - cy[l];
                        k_nb = k - cz[l];

                        if (mixture[i_nb][j_nb][k_nb].type==TYPE_O){
                            rho_loc += f_out[l];
                            rhou_loc += f_out[l]*cx[l];
                            rhov_loc += f_out[l]*cy[l];
                            rhow_loc += f_out[l]*cz[l];
                            rhoe_loc += g_out[l];
                        }
                        else{
                            rho_loc += mixture[i][j][k].f[l];
                            rhou_loc += mixture[i][j][k].f[l]*cx[l];
                            rhov_loc += mixture[i][j][k].f[l]*cy[l];
                            rhow_loc += mixture[i][j][k].f[l]*cz[l];
                            rhoe_loc += mixture[i][j][k].g[l];
                        }
                    }

                    double vel_loc[3];
                    vel_loc[0] = rhou_loc / rho_loc;
                    vel_loc[1] = rhov_loc / rho_loc;
                    vel_loc[2] = rhow_loc / rho_loc;

                    double internalEnergy=rhoe_loc/rho_loc - 0.5*v_sqr(vel_loc[0], vel_loc[1], vel_loc[2]);
                    
                    double cv = gas_const / (gamma - 1.0);
                    double T_loc = internalEnergy / cv;                        

                    calculate_feq_geq(f_loc, g_loc, rho_loc, vel_loc, T_loc);

                    for (int l=0; l < npop; ++l){
                        i_nb = i - cx[l];
                        j_nb = j - cy[l];
                        k_nb = k - cz[l];

                        if (mixture[i_nb][j_nb][k_nb].type==TYPE_O){
                            mixture[i][j][k].f[l] = 2*f_out[l] - f_loc[l];

                            #ifndef ISOTHERM
                            mixture[i][j][k].g[l] = 2*g_out[l] - g_loc[l];
                            #endif
                        }
                        else{
                            mixture[i][j][k].f[l] = f_out[l] + mixture[i][j][k].f[l] - f_loc[l];

                            #ifndef ISOTHERM
                            mixture[i][j][k].g[l] = g_out[l] + mixture[i][j][k].g[l] - g_loc[l];
                            #endif
                        }
                    }                  

                } 
            }
        }
    }
    #endif

    // Inflow =================================================================================================================================================
    #if defined MULTICOMP
        #ifdef PARALLEL 
        #pragma omp parallel for schedule(dynamic) 
        #endif
        for(int i=0; i<Nx; ++i){
            for(int j=0; j<Ny; ++j){
                for(int k = 0; k<Nz; ++k){
                    if(mixture[i][j][k].type==TYPE_F){
                        int i_nb, j_nb, k_nb;
                        double fa_in[nSpecies][npop];
                        double g_in[npop];
                        double fa_tgt[nSpecies][npop];
                        double g_tgt[npop];
                        double fa_loc[nSpecies][npop];
                        double g_loc[npop];
                        int l_interface = 999;

                        // check interface node                         
                        for (int l=0; l < npop; ++l){
                            i_nb = i - cx[l];
                            j_nb = j - cy[l];
                            k_nb = k - cz[l];

                            if( (mixture[i_nb][j_nb][k_nb].type==TYPE_I || mixture[i_nb][j_nb][k_nb].type==TYPE_I_C) && mixture[int (i+cx[l])][int (j+cy[l])][int (k+cz[l])].type==TYPE_F){    
                                l_interface = l;
                                break;
                            }                            
                        }

                        // check if there is the node is the interface with boundary
                        if (l_interface == 999 || !onlyOne_1(cx[l_interface], cy[l_interface], cz[l_interface]))
                            continue;

                        double vel_in[3] = {0.0};
                        double T_in = 0.0;
                        double rho_in = 0.0;
                        double rhoa_in[nSpecies] = {0.0};
                        double p_in = 0.0;

                        int i_bdr = i-cx[l_interface];
                        int j_bdr = j-cy[l_interface];
                        int k_bdr = k-cz[l_interface];

                        if (mixture[i_bdr][j_bdr][k_bdr].type == TYPE_I){
                            calculate_moment_point((int)(i+cx[l_interface]), (int)(j+cy[l_interface]), (int)(k+cz[l_interface]));

                            vel_in[0] = mixture[(int)(i-cx[l_interface])][(int)(j-cy[l_interface])][(int)(k-cz[l_interface])].u;
                            vel_in[1] = mixture[(int)(i-cx[l_interface])][(int)(j-cy[l_interface])][(int)(k-cz[l_interface])].v;
                            vel_in[2] = mixture[(int)(i-cx[l_interface])][(int)(j-cy[l_interface])][(int)(k-cz[l_interface])].w;
                            T_in = mixture[(int)(i-cx[l_interface])][(int)(j-cy[l_interface])][(int)(k-cz[l_interface])].temp;
                            rho_in = mixture[(int)(i-cx[l_interface])][(int)(j-cy[l_interface])][(int)(k-cz[l_interface])].rho;
                            for(size_t a = 0; a < nSpecies; ++a)
                                rhoa_in[a] = species[a][(int)(i-cx[l_interface])][(int)(j-cy[l_interface])][(int)(k-cz[l_interface])].rho;
                            
                            p_in = mixture[(int)(i+cx[l_interface])][(int)(j+cy[l_interface])][(int)(k+cz[l_interface])].p;
                        }
                        else if (mixture[i_bdr][j_bdr][k_bdr].type == TYPE_I_C){
                            double sigma = 5.0;

                            int rank = omp_get_thread_num();
                            auto gas = sols[rank]->thermo();   
                            std::vector <double> Y (gas->nSpecies());
                            for(size_t a = 0; a < nSpecies; ++a) Y[gas->speciesIndex(speciesName[a])] = species[a][i][j][k].rho / mixture[i][j][k].rho;
                            gas->setMassFractions(&Y[0]);
                            gas->setState_DP(units.si_rho(mixture[i][j][k].rho), units.si_p(mixture[i][j][k].p));

                            int i_1 = i + cx[l_interface];
                            int j_1 = j + cy[l_interface];
                            int k_1 = k + cz[l_interface];
                            int i_2 = i + 2*cx[l_interface];
                            int j_2 = j + 2*cy[l_interface];
                            int k_2 = k + 2*cz[l_interface];
                            double spd_sound = units.u(gas->soundSpeed());

                            if (cx[l_interface] != 0.0){
                                // double drho_dx  = fd_uw(mixture[i][j][k].rho   , mixture[i_1][j_1][k_1].rho, mixture[i_2][j_2][k_2].rho, dx, cx[l_interface]) ;
                                double du_dx    = fd_uw(mixture[i][j][k].u     , mixture[i_1][j_1][k_1].u  , mixture[i_2][j_2][k_2].u  , dx, cx[l_interface]) ;
                                // double dv_dx    = fd_uw(mixture[i][j][k].v     , mixture[i_1][j_1][k_1].v  , mixture[i_2][j_2][k_2].v  , dx, cx[l_interface]) ;
                                // double dw_dx    = fd_uw(mixture[i][j][k].w     , mixture[i_1][j_1][k_1].w  , mixture[i_2][j_2][k_2].w  , dx, cx[l_interface]) ;
                                double dp_dx    = fd_uw(mixture[i][j][k].p     , mixture[i_1][j_1][k_1].p  , mixture[i_2][j_2][k_2].p  , dx, cx[l_interface]) ;
                                double dYa_dx[nSpecies] = {};
                                for(size_t a = 0; a < nSpecies; ++a){
                                    dYa_dx[a] = fd_uw(species[a][i][j][k].rho/mixture[i][j][k].rho, species[a][i_1][j_1][k_1].rho/mixture[i_1][j_1][k_1].rho, species[a][i_2][j_2][k_2].rho/mixture[i_2][j_2][k_2].rho, dx, cx[l_interface]);                            
                                }

                                double lambda_1 = mixture[i][j][k].u - spd_sound;
                                // double lambda_2 = mixture[i][j][k].u;
                                // double lambda_3 = mixture[i][j][k].u;
                                // double lambda_4 = mixture[i][j][k].u;
                                double lambda_5 = mixture[i][j][k].u + spd_sound;

                                double L1 = lambda_1*(dp_dx - mixture[i][j][k].rho*spd_sound*du_dx);
                                double L2 = sigma*spd_sound/Nx* mixture[i][j][k].rho*units.cp(Cantera::GasConstant/gas->meanMolecularWeight())*(mixture[i][j][k].temp-mixture[i_bdr][j_bdr][k_bdr].temp);
                                double L3 = sigma*spd_sound/Nx* (mixture[i][j][k].v-mixture[i_bdr][j_bdr][k_bdr].v);
                                double L4 = sigma*spd_sound/Nx* (mixture[i][j][k].w-mixture[i_bdr][j_bdr][k_bdr].w);
                                double L5 = lambda_5*(dp_dx + mixture[i][j][k].rho*spd_sound*du_dx);

                                double relax_K = sigma*spd_sound*spd_sound/Nx *mixture[i][j][k].rho*(1.0 - (mixture[i][j][k].u/spd_sound)*(mixture[i][j][k].u/spd_sound));

                                if(cx[l_interface] > 0) // left boundary
                                    L5 = relax_K * ( mixture[i][j][k].u -  mixture[i_bdr][j_bdr][k_bdr].u);
                                else // right boundary
                                    L1 = relax_K * ( mixture[i][j][k].u -  mixture[i_bdr][j_bdr][k_bdr].u);

                                double d1 = 1.0/(spd_sound*spd_sound) * (L2+0.5*(L5+L1));
                                double d2 = 1.0/(2.0*mixture[i][j][k].rho*spd_sound) * (L5-L1);
                                double d3 = L3;
                                double d4 = L4;
                                double d5 = 0.5*(L5+L1);

                                rho_in    = mixture[i][j][k].rho - dt_sim*d1;
                                vel_in[0] = mixture[i][j][k].u - dt_sim*d2;
                                vel_in[1] = mixture[i][j][k].v - dt_sim*d3;
                                vel_in[2] = mixture[i][j][k].w - dt_sim*d4;
                                p_in      = mixture[i][j][k].p - dt_sim*d5;
                                for(size_t a = 0; a < nSpecies; ++a){
                                    rhoa_in[a] = (species[a][i][j][k].rho/mixture[i][j][k].rho - dt_sim*mixture[i][j][k].u*dYa_dx[a]) * rho_in;
                                }

                            }
                            else if (cy[l_interface] != 0.0){
                                // double drho_dx  = fd_uw(mixture[i][j][k].rho   , mixture[i_1][j_1][k_1].rho, mixture[i_2][j_2][k_2].rho, dy, cy[l_interface]) ;
                                // double du_dx    = fd_uw(mixture[i][j][k].u     , mixture[i_1][j_1][k_1].u  , mixture[i_2][j_2][k_2].u  , dy, cy[l_interface]) ;
                                double dv_dx    = fd_uw(mixture[i][j][k].v     , mixture[i_1][j_1][k_1].v  , mixture[i_2][j_2][k_2].v  , dy, cy[l_interface]) ;
                                // double dw_dx    = fd_uw(mixture[i][j][k].w     , mixture[i_1][j_1][k_1].w  , mixture[i_2][j_2][k_2].w  , dy, cy[l_interface]) ;
                                double dp_dx    = fd_uw(mixture[i][j][k].p     , mixture[i_1][j_1][k_1].p  , mixture[i_2][j_2][k_2].p  , dy, cy[l_interface]) ;
                                double dYa_dx[nSpecies] = {};
                                for(size_t a = 0; a < nSpecies; ++a){
                                    dYa_dx[a] = fd_uw(species[a][i][j][k].rho/mixture[i][j][k].rho, species[a][i_1][j_1][k_1].rho/mixture[i_1][j_1][k_1].rho, species[a][i_2][j_2][k_2].rho/mixture[i_2][j_2][k_2].rho, dy, cy[l_interface]);
                                }

                                double lambda_1 = mixture[i][j][k].v - spd_sound;
                                // double lambda_2 = mixture[i][j][k].v;
                                // double lambda_3 = mixture[i][j][k].v;
                                // double lambda_4 = mixture[i][j][k].v;
                                double lambda_5 = mixture[i][j][k].v + spd_sound;

                                double L1 = lambda_1*(dp_dx - mixture[i][j][k].rho*spd_sound*dv_dx);
                                double L2 = sigma*spd_sound/Ny* mixture[i][j][k].rho*units.cp(Cantera::GasConstant/gas->meanMolecularWeight())*(mixture[i][j][k].temp-mixture[i_bdr][j_bdr][k_bdr].temp);
                                double L3 = sigma*spd_sound/Ny* (mixture[i][j][k].w-mixture[i_bdr][j_bdr][k_bdr].w);
                                double L4 = sigma*spd_sound/Ny* (mixture[i][j][k].u-mixture[i_bdr][j_bdr][k_bdr].u);
                                double L5 = lambda_5*(dp_dx + mixture[i][j][k].rho*spd_sound*dv_dx);

                                double relax_K = sigma*spd_sound*spd_sound/Ny *mixture[i][j][k].rho*(1.0 - (mixture[i][j][k].v/spd_sound)*(mixture[i][j][k].v/spd_sound));

                                if(cy[l_interface] > 0) // left boundary
                                    L5 = relax_K * ( mixture[i][j][k].v -  mixture[i_bdr][j_bdr][k_bdr].v);
                                else // right boundary
                                    L1 = relax_K * ( mixture[i][j][k].v -  mixture[i_bdr][j_bdr][k_bdr].v);

                                double d1 = 1.0/(spd_sound*spd_sound) * (L2+0.5*(L5+L1));
                                double d2 = 1.0/(2.0*mixture[i][j][k].rho*spd_sound) * (L5-L1);
                                double d3 = L3;
                                double d4 = L4;
                                double d5 = 0.5*(L5+L1);

                                rho_in    = mixture[i][j][k].rho - dt_sim*d1;
                                vel_in[1] = mixture[i][j][k].v - dt_sim*d2;
                                vel_in[2] = mixture[i][j][k].w - dt_sim*d3;
                                vel_in[0] = mixture[i][j][k].u - dt_sim*d4;
                                p_in      = mixture[i][j][k].p - dt_sim*d5;
                                for(size_t a = 0; a < nSpecies; ++a)
                                    rhoa_in[a] = (species[a][i][j][k].rho/mixture[i][j][k].rho - dt_sim*mixture[i][j][k].v*dYa_dx[a]) * rho_in;
                            }
                            else if (cz[l_interface] != 0.0){
                                // double drho_dx  = fd_uw(mixture[i][j][k].rho   , mixture[i_1][j_1][k_1].rho, mixture[i_2][j_2][k_2].rho, dz, cz[l_interface]) ;
                                // double du_dx    = fd_uw(mixture[i][j][k].u     , mixture[i_1][j_1][k_1].u  , mixture[i_2][j_2][k_2].u  , dz, cz[l_interface]) ;
                                // double dv_dx    = fd_uw(mixture[i][j][k].v     , mixture[i_1][j_1][k_1].v  , mixture[i_2][j_2][k_2].v  , dz, cz[l_interface]) ;
                                double dw_dx    = fd_uw(mixture[i][j][k].w     , mixture[i_1][j_1][k_1].w  , mixture[i_2][j_2][k_2].w  , dz, cz[l_interface]) ;
                                double dp_dx    = fd_uw(mixture[i][j][k].p     , mixture[i_1][j_1][k_1].p  , mixture[i_2][j_2][k_2].p  , dz, cz[l_interface]) ;
                                double dYa_dx[nSpecies] = {};
                                for(size_t a = 0; a < nSpecies; ++a){
                                    dYa_dx[a] = fd_uw(species[a][i][j][k].rho/mixture[i][j][k].rho, species[a][i_1][j_1][k_1].rho/mixture[i_1][j_1][k_1].rho, species[a][i_2][j_2][k_2].rho/mixture[i_2][j_2][k_2].rho, dz, cz[l_interface]);
                                }

                                double lambda_1 = mixture[i][j][k].w - spd_sound;
                                // double lambda_2 = mixture[i][j][k].w;
                                // double lambda_3 = mixture[i][j][k].w;
                                // double lambda_4 = mixture[i][j][k].w;
                                double lambda_5 = mixture[i][j][k].w + spd_sound;

                                double L1 = lambda_1*(dp_dx - mixture[i][j][k].rho*spd_sound*dw_dx);
                                double L2 = sigma*spd_sound/Nz* mixture[i][j][k].rho*units.cp(Cantera::GasConstant/gas->meanMolecularWeight())*(mixture[i][j][k].temp-mixture[i_bdr][j_bdr][k_bdr].temp);
                                double L3 = sigma*spd_sound/Nz* (mixture[i][j][k].u-mixture[i_bdr][j_bdr][k_bdr].u);
                                double L4 = sigma*spd_sound/Nz* (mixture[i][j][k].v-mixture[i_bdr][j_bdr][k_bdr].v);
                                double L5 = lambda_5*(dp_dx + mixture[i][j][k].rho*spd_sound*dw_dx);

                                double relax_K = sigma*spd_sound*spd_sound/Nz *mixture[i][j][k].rho*(1.0 - (mixture[i][j][k].w/spd_sound)*(mixture[i][j][k].w/spd_sound));

                                if(cy[l_interface] > 0) // left boundary
                                    L5 = relax_K * ( mixture[i][j][k].w -  mixture[i_bdr][j_bdr][k_bdr].w);
                                else // right boundary
                                    L1 = relax_K * ( mixture[i][j][k].w -  mixture[i_bdr][j_bdr][k_bdr].w);

                                double d1 = 1.0/(spd_sound*spd_sound) * (L2+0.5*(L5+L1));
                                double d2 = 1.0/(2.0*mixture[i][j][k].rho*spd_sound) * (L5-L1);
                                double d3 = L3;
                                double d4 = L4;
                                double d5 = 0.5*(L5+L1);

                                rho_in    = mixture[i][j][k].rho - dt_sim*d1;
                                vel_in[2] = mixture[i][j][k].w - dt_sim*d2;
                                vel_in[0] = mixture[i][j][k].u - dt_sim*d3;
                                vel_in[1] = mixture[i][j][k].v - dt_sim*d4;
                                p_in      = mixture[i][j][k].p - dt_sim*d5;
                                for(size_t a = 0; a < nSpecies; ++a)
                                    rhoa_in[a] = (species[a][i][j][k].rho/mixture[i][j][k].rho - dt_sim*mixture[i][j][k].w*dYa_dx[a]) * rho_in;
                            }

                            std::fill(Y.begin(), Y.end(), 0.0);
                            for(size_t a = 0; a < nSpecies; ++a) Y[gas->speciesIndex(speciesName[a])] = rhoa_in[a] / rho_in;
                            gas->setMassFractions(&Y[0]);   
                            gas->setState_DP(units.si_rho(rho_in), units.si_p(p_in));
                            T_in = units.temp(gas->temperature());
                        }


                        int rank = omp_get_thread_num();
                        auto gas = sols[rank]->thermo();   
                        std::vector <double> Y (gas->nSpecies());
                        for(size_t a = 0; a < nSpecies; ++a) Y[gas->speciesIndex(speciesName[a])] = rhoa_in[a] / rho_in;
                        gas->setMassFractions(&Y[0]);   
                        gas->setState_TP(units.si_temp(T_in), units.si_p(p_in));

                        rho_in = units.rho(gas->density());
                        for(size_t a = 0; a < nSpecies; ++a) rhoa_in[a] = Y[gas->speciesIndex(speciesName[a])] * rho_in;


                        // std::cout << i << " | " << vel_in[0] << " | " << T_in << std::endl;
                        calculate_feq_geq(fa_in, g_in, rho_in, rhoa_in, vel_in, T_in);


                        // // step 2: calculate f_loc, g_loc, and fa_loc (local distribution function)
                        double rho_loc = 0.0;
                        double vel_loc[3] = {0.0};
                        double rhoe_loc = 0.0;
                        double rhoa_loc[nSpecies] = {0.0};
                        double vela_loc[nSpecies][3] = {0.0};
                        
                        for (int l=0; l < npop; ++l){
                            i_nb = i - cx[l];
                            j_nb = j - cy[l];
                            k_nb = k - cz[l];

                            if (mixture[i_nb][j_nb][k_nb].type==TYPE_I || mixture[i_nb][j_nb][k_nb].type==TYPE_I_C){
                                for (size_t a = 0; a < nSpecies; ++a){
                                    rhoa_loc[a] += fa_in[a][l];
                                    vela_loc[a][0] += fa_in[a][l]*cx[l];
                                    vela_loc[a][1] += fa_in[a][l]*cy[l];
                                    vela_loc[a][2] += fa_in[a][l]*cz[l];
                                }
                                #ifndef ISOTHERM
                                rhoe_loc += g_in[l];
                                #endif
                            }
                            else{
                                for (size_t a = 0; a < nSpecies; ++a){
                                    rhoa_loc[a] += species[a][i][j][k].f[l];
                                    vela_loc[a][0] += species[a][i][j][k].f[l]*cx[l];
                                    vela_loc[a][1] += species[a][i][j][k].f[l]*cy[l];
                                    vela_loc[a][2] += species[a][i][j][k].f[l]*cz[l];
                                }
                                #ifndef ISOTHERM
                                rhoe_loc += mixture[i][j][k].g[l];
                                #endif
                            }
                        }

                        for (size_t a = 0; a < nSpecies; ++a){
                            rho_loc += rhoa_loc[a];
                            vel_loc[0] += vela_loc[a][0];
                            vel_loc[1] += vela_loc[a][1];
                            vel_loc[2] += vela_loc[a][2];

                            if(rhoa_loc[a] > 0.0){
                                vela_loc[a][0] = vela_loc[a][0] / rhoa_loc[a];
                                vela_loc[a][1] = vela_loc[a][1] / rhoa_loc[a];
                                vela_loc[a][2] = vela_loc[a][2] / rhoa_loc[a];
                            }
                            else{
                                vela_loc[a][0] = 0.0;
                                vela_loc[a][1] = 0.0;
                                vela_loc[a][2] = 0.0;
                            }
                        }
                        if(rho_loc > 0){
                            vel_loc[0] = vel_loc[0] / rho_loc;
                            vel_loc[1] = vel_loc[1] / rho_loc;
                            vel_loc[2] = vel_loc[2] / rho_loc;
                        }
                        else{
                            vel_loc[0] = 0.0;
                            vel_loc[1] = 0.0;
                            vel_loc[2] = 0.0;
                        }

                        #ifndef ISOTHERM
                        double internalEnergy=rhoe_loc/rho_loc;// - 0.5*v_sqr(vel_loc[0], vel_loc[1], vel_loc[2]); 
                        for (size_t a = 0; a < nSpecies; ++a){
                            internalEnergy = internalEnergy - 0.5 * rhoa_loc[a]/rho_loc * v_sqr(vela_loc[a][0], vela_loc[a][1], vela_loc[a][2]);
                        }     
                        // double internalEnergy=rhoe_loc/rho_loc - 0.5*v_sqr(vel_loc[0], vel_loc[1], vel_loc[2]); 
                        double T_loc = calculate_temp(internalEnergy, rho_loc, rhoa_loc);
                        #else
                        double T_loc = T_in;
                        #endif

                        calculate_feq_geq(fa_loc, g_loc, rho_loc, rhoa_loc, vel_loc, vela_loc, T_loc);
                        calculate_feq_geq(fa_tgt, g_tgt, rho_loc, rhoa_loc, vel_in, T_in);

                        for (int l=0; l < npop; ++l){
                            i_nb = i - cx[l];
                            j_nb = j - cy[l];
                            k_nb = k - cz[l];

                            if (mixture[i_nb][j_nb][k_nb].type==TYPE_I || mixture[i_nb][j_nb][k_nb].type==TYPE_I_C){
                                for(size_t a = 0; a < nSpecies; ++a)
                                    // species[a][i][j][k].f[l] = 2.0*fa_in[a][l] - fa_loc[a][l];
                                    species[a][i][j][k].f[l] = fa_tgt[a][l] + fa_in[a][l] - fa_loc[a][l];
                                    // species[a][i][j][k].f[l] = fa_in[a][l];

                                #ifndef ISOTHERM
                                // mixture[i][j][k].g[l] = 2.0*g_in[l] - g_loc[l];
                                mixture[i][j][k].g[l] = g_tgt[l] + g_in[l] - g_loc[l];
                                // mixture[i][j][k].g[l] = g_in[l];
                                #endif
                            }
                            else{
                                for(size_t a = 0; a < nSpecies; ++a)
                                    // species[a][i][j][k].f[l] = fa_in[a][l] + species[a][i][j][k].f[l] - fa_loc[a][l];
                                    species[a][i][j][k].f[l] = fa_tgt[a][l] + species[a][i][j][k].f[l] - fa_loc[a][l];
                                    // species[a][i][j][k].f[l] = fa_in[a][l];

                                #ifndef ISOTHERM
                                // mixture[i][j][k].g[l] = g_in[l] + mixture[i][j][k].g[l] - g_loc[l];
                                mixture[i][j][k].g[l] = g_tgt[l] + mixture[i][j][k].g[l] - g_loc[l];
                                // mixture[i][j][k].g[l] = g_in[l];
                                #endif
                            }
                        }                  

                    } 
                }
            }
        }

    #else
    #ifdef PARALLEL 
    #pragma omp parallel for schedule(dynamic) 
    #endif
    for(int i=0; i<Nx; ++i){
        for(int j=0; j<Ny; ++j){
            for(int k = 0; k<Nz; ++k){
                if(mixture[i][j][k].type==TYPE_F){
                    int i_nb, j_nb, k_nb;
                    double f_in[npop];
                    double g_in[npop];
                    double f_tgt[npop];
                    double g_tgt[npop];
                    double f_loc[npop];
                    double g_loc[npop];
                    bool interface_nodes[npop] = {false}; 
                    int n_d = 0;

                    // check interface node                         
                    for (int l=0; l < npop; ++l){
                        i_nb = i - cx[l];
                        j_nb = j - cy[l];
                        k_nb = k - cz[l];

                        if(mixture[i_nb][j_nb][k_nb].type==TYPE_I && mixture[(int)(i+cx[l])][(int)(j+cy[l])][(int)(k+cz[l])].type==TYPE_F){    
                            interface_nodes[l] = true;
                            n_d = n_d + 1;
                            break;
                        }                            
                    }                        

                    // check if there is the node is the interface with boundary
                    if (n_d == 0)
                        continue;

                    double vel_in[3] = {0.0};
                    double T_in = 0.0;
                    double rho_in = 0.0;
                    double p_in = 0.0;
                    for (int l=0; l < npop; ++l){
                        if (interface_nodes[l] == true){
                            vel_in[0] += mixture[(int)(i-cx[l])][(int)(j-cy[l])][(int)(k-cz[l])].u;
                            vel_in[1] += mixture[(int)(i-cx[l])][(int)(j-cy[l])][(int)(k-cz[l])].v;
                            vel_in[2] += mixture[(int)(i-cx[l])][(int)(j-cy[l])][(int)(k-cz[l])].w;
                            T_in += mixture[(int)(i-cx[l])][(int)(j-cy[l])][(int)(k-cz[l])].temp;
                            rho_in += mixture[(int)(i-cx[l])][(int)(j-cy[l])][(int)(k-cz[l])].rho;
                            p_in += mixture[(int)(i+cx[l])][(int)(j+cy[l])][(int)(k+cz[l])].p;
                            break;
                        }
                    }

                    vel_in[0] = 1.0/n_d * vel_in[0];
                    vel_in[1] = 1.0/n_d * vel_in[1];
                    vel_in[2] = 1.0/n_d * vel_in[2];
                    T_in = 1.0/n_d * T_in;
                    p_in = 1.0/n_d * p_in;
                    // rho_in = 1.0/n_d * rho_in;
                    rho_in = p_in / (gas_const * T_in);

                    // std::cout << i << " | " << vel_in[0] << " | " << T_in << std::endl;
                    calculate_feq_geq(f_in, g_in, rho_in, vel_in, T_in);
                    
                    // step 2: calculate f_loc, g_loc, and fa_loc (local distribution function)
                    double rho_loc = 0.0;
                    double rhou_loc = 0.0;
                    double rhov_loc = 0.0;
                    double rhow_loc = 0.0;
                    double rhoe_loc = 0.0;

                    for (int l=0; l < npop; ++l){
                        i_nb = i - cx[l];
                        j_nb = j - cy[l];
                        k_nb = k - cz[l];

                        if (mixture[i_nb][j_nb][k_nb].type==TYPE_I){
                            rho_loc += f_in[l];
                            rhou_loc += f_in[l]*cx[l];
                            rhov_loc += f_in[l]*cy[l];
                            rhow_loc += f_in[l]*cz[l];
                            rhoe_loc += g_in[l];
                        }
                        else{
                            rho_loc += mixture[i][j][k].f[l];
                            rhou_loc += mixture[i][j][k].f[l]*cx[l];
                            rhov_loc += mixture[i][j][k].f[l]*cy[l];
                            rhow_loc += mixture[i][j][k].f[l]*cz[l];
                            rhoe_loc += mixture[i][j][k].g[l];
                        }
                    }

                    double vel_loc[3];
                    vel_loc[0] = rhou_loc / rho_loc;
                    vel_loc[1] = rhov_loc / rho_loc;
                    vel_loc[2] = rhow_loc / rho_loc;

                    double internalEnergy=rhoe_loc/rho_loc - 0.5*v_sqr(vel_loc[0], vel_loc[1], vel_loc[2]);
                    
                    double cv = gas_const / (gamma - 1.0);
                    double T_loc = internalEnergy / cv;                        

                    calculate_feq_geq(f_loc, g_loc, rho_loc, vel_loc, T_loc);
                    calculate_feq_geq(f_tgt, g_tgt, rho_loc, vel_in, T_in);

                    for (int l=0; l < npop; ++l){
                        i_nb = i - cx[l];
                        j_nb = j - cy[l];
                        k_nb = k - cz[l];

                        if (mixture[i_nb][j_nb][k_nb].type==TYPE_I){
                            mixture[i][j][k].f[l] = f_tgt[l] + f_in[l] - f_loc[l];

                            #ifndef ISOTHERM
                            mixture[i][j][k].g[l] = g_tgt[l] + g_in[l] - g_loc[l];
                            #endif
                        }
                        else{
                            mixture[i][j][k].f[l] = f_tgt[l] + mixture[i][j][k].f[l] - f_loc[l];

                            #ifndef ISOTHERM
                            mixture[i][j][k].g[l] = g_tgt[l] + mixture[i][j][k].g[l] - g_loc[l];
                            #endif
                        }

                    }                  

                } 
            }
        }
    }
    #endif

    // Solid =================================================================================================================================================
    #if defined MULTICOMP
     #ifdef PARALLEL 
    #pragma omp parallel for schedule(dynamic) 
    #endif
    for(int i=0; i<Nx; ++i){
        for(int j=0; j<Ny; ++j){
            for(int k = 0; k<Nz; ++k){
                if(mixture[i][j][k].type==TYPE_F){
                    int i_nb, j_nb, k_nb;
                    double fa_tgt[nSpecies][npop];
                    double g_tgt[npop];
                    double fa_loc[nSpecies][npop];
                    double g_loc[npop];
                    bool interface_nodes[npop] = {false}; 
                    int n_d = 0;

                    // check interface node                         
                    for (int l=0; l < npop; ++l){
                        i_nb = i - cx[l];
                        j_nb = j - cy[l];
                        k_nb = k - cz[l];

                        if(mixture[i_nb][j_nb][k_nb].type==TYPE_S && mixture[(int)(i+cx[l])][(int)(j+cy[l])][(int)(k+cz[l])].type==TYPE_F){    
                            interface_nodes[l] = true;
                            n_d = n_d + 1;
                        }                            
                    }

                    // check if there is the node is the interface with boundary
                    if (n_d == 0)
                        continue;

                    // step 1: calculate f_tgt, g_tft, and fa_tgt
                    double rho_bb = 0.0;
                    double rhoa_bb[nSpecies] = {0.0};
                    for (size_t a = 0; a < nSpecies; ++a){
                        for (int l=0; l < npop; ++l)
                            rhoa_bb[a] += species[a][i][j][k].f[l];
                        rho_bb += rhoa_bb[a];
                    }

                    double vel_tgt[3] = {0.0};
                    double vela_tgt[nSpecies][3] = {0.0};
                    double T_tgt = 0.0;
                    for (int l=0; l < npop; ++l){
                        if (interface_nodes[l] == true){
                            double q = 0.5;
                            vel_tgt[0] += (q*mixture[(int)(i+cx[l])][(int)(j+cy[l])][(int)(k+cz[l])].u+mixture[(int)(i-cx[l])][(int)(j-cy[l])][(int)(k-cz[l])].u)/(1.0+q);
                            vel_tgt[1] += (q*mixture[(int)(i+cx[l])][(int)(j+cy[l])][(int)(k+cz[l])].v+mixture[(int)(i-cx[l])][(int)(j-cy[l])][(int)(k-cz[l])].v)/(1.0+q);
                            vel_tgt[2] += (q*mixture[(int)(i+cx[l])][(int)(j+cy[l])][(int)(k+cz[l])].w+mixture[(int)(i-cx[l])][(int)(j-cy[l])][(int)(k-cz[l])].w)/(1.0+q);
                            T_tgt += (q*mixture[(int)(i+cx[l])][(int)(j+cy[l])][(int)(k+cz[l])].temp+mixture[(int)(i-cx[l])][(int)(j-cy[l])][(int)(k-cz[l])].temp)/(1.0+q);
                            for (size_t a = 0; a < nSpecies; ++a){
                                vela_tgt[a][0] += (q*species[a][(int)(i+cx[l])][(int)(j+cy[l])][(int)(k+cz[l])].u+mixture[(int)(i-cx[l])][(int)(j-cy[l])][(int)(k-cz[l])].u)/(1.0+q);
                                vela_tgt[a][1] += (q*species[a][(int)(i+cx[l])][(int)(j+cy[l])][(int)(k+cz[l])].v+mixture[(int)(i-cx[l])][(int)(j-cy[l])][(int)(k-cz[l])].v)/(1.0+q);
                                vela_tgt[a][2] += (q*species[a][(int)(i+cx[l])][(int)(j+cy[l])][(int)(k+cz[l])].w+mixture[(int)(i-cx[l])][(int)(j-cy[l])][(int)(k-cz[l])].w)/(1.0+q);
                            }
                        }
                    }
                    
                    vel_tgt[0] = 1.0/n_d * vel_tgt[0];
                    vel_tgt[1] = 1.0/n_d * vel_tgt[1];
                    vel_tgt[2] = 1.0/n_d * vel_tgt[2];
                    T_tgt = 1.0/n_d * T_tgt;
                    for (size_t a = 0; a < nSpecies; ++a){
                        vela_tgt[a][0] = 1.0/n_d * vela_tgt[a][0];
                        vela_tgt[a][1] = 1.0/n_d * vela_tgt[a][1];
                        vela_tgt[a][2] = 1.0/n_d * vela_tgt[a][2];
                    }

                    // std::cout << vel_tgt[0] << " | " << T_tgt << std::endl;
                    calculate_feq_geq(fa_tgt, g_tgt, rho_bb, rhoa_bb, vel_tgt, vela_tgt, T_tgt);

                    // // step 2: calculate f_loc, g_loc, and fa_loc (local distribution function)
                    double rho_loc = 0.0;
                    double vel_loc[3] = {0.0};
                    double rhoe_loc = 0.0;
                    double rhoa_loc[nSpecies] = {0.0};
                    double vela_loc[nSpecies][3] = {0.0};
                    
                    for (int l=0; l < npop; ++l){
                        i_nb = i - cx[l];
                        j_nb = j - cy[l];
                        k_nb = k - cz[l];

                        if (mixture[i_nb][j_nb][k_nb].type==TYPE_S){
                            for (size_t a = 0; a < nSpecies; ++a){
                                rhoa_loc[a] += fa_tgt[a][l];
                                vela_loc[a][0] += fa_tgt[a][l]*cx[l];
                                vela_loc[a][1] += fa_tgt[a][l]*cy[l];
                                vela_loc[a][2] += fa_tgt[a][l]*cz[l];
                            }
                            rhoe_loc += g_tgt[l];
                        }
                        else{
                            for (size_t a = 0; a < nSpecies; ++a){
                                rhoa_loc[a] += species[a][i][j][k].f[l];
                                vela_loc[a][0] += species[a][i][j][k].f[l]*cx[l];
                                vela_loc[a][1] += species[a][i][j][k].f[l]*cy[l];
                                vela_loc[a][2] += species[a][i][j][k].f[l]*cz[l];
                            }
                            #ifndef ISOTHERM
                            rhoe_loc += mixture[i][j][k].g[l];
                            #endif
                        }
                    }

                    for (size_t a = 0; a < nSpecies; ++a){
                        rho_loc += rhoa_loc[a];
                        vel_loc[0] += vela_loc[a][0];
                        vel_loc[1] += vela_loc[a][1];
                        vel_loc[2] += vela_loc[a][2];

                        if(rhoa_loc[a] > 0){
                            vela_loc[a][0] = vela_loc[a][0] / rhoa_loc[a];
                            vela_loc[a][1] = vela_loc[a][1] / rhoa_loc[a];
                            vela_loc[a][2] = vela_loc[a][2] / rhoa_loc[a];
                        }
                        else {
                            vela_loc[a][0] = 0.0;
                            vela_loc[a][1] = 0.0;
                            vela_loc[a][2] = 0.0;
                        }
                    }
                    if(rho_loc > 0){
                        vel_loc[0] = vel_loc[0] / rho_loc;
                        vel_loc[1] = vel_loc[1] / rho_loc;
                        vel_loc[2] = vel_loc[2] / rho_loc;
                    }
                    else{
                        vel_loc[0] = 0.0; 
                        vel_loc[1] = 0.0;
                        vel_loc[2] = 0.0;
                    }

                    double internalEnergy=rhoe_loc/rho_loc;// - 0.5*v_sqr(vel_loc[0], vel_loc[1], vel_loc[2]); 
                    for (size_t a = 0; a < nSpecies; ++a){
                        internalEnergy = internalEnergy - 0.5 * rhoa_loc[a]/rho_loc * v_sqr(vela_loc[a][0], vela_loc[a][1], vela_loc[a][2]);
                    }                         
                    double T_loc = calculate_temp(internalEnergy, rho_loc, rhoa_loc);

                    calculate_feq_geq(fa_loc, g_loc, rho_loc, rhoa_loc, vel_loc, vela_loc, T_loc);

                    for (int l=0; l < npop; ++l){
                        i_nb = i - cx[l];
                        j_nb = j - cy[l];
                        k_nb = k - cz[l];

                        if (mixture[i_nb][j_nb][k_nb].type==TYPE_S){
                            for(size_t a = 0; a < nSpecies; ++a)
                                species[a][i][j][k].f[l] = 2*fa_tgt[a][l] - fa_loc[a][l];

                            #ifndef ISOTHERM
                            mixture[i][j][k].g[l] = 2*g_tgt[l] - g_loc[l];
                            #endif
                        }
                        else{
                            for(size_t a = 0; a < nSpecies; ++a)
                                species[a][i][j][k].f[l] = fa_tgt[a][l] + species[a][i][j][k].f[l] - fa_loc[a][l];

                            #ifndef ISOTHERM
                            mixture[i][j][k].g[l] = g_tgt[l] + mixture[i][j][k].g[l] - g_loc[l];
                            #endif
                        }
                    }                   

                } 
            }
        }
    }

    #else
    #ifdef PARALLEL 
    #pragma omp parallel for schedule(dynamic) 
    #endif
    for(int i=0; i<Nx; ++i){
        for(int j=0; j<Ny; ++j){
            for(int k = 0; k<Nz; ++k){
                if(mixture[i][j][k].type==TYPE_F){
                    int i_nb, j_nb, k_nb;
                    double f_tgt[npop];
                    double g_tgt[npop];
                    double f_loc[npop];
                    double g_loc[npop];
                    bool interface_nodes[npop] = {false}; 
                    int n_d = 0;

                    // check interface node                         
                    for (int l=0; l < npop; ++l){
                        i_nb = i - cx[l];
                        j_nb = j - cy[l];
                        k_nb = k - cz[l];

                        if(mixture[i_nb][j_nb][k_nb].type==TYPE_S && mixture[(int)(i+cx[l])][(int)(j+cy[l])][(int)(k+cz[l])].type==TYPE_F){    
                            interface_nodes[l] = true;
                            n_d = n_d + 1;
                        }                            
                    }

                    // check if there is the node is the interface with boundary
                    if (n_d == 0)
                        continue;

                    // step 1: calculate f_tgt, g_tft
                    double rho_bb = 0.0;
                    for (int l=0; l < npop; ++l){
                        rho_bb += mixture[i][j][k].f[l];
                    }

                    double vel_tgt[3] = {0.0};
                    double T_tgt = 0.0;
                    for (int l=0; l < npop; ++l){
                        double q = 0.5;
                        if (interface_nodes[l] == true){
                            vel_tgt[0] += (q*mixture[(int)(i+cx[l])][(int)(j+cy[l])][(int)(k+cz[l])].u+mixture[(int)(i-cx[l])][(int)(j-cy[l])][(int)(k-cz[l])].u)/(1.0+q);
                            vel_tgt[1] += (q*mixture[(int)(i+cx[l])][(int)(j+cy[l])][(int)(k+cz[l])].v+mixture[(int)(i-cx[l])][(int)(j-cy[l])][(int)(k-cz[l])].v)/(1.0+q);
                            vel_tgt[2] += (q*mixture[(int)(i+cx[l])][(int)(j+cy[l])][(int)(k+cz[l])].w+mixture[(int)(i-cx[l])][(int)(j-cy[l])][(int)(k-cz[l])].w)/(1.0+q);
                            T_tgt += (q*mixture[(int)(i+cx[l])][(int)(j+cy[l])][(int)(k+cz[l])].temp+mixture[(int)(i-cx[l])][(int)(j-cy[l])][(int)(k-cz[l])].temp)/(1.0+q);
                            // std::cout << j-cy[l] << " | " << mixture[(int)(i-cx[l])][(int)(j-cy[l])][(int)(k-cz[l])].temp << std::endl;
                        }
                    }

                    vel_tgt[0] = 1.0/n_d * vel_tgt[0];
                    vel_tgt[1] = 1.0/n_d * vel_tgt[1];
                    vel_tgt[2] = 1.0/n_d * vel_tgt[2];
                    T_tgt = 1.0/n_d * T_tgt;

                    // std::cout << vel_tgt[0] << " | " << T_tgt << std::endl;
                    calculate_feq_geq(f_tgt, g_tgt, rho_bb, vel_tgt, T_tgt);

                    // step 2: calculate f_loc, g_loc, and fa_loc (local distribution function)
                    double rho_loc = 0.0;
                    double rhou_loc = 0.0;
                    double rhov_loc = 0.0;
                    double rhow_loc = 0.0;
                    double rhoe_loc = 0.0;

                    for (int l=0; l < npop; ++l){
                        i_nb = i - cx[l];
                        j_nb = j - cy[l];
                        k_nb = k - cz[l];

                        if (mixture[i_nb][j_nb][k_nb].type==TYPE_S){
                            rho_loc += f_tgt[l];
                            rhou_loc += f_tgt[l]*cx[l];
                            rhov_loc += f_tgt[l]*cy[l];
                            rhow_loc += f_tgt[l]*cz[l];
                            rhoe_loc += g_tgt[l];
                        }
                        else{
                            rho_loc += mixture[i][j][k].f[l];
                            rhou_loc += mixture[i][j][k].f[l]*cx[l];
                            rhov_loc += mixture[i][j][k].f[l]*cy[l];
                            rhow_loc += mixture[i][j][k].f[l]*cz[l];
                            rhoe_loc += mixture[i][j][k].g[l];
                        }
                    }

                    double vel_loc[3];
                    vel_loc[0] = rhou_loc / rho_loc;
                    vel_loc[1] = rhov_loc / rho_loc;
                    vel_loc[2] = rhow_loc / rho_loc;

                    double internalEnergy=rhoe_loc/rho_loc - 0.5*v_sqr(vel_loc[0], vel_loc[1], vel_loc[2]);
                    
                    double cv = gas_const / (gamma - 1.0);
                    double T_loc = internalEnergy / cv;                        

                    calculate_feq_geq(f_loc, g_loc, rho_loc, vel_loc, T_loc);

                    for (int l=0; l < npop; ++l){
                        i_nb = i - cx[l];
                        j_nb = j - cy[l];
                        k_nb = k - cz[l];

                        if (mixture[i_nb][j_nb][k_nb].type==TYPE_S){
                            mixture[i][j][k].f[l] = 2*f_tgt[l] - f_loc[l];
                            mixture[i][j][k].g[l] = 2*g_tgt[l] - g_loc[l];
                        }
                        else{
                            mixture[i][j][k].f[l] = f_tgt[l] + mixture[i][j][k].f[l] - f_loc[l];
                            mixture[i][j][k].g[l] = g_tgt[l] + mixture[i][j][k].g[l] - g_loc[l];
                        }
                    }                   

                } 
            }
        }
    }

        #ifdef CONJUGATE
        #ifdef PARALLEL 
        #pragma omp parallel for schedule(dynamic) 
        #endif
        for(int i=0; i<Nx; ++i){
            for(int j=0; j<Ny; ++j){
                for(int k = 0; k<Nz; ++k){
                    if(mixture[i][j][k].type==TYPE_S){
                        int i_nb, j_nb, k_nb;
                        double f_tgt[npop];
                        double g_tgt[npop];
                        double f_loc[npop];
                        double g_loc[npop];
                        bool interface_nodes[npop] = {false}; 
                        int n_d = 0;

                        // check interface node                         
                        for (int l=0; l < npop; ++l){
                            i_nb = i - cx[l];
                            j_nb = j - cy[l];
                            k_nb = k - cz[l];

                            if((mixture[i_nb][j_nb][k_nb].type==TYPE_F || mixture[i_nb][j_nb][k_nb].type==TYPE_Q ) && mixture[(int)(i+cx[l])][(int)(j+cy[l])][(int)(k+cz[l])].type==TYPE_S){    //|| mixture[i_nb][j_nb][k_nb].type==TYPE_Q
                                interface_nodes[l] = true;
                                n_d = n_d + 1;
                            }                            
                        }

                        // check if there is the node is the interface with boundary
                        if (n_d == 0)
                            continue;

                        // step 1: calculate f_tgt, g_tft
                        double rho_bb = 0.0;
                        for (int l=0; l < npop; ++l){
                            rho_bb += mixture[i][j][k].f[l];
                        }

                        double vel_tgt[3] = {0.0};
                        double T_tgt = 0.0;
                        for (int l=0; l < npop; ++l){
                            double q = 0.5;
                            if (interface_nodes[l] == true){
                                vel_tgt[0] += 0.0;//(q*mixture[(int)(i+cx[l])][(int)(j+cy[l])][(int)(k+cz[l])].u+mixture[(int)(i-cx[l])][(int)(j-cy[l])][(int)(k-cz[l])].u)/(1.0+q);
                                vel_tgt[1] += 0.0;//(q*mixture[(int)(i+cx[l])][(int)(j+cy[l])][(int)(k+cz[l])].v+mixture[(int)(i-cx[l])][(int)(j-cy[l])][(int)(k-cz[l])].v)/(1.0+q);
                                vel_tgt[2] += 0.0;//(q*mixture[(int)(i+cx[l])][(int)(j+cy[l])][(int)(k+cz[l])].w+mixture[(int)(i-cx[l])][(int)(j-cy[l])][(int)(k-cz[l])].w)/(1.0+q);
                                T_tgt += (q*mixture[(int)(i+cx[l])][(int)(j+cy[l])][(int)(k+cz[l])].temp+mixture[(int)(i-cx[l])][(int)(j-cy[l])][(int)(k-cz[l])].temp)/(1.0+q);
                                // std::cout << j-cy[l] << " | " << mixture[(int)(i-cx[l])][(int)(j-cy[l])][(int)(k-cz[l])].temp << std::endl;
                            }
                        }

                        vel_tgt[0] = 1.0/n_d * vel_tgt[0];
                        vel_tgt[1] = 1.0/n_d * vel_tgt[1];
                        vel_tgt[2] = 1.0/n_d * vel_tgt[2];
                        T_tgt = 1.0/n_d * T_tgt;

                        // std::cout << vel_tgt[0] << " | " << T_tgt << std::endl;
                        calculate_feq_geq(f_tgt, g_tgt, rho_bb, vel_tgt, T_tgt);

                        // step 2: calculate f_loc, g_loc, and fa_loc (local distribution function)
                        double rho_loc = 0.0;
                        double rhou_loc = 0.0;
                        double rhov_loc = 0.0;
                        double rhow_loc = 0.0;
                        double rhoe_loc = 0.0;

                        for (int l=0; l < npop; ++l){
                            i_nb = i - cx[l];
                            j_nb = j - cy[l];
                            k_nb = k - cz[l];

                            if (mixture[i_nb][j_nb][k_nb].type==TYPE_F || mixture[i_nb][j_nb][k_nb].type==TYPE_Q){ //|| mixture[i_nb][j_nb][k_nb].type==TYPE_Q
                                rho_loc += f_tgt[l];
                                rhou_loc += f_tgt[l]*cx[l];
                                rhov_loc += f_tgt[l]*cy[l];
                                rhow_loc += f_tgt[l]*cz[l];
                                rhoe_loc += g_tgt[l];
                            }
                            else{
                                rho_loc += mixture[i][j][k].f[l];
                                rhou_loc += mixture[i][j][k].f[l]*cx[l];
                                rhov_loc += mixture[i][j][k].f[l]*cy[l];
                                rhow_loc += mixture[i][j][k].f[l]*cz[l];
                                rhoe_loc += mixture[i][j][k].g[l];
                            }
                        }

                        double vel_loc[3];
                        vel_loc[0] = rhou_loc / rho_loc;
                        vel_loc[1] = rhov_loc / rho_loc;
                        vel_loc[2] = rhow_loc / rho_loc;

                        double internalEnergy=rhoe_loc/rho_loc - 0.5*v_sqr(vel_loc[0], vel_loc[1], vel_loc[2]);
                        
                        double cv = gas_const / (gamma - 1.0);
                        double T_loc = internalEnergy / cv;                        

                        calculate_feq_geq(f_loc, g_loc, rho_loc, vel_loc, T_loc);

                        for (int l=0; l < npop; ++l){
                            i_nb = i - cx[l];
                            j_nb = j - cy[l];
                            k_nb = k - cz[l];

                            if (mixture[i_nb][j_nb][k_nb].type==TYPE_F || mixture[i_nb][j_nb][k_nb].type==TYPE_Q){ //|| mixture[i_nb][j_nb][k_nb].type==TYPE_Q
                                // mixture[i][j][k].f[l] = 2*f_tgt[l] - f_loc[l];
                                mixture[i][j][k].g[l] = 2*g_tgt[l] - g_loc[l];
                            }
                            else{
                                // mixture[i][j][k].f[l] = f_tgt[l] + mixture[i][j][k].f[l] - f_loc[l];
                                mixture[i][j][k].g[l] = g_tgt[l] + mixture[i][j][k].g[l] - g_loc[l];
                            }
                        }                   

                    } 
                }
            }
        }
        #endif
    #endif
}

void LBM::dirSlip(int l, int i, int j, int k, int &lp, int &ip, int &jp, int &kp)
{
    int i_nb, j_nb, k_nb;
    i_nb = i - cx[l];
    j_nb = j - cy[l];
    k_nb = k - cz[l];

    bool i_chk = false;
    bool j_chk = false;
    bool k_chk = false;
    if (mixture[i_nb][j][k].type==TYPE_FS)
        i_chk = true;
    if (mixture[i][j_nb][k].type==TYPE_FS)
        j_chk = true;
    if (mixture[i][j][k_nb].type==TYPE_FS)
        k_chk = true;

    double cxp = cx[l];
    double cyp = cy[l];
    double czp = cz[l];
    ip = i_nb;
    jp = j_nb;
    kp = k_nb;
    if (i_chk == true){
        cxp = -1.0*cx[l];
        ip = i;
    }
    else if (j_chk == true){
        cyp = -1.0*cy[l];
        jp = j;
    }
    else if (k_chk == true){
        czp = -1.0*cz[l];
        kp = k;
    } 

    ip = ((ip - 1 + (Nx-2)) % (Nx-2)) + 1;
    jp = ((jp - 1 + (Ny-2)) % (Ny-2)) + 1;
    kp = ((kp - 1 + (Nz-2)) % (Nz-2)) + 1;
    
    for(size_t a = 0; a < npop; ++a){
        if (cxp == cx[a] && cyp == cy[a] && czp == cz[a]){
            lp = a;
            break;
        }
    }

}


//==================================================================================================================================================
