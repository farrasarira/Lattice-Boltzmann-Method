// #include "lbm.hpp"
// #include "math_util.hpp"
// #include "FD.hpp"
// #include "units.hpp"
// #include <omp.h>

// void LBM::fill_BC()
// {
//     #ifdef PARALLEL 
//     #pragma omp parallel for schedule(static, 1) 
//     #endif
//     for(int i = 0; i < Nx; ++i){
//         for(int j = 0; j < Ny; ++j){
//             for(int k = 0; k < Nz; ++k){

//                 // Inlet --------------------------------------------------------------------------------------------------------------------
//                 if(mixture[i][j][k].type == TYPE_I_C){
//                     int l_fluid = 999;

//                     for(int l = 0; l < npop; ++l){               
//                         if( i+cx[l] < 0 || i+cx[l] >= Nx || j+cy[l] < 0 || j+cy[l] >= Ny || k+cz[l] < 0 || k+cz[l] >= Nz)     
//                             continue;

//                         if(mixture[(int)(i+cx[l])][(int)(j+cy[l])][(int)(k+cz[l])].type==TYPE_F){   
//                             l_fluid = l;
//                             break;
//                         }
//                     }
                    
//                     if(l_fluid == 999 || !onlyOne_1(cx[l_fluid], cy[l_fluid], cz[l_fluid]))
//                         continue;
                    
//                     if (cx[l_fluid] != 0) {
//                         double sign_side = -1*cx[l_fluid];
//                         double spd_sound = get_soundspeed(mixture[i][j][k].temp);
//                         double lambda_1 = mixture[i][j][k].u + sign_side*spd_sound;

//                         int i_1 = i + cx[l_fluid];
//                         int j_1 = j + cy[l_fluid];
//                         int k_1 = k + cz[l_fluid];
//                         int i_2 = i + 2*cx[l_fluid];
//                         int j_2 = j + 2*cy[l_fluid];
//                         int k_2 = k + 2*cz[l_fluid];

//                         double dp_dxi = fd_uw(mixture[i][j][k].p, mixture[i_1][j_1][k_1].p, mixture[i_2][j_2][k_2].p, dx, sign_side);
//                         double du_dxi = fd_uw(mixture[i][j][k].u, mixture[i_1][j_1][k_1].u, mixture[i_2][j_2][k_2].u, dx, sign_side);

//                         double rho = mixture[i][j][k].rho;
//                         double p = mixture[i][j][k].p;
//                         double u = mixture[i][j][k].u;
//                         double v = mixture[i][j][k].v;
//                         double w = mixture[i][j][k].w;
//                         double temp = mixture[i][j][k].temp;

                        
//                         double L1 = lambda_1 * (dp_dxi - rho*spd_sound*du_dxi);
//                         double L5 = L1;
//                         double L2 = 0.5*(gamma-1.0)*(L5+L1);

//                         double d1 = 1.0/(spd_sound*spd_sound) * (L2 + 0.5*(L5+L1));
//                         double d2 = 0.5*(L5+L1);

//                         mixture[i][j][k].rho = rho + dt_sim*(-d1); 
//                         mixture[i][j][k].temp = temp ;
//                         mixture[i][j][k].u = u;
//                         mixture[i][j][k].v = v;
//                         mixture[i][j][k].w = w;
//                         mixture[i][j][k].p = p + dt_sim*(-d2);
//                     }
//                     else if (cy[l_fluid] != 0){
//                         double sign_side = -1*cy[l_fluid];
//                         double spd_sound = get_soundspeed(mixture[i][j][k].temp);
//                         double lambda_1 = mixture[i][j][k].v + sign_side*spd_sound;

//                         int i_1 = i + cx[l_fluid];
//                         int j_1 = j + cy[l_fluid];
//                         int k_1 = k + cz[l_fluid];
//                         int i_2 = i + 2*cx[l_fluid];
//                         int j_2 = j + 2*cy[l_fluid];
//                         int k_2 = k + 2*cz[l_fluid];

//                         double dp_dxi = fd_uw(mixture[i][j][k].p, mixture[i_1][j_1][k_1].p, mixture[i_2][j_2][k_2].p, dx, sign_side);
//                         double dv_dxi = fd_uw(mixture[i][j][k].v, mixture[i_1][j_1][k_1].v, mixture[i_2][j_2][k_2].v, dx, sign_side);

//                         double rho = mixture[i][j][k].rho;
//                         double p = mixture[i][j][k].p;
//                         double u = mixture[i][j][k].u;
//                         double v = mixture[i][j][k].v;
//                         double w = mixture[i][j][k].w;
//                         double temp = mixture[i][j][k].temp;
                        
//                         double L1 = lambda_1 * (dp_dxi - rho*spd_sound*dv_dxi);
//                         double L5 = L1;
//                         double L2 = 0.5*(gamma-1.0)*(L5+L1);

//                         double d1 = 1.0/(spd_sound*spd_sound) * (L2 + 0.5*(L5+L1));
//                         double d2 = 0.5*(L5+L1);

//                         mixture[i][j][k].rho = rho + dt_sim*(-d1); 
//                         mixture[i][j][k].temp = temp ;
//                         mixture[i][j][k].u = u;
//                         mixture[i][j][k].v = v;
//                         mixture[i][j][k].w = w;
//                         mixture[i][j][k].p = p + dt_sim*(-d2);
//                     }
//                     else{
//                         double sign_side = -1*cz[l_fluid];
//                         double spd_sound = get_soundspeed(mixture[i][j][k].temp);
//                         double lambda_1 = mixture[i][j][k].w - sign_side*spd_sound;

//                         int i_1 = i + cx[l_fluid];
//                         int j_1 = j + cy[l_fluid];
//                         int k_1 = k + cz[l_fluid];
//                         int i_2 = i + 2*cx[l_fluid];
//                         int j_2 = j + 2*cy[l_fluid];
//                         int k_2 = k + 2*cz[l_fluid];

//                         double dp_dxi = fd_uw(mixture[i][j][k].p, mixture[i_1][j_1][k_1].p, mixture[i_2][j_2][k_2].p, dx, sign_side);
//                         double dw_dxi = fd_uw(mixture[i][j][k].w, mixture[i_1][j_1][k_1].w, mixture[i_2][j_2][k_2].w, dx, sign_side);

//                         double rho = mixture[i][j][k].rho;
//                         double p = mixture[i][j][k].p;
//                         double u = mixture[i][j][k].u;
//                         double v = mixture[i][j][k].v;
//                         double w = mixture[i][j][k].w;
//                         double temp = mixture[i][j][k].temp;
                        
//                         double L1 = lambda_1 * (dp_dxi - rho*spd_sound*dw_dxi);
//                         double L5 = L1;
//                         double L2 = 0.5*(gamma-1.0)*(L5+L1);

//                         double d1 = 1.0/(spd_sound*spd_sound) * (L2 + 0.5*(L5+L1));
//                         double d2 = 0.5*(L5+L1);

//                         mixture[i][j][k].rho = rho + dt_sim*(-d1); 
//                         mixture[i][j][k].temp = temp ;
//                         mixture[i][j][k].u = u;
//                         mixture[i][j][k].v = v;
//                         mixture[i][j][k].w = w;
//                         mixture[i][j][k].p = p + dt_sim*(-d2);
//                     }
//                 }
                
//                 // Outlet -------------------------------------------------------------------------------------------------------
//                 if(mixture[i][j][k].type == TYPE_O_C){
//                     double sigma = 0.7;
//                     double Mach_ref = 1E-1;
//                     double L_ref = 50;

//                     int l_fluid = 999;

//                     for(int l = 0; l < npop; ++l){               
//                         if( i+cx[l] < 0 || i+cx[l] >= Nx || j+cy[l] < 0 || j+cy[l] >= Ny || k+cz[l] < 0 || k+cz[l] >= Nz)     
//                             continue;

//                         if(mixture[(int)(i+cx[l])][(int)(j+cy[l])][(int)(k+cz[l])].type==TYPE_F){   
//                             l_fluid = l;
//                             break;
//                         }
//                     }
                    
//                     if(l_fluid == 999 || !onlyOne_1(cx[l_fluid], cy[l_fluid], cz[l_fluid]))
//                         continue;
                    
//                     if (cx[l_fluid] != 0) {
                        
//                         int i_0 = i + cx[l_fluid];
//                         int j_0 = j + cy[l_fluid];
//                         int k_0 = k + cz[l_fluid];
//                         int i_1 = i + 2*cx[l_fluid];
//                         int j_1 = j + 2*cy[l_fluid];
//                         int k_1 = k + 2*cz[l_fluid];
//                         int i_2 = i + 3*cx[l_fluid];
//                         int j_2 = j + 3*cy[l_fluid];
//                         int k_2 = k + 3*cz[l_fluid];

//                         double sign_side = -1*cx[l_fluid];
//                         double spd_sound = get_soundspeed(mixture[i_0][j_0][k_0].temp);
//                         // double lambda_1 = mixture[i][j][k].u - spd_sound;
//                         double lambda_2 = mixture[i_0][j_0][k_0].u;
//                         double lambda_3 = mixture[i_0][j_0][k_0].u;
//                         double lambda_4 = mixture[i_0][j_0][k_0].u;
//                         double lambda_5 = mixture[i_0][j_0][k_0].u + sign_side*spd_sound;

//                         double dp_dxi   = fd_uw(mixture[i_0][j_0][k_0].p, mixture[i_1][j_1][k_1].p, mixture[i_2][j_2][k_2].p, dx, sign_side);
//                         double drho_dxi = fd_uw(mixture[i_0][j_0][k_0].rho, mixture[i_1][j_1][k_1].rho, mixture[i_2][j_2][k_2].rho, dx, sign_side);
//                         double du_dxi   = fd_uw(mixture[i_0][j_0][k_0].u, mixture[i_1][j_1][k_1].u, mixture[i_2][j_2][k_2].u, dx, sign_side);
//                         double dv_dxi   = fd_uw(mixture[i_0][j_0][k_0].v, mixture[i_1][j_1][k_1].v, mixture[i_2][j_2][k_2].v, dx, sign_side);
//                         double dw_dxi   = fd_uw(mixture[i_0][j_0][k_0].w, mixture[i_1][j_1][k_1].w, mixture[i_2][j_2][k_2].w, dx, sign_side);

//                         double K    = sigma*(1.0-Mach_ref*Mach_ref)*spd_sound/L_ref;
//                         double rho  = mixture[i_0][j_0][k_0].rho;
//                         double p    = mixture[i_0][j_0][k_0].p;
//                         double u    = mixture[i_0][j_0][k_0].u;
//                         double v    = mixture[i_0][j_0][k_0].v;
//                         double w    = mixture[i_0][j_0][k_0].w;
//                         double temp = mixture[i_0][j_0][k_0].temp;

                        
//                         double L1 = K*(mixture[i_0][j_0][k_0].p - mixture[i][j][k].p);
//                         double L2 = lambda_2 * (spd_sound*spd_sound*drho_dxi - dp_dxi);
//                         double L3 = lambda_3 * dv_dxi;
//                         double L4 = lambda_4 * dw_dxi;
//                         double L5 = lambda_5 * (dp_dxi + rho*spd_sound*du_dxi);

//                         double d1 = 1.0/(spd_sound*spd_sound) * (L2 + 0.5*(L5+L1));
//                         double d2 = temp/(rho*spd_sound*spd_sound) * (-L2 + 0.5*(gamma-1.0)*(L5+L1));
//                         double d3 = 1.0/(2.0*rho*spd_sound) * (L5-L1);
//                         double d4 = L3;
//                         double d5 = L4;

//                         mixture[i][j][k].rho = rho - dt_sim*(d1); 
//                         mixture[i][j][k].temp = temp - dt_sim*(d2);
//                         mixture[i][j][k].u = u - dt_sim*(d3);
//                         mixture[i][j][k].v = v - dt_sim*(d4);
//                         mixture[i][j][k].w = w - dt_sim*(d5);
//                         mixture[i][j][k].p = p;
//                     }
//                     else if (cy[l_fluid] != 0){
                        
//                         int i_0 = i + cx[l_fluid];
//                         int j_0 = j + cy[l_fluid];
//                         int k_0 = k + cz[l_fluid];
//                         int i_1 = i + 2*cx[l_fluid];
//                         int j_1 = j + 2*cy[l_fluid];
//                         int k_1 = k + 2*cz[l_fluid];
//                         int i_2 = i + 3*cx[l_fluid];
//                         int j_2 = j + 3*cy[l_fluid];
//                         int k_2 = k + 3*cz[l_fluid];

//                         double sign_side = -1*cy[l_fluid];
//                         double spd_sound = get_soundspeed(mixture[i_0][j_0][k_0].temp);
//                         // double lambda_1 = mixture[i][j][k].v - spd_sound;
//                         double lambda_2 = mixture[i_0][j_0][k_0].v;
//                         double lambda_3 = mixture[i_0][j_0][k_0].v;
//                         double lambda_4 = mixture[i_0][j_0][k_0].v;
//                         double lambda_5 = mixture[i_0][j_0][k_0].v + sign_side*spd_sound;

//                         double dp_dxi   = fd_uw(mixture[i_0][j_0][k_0].p    , mixture[i_1][j_1][k_1].p  , mixture[i_2][j_2][k_2].p  , dy, sign_side);
//                         double drho_dxi = fd_uw(mixture[i_0][j_0][k_0].rho  , mixture[i_1][j_1][k_1].rho, mixture[i_2][j_2][k_2].rho, dy, sign_side);
//                         double du_dxi   = fd_uw(mixture[i_0][j_0][k_0].u    , mixture[i_1][j_1][k_1].u  , mixture[i_2][j_2][k_2].u  , dy, sign_side);
//                         double dv_dxi   = fd_uw(mixture[i_0][j_0][k_0].v    , mixture[i_1][j_1][k_1].v  , mixture[i_2][j_2][k_2].v  , dy, sign_side);
//                         double dw_dxi   = fd_uw(mixture[i_0][j_0][k_0].w    , mixture[i_1][j_1][k_1].w  , mixture[i_2][j_2][k_2].w  , dy, sign_side);

//                         double K    = sigma*(1.0-Mach_ref*Mach_ref)*spd_sound/L_ref;
//                         double rho  = mixture[i_0][j_0][k_0].rho;
//                         double p    = mixture[i_0][j_0][k_0].p;
//                         double u    = mixture[i_0][j_0][k_0].u;
//                         double v    = mixture[i_0][j_0][k_0].v;
//                         double w    = mixture[i_0][j_0][k_0].w;
//                         double temp = mixture[i_0][j_0][k_0].temp;

//                         double L1 = K*(mixture[i_0][j_0][k_0].p - mixture[i][j][k].p);
//                         double L2 = lambda_2 * (spd_sound*spd_sound*drho_dxi - dp_dxi);
//                         double L3 = lambda_3 * du_dxi;
//                         double L4 = lambda_4 * dw_dxi;
//                         double L5 = lambda_5 * (dp_dxi + rho*spd_sound*dv_dxi);

//                         double d1 = 1.0/(spd_sound*spd_sound) * (L2 + 0.5*(L5+L1));
//                         double d2 = temp/(rho*spd_sound*spd_sound) * (-L2 + 0.5*(gamma-1.0)*(L5+L1));
//                         double d3 = 1.0/(2.0*rho*spd_sound) * (L5-L1);
//                         double d4 = L3;
//                         double d5 = L4;

//                         mixture[i][j][k].rho = rho + dt_sim*(-d1); 
//                         mixture[i][j][k].temp = temp + dt_sim*(-d2);
//                         mixture[i][j][k].u = u + dt_sim*(-d4);
//                         mixture[i][j][k].v = v + dt_sim*(-d3);
//                         mixture[i][j][k].w = w + dt_sim*(-d5);
//                         mixture[i][j][k].p = p;

//                         // std::cout << mixture[i][j][k].v << std::endl;
//                     }
//                     else{
//                         double sign_side = -1*cz[l_fluid];
//                         double spd_sound = get_soundspeed(mixture[i][j][k].temp);
//                         // double lambda_1 = mixture[i][j][k].w - spd_sound;
//                         double lambda_2 = mixture[i][j][k].w;
//                         double lambda_3 = mixture[i][j][k].w;
//                         double lambda_4 = mixture[i][j][k].w;
//                         double lambda_5 = mixture[i][j][k].w + sign_side*spd_sound;

//                         int i_1 = i + cx[l_fluid];
//                         int j_1 = j + cy[l_fluid];
//                         int k_1 = k + cz[l_fluid];
//                         int i_2 = i + 2*cx[l_fluid];
//                         int j_2 = j + 2*cy[l_fluid];
//                         int k_2 = k + 2*cz[l_fluid];

//                         double dp_dxi = fd_uw(mixture[i][j][k].p, mixture[i_1][j_1][k_1].p, mixture[i_2][j_2][k_2].p, dz, sign_side);
//                         double drho_dxi = fd_uw(mixture[i][j][k].rho, mixture[i_1][j_1][k_1].rho, mixture[i_2][j_2][k_2].rho, dz, sign_side);
//                         double du_dxi = fd_uw(mixture[i][j][k].u, mixture[i_1][j_1][k_1].u, mixture[i_2][j_2][k_2].u, dz, sign_side);
//                         double dv_dxi = fd_uw(mixture[i][j][k].v, mixture[i_1][j_1][k_1].v, mixture[i_2][j_2][k_2].v, dz, sign_side);
//                         double dw_dxi = fd_uw(mixture[i][j][k].w, mixture[i_1][j_1][k_1].w, mixture[i_2][j_2][k_2].w, dz, sign_side);

//                         double K = sigma*(1.0-Mach_ref*Mach_ref)*spd_sound/L_ref;
//                         double rho = mixture[i][j][k].rho;
//                         double p = mixture[i][j][k].p;
//                         double u = mixture[i][j][k].u;
//                         double v = mixture[i][j][k].v;
//                         double w = mixture[i][j][k].w;
//                         double temp = mixture[i][j][k].temp;

                        
//                         double L1 = K*(mixture[i_1][j_1][k_1].p - mixture[i][j][k].p);
//                         double L2 = lambda_2 * (spd_sound*spd_sound*drho_dxi - dp_dxi);
//                         double L3 = lambda_3 * du_dxi;
//                         double L4 = lambda_4 * dv_dxi;
//                         double L5 = lambda_5 * (dp_dxi + rho*spd_sound*dw_dxi);

//                         double d1 = 1.0/(spd_sound*spd_sound) * (L2 + 0.5*(L5+L1));
//                         double d2 = temp/(rho*spd_sound*spd_sound) * (-L2 + 0.5*(gamma-1.0)*(L5+L1));
//                         double d3 = 1.0/(2.0*rho*spd_sound) * (L5-L1);
//                         double d4 = L3;
//                         double d5 = L4;

//                         mixture[i][j][k].rho = rho + dt_sim*(-d1); 
//                         mixture[i][j][k].temp = temp + dt_sim*(-d2);
//                         mixture[i][j][k].u = u + dt_sim*(-d4);
//                         mixture[i][j][k].v = v + dt_sim*(-d5);
//                         mixture[i][j][k].w = w + dt_sim*(-d3);
//                         mixture[i][j][k].p = p;
//                     }


//                 }
//             }
//         }
//     }
// }

// void LBM::TMS_BC()
// {
//     // TMS Boundary Conditions ---------------------------------------------------------------------------------------------
    
//     // SOLID ===============================================================================================================
//     #ifndef MULTICOMP
//     #ifdef PARALLEL 
//     #pragma omp parallel for schedule(static, 1) 
//     #endif
//     for(int i=0; i<Nx; ++i){
//         for(int j=0; j<Ny; ++j){
//             for(int k = 0; k<Nz; ++k){
//                 if(mixture[i][j][k].type==TYPE_F){
//                     int i_nb, j_nb, k_nb;
//                     double f_tgt[npop];
//                     double g_tgt[npop];
//                     double f_loc[npop];
//                     double g_loc[npop];
//                     bool interface_nodes[npop] = {false}; 
//                     int n_d = 0;

//                     // check interface node                         
//                     for (int l=0; l < npop; ++l){
//                         i_nb = i - cx[l];
//                         j_nb = j - cy[l];
//                         k_nb = k - cz[l];

//                         if(mixture[i_nb][j_nb][k_nb].type==TYPE_S && mixture[(int)(i+cx[l])][(int)(j+cy[l])][(int)(k+cz[l])].type==TYPE_F){    
//                             interface_nodes[l] = true;
//                             n_d = n_d + 1;
//                         }                            
//                     }

//                     // check if there is the node is the interface with boundary
//                     if (n_d == 0)
//                         continue;

//                     // step 1: calculate f_tgt, g_tft
//                     double rho_bb = 0.0;
//                     for (int l=0; l < npop; ++l){
//                         rho_bb += mixture[i][j][k].f[l];
//                     }

//                     double vel_tgt[3] = {0.0};
//                     double T_tgt = 0.0;
//                     for (int l=0; l < npop; ++l){
//                         double q = 0.5;
//                         if (interface_nodes[l] == true){
//                             vel_tgt[0] += (q*mixture[(int)(i+cx[l])][(int)(j+cy[l])][(int)(k+cz[l])].u+mixture[(int)(i-cx[l])][(int)(j-cy[l])][(int)(k-cz[l])].u)/(1.0+q);
//                             vel_tgt[1] += (q*mixture[(int)(i+cx[l])][(int)(j+cy[l])][(int)(k+cz[l])].v+mixture[(int)(i-cx[l])][(int)(j-cy[l])][(int)(k-cz[l])].v)/(1.0+q);
//                             vel_tgt[2] += (q*mixture[(int)(i+cx[l])][(int)(j+cy[l])][(int)(k+cz[l])].w+mixture[(int)(i-cx[l])][(int)(j-cy[l])][(int)(k-cz[l])].w)/(1.0+q);
//                             T_tgt += (q*mixture[(int)(i+cx[l])][(int)(j+cy[l])][(int)(k+cz[l])].temp+mixture[(int)(i-cx[l])][(int)(j-cy[l])][(int)(k-cz[l])].temp)/(1.0+q);
//                             // std::cout << j-cy[l] << " | " << mixture[(int)(i-cx[l])][(int)(j-cy[l])][(int)(k-cz[l])].temp << std::endl;
//                         }
//                     }

//                     vel_tgt[0] = 1.0/n_d * vel_tgt[0];
//                     vel_tgt[1] = 1.0/n_d * vel_tgt[1];
//                     vel_tgt[2] = 1.0/n_d * vel_tgt[2];
//                     T_tgt = 1.0/n_d * T_tgt;

//                     // std::cout << vel_tgt[0] << " | " << T_tgt << std::endl;
//                     calculate_feq_geq(f_tgt, g_tgt, rho_bb, vel_tgt, T_tgt);

//                     // step 2: calculate f_loc, g_loc, and fa_loc (local distribution function)
//                     double rho_loc = 0.0;
//                     double rhou_loc = 0.0;
//                     double rhov_loc = 0.0;
//                     double rhow_loc = 0.0;
//                     double rhoe_loc = 0.0;

//                     for (int l=0; l < npop; ++l){
//                         i_nb = i - cx[l];
//                         j_nb = j - cy[l];
//                         k_nb = k - cz[l];

//                         if (mixture[i_nb][j_nb][k_nb].type==TYPE_S){
//                             rho_loc += f_tgt[l];
//                             rhou_loc += f_tgt[l]*cx[l];
//                             rhov_loc += f_tgt[l]*cy[l];
//                             rhow_loc += f_tgt[l]*cz[l];
//                             #ifndef ISOTHERM
//                             rhoe_loc += g_tgt[l];
//                             #endif
//                         }
//                         else{
//                             rho_loc += mixture[i][j][k].f[l];
//                             rhou_loc += mixture[i][j][k].f[l]*cx[l];
//                             rhov_loc += mixture[i][j][k].f[l]*cy[l];
//                             rhow_loc += mixture[i][j][k].f[l]*cz[l];
//                             #ifndef ISOTHERM
//                             rhoe_loc += mixture[i][j][k].g[l];
//                             #endif
//                         }
//                     }

//                     double vel_loc[3];
//                     vel_loc[0] = rhou_loc / rho_loc;
//                     vel_loc[1] = rhov_loc / rho_loc;
//                     vel_loc[2] = rhow_loc / rho_loc;

//                     double internalEnergy=rhoe_loc/rho_loc - 0.5*v_sqr(vel_loc[0], vel_loc[1], vel_loc[2]);
                    
//                     double cv = gas_const / (gamma - 1.0);
//                     double T_loc = internalEnergy / cv;                        

//                     calculate_feq_geq(f_loc, g_loc, rho_loc, vel_loc, T_loc);

//                     for (int l=0; l < npop; ++l){
//                         i_nb = i - cx[l];
//                         j_nb = j - cy[l];
//                         k_nb = k - cz[l];

//                         if (mixture[i_nb][j_nb][k_nb].type==TYPE_S){
//                             mixture[i][j][k].f[l] = 2*f_tgt[l] - f_loc[l];
//                             #ifndef ISOTHERM
//                             mixture[i][j][k].g[l] = 2*g_tgt[l] - g_loc[l];
//                             #endif
//                         }
//                         else{
//                             mixture[i][j][k].f[l] = f_tgt[l] + mixture[i][j][k].f[l] - f_loc[l];
//                             #ifndef ISOTHERM
//                             mixture[i][j][k].g[l] = g_tgt[l] + mixture[i][j][k].g[l] - g_loc[l];
//                             #endif
//                         }
//                     }                   

//                 } 
//             }
//         }
//     }

//     #elif defined MULTICOMP

//     #ifdef PARALLEL 
//     #pragma omp parallel for schedule(static, 1) 
//     #endif
//     for(int i=0; i<Nx; ++i){
//         for(int j=0; j<Ny; ++j){
//             for(int k = 0; k<Nz; ++k){
//                 if(mixture[i][j][k].type==TYPE_F){
//                     int i_nb, j_nb, k_nb;
//                     double fa_tgt[nSpecies][npop];
//                     double g_tgt[npop];
//                     double fa_loc[nSpecies][npop];
//                     double g_loc[npop];
//                     bool interface_nodes[npop] = {false}; 
//                     int n_d = 0;

//                     // check interface node                         
//                     for (int l=0; l < npop; ++l){
//                         i_nb = i - cx[l];
//                         j_nb = j - cy[l];
//                         k_nb = k - cz[l];

//                         if(mixture[i_nb][j_nb][k_nb].type==TYPE_S && mixture[(int)(i+cx[l])][(int)(j+cy[l])][(int)(k+cz[l])].type==TYPE_F){    
//                             interface_nodes[l] = true;
//                             n_d = n_d + 1;
//                         }                            
//                     }

//                     // check if there is the node is the interface with boundary
//                     if (n_d == 0)
//                         continue;

//                     // step 1: calculate f_tgt, g_tft, and fa_tgt
//                     double rho_bb = 0.0;
//                     double rhoa_bb[nSpecies] = {0.0};
//                     for (size_t a = 0; a < nSpecies; ++a){
//                         for (int l=0; l < npop; ++l)
//                             rhoa_bb[a] += species[a][i][j][k].f[l];
//                         rho_bb += rhoa_bb[a];
//                     }

//                     double vel_tgt[3] = {0.0};
//                     double vela_tgt[nSpecies][3] = {0.0};
//                     double T_tgt = 0.0;
//                     for (int l=0; l < npop; ++l){
//                         if (interface_nodes[l] == true){
//                             double q = 0.5;
//                             vel_tgt[0] += (q*mixture[(int)(i+cx[l])][(int)(j+cy[l])][(int)(k+cz[l])].u+mixture[(int)(i-cx[l])][(int)(j-cy[l])][(int)(k-cz[l])].u)/(1.0+q);
//                             vel_tgt[1] += (q*mixture[(int)(i+cx[l])][(int)(j+cy[l])][(int)(k+cz[l])].v+mixture[(int)(i-cx[l])][(int)(j-cy[l])][(int)(k-cz[l])].v)/(1.0+q);
//                             vel_tgt[2] += (q*mixture[(int)(i+cx[l])][(int)(j+cy[l])][(int)(k+cz[l])].w+mixture[(int)(i-cx[l])][(int)(j-cy[l])][(int)(k-cz[l])].w)/(1.0+q);
//                             T_tgt += (q*mixture[(int)(i+cx[l])][(int)(j+cy[l])][(int)(k+cz[l])].temp+mixture[(int)(i-cx[l])][(int)(j-cy[l])][(int)(k-cz[l])].temp)/(1.0+q);
//                             for (size_t a = 0; a < nSpecies; ++a){
//                                 vela_tgt[a][0] += (q*species[a][(int)(i+cx[l])][(int)(j+cy[l])][(int)(k+cz[l])].u+mixture[(int)(i-cx[l])][(int)(j-cy[l])][(int)(k-cz[l])].u)/(1.0+q);
//                                 vela_tgt[a][1] += (q*species[a][(int)(i+cx[l])][(int)(j+cy[l])][(int)(k+cz[l])].v+mixture[(int)(i-cx[l])][(int)(j-cy[l])][(int)(k-cz[l])].v)/(1.0+q);
//                                 vela_tgt[a][2] += (q*species[a][(int)(i+cx[l])][(int)(j+cy[l])][(int)(k+cz[l])].w+mixture[(int)(i-cx[l])][(int)(j-cy[l])][(int)(k-cz[l])].w)/(1.0+q);
//                             }
//                         }
//                     }
                    
//                     vel_tgt[0] = 1.0/n_d * vel_tgt[0];
//                     vel_tgt[1] = 1.0/n_d * vel_tgt[1];
//                     vel_tgt[2] = 1.0/n_d * vel_tgt[2];
//                     T_tgt = 1.0/n_d * T_tgt;
//                     for (size_t a = 0; a < nSpecies; ++a){
//                         vela_tgt[a][0] = 1.0/n_d * vela_tgt[a][0];
//                         vela_tgt[a][1] = 1.0/n_d * vela_tgt[a][1];
//                         vela_tgt[a][2] = 1.0/n_d * vela_tgt[a][2];
//                     }

//                     // std::cout << vel_tgt[0] << " | " << T_tgt << std::endl;
//                     calculate_feq_geq(fa_tgt, g_tgt, rho_bb, rhoa_bb, vel_tgt, vela_tgt, T_tgt);

//                     // // step 2: calculate f_loc, g_loc, and fa_loc (local distribution function)
//                     double rho_loc = 0.0;
//                     double vel_loc[3] = {0.0};
//                     double rhoe_loc = 0.0;
//                     double rhoa_loc[nSpecies] = {0.0};
//                     double vela_loc[nSpecies][3] = {0.0};
                    
//                     for (int l=0; l < npop; ++l){
//                         i_nb = i - cx[l];
//                         j_nb = j - cy[l];
//                         k_nb = k - cz[l];

//                         if (mixture[i_nb][j_nb][k_nb].type==TYPE_S){
//                             for (size_t a = 0; a < nSpecies; ++a){
//                                 rhoa_loc[a] += fa_tgt[a][l];
//                                 vela_loc[a][0] += fa_tgt[a][l]*cx[l];
//                                 vela_loc[a][1] += fa_tgt[a][l]*cy[l];
//                                 vela_loc[a][2] += fa_tgt[a][l]*cz[l];
//                             }
//                             rhoe_loc += g_tgt[l];
//                         }
//                         else{
//                             for (size_t a = 0; a < nSpecies; ++a){
//                                 rhoa_loc[a] += species[a][i][j][k].f[l];
//                                 vela_loc[a][0] += species[a][i][j][k].f[l]*cx[l];
//                                 vela_loc[a][1] += species[a][i][j][k].f[l]*cy[l];
//                                 vela_loc[a][2] += species[a][i][j][k].f[l]*cz[l];
//                             }
//                             #ifndef ISOTHERM
//                             rhoe_loc += mixture[i][j][k].g[l];
//                             #endif
//                         }
//                     }

//                     for (size_t a = 0; a < nSpecies; ++a){
//                         rho_loc += rhoa_loc[a];
//                         vel_loc[0] += vela_loc[a][0];
//                         vel_loc[1] += vela_loc[a][1];
//                         vel_loc[2] += vela_loc[a][2];

//                         vela_loc[a][0] = vela_loc[a][0] / rhoa_loc[a];
//                         vela_loc[a][1] = vela_loc[a][1] / rhoa_loc[a];
//                         vela_loc[a][2] = vela_loc[a][2] / rhoa_loc[a];
//                     }
//                     vel_loc[0] = vel_loc[0] / rho_loc;
//                     vel_loc[1] = vel_loc[1] / rho_loc;
//                     vel_loc[2] = vel_loc[2] / rho_loc;

//                     double internalEnergy=rhoe_loc/rho_loc - 0.5*v_sqr(vel_loc[0], vel_loc[1], vel_loc[2]);                     
//                     double T_loc = calculate_temp(internalEnergy, rho_loc, rhoa_loc);

//                     calculate_feq_geq(fa_loc, g_loc, rho_loc, rhoa_loc, vel_loc, vela_loc, T_loc);

//                     for (int l=0; l < npop; ++l){
//                         i_nb = i - cx[l];
//                         j_nb = j - cy[l];
//                         k_nb = k - cz[l];

//                         if (mixture[i_nb][j_nb][k_nb].type==TYPE_S){
//                             for(size_t a = 0; a < nSpecies; ++a)
//                                 species[a][i][j][k].f[l] = 2*fa_tgt[a][l] - fa_loc[a][l];

//                             #ifndef ISOTHERM
//                             mixture[i][j][k].g[l] = 2*g_tgt[l] - g_loc[l];
//                             #endif
//                         }
//                         else{
//                             for(size_t a = 0; a < nSpecies; ++a)
//                                 species[a][i][j][k].f[l] = fa_tgt[a][l] + species[a][i][j][k].f[l] - fa_loc[a][l];

//                             #ifndef ISOTHERM
//                             mixture[i][j][k].g[l] = g_tgt[l] + mixture[i][j][k].g[l] - g_loc[l];
//                             #endif
//                         }
//                     }                   

//                 } 
//             }
//         }
//     }
//     #endif

//     // INFLOW ============================================================================================================= 
//     #ifndef MULTICOMP
//         #ifdef PARALLEL 
//         #pragma omp parallel for schedule(static, 1) 
//         #endif
//         for(int i=0; i<Nx; ++i){
//             for(int j=0; j<Ny; ++j){
//                 for(int k = 0; k<Nz; ++k){
//                     if(mixture[i][j][k].type==TYPE_F){
//                         int i_nb, j_nb, k_nb;
//                         double f_in[npop];
//                         double g_in[npop];
//                         double f_tgt[npop];
//                         double g_tgt[npop];
//                         double f_loc[npop];
//                         double g_loc[npop];
//                         int l_fluid = 999;

//                         // check interface node                         
//                         for (int l=0; l < npop; ++l){
//                             i_nb = i - cx[l];
//                             j_nb = j - cy[l];
//                             k_nb = k - cz[l];

//                             if((mixture[i_nb][j_nb][k_nb].type==TYPE_I_E || mixture[i_nb][j_nb][k_nb].type==TYPE_I_C)&& mixture[(int)(i+cx[l])][(int)(j+cy[l])][(int)(k+cz[l])].type==TYPE_F){    
//                                 l_fluid = l;
//                                 break;
//                             }                            
//                         }       

//                         if(l_fluid == 999 || !onlyOne_1(cx[l_fluid], cy[l_fluid], cz[l_fluid]))
//                             continue;            
                        
//                         double vel_in[3] = {0.0};
//                         double T_in = 0.0;
//                         double rho_in = 0.0;
//                         double p_in = 0.0;
//                         double sign_side = -1*(cx[l_fluid] + cy[l_fluid] + cz[l_fluid]);

//                         if (mixture[(int)(i-cx[l_fluid])][(int)(j-cy[l_fluid])][(int)(k-cz[l_fluid])].type==TYPE_I_E){
//                             vel_in[0]   = mixture[(int)(i-cx[l_fluid])][(int)(j-cy[l_fluid])][(int)(k-cz[l_fluid])].u;
//                             vel_in[1]   = mixture[(int)(i-cx[l_fluid])][(int)(j-cy[l_fluid])][(int)(k-cz[l_fluid])].v;
//                             vel_in[2]   = mixture[(int)(i-cx[l_fluid])][(int)(j-cy[l_fluid])][(int)(k-cz[l_fluid])].w;
//                             T_in        = mixture[(int)(i-cx[l_fluid])][(int)(j-cy[l_fluid])][(int)(k-cz[l_fluid])].temp;
//                             rho_in      = mixture[(int)(i-cx[l_fluid])][(int)(j-cy[l_fluid])][(int)(k-cz[l_fluid])].rho;
//                             p_in        = mixture[(int)(i+cx[l_fluid])][(int)(j+cy[l_fluid])][(int)(k+cz[l_fluid])].p;
//                             rho_in      = p_in / (gas_const * T_in);
//                         }
//                         else if (mixture[i_nb][j_nb][k_nb].type==TYPE_I_C){
                           
//                             int i_1 = i + cx[l_fluid];
//                             int j_1 = j + cy[l_fluid];
//                             int k_1 = k + cz[l_fluid];
//                             int i_2 = i + 2*cx[l_fluid];
//                             int j_2 = j + 2*cy[l_fluid];
//                             int k_2 = k + 2*cz[l_fluid];
                            
//                             double dp_dx    = fd_uw(mixture[i][j][k].p  , mixture[i_1][j_1][k_1].p  , mixture[i_2][j_2][k_2].p  , dx, sign_side);
//                             // double drho_dx  = fd_uw(mixture[i][j][k].rho, mixture[i_1][j_1][k_1].rho, mixture[i_2][j_2][k_2].rho, dx, sign_side);
//                             double du_dx    = fd_uw(mixture[i][j][k].u  , mixture[i_1][j_1][k_1].u  , mixture[i_2][j_2][k_2].u  , dx, sign_side);
//                             // double dv_dx    = fd_uw(mixture[i][j][k].v  , mixture[i_1][j_1][k_1].v  , mixture[i_2][j_2][k_2].v  , dx, sign_side);
//                             // double dw_dx    = fd_uw(mixture[i][j][k].w  , mixture[i_1][j_1][k_1].w  , mixture[i_2][j_2][k_2].w  , dx, sign_side);

//                             double u    = mixture[i][j][k].u*abs(cx[l_fluid]) + mixture[i][j][k].v*abs(cy[l_fluid]) + mixture[i][j][k].w*abs(cz[l_fluid]);
//                             double rho  = mixture[i][j][k].rho;
//                             double cs   = 1.0*sign_side*get_soundspeed(mixture[i][j][k].temp);

//                             // double L1 = 0.0;//u*(drho_dx - 1.0/(cs*cs) * dp_dx);
//                             // double L2 = 0.0;//u*(abs(cz[l_fluid])*du_dx + abs(cx[l_fluid])*dv_dx + abs(cy[l_fluid])*dw_dx);
//                             // double L3 = 0.0;//u*(abs(cy[l_fluid])*du_dx + abs(cz[l_fluid])*dv_dx + abs(cx[l_fluid])*dw_dx);
//                             // double L4 = 0.0;//(u + cs) * (1.0/(rho*cs) * dp_dx + abs(cx[l_fluid])*du_dx + abs(cy[l_fluid])*dv_dx + abs(cz[l_fluid])*dw_dx);
//                             // double L5 = (u - cs) * (1.0/(rho*cs) * dp_dx - abs(cx[l_fluid])*du_dx - abs(cy[l_fluid])*dv_dx - abs(cz[l_fluid])*dw_dx);

//                             // rho_in     = mixture[i_nb][j_nb][k_nb].rho  - dt_sim*(L1 + rho/(2.0*cs) * (L4+L5));
//                             // p_in       = mixture[i_nb][j_nb][k_nb].p    - dt_sim*(rho*cs/2.0 * (L4+L5));
//                             // vel_in[0]  = mixture[i_nb][j_nb][k_nb].u    - dt_sim*(abs(cx[l_fluid])*0.5*(L4-L5)   + abs(cy[l_fluid])*L3            + abs(cz[l_fluid])*L2            );
//                             // vel_in[1]  = mixture[i_nb][j_nb][k_nb].v    - dt_sim*(abs(cx[l_fluid])*L2            + abs(cy[l_fluid])*0.5*(L4-L5)   + abs(cz[l_fluid])*L3            );
//                             // vel_in[2]  = mixture[i_nb][j_nb][k_nb].w    - dt_sim*(abs(cx[l_fluid])*L3            + abs(cy[l_fluid])*L2            + abs(cz[l_fluid])*0.5*(L4-L5)   );
//                             // T_in       = p_in / (rho_in * gas_const);

//                             // mixture[i_nb][j_nb][k_nb].rho = rho_in;
//                             // mixture[i_nb][j_nb][k_nb].p = p_in;
//                             // mixture[i_nb][j_nb][k_nb].u = vel_in[0];
//                             // mixture[i_nb][j_nb][k_nb].v = vel_in[1]; 
//                             // mixture[i_nb][j_nb][k_nb].w = vel_in[2]; 

//                             double L1 = (u-cs) * (dp_dx - rho*cs*du_dx);
//                             double L5 = L1;//(u+cs) * (dp_dx + rho*cs*du_dx);
//                             double L2 = 0.5*(gamma-1.0)*(L5+L1);//(u   ) * (cs*cs*drho_dx - dp_dx);
//                             double L3 = 0.0;//(u   ) * (abs(cz[l_fluid])*du_dx + abs(cx[l_fluid])*dv_dx + abs(cy[l_fluid])*dw_dx);
//                             double L4 = 0.0;//(u   ) * (abs(cy[l_fluid])*du_dx + abs(cz[l_fluid])*dv_dx + abs(cx[l_fluid])*dw_dx);

//                             rho_in     = mixture[i_nb][j_nb][k_nb].rho  - dt_sim*(1.0/(cs*cs) * (L2+0.5*(L5+L1)));
//                             p_in       = mixture[i_nb][j_nb][k_nb].p    - dt_sim*(0.5*(L5+L1));
//                             vel_in[0]  = mixture[i_nb][j_nb][k_nb].u    - dt_sim*(abs(cx[l_fluid])*0.5/(rho*cs)*(L5-L1) + abs(cy[l_fluid])*L4                   + abs(cz[l_fluid])*L3                   );
//                             vel_in[1]  = mixture[i_nb][j_nb][k_nb].v    - dt_sim*(abs(cx[l_fluid])*L3                   + abs(cy[l_fluid])*0.5/(rho*cs)*(L5-L1) + abs(cz[l_fluid])*L4                   );
//                             vel_in[2]  = mixture[i_nb][j_nb][k_nb].w    - dt_sim*(abs(cx[l_fluid])*L4                   + abs(cy[l_fluid])*L3                   + abs(cz[l_fluid])*0.5/(rho*cs)*(L5-L1) );
//                             // T_in       = mixture[i_nb][j_nb][k_nb].temp - dt_sim*(mixture[i_nb][j_nb][k_nb].temp/(rho*cs*cs) * (-L2+0.5*(gamma-1.0)*(L5+L1))     );
//                             T_in =  p_in / (rho_in * gas_const);

//                             mixture[i_nb][j_nb][k_nb].rho   = rho_in;
//                             mixture[i_nb][j_nb][k_nb].p     = p_in;
//                             mixture[i_nb][j_nb][k_nb].temp  = T_in;
//                             mixture[i_nb][j_nb][k_nb].u     = vel_in[0];
//                             mixture[i_nb][j_nb][k_nb].v     = vel_in[1]; 
//                             mixture[i_nb][j_nb][k_nb].w     = vel_in[2]; 
//                         }


//                         // std::cout << i << " | " << vel_in[0] << " | " << T_in << std::endl;
//                         calculate_feq_geq(f_in, g_in, rho_in, vel_in, T_in);
                        
//                         // step 2: calculate f_loc, g_loc, and fa_loc (local distribution function)
//                         double rho_loc = 0.0;
//                         double rhou_loc = 0.0;
//                         double rhov_loc = 0.0;
//                         double rhow_loc = 0.0;
//                         double rhoe_loc = 0.0;

//                         for (int l=0; l < npop; ++l){
//                             i_nb = i - cx[l];
//                             j_nb = j - cy[l];
//                             k_nb = k - cz[l];

//                             if (mixture[i_nb][j_nb][k_nb].type==TYPE_I_E || mixture[i_nb][j_nb][k_nb].type==TYPE_I_C){
//                                 rho_loc += f_in[l];
//                                 rhou_loc += f_in[l]*cx[l];
//                                 rhov_loc += f_in[l]*cy[l];
//                                 rhow_loc += f_in[l]*cz[l];
//                                 rhoe_loc += g_in[l];
//                             }
//                             else{
//                                 rho_loc += mixture[i][j][k].f[l];
//                                 rhou_loc += mixture[i][j][k].f[l]*cx[l];
//                                 rhov_loc += mixture[i][j][k].f[l]*cy[l];
//                                 rhow_loc += mixture[i][j][k].f[l]*cz[l];
//                                 rhoe_loc += mixture[i][j][k].g[l];
//                             }
//                         }

//                         double vel_loc[3];
//                         vel_loc[0] = rhou_loc / rho_loc;
//                         vel_loc[1] = rhov_loc / rho_loc;
//                         vel_loc[2] = rhow_loc / rho_loc;

//                         double internalEnergy=rhoe_loc/rho_loc - 0.5*v_sqr(vel_loc[0], vel_loc[1], vel_loc[2]);
                        
//                         double cv = gas_const / (gamma - 1.0);
//                         double T_loc = internalEnergy / cv;                        

//                         calculate_feq_geq(f_loc, g_loc, rho_loc, vel_loc, T_loc);
//                         calculate_feq_geq(f_tgt, g_tgt, rho_loc, vel_in, T_in);

//                         for (int l=0; l < npop; ++l){
//                             i_nb = i - cx[l];
//                             j_nb = j - cy[l];
//                             k_nb = k - cz[l];

//                             if (mixture[i_nb][j_nb][k_nb].type==TYPE_I_E || mixture[i_nb][j_nb][k_nb].type==TYPE_I_C){
//                                 mixture[i][j][k].f[l] = f_tgt[l] + f_in[l] - f_loc[l];

//                                 #ifndef ISOTHERM
//                                 mixture[i][j][k].g[l] = g_tgt[l] + g_in[l] - g_loc[l];
//                                 #endif
//                             }
//                             else{
//                                 mixture[i][j][k].f[l] = f_tgt[l] + mixture[i][j][k].f[l] - f_loc[l];

//                                 #ifndef ISOTHERM
//                                 mixture[i][j][k].g[l] = g_tgt[l] + mixture[i][j][k].g[l] - g_loc[l];
//                                 #endif
//                             }

//                         }                  

//                     } 
//                 }
//             }
//         }

//     #elif defined MULTICOMP
//         #ifdef PARALLEL 
//         #pragma omp parallel for schedule(static, 1) 
//         #endif
//         for(int i=0; i<Nx; ++i){
//             for(int j=0; j<Ny; ++j){
//                 for(int k = 0; k<Nz; ++k){
//                     if(mixture[i][j][k].type==TYPE_F){
//                         int i_nb, j_nb, k_nb;
//                         double fa_in[nSpecies][npop];
//                         double g_in[npop];
//                         double fa_tgt[nSpecies][npop];
//                         double g_tgt[npop];
//                         double fa_loc[nSpecies][npop];
//                         double g_loc[npop];
//                         int l_fluid = 999;

//                         // check interface node                         
//                         for (int l=0; l < npop; ++l){
//                             i_nb = i - cx[l];
//                             j_nb = j - cy[l];
//                             k_nb = k - cz[l];

//                             if(mixture[i_nb][j_nb][k_nb].type==TYPE_I_E&& mixture[int (i+cx[l])][int (j+cy[l])][int (k+cz[l])].type==TYPE_F){    
//                                 l_fluid = l;
//                                 break;
//                             }                            
//                         }

//                         // check if there is the node is the interface with boundary
//                         if ( l_fluid == 999 || !onlyOne_1(cx[l_fluid], cy[l_fluid], cz[l_fluid]) )
//                             continue;

//                         double vel_in[3] = {0.0};
//                         double T_in = 0.0;
//                         double rho_in = 0.0;
//                         double rhoa_in[nSpecies] = {0.0};

//                         vel_in[0] = mixture[(int)(i-cx[l_fluid])][(int)(j-cy[l_fluid])][(int)(k-cz[l_fluid])].u;
//                         vel_in[1] = mixture[(int)(i-cx[l_fluid])][(int)(j-cy[l_fluid])][(int)(k-cz[l_fluid])].v;
//                         vel_in[2] = mixture[(int)(i-cx[l_fluid])][(int)(j-cy[l_fluid])][(int)(k-cz[l_fluid])].w;
//                         T_in = mixture[(int)(i-cx[l_fluid])][(int)(j-cy[l_fluid])][(int)(k-cz[l_fluid])].temp;
//                         rho_in = mixture[(int)(i-cx[l_fluid])][(int)(j-cy[l_fluid])][(int)(k-cz[l_fluid])].rho;
//                         for(size_t a = 0; a < nSpecies; ++a)
//                             rhoa_in[a] = species[a][(int)(i-cx[l_fluid])][(int)(j-cy[l_fluid])][(int)(k-cz[l_fluid])].rho;

//                         double pres_field = mixture[(int)(i+cx[l_fluid])][(int)(j+cy[l_fluid])][(int)(k+cz[l_fluid])].p;

//                         int rank = omp_get_thread_num();
//                         auto gas = sols[rank]->thermo();   
//                         std::vector <double> Y (gas->nSpecies());
//                         for(size_t a = 0; a < nSpecies; ++a) Y[gas->speciesIndex(speciesName[a])] = rhoa_in[a] / rho_in;
//                         gas->setState_TP(units.si_temp(T_in), units.si_p(pres_field));
//                         gas->setMassFractions(&Y[0]);  

//                         rho_in = units.rho(gas->density());   
//                         for(size_t a = 0; a < nSpecies; ++a) rhoa_in[a] = Y[gas->speciesIndex(speciesName[a])] * rho_in;

                                                 
//                         // std::cout << i << " | " << vel_in[0] << " | " << T_in << std::endl;
//                         calculate_feq_geq(fa_in, g_in, rho_in, rhoa_in, vel_in, T_in);

//                         // // step 2: calculate f_loc, g_loc, and fa_loc (local distribution function)
//                         double rho_loc = 0.0;
//                         double vel_loc[3] = {0.0};
//                         double rhoe_loc = 0.0;
//                         double rhoa_loc[nSpecies] = {0.0};
//                         double vela_loc[nSpecies][3] = {0.0};
                        
//                         for (int l=0; l < npop; ++l){
//                             i_nb = i - cx[l];
//                             j_nb = j - cy[l];
//                             k_nb = k - cz[l];

//                             if (mixture[i_nb][j_nb][k_nb].type==TYPE_I_E){
//                                 for (size_t a = 0; a < nSpecies; ++a){
//                                     rhoa_loc[a] += fa_in[a][l];
//                                     vela_loc[a][0] += fa_in[a][l]*cx[l];
//                                     vela_loc[a][1] += fa_in[a][l]*cy[l];
//                                     vela_loc[a][2] += fa_in[a][l]*cz[l];
//                                 }
//                                 #ifndef ISOTHERM
//                                 rhoe_loc += g_in[l];
//                                 #endif
//                             }
//                             else{
//                                 for (size_t a = 0; a < nSpecies; ++a){
//                                     rhoa_loc[a] += species[a][i][j][k].f[l];
//                                     vela_loc[a][0] += species[a][i][j][k].f[l]*cx[l];
//                                     vela_loc[a][1] += species[a][i][j][k].f[l]*cy[l];
//                                     vela_loc[a][2] += species[a][i][j][k].f[l]*cz[l];
//                                 }
//                                 #ifndef ISOTHERM
//                                 rhoe_loc += mixture[i][j][k].g[l];
//                                 #endif
//                             }
//                         }

//                         for (size_t a = 0; a < nSpecies; ++a){
//                             rho_loc += rhoa_loc[a];
//                             vel_loc[0] += vela_loc[a][0];
//                             vel_loc[1] += vela_loc[a][1];
//                             vel_loc[2] += vela_loc[a][2];

//                             if(rhoa_loc[a] > 0.0){
//                                 vela_loc[a][0] = vela_loc[a][0] / rhoa_loc[a];
//                                 vela_loc[a][1] = vela_loc[a][1] / rhoa_loc[a];
//                                 vela_loc[a][2] = vela_loc[a][2] / rhoa_loc[a];
//                             }
//                             else{
//                                 vela_loc[a][0] = 0.0;
//                                 vela_loc[a][1] = 0.0;
//                                 vela_loc[a][2] = 0.0;
//                             }
//                         }
//                         if(rho_loc > 0){
//                             vel_loc[0] = vel_loc[0] / rho_loc;
//                             vel_loc[1] = vel_loc[1] / rho_loc;
//                             vel_loc[2] = vel_loc[2] / rho_loc;
//                         }
//                         else{
//                             vel_loc[0] = 0.0;
//                             vel_loc[1] = 0.0;
//                             vel_loc[2] = 0.0;
//                         }

//                         #ifndef ISOTHERM
//                         double internalEnergy=rhoe_loc/rho_loc - 0.5*v_sqr(vel_loc[0], vel_loc[1], vel_loc[2]);                     
//                         double T_loc = calculate_temp(internalEnergy, rho_loc, rhoa_loc);
//                         #else
//                         double T_loc = T_in;
//                         #endif

//                         calculate_feq_geq(fa_loc, g_loc, rho_loc, rhoa_loc, vel_loc, vela_loc, T_loc);
//                         calculate_feq_geq(fa_tgt, g_tgt, rho_loc, rhoa_loc, vel_in, T_in);

//                         for (int l=0; l < npop; ++l){
//                             i_nb = i - cx[l];
//                             j_nb = j - cy[l];
//                             k_nb = k - cz[l];

//                             if (mixture[i_nb][j_nb][k_nb].type==TYPE_I_E){
//                                 for(size_t a = 0; a < nSpecies; ++a)
//                                     species[a][i][j][k].f[l] = fa_tgt[a][l] + fa_in[a][l] - fa_loc[a][l];

//                                 #ifndef ISOTHERM
//                                 mixture[i][j][k].g[l] = g_tgt[l] + g_in[l] - g_loc[l];
//                                 #endif
//                             }
//                             else{
//                                 for(size_t a = 0; a < nSpecies; ++a)
//                                     species[a][i][j][k].f[l] = fa_tgt[a][l] + species[a][i][j][k].f[l] - fa_loc[a][l];

//                                 #ifndef ISOTHERM
//                                 mixture[i][j][k].g[l] = g_tgt[l] + mixture[i][j][k].g[l] - g_loc[l];
//                                 #endif
//                             }
//                         }                  

//                     } 
//                 }
//             }
//         }
//     #endif

//     // OUTFLOW
//     #ifndef MULTICOMP
//         #ifdef PARALLEL 
//         #pragma omp parallel for schedule(static, 1) 
//         #endif
//         for(int i=0; i<Nx; ++i){
//             for(int j=0; j<Ny; ++j){
//                 for(int k = 0; k<Nz; ++k){
//                     if(mixture[i][j][k].type==TYPE_F){
//                         int i_nb, j_nb, k_nb;
//                         double f_out[npop];
//                         double g_out[npop];
//                         double f_loc[npop];
//                         double g_loc[npop];
//                         int l_fluid = 999;

//                         // check interface node                         
//                         for (int l=0; l < npop; ++l){
//                             i_nb = i - cx[l];
//                             j_nb = j - cy[l];
//                             k_nb = k - cz[l];

//                             if((mixture[i_nb][j_nb][k_nb].type==TYPE_O_E||mixture[i_nb][j_nb][k_nb].type==TYPE_O_C) && mixture[int (i+cx[l])][int (j+cy[l])][int (k+cz[l])].type==TYPE_F){    
//                                 l_fluid = l;
//                                 break;
//                             }                            
//                         }

//                         if(l_fluid == 999 || !onlyOne_1(cx[l_fluid], cy[l_fluid], cz[l_fluid]))
//                             continue;

//                         double vel_out[3] = {0.0};
//                         double T_out = 0.0;
//                         double rho_out = 0.0;
//                         double p_out = 0.0;
//                         double sign_side = -1*(cx[l_fluid] + cy[l_fluid] + cz[l_fluid]);

//                         if (mixture[i_nb][j_nb][k_nb].type==TYPE_O_E){
//                             double p_a      = mixture[(int)(i-cx[l_fluid])][(int)(j-cy[l_fluid])][(int)(k-cz[l_fluid])].p;
//                             double p_d      = mixture[(int)(i+cx[l_fluid])][(int)(j+cy[l_fluid])][(int)(k+cz[l_fluid])].p;
//                             double rho_d    = mixture[(int)(i+cx[l_fluid])][(int)(j+cy[l_fluid])][(int)(k+cz[l_fluid])].rho;
//                             double u_d      = mixture[(int)(i+cx[l_fluid])][(int)(j+cy[l_fluid])][(int)(k+cz[l_fluid])].u;
//                             double v_d      = mixture[(int)(i+cx[l_fluid])][(int)(j+cy[l_fluid])][(int)(k+cz[l_fluid])].v;
//                             double w_d      = mixture[(int)(i+cx[l_fluid])][(int)(j+cy[l_fluid])][(int)(k+cz[l_fluid])].w;
//                             double temp_d   = mixture[(int)(i+cx[l_fluid])][(int)(j+cy[l_fluid])][(int)(k+cz[l_fluid])].temp;

//                             double cs = get_soundspeed(temp_d);

//                             p_out       = p_a;
//                             rho_out     = rho_d ;//+ (p_out - p_d) / (cs*cs);
//                             vel_out[0]  = u_d   ;//+ sign_side*(p_d - p_out) / (rho_d*cs);
//                             vel_out[1]  = v_d   ;//+ sign_side*(p_d - p_out) / (rho_d*cs);
//                             vel_out[2]  = w_d   ;//+ sign_side*(p_d - p_out) / (rho_d*cs);
//                             T_out       = temp_d;
//                             rho_out = p_out / (gas_const * T_out);
//                         }
//                         else if (mixture[i_nb][j_nb][k_nb].type==TYPE_O_C){   
//                             int i_1 = i + cx[l_fluid];
//                             int j_1 = j + cy[l_fluid];
//                             int k_1 = k + cz[l_fluid];
//                             int i_2 = i + 2*cx[l_fluid];
//                             int j_2 = j + 2*cy[l_fluid];
//                             int k_2 = k + 2*cz[l_fluid];
                            
//                             double drho_dx  = fd_uw(mixture[i][j][k].rho, mixture[i_1][j_1][k_1].rho, mixture[i_2][j_2][k_2].rho, dx, sign_side);
//                             double du_dx    = fd_uw(mixture[i][j][k].u  , mixture[i_1][j_1][k_1].u  , mixture[i_2][j_2][k_2].u  , dx, sign_side);
//                             double dv_dx    = fd_uw(mixture[i][j][k].v  , mixture[i_1][j_1][k_1].v  , mixture[i_2][j_2][k_2].v  , dx, sign_side);
//                             double dw_dx    = fd_uw(mixture[i][j][k].w  , mixture[i_1][j_1][k_1].w  , mixture[i_2][j_2][k_2].w  , dx, sign_side);
//                             double dp_dx    = fd_uw(mixture[i][j][k].p  , mixture[i_1][j_1][k_1].p  , mixture[i_2][j_2][k_2].p  , dx, sign_side);

//                             // std::cout << mixture[i][j][k].p << "| " << mixture[i_1][j_1][k_1].p << " | " << mixture[i_2][j_2][k_2].p << " | " << dp_dx << std::endl;                         

//                             double u    = mixture[i][j][k].u*abs(cx[l_fluid]) + mixture[i][j][k].v*abs(cy[l_fluid]) + mixture[i][j][k].w*abs(cz[l_fluid]);
//                             double rho  = mixture[i][j][k].rho;
//                             double temp = mixture[i][j][k].temp;
//                             double cs   = sign_side*get_soundspeed(mixture[i][j][k].temp);
//                             double sigma = 0.7;
//                             double Mach_ref = 1E-1;
//                             double L_ref = 50; 
//                             double K    = 0.0;//sigma*(1.0-Mach_ref*Mach_ref)*abs(cs)/L_ref;

//                             // double L1 = K*(mixture[i][j][k].p - mixture[i_nb][j_nb][k_nb].p);//u*(drho_dx - 1.0/(cs*cs) * dp_dx);
//                             // double L2 = u*(abs(cz[l_fluid])*du_dx + abs(cx[l_fluid])*dv_dx + abs(cy[l_fluid])*dw_dx);
//                             // double L3 = u*(abs(cy[l_fluid])*du_dx + abs(cz[l_fluid])*dv_dx + abs(cx[l_fluid])*dw_dx);
//                             // double L4 = (u + cs) * (1.0/(rho*cs) * dp_dx + abs(cx[l_fluid])*du_dx + abs(cy[l_fluid])*dv_dx + abs(cz[l_fluid])*dw_dx);
//                             // double L5 = 0.0;

//                             double L1 = K*(mixture[i][j][k].p - mixture[i_nb][j_nb][k_nb].p);//(u-cs) * (dp_dx - rho*cs*du_dx);
//                             double L2 = (u   ) * (cs*cs*drho_dx - dp_dx);
//                             double L3 = (u   ) * (abs(cz[l_fluid])*du_dx + abs(cx[l_fluid])*dv_dx + abs(cy[l_fluid])*dw_dx);
//                             double L4 = (u   ) * (abs(cy[l_fluid])*du_dx + abs(cz[l_fluid])*dv_dx + abs(cx[l_fluid])*dw_dx);
//                             double L5 = (u+cs) * (dp_dx + rho*cs*du_dx);

//                             rho_out     = mixture[i][j][k].rho  - dt_sim*(1.0/(cs*cs) * (L2+0.5*(L5+L1)));
//                             p_out       = mixture[i][j][k].p    - dt_sim*(0.5*(L5+L1));
//                             vel_out[0]  = mixture[i][j][k].u    - dt_sim*(abs(cx[l_fluid])*0.5/(rho*cs)*(L5-L1) + abs(cy[l_fluid])*L4                   + abs(cz[l_fluid])*L3                   );
//                             vel_out[1]  = mixture[i][j][k].v    - dt_sim*(abs(cx[l_fluid])*L3                   + abs(cy[l_fluid])*0.5/(rho*cs)*(L5-L1) + abs(cz[l_fluid])*L4                   );
//                             vel_out[2]  = mixture[i][j][k].w    - dt_sim*(abs(cx[l_fluid])*L4                   + abs(cy[l_fluid])*L3                   + abs(cz[l_fluid])*0.5/(rho*cs)*(L5-L1) );
//                             T_out       = mixture[i][j][k].temp - dt_sim*(temp/(rho*cs*cs) * (-L2+0.5*(gamma-1.0)*(L5+L1))     );
//                             // rho_out = p_out / (gas_const * T_out);
//                             T_out     = p_out / (gas_const * rho_out);
//                         }

//                         // std::cout << i << " | " << vel_in[0] << " | " << T_in << std::endl;
//                         calculate_feq_geq(f_out, g_out, rho_out, vel_out, T_out);
                        
//                         // step 2: calculate f_loc, g_loc, and fa_loc (local distribution function)
//                         double rho_loc = 0.0;
//                         double rhou_loc = 0.0;
//                         double rhov_loc = 0.0;
//                         double rhow_loc = 0.0;
//                         double rhoe_loc = 0.0;

//                         for (int l=0; l < npop; ++l){
//                             i_nb = i - cx[l];
//                             j_nb = j - cy[l];
//                             k_nb = k - cz[l];

//                             if (mixture[i_nb][j_nb][k_nb].type==TYPE_O_E||mixture[i_nb][j_nb][k_nb].type==TYPE_O_C){
//                                 rho_loc += f_out[l];
//                                 rhou_loc += f_out[l]*cx[l];
//                                 rhov_loc += f_out[l]*cy[l];
//                                 rhow_loc += f_out[l]*cz[l];

//                                 #ifndef ISOTHERM
//                                 rhoe_loc += g_out[l];
//                                 #endif
//                             }
//                             else{
//                                 rho_loc += mixture[i][j][k].f[l];
//                                 rhou_loc += mixture[i][j][k].f[l]*cx[l];
//                                 rhov_loc += mixture[i][j][k].f[l]*cy[l];
//                                 rhow_loc += mixture[i][j][k].f[l]*cz[l];

//                                 #ifndef ISOTHERM
//                                 rhoe_loc += mixture[i][j][k].g[l];
//                                 #endif
//                             }
//                         }

//                         double vel_loc[3];
//                         vel_loc[0] = rhou_loc / rho_loc;
//                         vel_loc[1] = rhov_loc / rho_loc;
//                         vel_loc[2] = rhow_loc / rho_loc;

//                         double internalEnergy=rhoe_loc/rho_loc - 0.5*v_sqr(vel_loc[0], vel_loc[1], vel_loc[2]);
                        
//                         double cv = gas_const / (gamma - 1.0);
//                         double T_loc = internalEnergy / cv;                        

//                         calculate_feq_geq(f_loc, g_loc, rho_loc, vel_loc, T_loc);

//                         for (int l=0; l < npop; ++l){
//                             i_nb = i - cx[l];
//                             j_nb = j - cy[l];
//                             k_nb = k - cz[l];

//                             if (mixture[i_nb][j_nb][k_nb].type==TYPE_O_E||mixture[i_nb][j_nb][k_nb].type==TYPE_O_C){
//                                 mixture[i][j][k].f[l] = 2*f_out[l] - f_loc[l];

//                                 #ifndef ISOTHERM
//                                 mixture[i][j][k].g[l] = 2*g_out[l] - g_loc[l];
//                                 #endif
//                             }
//                             else{
//                                 mixture[i][j][k].f[l] = f_out[l] + mixture[i][j][k].f[l] - f_loc[l];

//                                 #ifndef ISOTHERM
//                                 mixture[i][j][k].g[l] = g_out[l] + mixture[i][j][k].g[l] - g_loc[l];
//                                 #endif
//                             }
//                         }                  

//                     } 
//                 }
//             }
//         }

//     #elif defined MULTICOMP
//         #ifdef PARALLEL 
//         #pragma omp parallel for schedule(static, 1) 
//         #endif
//         for(int i=0; i<Nx; ++i){
//             for(int j=0; j<Ny; ++j){
//                 for(int k = 0; k<Nz; ++k){
//                     if(mixture[i][j][k].type==TYPE_F){
//                         int i_nb, j_nb, k_nb;
//                         double fa_out[nSpecies][npop];
//                         double g_out[npop];
//                         double fa_loc[nSpecies][npop];
//                         double g_loc[npop];
//                         int l_fluid = 999;

//                         // check interface node                         
//                         for (int l=0; l < npop; ++l){
//                             i_nb = i - cx[l];
//                             j_nb = j - cy[l];
//                             k_nb = k - cz[l];

//                             if(mixture[i_nb][j_nb][k_nb].type==TYPE_O_E&& mixture[int (i+cx[l])][int (j+cy[l])][int (k+cz[l])].type==TYPE_F){    
//                                 l_fluid = l;
//                                 break;
//                             }                            
//                         }

//                         // check if there is the node is the interface with boundary
//                         if ( l_fluid == 999 || !onlyOne_1(cx[l_fluid], cy[l_fluid], cz[l_fluid]) )
//                             continue;

//                         double rho_out = 0.0;
//                         double vel_out[3] = {0.0};
//                         double T_out = 0.0;
//                         double rhoa_out[nSpecies] = {0.0};

//                         vel_out[0] = mixture[(int)(i+cx[l_fluid])][(int)(j+cy[l_fluid])][(int)(k+cz[l_fluid])].u;
//                         vel_out[1] = mixture[(int)(i+cx[l_fluid])][(int)(j+cy[l_fluid])][(int)(k+cz[l_fluid])].v;
//                         vel_out[2] = mixture[(int)(i+cx[l_fluid])][(int)(j+cy[l_fluid])][(int)(k+cz[l_fluid])].w;
//                         T_out      = mixture[(int)(i+cx[l_fluid])][(int)(j+cy[l_fluid])][(int)(k+cz[l_fluid])].temp;
//                         rho_out    = mixture[(int)(i+cx[l_fluid])][(int)(j+cy[l_fluid])][(int)(k+cz[l_fluid])].rho;
//                         for(size_t a = 0; a < nSpecies; ++a)
//                             rhoa_out[a]    = species[a][(int)(i+cx[l_fluid])][(int)(j+cy[l_fluid])][(int)(k+cz[l_fluid])].rho;

//                         // double pres_far = mixture[(int)(i-cx[l_fluid])][(int)(j-cy[l_fluid])][(int)(k-cz[l_fluid])].p;

//                         // int rank = omp_get_thread_num();
//                         // auto gas = sols[rank]->thermo();   
//                         // std::vector <double> Y (gas->nSpecies());
//                         // for(size_t a = 0; a < nSpecies; ++a) Y[gas->speciesIndex(speciesName[a])] = rhoa_out[a] / rho_out;
//                         // gas->setState_TP(units.si_temp(T_out), units.si_p(pres_far));
//                         // gas->setMassFractions(&Y[0]);  

//                         // rho_out = units.rho(gas->density());  
//                         // for(size_t a = 0; a < nSpecies; ++a) rhoa_out[a] = Y[gas->speciesIndex(speciesName[a])] * rho_out;



//                         // std::cout << i << " | " << vel_in[0] << " | " << T_in << std::endl;
//                         calculate_feq_geq(fa_out, g_out, rho_out, rhoa_out, vel_out, T_out);

//                         // // step 2: calculate f_loc, g_loc, and fa_loc (local distribution function)
//                         double rho_loc = 0.0;
//                         double vel_loc[3] = {0.0};
//                         double rhoe_loc = 0.0;
//                         double rhoa_loc[nSpecies] = {0.0};
//                         double vela_loc[nSpecies][3] = {0.0};
                        
//                         for (int l=0; l < npop; ++l){
//                             i_nb = i - cx[l];
//                             j_nb = j - cy[l];
//                             k_nb = k - cz[l];

//                             if (mixture[i_nb][j_nb][k_nb].type==TYPE_O_E){
//                                 for (size_t a = 0; a < nSpecies; ++a){
//                                     rhoa_loc[a] += fa_out[a][l];
//                                     vela_loc[a][0] += fa_out[a][l]*cx[l];
//                                     vela_loc[a][1] += fa_out[a][l]*cy[l];
//                                     vela_loc[a][2] += fa_out[a][l]*cz[l];
//                                 }
//                                 #ifndef ISOTHERM
//                                 rhoe_loc += g_out[l];
//                                 #endif
//                             }
//                             else{
//                                 for (size_t a = 0; a < nSpecies; ++a){
//                                     rhoa_loc[a] += species[a][i][j][k].f[l];
//                                     vela_loc[a][0] += species[a][i][j][k].f[l]*cx[l];
//                                     vela_loc[a][1] += species[a][i][j][k].f[l]*cy[l];
//                                     vela_loc[a][2] += species[a][i][j][k].f[l]*cz[l];
//                                 }
//                                 #ifndef ISOTHERM
//                                 rhoe_loc += mixture[i][j][k].g[l];
//                                 #endif
//                             }
//                         }

//                         for (size_t a = 0; a < nSpecies; ++a){
//                             rho_loc += rhoa_loc[a];
//                             vel_loc[0] += vela_loc[a][0];
//                             vel_loc[1] += vela_loc[a][1];
//                             vel_loc[2] += vela_loc[a][2];

//                             if(rhoa_loc[a] > 0.0){
//                                 vela_loc[a][0] = vela_loc[a][0] / rhoa_loc[a];
//                                 vela_loc[a][1] = vela_loc[a][1] / rhoa_loc[a];
//                                 vela_loc[a][2] = vela_loc[a][2] / rhoa_loc[a];
//                             }
//                             else{
//                                 vela_loc[a][0] = 0.0;
//                                 vela_loc[a][1] = 0.0;
//                                 vela_loc[a][2] = 0.0;
//                             }
//                         }
//                         if(rho_loc > 0){
//                             vel_loc[0] = vel_loc[0] / rho_loc;
//                             vel_loc[1] = vel_loc[1] / rho_loc;
//                             vel_loc[2] = vel_loc[2] / rho_loc;
//                         }
//                         else{
//                             vel_loc[0] = 0.0;
//                             vel_loc[1] = 0.0;
//                             vel_loc[2] = 0.0;
//                         }

//                         #ifndef ISOTHERM
//                         double internalEnergy=rhoe_loc/rho_loc - 0.5*v_sqr(vel_loc[0], vel_loc[1], vel_loc[2]);      
//                         double T_loc = calculate_temp(internalEnergy, rho_loc, rhoa_loc);
//                         #else
//                         double T_loc = T_out;
//                         #endif
//                         calculate_feq_geq(fa_loc, g_loc, rho_loc, rhoa_loc, vel_loc, vela_loc, T_out);

//                         for (int l=0; l < npop; ++l){
//                             i_nb = i - cx[l];
//                             j_nb = j - cy[l];
//                             k_nb = k - cz[l];

//                             // for(size_t a = 0; a < nSpecies; ++a)
//                             //     std::cout << a << " | " << i << " | " << j << " | " << k << " | " << fa_loc[a][l] << std::endl;

//                             if (mixture[i_nb][j_nb][k_nb].type==TYPE_O_E){
//                                 for(size_t a = 0; a < nSpecies; ++a)
//                                     species[a][i][j][k].f[l] = 2*fa_out[a][l] - fa_loc[a][l];

//                                 #ifndef ISOTHERM
//                                 mixture[i][j][k].g[l] = 2*g_out[l] - g_loc[l];
//                                 #endif
//                             }
//                             else{
//                                 for(size_t a = 0; a < nSpecies; ++a)
//                                     species[a][i][j][k].f[l] = fa_out[a][l] + species[a][i][j][k].f[l] - fa_loc[a][l];

//                                 #ifndef ISOTHERM
//                                 mixture[i][j][k].g[l] = g_out[l] + mixture[i][j][k].g[l] - g_loc[l];
//                                 #endif
//                             }
//                         }                  

//                     } 
//                 }
//             }
//         }

//     #endif
// }

// void LBM::dirSlip(int l, int i, int j, int k, int &lp, int &ip, int &jp, int &kp)
// {
//     int i_nb, j_nb, k_nb;
//     i_nb = i - cx[l];
//     j_nb = j - cy[l];
//     k_nb = k - cz[l];

//     bool i_chk = false;
//     bool j_chk = false;
//     bool k_chk = false;
//     if (mixture[i_nb][j][k].type==TYPE_FS)
//         i_chk = true;
//     if (mixture[i][j_nb][k].type==TYPE_FS)
//         j_chk = true;
//     if (mixture[i][j][k_nb].type==TYPE_FS)
//         k_chk = true;

//     double cxp = cx[l];
//     double cyp = cy[l];
//     double czp = cz[l];
//     ip = i_nb;
//     jp = j_nb;
//     kp = k_nb;
//     if (i_chk == true){
//         cxp = -1.0*cx[l];
//         ip = i;
//     }
//     if (j_chk == true){
//         cyp = -1.0*cy[l];
//         jp = j;
//     }
//     if (k_chk == true){
//         czp = -1.0*cz[l];
//         kp = k;
//     } 
    
//     for(size_t a = 0; a < npop; ++a){
//         if (cxp == cx[a] && cyp == cy[a] && czp == cz[a]){
//             lp = a;
//             break;
//         }
//     }

// }


// ===================================================================================================================================================

#include "lbm.hpp"
#include "math_util.hpp"
#include <omp.h>
#include "units.hpp"

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
    // TMS Boundary Conditions ---------------------------------------------------------------------------------------------
    
    // SOLID ===============================================================================================================
    #ifndef MULTICOMP
    #ifdef PARALLEL 
    #pragma omp parallel for schedule(static, 1) 
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

    #elif defined MULTICOMP

    #ifdef PARALLEL 
    #pragma omp parallel for schedule(static, 1) 
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

                        vela_loc[a][0] = vela_loc[a][0] / rhoa_loc[a];
                        vela_loc[a][1] = vela_loc[a][1] / rhoa_loc[a];
                        vela_loc[a][2] = vela_loc[a][2] / rhoa_loc[a];
                    }
                    vel_loc[0] = vel_loc[0] / rho_loc;
                    vel_loc[1] = vel_loc[1] / rho_loc;
                    vel_loc[2] = vel_loc[2] / rho_loc;

                    double internalEnergy=rhoe_loc/rho_loc - 0.5*v_sqr(vel_loc[0], vel_loc[1], vel_loc[2]);                     
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
    #endif

    // INFLOW ============================================================================================================= 
    #ifndef MULTICOMP
        #ifdef PARALLEL 
        #pragma omp parallel for schedule(static, 1) 
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
                        int l_fluid = 999;

                        // check interface node                         
                        for (int l=0; l < npop; ++l){
                            i_nb = i - cx[l];
                            j_nb = j - cy[l];
                            k_nb = k - cz[l];

                            if(mixture[i_nb][j_nb][k_nb].type==TYPE_I_E && mixture[(int)(i+cx[l])][(int)(j+cy[l])][(int)(k+cz[l])].type==TYPE_F){    
                                l_fluid = l;
                                break;
                            }                            
                        }                        

                        // check if there is the node is the interface with boundary
                        if(l_fluid == 999 || !onlyOne_1(cx[l_fluid], cy[l_fluid], cz[l_fluid]))
                            continue;

                        double vel_in[3] = {0.0};
                        double T_in = 0.0;
                        double rho_in = 0.0;
                        double p_in = 0.0;
 
                        vel_in[0]   = mixture[(int)(i-cx[l_fluid])][(int)(j-cy[l_fluid])][(int)(k-cz[l_fluid])].u;
                        vel_in[1]   = mixture[(int)(i-cx[l_fluid])][(int)(j-cy[l_fluid])][(int)(k-cz[l_fluid])].v;
                        vel_in[2]   = mixture[(int)(i-cx[l_fluid])][(int)(j-cy[l_fluid])][(int)(k-cz[l_fluid])].w;
                        T_in        = mixture[(int)(i-cx[l_fluid])][(int)(j-cy[l_fluid])][(int)(k-cz[l_fluid])].temp;
                        rho_in      = mixture[(int)(i-cx[l_fluid])][(int)(j-cy[l_fluid])][(int)(k-cz[l_fluid])].rho;
                        p_in        = mixture[(int)(i+cx[l_fluid])][(int)(j+cy[l_fluid])][(int)(k+cz[l_fluid])].p;

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

                            if (mixture[i_nb][j_nb][k_nb].type==TYPE_I_E){
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

                            if (mixture[i_nb][j_nb][k_nb].type==TYPE_I_E){
                                // mixture[i][j][k].f[l] = f_tgt[l] + f_in[l] - f_loc[l];
                                mixture[i][j][k].f[l] = 2.0*f_in[l] - f_loc[l];

                                #ifndef ISOTHERM
                                // mixture[i][j][k].g[l] = g_tgt[l] + g_in[l] - g_loc[l];
                                mixture[i][j][k].g[l] = 2.0*g_in[l] - g_loc[l];
                                #endif
                            }
                            else{
                                // mixture[i][j][k].f[l] = f_tgt[l] + mixture[i][j][k].f[l] - f_loc[l];
                                mixture[i][j][k].f[l] = f_in[l] + mixture[i][j][k].f[l] - f_loc[l];

                                #ifndef ISOTHERM
                                // mixture[i][j][k].g[l] = g_tgt[l] + mixture[i][j][k].g[l] - g_loc[l];
                                mixture[i][j][k].g[l] = g_in[l] + mixture[i][j][k].g[l] - g_loc[l];
                                #endif
                            }
                        }                  

                    } 
                }
            }
        }

    #elif defined MULTICOMP
        #ifdef PARALLEL 
        #pragma omp parallel for schedule(static, 1) 
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
                        int l_fluid = 999;

                        // check interface node                         
                        for (int l=0; l < npop; ++l){
                            i_nb = i - cx[l];
                            j_nb = j - cy[l];
                            k_nb = k - cz[l];

                            if(mixture[i_nb][j_nb][k_nb].type==TYPE_I_E && mixture[int (i+cx[l])][int (j+cy[l])][int (k+cz[l])].type==TYPE_F){    
                                l_fluid = l;
                                break;
                            }                            
                        }

                        // check if there is the node is the interface with boundary
                        if (l_fluid == 999 || !onlyOne_1(cx[l_fluid], cy[l_fluid], cz[l_fluid]) )
                            continue;

                        double vel_in[3] = {0.0};
                        double T_in = 0.0;
                        double rho_in = 0.0;
                        double rhoa_in[nSpecies] = {0.0};
                        double p_in = 0.0;

                        vel_in[0]   = mixture[(int)(i-cx[l_fluid])][(int)(j-cy[l_fluid])][(int)(k-cz[l_fluid])].u;
                        vel_in[1]   = mixture[(int)(i-cx[l_fluid])][(int)(j-cy[l_fluid])][(int)(k-cz[l_fluid])].v;
                        vel_in[2]   = mixture[(int)(i-cx[l_fluid])][(int)(j-cy[l_fluid])][(int)(k-cz[l_fluid])].w;
                        T_in        = mixture[(int)(i-cx[l_fluid])][(int)(j-cy[l_fluid])][(int)(k-cz[l_fluid])].temp;
                        rho_in      = mixture[(int)(i-cx[l_fluid])][(int)(j-cy[l_fluid])][(int)(k-cz[l_fluid])].rho;
                        for(size_t a = 0; a < nSpecies; ++a)
                            rhoa_in[a] = species[a][(int)(i-cx[l_fluid])][(int)(j-cy[l_fluid])][(int)(k-cz[l_fluid])].rho;
                        p_in        = mixture[(int)(i+cx[l_fluid])][(int)(j+cy[l_fluid])][(int)(k+cz[l_fluid])].p;

                        int rank = omp_get_thread_num();
                        auto gas = sols[rank]->thermo();   
                        std::vector <double> Y (gas->nSpecies());
                        for(size_t a = 0; a < nSpecies; ++a) Y[gas->speciesIndex(speciesName[a])] = rhoa_in[a] / rho_in;
                        gas->setState_TP(units.si_temp(T_in), units.si_p(p_in));
                        gas->setMassFractions(&Y[0]);  

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

                            if (mixture[i_nb][j_nb][k_nb].type==TYPE_I_E){
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
                        double internalEnergy=rhoe_loc/rho_loc - 0.5*v_sqr(vel_loc[0], vel_loc[1], vel_loc[2]);                     
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

                            if (mixture[i_nb][j_nb][k_nb].type==TYPE_I_E){
                                for(size_t a = 0; a < nSpecies; ++a)
                                    species[a][i][j][k].f[l] = 2.0*fa_in[a][l] - fa_loc[a][l];

                                #ifndef ISOTHERM
                                mixture[i][j][k].g[l] = 2.0*g_in[l] - g_loc[l];
                                #endif
                            }
                            else{
                                for(size_t a = 0; a < nSpecies; ++a)
                                    species[a][i][j][k].f[l] = fa_in[a][l] + species[a][i][j][k].f[l] - fa_loc[a][l];

                                #ifndef ISOTHERM
                                mixture[i][j][k].g[l] = g_in[l] + mixture[i][j][k].g[l] - g_loc[l];
                                #endif
                            }
                        }                  

                    } 
                }
            }
        }
    #endif

    // OUTFLOW
    #ifndef MULTICOMP
        #ifdef PARALLEL 
        #pragma omp parallel for schedule(static, 1) 
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
                        int l_fluid = 999;

                        // check interface node                         
                        for (int l=0; l < npop; ++l){
                            i_nb = i - cx[l];
                            j_nb = j - cy[l];
                            k_nb = k - cz[l];

                            if(mixture[i_nb][j_nb][k_nb].type==TYPE_O_E && mixture[int (i+cx[l])][int (j+cy[l])][int (k+cz[l])].type==TYPE_F){    
                                l_fluid = l;
                                break;
                            }                            
                        }

                        // check if there is the node is the interface with boundary
                        if (l_fluid == 999 )
                            continue;

                        double vel_out[3] = {0.0};
                        double T_out = 0.0;
                        double rho_out = 0.0;
                        double p_out = 0.0;

                        vel_out[0]  = mixture[(int)(i+cx[l_fluid])][(int)(j+cy[l_fluid])][(int)(k+cz[l_fluid])].u;
                        vel_out[1]  = mixture[(int)(i+cx[l_fluid])][(int)(j+cy[l_fluid])][(int)(k+cz[l_fluid])].v;
                        vel_out[2]  = mixture[(int)(i+cx[l_fluid])][(int)(j+cy[l_fluid])][(int)(k+cz[l_fluid])].w;
                        T_out       = mixture[(int)(i+cx[l_fluid])][(int)(j+cy[l_fluid])][(int)(k+cz[l_fluid])].temp;
                        rho_out     = mixture[(int)(i+cx[l_fluid])][(int)(j+cy[l_fluid])][(int)(k+cz[l_fluid])].rho;
                        p_out       = mixture[(int)(i-cx[l_fluid])][(int)(j-cy[l_fluid])][(int)(k-cz[l_fluid])].p;

                        // rho_out     = p_out/(gas_const * T_out);

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

                            if (mixture[i_nb][j_nb][k_nb].type==TYPE_O_E){
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

                            if (mixture[i_nb][j_nb][k_nb].type==TYPE_O_E){
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

    #elif defined MULTICOMP
        #ifdef PARALLEL 
        #pragma omp parallel for schedule(static, 1) 
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
                        int l_fluid = 999;

                        // check interface node                         
                        for (int l=0; l < npop; ++l){
                            i_nb = i - cx[l];
                            j_nb = j - cy[l];
                            k_nb = k - cz[l];

                            if(mixture[i_nb][j_nb][k_nb].type==TYPE_O_E && mixture[int (i+cx[l])][int (j+cy[l])][int (k+cz[l])].type==TYPE_F){    
                                l_fluid = l;
                                break;
                            }                            
                        }

                        // check if there is the node is the interface with boundary
                        if ( l_fluid == 999 )
                            continue;

                        double rho_out = 0.0;
                        double vel_out[3] = {0.0};
                        double T_out = 0.0;
                        double rhoa_out[nSpecies] = {0.0};
                        double vela_out [nSpecies][3] = {0.0};

                        vel_out[0] = mixture[(int)(i+cx[l_fluid])][(int)(j+cy[l_fluid])][(int)(k+cz[l_fluid])].u;
                        vel_out[1] = mixture[(int)(i+cx[l_fluid])][(int)(j+cy[l_fluid])][(int)(k+cz[l_fluid])].v;
                        vel_out[2] = mixture[(int)(i+cx[l_fluid])][(int)(j+cy[l_fluid])][(int)(k+cz[l_fluid])].w;
                        T_out      = mixture[(int)(i+cx[l_fluid])][(int)(j+cy[l_fluid])][(int)(k+cz[l_fluid])].temp;
                        rho_out    = mixture[(int)(i+cx[l_fluid])][(int)(j+cy[l_fluid])][(int)(k+cz[l_fluid])].rho;
                        for(size_t a = 0; a < nSpecies; ++a){
                            rhoa_out[a]     = species[a][(int)(i+cx[l_fluid])][(int)(j+cy[l_fluid])][(int)(k+cz[l_fluid])].rho;
                            vela_out[a][0]  = species[a][(int)(i+cx[l_fluid])][(int)(j+cy[l_fluid])][(int)(k+cz[l_fluid])].u;
                            vela_out[a][1]  = species[a][(int)(i+cx[l_fluid])][(int)(j+cy[l_fluid])][(int)(k+cz[l_fluid])].v;
                            vela_out[a][2]  = species[a][(int)(i+cx[l_fluid])][(int)(j+cy[l_fluid])][(int)(k+cz[l_fluid])].w;
                        }

                        // std::cout << i << " | " << vel_in[0] << " | " << T_in << std::endl;
                        calculate_feq_geq(fa_out, g_out, rho_out, rhoa_out, vel_out, vela_out, T_out);

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

                            if (mixture[i_nb][j_nb][k_nb].type==TYPE_O_E){
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
                        double internalEnergy=rhoe_loc/rho_loc - 0.5*v_sqr(vel_loc[0], vel_loc[1], vel_loc[2]);      
                        double T_loc = calculate_temp(internalEnergy, rho_loc, rhoa_loc);
                        #else
                        double T_loc = T_out;
                        #endif
                        calculate_feq_geq(fa_loc, g_loc, rho_loc, rhoa_loc, vel_loc, vela_loc, T_out);

                        for (int l=0; l < npop; ++l){
                            i_nb = i - cx[l];
                            j_nb = j - cy[l];
                            k_nb = k - cz[l];

                            // for(size_t a = 0; a < nSpecies; ++a)
                            //     std::cout << a << " | " << i << " | " << j << " | " << k << " | " << fa_loc[a][l] << std::endl;

                            if (mixture[i_nb][j_nb][k_nb].type==TYPE_O_E){
                                for(size_t a = 0; a < nSpecies; ++a)
                                    species[a][i][j][k].f[l] = 2*fa_out[a][l] - fa_loc[a][l];

                                #ifndef ISOTHERM
                                mixture[i][j][k].g[l] = 2*g_out[l] - g_loc[l];
                                #endif
                            }
                            else{
                                for(size_t a = 0; a < nSpecies; ++a)
                                    species[a][i][j][k].f[l] = fa_out[a][l] + species[a][i][j][k].f[l] - fa_loc[a][l];

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
    if (j_chk == true){
        cyp = -1.0*cy[l];
        jp = j;
    }
    if (k_chk == true){
        czp = -1.0*cz[l];
        kp = k;
    } 
    
    for(size_t a = 0; a < npop; ++a){
        if (cxp == cx[a] && cyp == cy[a] && czp == cz[a]){
            lp = a;
            break;
        }
    }

}