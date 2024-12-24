// #include <iostream>
// #include <omp.h>
// #include "lbm.hpp"
// #include "units.hpp"

// void impose_NSCBC(LBM &lb, int i, int j, int k, int l_interface)
// {
// #if defined MULTICOMP

//     // Initialize local variables for derivatives
//     double dpdx = 0.0, dudx = 0.0, dvdx = 0.0, dwdx = 0.0, drhodx = 0.0;
//     double dpdy = 0.0, dudy = 0.0, dvdy = 0.0, dwdy = 0.0, drhody = 0.0; 
//     double dpdz = 0.0, dudz = 0.0, dvdz = 0.0, dwdz = 0.0, drhodz = 0.0; 
//     int dx = lb.get_dx(); int dy = lb.get_dy(); int dz = lb.get_dz();

//     double T1 = 0.0, T2 = 0.0, T3 = 0.0, T4 = 0.0, T5 = 0.0;
    
//     if (cx[l_interface] != 0.0){
//         normal_derivative(lb, i, j, k, 1, cx[l_interface], dx, dpdx, dudx, dvdx, dwdx, drhodx);
//         tangential_derivative(lb, i, j, k, 2, dy, dpdy, dudy, dvdy, dwdy, drhody);
//         tangential_derivative(lb, i, j, k, 3, dz, dpdz, dudz, dvdz, dwdz, drhodz);
//         compute_tranverse_terms(lb, i, j, k, 1, T1, T2, T3, T4, T5, dpdx, dudx, dvdx, dwdx, drhodx, dpdy, dudy, dvdy, dwdy, drhody, dpdz, dudz, dvdz, dwdz, drhodz);


//     }
//     else if (cy[l_interface] != 0.0){

//     }
//     else if (cz[l_interface] != 0.0){

//     }
    
    

    
// #endif

// }


// void normal_derivative(LBM &lb, int i, int j, int k, int idir, int isign, double delta, double &dp, double &du, double &dv, double &dw, double &drho)
// {
//     if (idir == 1){
//         if (isign == 1){
//             dp = (-1.5*lb.mixture[i][j][k].p + 2.0*lb.mixture[i+1][j][k].p - 0.5*lb.mixture[i+2][j][k].p ) / delta;
//             du = (-1.5*lb.mixture[i][j][k].u + 2.0*lb.mixture[i+1][j][k].u - 0.5*lb.mixture[i+2][j][k].u ) / delta;
//             dv = (-1.5*lb.mixture[i][j][k].v + 2.0*lb.mixture[i+1][j][k].v - 0.5*lb.mixture[i+2][j][k].v ) / delta;
//             dw = (-1.5*lb.mixture[i][j][k].w + 2.0*lb.mixture[i+1][j][k].w - 0.5*lb.mixture[i+2][j][k].w ) / delta;
//             drho = (-1.5*lb.mixture[i][j][k].rho + 2.0*lb.mixture[i+1][j][k].rho - 0.5*lb.mixture[i+2][j][k].rho ) / delta;
//         }
//         else if (isign == -1){
//             dp = (1.5*lb.mixture[i][j][k].p - 2.0*lb.mixture[i-1][j][k].p + 0.5*lb.mixture[i-2][j][k].p ) / delta;
//             du = (1.5*lb.mixture[i][j][k].u - 2.0*lb.mixture[i-1][j][k].u + 0.5*lb.mixture[i-2][j][k].u ) / delta;
//             dv = (1.5*lb.mixture[i][j][k].v - 2.0*lb.mixture[i-1][j][k].v + 0.5*lb.mixture[i-2][j][k].v ) / delta;
//             dw = (1.5*lb.mixture[i][j][k].w - 2.0*lb.mixture[i-1][j][k].w + 0.5*lb.mixture[i-2][j][k].w ) / delta;
//             drho = (1.5*lb.mixture[i][j][k].rho - 2.0*lb.mixture[i-1][j][k].rho + 0.5*lb.mixture[i-2][j][k].rho ) / delta;
//         }
//     }
//     else if (idir == 2){
//         if (isign == 1){
//             dp = (-1.5*lb.mixture[i][j][k].p + 2.0*lb.mixture[i][j+1][k].p - 0.5*lb.mixture[i][j+2][k].p ) / delta;
//             du = (-1.5*lb.mixture[i][j][k].u + 2.0*lb.mixture[i][j+1][k].u - 0.5*lb.mixture[i][j+2][k].u ) / delta;
//             dv = (-1.5*lb.mixture[i][j][k].v + 2.0*lb.mixture[i][j+1][k].v - 0.5*lb.mixture[i][j+2][k].v ) / delta;
//             dw = (-1.5*lb.mixture[i][j][k].w + 2.0*lb.mixture[i][j+1][k].w - 0.5*lb.mixture[i][j+2][k].w ) / delta;
//             drho = (-1.5*lb.mixture[i][j][k].rho + 2.0*lb.mixture[i][j+1][k].rho - 0.5*lb.mixture[i][j+2][k].rho ) / delta;
//         }
//         else if (isign == -1){
//             dp = (1.5*lb.mixture[i][j][k].p - 2.0*lb.mixture[i][j-1][k].p + 0.5*lb.mixture[i][j-2][k].p ) / delta;
//             du = (1.5*lb.mixture[i][j][k].u - 2.0*lb.mixture[i][j-1][k].u + 0.5*lb.mixture[i][j-2][k].u ) / delta;
//             dv = (1.5*lb.mixture[i][j][k].v - 2.0*lb.mixture[i][j-1][k].v + 0.5*lb.mixture[i][j-2][k].v ) / delta;
//             dw = (1.5*lb.mixture[i][j][k].w - 2.0*lb.mixture[i][j-1][k].w + 0.5*lb.mixture[i][j-2][k].w ) / delta;
//             drho = (1.5*lb.mixture[i][j][k].rho - 2.0*lb.mixture[i][j-1][k].rho + 0.5*lb.mixture[i][j-2][k].rho ) / delta;
//         }
//     }
//     else if (idir == 3){
//         if (isign == 1){
//             dp = (-1.5*lb.mixture[i][j][k].p + 2.0*lb.mixture[i][j][k+1].p - 0.5*lb.mixture[i][j][k+2].p ) / delta;
//             du = (-1.5*lb.mixture[i][j][k].u + 2.0*lb.mixture[i][j][k+1].u - 0.5*lb.mixture[i][j][k+2].u ) / delta;
//             dv = (-1.5*lb.mixture[i][j][k].v + 2.0*lb.mixture[i][j][k+1].v - 0.5*lb.mixture[i][j][k+2].v ) / delta;
//             dw = (-1.5*lb.mixture[i][j][k].w + 2.0*lb.mixture[i][j][k+1].w - 0.5*lb.mixture[i][j][k+2].w ) / delta;
//             drho = (-1.5*lb.mixture[i][j][k].rho + 2.0*lb.mixture[i][j][k+1].rho - 0.5*lb.mixture[i][j][k+2].rho ) / delta;
//         }
//         else if (isign == -1){
//             dp = (1.5*lb.mixture[i][j][k].p - 2.0*lb.mixture[i][j][k-1].p + 0.5*lb.mixture[i][j][k-2].p ) / delta;
//             du = (1.5*lb.mixture[i][j][k].u - 2.0*lb.mixture[i][j][k-1].u + 0.5*lb.mixture[i][j][k-2].u ) / delta;
//             dv = (1.5*lb.mixture[i][j][k].v - 2.0*lb.mixture[i][j][k-1].v + 0.5*lb.mixture[i][j][k-2].v ) / delta;
//             dw = (1.5*lb.mixture[i][j][k].w - 2.0*lb.mixture[i][j][k-1].w + 0.5*lb.mixture[i][j][k-2].w ) / delta;
//             drho = (1.5*lb.mixture[i][j][k].rho - 2.0*lb.mixture[i][j][k-1].rho + 0.5*lb.mixture[i][j][k-2].rho ) / delta;
//         }
//     }

// }

// void tangential_derivative(LBM &lb, int i, int j, int k, int idir, double delta, double &dp, double &du, double &dv, double &dw, double &drho)
// {
//     if (idir == 1){
//         if (lb.mixture[i-1][j][k].type != TYPE_F){
//             dp = (-1.5*lb.mixture[i][j][k].p + 2.0*lb.mixture[i+1][j][k].p - 0.5*lb.mixture[i+2][j][k].p ) / delta;
//             du = (-1.5*lb.mixture[i][j][k].u + 2.0*lb.mixture[i+1][j][k].u - 0.5*lb.mixture[i+2][j][k].u ) / delta;
//             dv = (-1.5*lb.mixture[i][j][k].v + 2.0*lb.mixture[i+1][j][k].v - 0.5*lb.mixture[i+2][j][k].v ) / delta;
//             dw = (-1.5*lb.mixture[i][j][k].w + 2.0*lb.mixture[i+1][j][k].w - 0.5*lb.mixture[i+2][j][k].w ) / delta;
//             drho = (-1.5*lb.mixture[i][j][k].rho + 2.0*lb.mixture[i+1][j][k].rho - 0.5*lb.mixture[i+2][j][k].rho ) / delta;
//         }
//         else if (lb.mixture[i+1][j][k].type != TYPE_F){
//             dp = (1.5*lb.mixture[i][j][k].p - 2.0*lb.mixture[i-1][j][k].p + 0.5*lb.mixture[i-2][j][k].p ) / delta;
//             du = (1.5*lb.mixture[i][j][k].u - 2.0*lb.mixture[i-1][j][k].u + 0.5*lb.mixture[i-2][j][k].u ) / delta;
//             dv = (1.5*lb.mixture[i][j][k].v - 2.0*lb.mixture[i-1][j][k].v + 0.5*lb.mixture[i-2][j][k].v ) / delta;
//             dw = (1.5*lb.mixture[i][j][k].w - 2.0*lb.mixture[i-1][j][k].w + 0.5*lb.mixture[i-2][j][k].w ) / delta;
//             drho = (1.5*lb.mixture[i][j][k].rho - 2.0*lb.mixture[i-1][j][k].rho + 0.5*lb.mixture[i-2][j][k].rho ) / delta;
//         }
//         else{
//             dp = (lb.mixture[i+1][j][k].p - lb.mixture[i-1][j][k].p ) / (2.0 * delta);
//             du = (lb.mixture[i+1][j][k].u - lb.mixture[i-1][j][k].u ) / (2.0 * delta);
//             dv = (lb.mixture[i+1][j][k].v - lb.mixture[i-1][j][k].v ) / (2.0 * delta);
//             dw = (lb.mixture[i+1][j][k].w - lb.mixture[i-1][j][k].w ) / (2.0 * delta);
//             drho = (lb.mixture[i+1][j][k].rho - lb.mixture[i-1][j][k].rho ) / (2.0 * delta);
//         }
//     }
//     else if (idir == 2){
//         if (lb.mixture[i][j-1][k].type != TYPE_F){
//             dp = (-1.5*lb.mixture[i][j][k].p + 2.0*lb.mixture[i][j+1][k].p - 0.5*lb.mixture[i][j+2][k].p ) / delta;
//             du = (-1.5*lb.mixture[i][j][k].u + 2.0*lb.mixture[i][j+1][k].u - 0.5*lb.mixture[i][j+2][k].u ) / delta;
//             dv = (-1.5*lb.mixture[i][j][k].v + 2.0*lb.mixture[i][j+1][k].v - 0.5*lb.mixture[i][j+2][k].v ) / delta;
//             dw = (-1.5*lb.mixture[i][j][k].w + 2.0*lb.mixture[i][j+1][k].w - 0.5*lb.mixture[i][j+2][k].w ) / delta;
//             drho = (-1.5*lb.mixture[i][j][k].rho + 2.0*lb.mixture[i][j+1][k].rho - 0.5*lb.mixture[i][j+2][k].rho ) / delta;
//         }
//         else if (lb.mixture[i][j+1][k].type != TYPE_F){
//             dp = (1.5*lb.mixture[i][j][k].p - 2.0*lb.mixture[i][j-1][k].p + 0.5*lb.mixture[i][j-2][k].p ) / delta;
//             du = (1.5*lb.mixture[i][j][k].u - 2.0*lb.mixture[i][j-1][k].u + 0.5*lb.mixture[i][j-2][k].u ) / delta;
//             dv = (1.5*lb.mixture[i][j][k].v - 2.0*lb.mixture[i][j-1][k].v + 0.5*lb.mixture[i][j-2][k].v ) / delta;
//             dw = (1.5*lb.mixture[i][j][k].w - 2.0*lb.mixture[i][j-1][k].w + 0.5*lb.mixture[i][j-2][k].w ) / delta;
//             drho = (1.5*lb.mixture[i][j][k].rho - 2.0*lb.mixture[i][j-1][k].rho + 0.5*lb.mixture[i][j-2][k].rho ) / delta;
//         }
//         else {
//             dp = (lb.mixture[i][j+1][k].p - lb.mixture[i][j-1][k].p ) / (2.0 * delta);
//             du = (lb.mixture[i][j+1][k].u - lb.mixture[i][j-1][k].u ) / (2.0 * delta);
//             dv = (lb.mixture[i][j+1][k].v - lb.mixture[i][j-1][k].v ) / (2.0 * delta);
//             dw = (lb.mixture[i][j+1][k].w - lb.mixture[i][j-1][k].w ) / (2.0 * delta);
//             drho = (lb.mixture[i][j+1][k].rho - lb.mixture[i][j-1][k].rho ) / (2.0 * delta);
//         }
//     }
//     else if (idir == 3){
//         if (lb.mixture[i][j][k-1].type != TYPE_F){
//             dp = (-1.5*lb.mixture[i][j][k].p + 2.0*lb.mixture[i][j][k+1].p - 0.5*lb.mixture[i][j][k+2].p ) / delta;
//             du = (-1.5*lb.mixture[i][j][k].u + 2.0*lb.mixture[i][j][k+1].u - 0.5*lb.mixture[i][j][k+2].u ) / delta;
//             dv = (-1.5*lb.mixture[i][j][k].v + 2.0*lb.mixture[i][j][k+1].v - 0.5*lb.mixture[i][j][k+2].v ) / delta;
//             dw = (-1.5*lb.mixture[i][j][k].w + 2.0*lb.mixture[i][j][k+1].w - 0.5*lb.mixture[i][j][k+2].w ) / delta;
//             drho = (-1.5*lb.mixture[i][j][k].rho + 2.0*lb.mixture[i][j][k+1].rho - 0.5*lb.mixture[i][j][k+2].rho ) / delta;
//         }
//         else if (lb.mixture[i][j][k-1].type != TYPE_F){
//             dp = (1.5*lb.mixture[i][j][k].p - 2.0*lb.mixture[i][j][k-1].p + 0.5*lb.mixture[i][j][k-2].p ) / delta;
//             du = (1.5*lb.mixture[i][j][k].u - 2.0*lb.mixture[i][j][k-1].u + 0.5*lb.mixture[i][j][k-2].u ) / delta;
//             dv = (1.5*lb.mixture[i][j][k].v - 2.0*lb.mixture[i][j][k-1].v + 0.5*lb.mixture[i][j][k-2].v ) / delta;
//             dw = (1.5*lb.mixture[i][j][k].w - 2.0*lb.mixture[i][j][k-1].w + 0.5*lb.mixture[i][j][k-2].w ) / delta;
//             drho = (1.5*lb.mixture[i][j][k].rho - 2.0*lb.mixture[i][j][k-1].rho + 0.5*lb.mixture[i][j][k-2].rho ) / delta;
//         }
//         else {
//             dp = (lb.mixture[i][j][k+1].p - lb.mixture[i][j][k-1].p ) / (2.0 * delta);
//             du = (lb.mixture[i][j][k+1].u - lb.mixture[i][j][k-1].u ) / (2.0 * delta);
//             dv = (lb.mixture[i][j][k+1].v - lb.mixture[i][j][k-1].v ) / (2.0 * delta);
//             dw = (lb.mixture[i][j][k+1].w - lb.mixture[i][j][k-1].w ) / (2.0 * delta);
//             drho = (lb.mixture[i][j][k+1].rho - lb.mixture[i][j][k-1].rho ) / (2.0 * delta);
//         }
//     }

// }

// void compute_tranverse_terms(LBM &lb, int i, int j, int k, int idir, double& T1, double& T2, double& T3, double& T4, double& T5,
// double dpdx, double dudx, double dvdx, double dwdx, double drhodx,
//     double dpdy, double dudy, double dvdy, double dwdy, double drhody,
//     double dpdz, double dudz, double dvdz, double dwdz, double drhodz)
// {
//     int nSpecies = lb.get_nSpecies();
//     std::vector<std::string> speciesName = lb.get_speciesName();

//     int rank = omp_get_thread_num();
//     auto gas = sols[rank]->thermo();   
//     std::vector <double> Y (gas->nSpecies());
//     for(size_t a = 0; a < nSpecies; ++a) Y[gas->speciesIndex(speciesName[a])] = lb.species[a][i][j][k].rho / lb.mixture[i][j][k].rho;
//     gas->setMassFractions(&Y[0]);
//     gas->setState_DP(units.si_rho(lb.mixture[i][j][k].rho), units.si_p(lb.mixture[i][j][k].p));

//     double soundSpeed = units.u(gas->soundSpeed());
//     double gamma = gas->cp_mass() / gas->cv_mass();

//     double inv_rho = 1.0 / lb.mixture[i][j][k].rho;

//     if (idir == 1) {
//         T1 = (lb.mixture[i][j][k].v * (dpdy - lb.mixture[i][j][k].rho * soundSpeed * dudy)) +
//              (lb.mixture[i][j][k].w * (dpdz - lb.mixture[i][j][k].rho * soundSpeed * dudz)) +
//              (gamma * lb.mixture[i][j][k].p * (dvdy + dwdz));

//         T2 = (lb.mixture[i][j][k].v * ((soundSpeed * soundSpeed * drhody) - dpdy)) +
//              (lb.mixture[i][j][k].w * ((soundSpeed * soundSpeed * drhodz) - dpdz));

//         T3 = lb.mixture[i][j][k].v * dvdy + lb.mixture[i][j][k].w * dvdz + dpdy * inv_rho;

//         T4 = lb.mixture[i][j][k].v * dwdy + lb.mixture[i][j][k].w * dwdz + dpdz * inv_rho;

//         T5 = (lb.mixture[i][j][k].v * (dpdy + lb.mixture[i][j][k].rho * soundSpeed * dudy)) +
//              (lb.mixture[i][j][k].w * (dpdz + lb.mixture[i][j][k].rho * soundSpeed * dudz)) +
//              (gamma * lb.mixture[i][j][k].p * (dvdy + dwdz));

//     } else if (idir == 2) {
//         T1 = (lb.mixture[i][j][k].u * (dpdx - lb.mixture[i][j][k].rho * soundSpeed * dvdx)) +
//              (lb.mixture[i][j][k].w * (dpdz - lb.mixture[i][j][k].rho * soundSpeed * dvdz)) +
//              (gamma * lb.mixture[i][j][k].p * (dudx + dwdz));

//         T2 = lb.mixture[i][j][k].u * dudx + lb.mixture[i][j][k].w * dudz + dpdx * inv_rho;

//         T3 = (lb.mixture[i][j][k].u * ((soundSpeed * soundSpeed * drhodx) - dpdx)) +
//              (lb.mixture[i][j][k].w * ((soundSpeed * soundSpeed * drhodz) - dpdz));

//         T4 = lb.mixture[i][j][k].u * dwdx + lb.mixture[i][j][k].w * dwdz + dpdz * inv_rho;

//         T5 = (lb.mixture[i][j][k].u * (dpdx + lb.mixture[i][j][k].rho * soundSpeed * dvdx)) +
//              (lb.mixture[i][j][k].w * (dpdz + lb.mixture[i][j][k].rho * soundSpeed * dvdz)) +
//              (gamma * lb.mixture[i][j][k].p * (dudx + dwdz));

//     } else if (idir == 3) {
//         T1 = (lb.mixture[i][j][k].u * (dpdx - lb.mixture[i][j][k].rho * soundSpeed * dwdx)) +
//              (lb.mixture[i][j][k].v * (dpdy - lb.mixture[i][j][k].rho * soundSpeed * dwdy)) +
//              (gamma * lb.mixture[i][j][k].p * (dudx + dvdy));

//         T2 = lb.mixture[i][j][k].u * dudx + lb.mixture[i][j][k].v * dudy + dpdx * inv_rho;

//         T3 = lb.mixture[i][j][k].u * dvdx + lb.mixture[i][j][k].v * dvdy + dpdy * inv_rho;

//         T4 = (lb.mixture[i][j][k].u * ((soundSpeed * soundSpeed * drhodx) - dpdx)) +
//              (lb.mixture[i][j][k].v * ((soundSpeed * soundSpeed * drhody) - dpdy));

//         T5 = (lb.mixture[i][j][k].u * (dpdx + lb.mixture[i][j][k].rho * soundSpeed * dwdx)) +
//              (lb.mixture[i][j][k].v * (dpdy + lb.mixture[i][j][k].rho * soundSpeed * dwdy)) +
//              (gamma * lb.mixture[i][j][k].p * (dudx + dvdy));

//     } else {
//         throw std::invalid_argument("Invalid idir in compute_transverse_terms");
//     }
// }

// void compute_waves(LBM &lb,
//     int i, int j, int idir, int isign,
//     double T1, double T2, double T3, double T4,
//     double &L1, double &L2, double &L3, double &L4,
//     double dp, double du, double dv, double drho,
// ) {
//     if (idir != 1 && idir != 2) {
//         throw std::runtime_error("Problem of idir in compute_waves");
//     }

//     if (isign != 1 && isign != -1) {
//         throw std::runtime_error("Problem of isign in compute_waves");
//     }

//     int nSpecies = lb.get_nSpecies();
//     std::vector<std::string> speciesName = lb.get_speciesName();

//     int rank = omp_get_thread_num();
//     auto gas = sols[rank]->thermo();   
//     std::vector <double> Y (gas->nSpecies());
//     for(size_t a = 0; a < nSpecies; ++a) Y[gas->speciesIndex(speciesName[a])] = lb.species[a][i][j][k].rho / lb.mixture[i][j][k].rho;
//     gas->setMassFractions(&Y[0]);
//     gas->setState_DP(units.si_rho(lb.mixture[i][j][k].rho), units.si_p(lb.mixture[i][j][k].p));

//     double mach = v_mag(lb.mixture[i][j][k].u, lb.mixture[i][j][k].v, lb.mixture[i][j][k].w) / units.u(gas->soundSpeed());



//     double TARGET_VX = bc_target[0];
//     double TARGET_VY = bc_target[1];
//     double TARGET_TEMPERATURE = bc_target[3];
//     double TARGET_PRESSURE = bc_target[4];

//     double relax_T = bc_params[0];
//     double relax_U = bc_params[1];
//     double relax_V = bc_params[2];

//     double beta = (bc_params[4] < 0.0) ? mach_local : bc_params[4];
//     double sigma_out = bc_params[5];

//     if (idir == 1) {
//         L1 = (q[(i * q_h2 + j) * QVAR + 1] - qaux[(i * qa_h2 + j) * NQAUX + 0]) *
//              (dp - q[(i * q_h2 + j) * QVAR + 0] * qaux[(i * qa_h2 + j) * NQAUX + 0] * du);
//         L2 = q[(i * q_h2 + j) * QVAR + 1] * (pow(qaux[(i * qa_h2 + j) * NQAUX + 0], 2.0) * drho - dp);
//         L3 = q[(i * q_h2 + j) * QVAR + 1] * dv;
//         L4 = (q[(i * q_h2 + j) * QVAR + 1] + qaux[(i * qa_h2 + j) * NQAUX + 0]) *
//              (dp + q[(i * q_h2 + j) * QVAR + 0] * qaux[(i * qa_h2 + j) * NQAUX + 0] * du);
//     } else if (idir == 2) {
//         L1 = (q[(i * q_h2 + j) * QVAR + 2] - qaux[(i * qa_h2 + j) * NQAUX + 0]) *
//              (dp - q[(i * q_h2 + j) * QVAR + 0] * qaux[(i * qa_h2 + j) * NQAUX + 0] * dv);
//         L2 = q[(i * q_h2 + j) * QVAR + 2] * du;
//         L3 = q[(i * q_h2 + j) * QVAR + 2] * (pow(qaux[(i * qa_h2 + j) * NQAUX + 0], 2.0) * drho - dp);
//         L4 = (q[(i * q_h2 + j) * QVAR + 2] + qaux[(i * qa_h2 + j) * NQAUX + 0]) *
//              (dp + q[(i * q_h2 + j) * QVAR + 0] * qaux[(i * qa_h2 + j) * NQAUX + 0] * dv);
//     }

//     if (bc_type == Inflow) {
//         if (idir == 1) {
//             if (isign == 1) {
//                 L4 = relax_U *
//                      ((q[(i * q_h2 + j) * QVAR + 0] * pow(qaux[(i * qa_h2 + j) * NQAUX + 0], 2.0)) *
//                       (1.0 - pow(mach_local, 2.0)) / probhi[idir - 1]) *
//                      (q[(i * q_h2 + j) * QVAR + 1] - TARGET_VX) -
//                      ((1.0 - beta) * T4);
//             } else if (isign == -1) {
//                 L1 = relax_U *
//                      ((q[(i * q_h2 + j) * QVAR + 0] * pow(qaux[(i * qa_h2 + j) * NQAUX + 0], 2.0)) *
//                       (1.0 - pow(mach_local, 2.0)) / probhi[idir - 1]) *
//                      (q[(i * q_h2 + j) * QVAR + 1] - TARGET_VX) -
//                      ((1.0 - beta) * T1);
//             }

//             L2 = relax_T *
//                  (q[(i * q_h2 + j) * QVAR + 0] * qaux[(i * qa_h2 + j) * NQAUX + 0] *
//                   qaux[(i * qa_h2 + j) * NQAUX + 3] / probhi[idir - 1]) *
//                  (q[(i * q_h2 + j) * QVAR + 4] - TARGET_TEMPERATURE) -
//                  ((1.0 - beta) * T2);

//             L3 = relax_V *
//                  (qaux[(i * qa_h2 + j) * NQAUX + 0] / probhi[idir - 1]) *
//                  (q[(i * q_h2 + j) * QVAR + 2] - TARGET_VY) -
//                  ((1.0 - beta) * T3);
//         } else if (idir == 2) {
//             // Similar computations for idir == 2...
//         }
//     } 
//     else if (bc_type == Outflow) {
//         double Kout = sigma_out * (1.0 - pow(mach_local, 2.0)) *
//                       (qaux[(i * qa_h2 + j) * NQAUX + 0] / probhi[idir - 1]);
//         if (isign == 1) {
//             L4 = (Kout * (q[(i * q_h2 + j) * QVAR + 3] - TARGET_PRESSURE)) -
//                  ((1.0 - beta) * T4);
//         } else if (isign == -1) {
//             L1 = (Kout * (q[(i * q_h2 + j) * QVAR + 3] - TARGET_PRESSURE)) -
//                  ((1.0 - beta) * T1);
//         }
//     } else {
//         throw std::runtime_error("Error: This BC is not yet implemented in characteristic form");
//     }

//     if (idir == 1) {
//         L1 /= (q[(i * q_h2 + j) * QVAR + 1] - qaux[(i * qa_h2 + j) * NQAUX + 0]);
//         L4 /= (q[(i * q_h2 + j) * QVAR + 1] + qaux[(i * qa_h2 + j) * NQAUX + 0]);
//         if (q[(i * q_h2 + j) * QVAR + 1] == 0.0) {
//             L2 = 0.0;
//             L3 = 0.0;
//         } else {
//             L2 /= q[(i * q_h2 + j) * QVAR + 1];
//             L3 /= q[(i * q_h2 + j) * QVAR + 1];
//         }
//     } else if (idir == 2) {
//         L1 /= (q[(i * q_h2 + j) * QVAR + 2] - qaux[(i * qa_h2 + j) * NQAUX + 0]);
//         L4 /= (q[(i * q_h2 + j) * QVAR + 2] + qaux[(i * qa_h2 + j) * NQAUX + 0]);
//         if (q[(i * q_h2 + j) * QVAR + 2] == 0.0) {
//             L2 = 0.0;
//             L3 = 0.0;
//         } else {
//             L2 /= q[(i * q_h2 + j) * QVAR + 2];
//             L3 /= q[(i * q_h2 + j) * QVAR + 2];
//         }
//     }
// }

// void compute_waves(lbm
//     int i, int j, int k, int idir, int isign,
//     int bc_type, const std::vector<double>& bc_params, const std::vector<double>& bc_target,
//     double T1, double T2, double T3, double T4, double T5,
//     double& L1, double& L2, double& L3, double& L4, double& L5,
//     double dp, double du, double dv, double dw, double drho,
//     const std::vector<std::vector<std::vector<std::vector<double>>>>& q,
//     const std::vector<std::vector<std::vector<std::vector<double>>>>& qaux,
//     int q_l1, int q_l2, int q_l3, int q_h1, int q_h2, int q_h3,
//     int qa_l1, int qa_l2, int qa_l3, int qa_h1, int qa_h2, int qa_h3)
// {
//     const int QVAR = 5;  // Placeholder for variable count
//     const int NQAUX = 5; // Placeholder for auxiliary variable count
//     const int QU = 0, QV = 1, QW = 2, QRHO = 3, QC = 4, QPRES = 5, QTEMP = 6, QRSPEC = 7;

//     // Validate idir and isign
//     if (idir < 1 || idir > 3) {
//         throw std::invalid_argument("Problem with idir in compute_waves");
//     }
//     if (isign != 1 && isign != -1) {
//         throw std::invalid_argument("Problem with isign in compute_waves");
//     }

//     double mach_local = std::sqrt(
//         std::pow(q[i][j][k][QU], 2.0) + std::pow(q[i][j][k][QV], 2.0) + std::pow(q[i][j][k][QW], 2.0)) /
//         qaux[i][j][k][QC];

//     // Recasting target values and numerical parameters
//     double TARGET_VX = bc_target[0];
//     double TARGET_VY = bc_target[1];
//     double TARGET_VZ = bc_target[2];
//     double TARGET_TEMPERATURE = bc_target[3];
//     double TARGET_PRESSURE = bc_target[4];

//     double relax_T = bc_params[0];
//     double relax_U = bc_params[1];
//     double relax_V = bc_params[2];
//     double relax_W = bc_params[3];
//     double beta = (bc_params[4] < 0.0) ? mach_local : bc_params[4];
//     double sigma_out = bc_params[5];

//     // Computing known numerical LODI waves
//     if (idir == 1) {
//         L1 = (q[i][j][k][QU] - qaux[i][j][k][QC]) * (dp - (q[i][j][k][QRHO] * qaux[i][j][k][QC]) * du);
//         L2 = q[i][j][k][QU] * (std::pow(qaux[i][j][k][QC], 2.0) * drho - dp);
//         L3 = q[i][j][k][QU] * dv;
//         L4 = q[i][j][k][QU] * dw;
//         L5 = (q[i][j][k][QU] + qaux[i][j][k][QC]) * (dp + (q[i][j][k][QRHO] * qaux[i][j][k][QC]) * du);
//     } else if (idir == 2) {
//         L1 = (q[i][j][k][QV] - qaux[i][j][k][QC]) * (dp - (q[i][j][k][QRHO] * qaux[i][j][k][QC]) * dv);
//         L2 = q[i][j][k][QV] * du;
//         L3 = q[i][j][k][QV] * (std::pow(qaux[i][j][k][QC], 2.0) * drho - dp);
//         L4 = q[i][j][k][QV] * dw;
//         L5 = (q[i][j][k][QV] + qaux[i][j][k][QC]) * (dp + (q[i][j][k][QRHO] * qaux[i][j][k][QC]) * dv);
//     } else if (idir == 3) {
//         L1 = (q[i][j][k][QW] - qaux[i][j][k][QC]) * (dp - (q[i][j][k][QRHO] * qaux[i][j][k][QC]) * dw);
//         L2 = q[i][j][k][QW] * du;
//         L3 = q[i][j][k][QW] * dv;
//         L4 = q[i][j][k][QW] * (std::pow(qaux[i][j][k][QC], 2.0) * drho - dp);
//         L5 = (q[i][j][k][QW] + qaux[i][j][k][QC]) * (dp + (q[i][j][k][QRHO] * qaux[i][j][k][QC]) * dw);
//     }

//     // Additional computations for specific boundary conditions (Inflow, SlipWall, Outflow, etc.)
//     if (bc_type == 1 /* Inflow */) {
//         // Add Inflow boundary-specific logic here
//     } else if (bc_type == 2 /* SlipWall */) {
//         L1 = L2 = L3 = L4 = L5 = 0.0;
//     } else if (bc_type == 3 /* Outflow */) {
//         double Kout = sigma_out * (1.0 - std::pow(mach_local, 2.0)) * (qaux[i][j][k][QC]);
//         if (isign == 1) {
//             L5 = Kout * (q[i][j][k][QPRES] - TARGET_PRESSURE) - (1.0 - beta) * T5;
//         } else {
//             L1 = Kout * (q[i][j][k][QPRES] - TARGET_PRESSURE) - (1.0 - beta) * T1;
//         }
//     } else {
//         throw std::runtime_error("Error: Unsupported boundary condition");
//     }

//     // Shaping the waves
//     if (idir == 1) {
//         L1 /= (q[i][j][k][QU] - qaux[i][j][k][QC]);
//         L5 /= (q[i][j][k][QU] + qaux[i][j][k][QC]);
//         if (q[i][j][k][QU] == 0.0) {
//             L2 = L3 = L4 = 0.0;
//         } else {
//             L2 /= q[i][j][k][QU];
//             L3 /= q[i][j][k][QU];
//             L4 /= q[i][j][k][QU];
//         }
//     }
//     // Similarly handle idir == 2 and idir == 3
// }