#ifndef IMPOSE_NSCBC
#define IMPOASE_NSCBC

#include "lbm.hpp"

void impose_NSCBC(LBM lb, int i, int j, int k, int l_interface, double &rho_out, double rhoa_out[], double vel_out[], double &T_out );
void normal_derivative(LBM lb, int i, int j, int k, int idir, int isign, double delta, double &dp, double &du, double &dv, double &dw, double &drho, double *dY, size_t nSpecies);
void tangential_derivative(LBM lb, int i, int j, int k, int idir, double delta, double &dp, double &du, double &dv, double &dw, double &drho, double *dY, size_t nSpecies);
void compute_tranverse_terms(LBM lb, int i, int j, int k, int idir, double& T1, double& T2, double& T3, double& T4, double& T5,
double dpdx, double dudx, double dvdx, double dwdx, double drhodx,
    double dpdy, double dudy, double dvdy, double dwdy, double drhody,
    double dpdz, double dudz, double dvdz, double dwdz, double drhodz);
void compute_waves(LBM lb,
    int i, int j, int k, int idir, int isign,
    double T1, double T2, double T3, double T4, double T5,
    double& L1, double& L2, double& L3, double& L4, double& L5, double *L6,
    double dp, double du, double dv, double dw, double drho, double dYdx[], double dYdy[], double dYdz[]);
void update_bc_cells(LBM lb,
    int i, int j, int k, int idir, int isign,
    double L1, double L2, double L3, double L4, double L5, double L6[],
    double &rho_out, double rhoa_out[], double vel_out[], double &T_out);

#endif