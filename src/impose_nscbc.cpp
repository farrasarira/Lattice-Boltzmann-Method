#include <iostream>
#include <omp.h>
#include "lbm.hpp"
#include "units.hpp"
#include "impose_nscbc.hpp"

void impose_NSCBC(LBM lb, int i, int j, int k, int l_interface, double &rho_out, double rhoa_out[], double vel_out[], double &T_out )
{
#if defined MULTICOMP
    size_t nSpecies = lb.get_nSpecies();

    // Initialize local variables for derivatives
    double dpdx = 0.0, dudx = 0.0, dvdx = 0.0, dwdx = 0.0, drhodx = 0.0, dYdx[nSpecies] = {};
    double dpdy = 0.0, dudy = 0.0, dvdy = 0.0, dwdy = 0.0, drhody = 0.0, dYdy[nSpecies] = {}; 
    double dpdz = 0.0, dudz = 0.0, dvdz = 0.0, dwdz = 0.0, drhodz = 0.0, dYdz[nSpecies] = {}; 
    int dx = lb.get_dx(); int dy = lb.get_dy(); int dz = lb.get_dz();

    double T1 = 0.0, T2 = 0.0, T3 = 0.0, T4 = 0.0, T5 = 0.0;
    double L1 = 0.0, L2 = 0.0, L3 = 0.0, L4 = 0.0, L5 = 0.0, L6[nSpecies] = {};
    
    if (cx[l_interface] != 0.0){
        normal_derivative(lb, i, j, k, 1, cx[l_interface], dx, dpdx, dudx, dvdx, dwdx, drhodx, dYdx, nSpecies);
        // std::cout << dudx << std::endl;

        // if (lb.get_NY() > 3)    
            // tangential_derivative(lb, i, j, k, 2, dy, dpdy, dudy, dvdy, dwdy, drhody, dYdy, nSpecies);
        // if (lb.get_NZ() > 3)    
            // tangential_derivative(lb, i, j, k, 3, dz, dpdz, dudz, dvdz, dwdz, drhodz, dYdz, nSpecies);
            
        // compute_tranverse_terms(lb, i, j, k, 1, T1, T2, T3, T4, T5, dpdx, dudx, dvdx, dwdx, drhodx, dpdy, dudy, dvdy, dwdy, drhody, dpdz, dudz, dvdz, dwdz, drhodz);
        compute_waves(lb, i, j, k, 1, cx[l_interface], T1, T2, T3, T4, T5, L1, L2, L3, L4, L5, L6, dpdx, dudx, dvdx, dwdx, drhodx, dYdx, dYdy, dYdz);
        update_bc_cells(lb, i, j, k, 1, cx[l_interface], L1, L2, L3, L4, L5, L6, rho_out, rhoa_out, vel_out, T_out);

    }
    else if (cy[l_interface] != 0.0){

    }
    else if (cz[l_interface] != 0.0){

    }
    
    

    
#endif

}


void normal_derivative(LBM lb, int i, int j, int k, int idir, int isign, double delta, double &dp, double &du, double &dv, double &dw, double &drho, double *dY, size_t nSpecies)
{
    if (idir == 1){
        if (isign == 1){
            dp = (-1.5*lb.mixture[i][j][k].p + 2.0*lb.mixture[i+1][j][k].p - 0.5*lb.mixture[i+2][j][k].p ) / delta;
            du = (-1.5*lb.mixture[i][j][k].u + 2.0*lb.mixture[i+1][j][k].u - 0.5*lb.mixture[i+2][j][k].u ) / delta;
            dv = (-1.5*lb.mixture[i][j][k].v + 2.0*lb.mixture[i+1][j][k].v - 0.5*lb.mixture[i+2][j][k].v ) / delta;
            dw = (-1.5*lb.mixture[i][j][k].w + 2.0*lb.mixture[i+1][j][k].w - 0.5*lb.mixture[i+2][j][k].w ) / delta;
            drho = (-1.5*lb.mixture[i][j][k].rho + 2.0*lb.mixture[i+1][j][k].rho - 0.5*lb.mixture[i+2][j][k].rho ) / delta;
            for (size_t a = 0; a < nSpecies; ++a)
                dY[a] = (-1.5*lb.species[a][i][j][k].rho/lb.mixture[i][j][k].rho + 2.0*lb.species[a][i+1][j][k].rho/lb.mixture[i+1][j][k].rho - 0.5*lb.species[a][i+2][j][k].rho/lb.mixture[i+2][j][k].rho ) / delta;
        }
        else if (isign == -1){
            dp = (1.5*lb.mixture[i][j][k].p - 2.0*lb.mixture[i-1][j][k].p + 0.5*lb.mixture[i-2][j][k].p ) / delta;
            du = (1.5*lb.mixture[i][j][k].u - 2.0*lb.mixture[i-1][j][k].u + 0.5*lb.mixture[i-2][j][k].u ) / delta;
            dv = (1.5*lb.mixture[i][j][k].v - 2.0*lb.mixture[i-1][j][k].v + 0.5*lb.mixture[i-2][j][k].v ) / delta;
            dw = (1.5*lb.mixture[i][j][k].w - 2.0*lb.mixture[i-1][j][k].w + 0.5*lb.mixture[i-2][j][k].w ) / delta;
            drho = (1.5*lb.mixture[i][j][k].rho - 2.0*lb.mixture[i-1][j][k].rho + 0.5*lb.mixture[i-2][j][k].rho ) / delta;
            for (size_t a = 0; a < nSpecies; ++a)
                dY[a] = (1.5*lb.species[a][i][j][k].rho/lb.mixture[i][j][k].rho - 2.0*lb.species[a][i-1][j][k].rho/lb.mixture[i-1][j][k].rho + 0.5*lb.species[a][i-2][j][k].rho/lb.mixture[i-2][j][k].rho ) / delta;
        }
    }
    else if (idir == 2){
        if (isign == 1){
            dp = (-1.5*lb.mixture[i][j][k].p + 2.0*lb.mixture[i][j+1][k].p - 0.5*lb.mixture[i][j+2][k].p ) / delta;
            du = (-1.5*lb.mixture[i][j][k].u + 2.0*lb.mixture[i][j+1][k].u - 0.5*lb.mixture[i][j+2][k].u ) / delta;
            dv = (-1.5*lb.mixture[i][j][k].v + 2.0*lb.mixture[i][j+1][k].v - 0.5*lb.mixture[i][j+2][k].v ) / delta;
            dw = (-1.5*lb.mixture[i][j][k].w + 2.0*lb.mixture[i][j+1][k].w - 0.5*lb.mixture[i][j+2][k].w ) / delta;
            drho = (-1.5*lb.mixture[i][j][k].rho + 2.0*lb.mixture[i][j+1][k].rho - 0.5*lb.mixture[i][j+2][k].rho ) / delta;
            for (size_t a = 0; a < nSpecies; ++a)
                dY[a] = (-1.5*lb.species[a][i][j][k].rho/lb.mixture[i][j][k].rho + 2.0*lb.species[a][i][j+1][k].rho/lb.mixture[i][j+1][k].rho - 0.5*lb.species[a][i][j+2][k].rho/lb.mixture[i][j+2][k].rho ) / delta;
        }
        else if (isign == -1){
            dp = (1.5*lb.mixture[i][j][k].p - 2.0*lb.mixture[i][j-1][k].p + 0.5*lb.mixture[i][j-2][k].p ) / delta;
            du = (1.5*lb.mixture[i][j][k].u - 2.0*lb.mixture[i][j-1][k].u + 0.5*lb.mixture[i][j-2][k].u ) / delta;
            dv = (1.5*lb.mixture[i][j][k].v - 2.0*lb.mixture[i][j-1][k].v + 0.5*lb.mixture[i][j-2][k].v ) / delta;
            dw = (1.5*lb.mixture[i][j][k].w - 2.0*lb.mixture[i][j-1][k].w + 0.5*lb.mixture[i][j-2][k].w ) / delta;
            drho = (1.5*lb.mixture[i][j][k].rho - 2.0*lb.mixture[i][j-1][k].rho + 0.5*lb.mixture[i][j-2][k].rho ) / delta;
            for (size_t a = 0; a < nSpecies; ++a)
                dY[a] = (1.5*lb.species[a][i][j][k].rho/lb.mixture[i][j][k].rho - 2.0*lb.species[a][i][j-1][k].rho/lb.mixture[i][j-1][k].rho + 0.5*lb.species[a][i][j-2][k].rho/lb.mixture[i][j-2][k].rho ) / delta;
        }
    }
    else if (idir == 3){
        if (isign == 1){
            dp = (-1.5*lb.mixture[i][j][k].p + 2.0*lb.mixture[i][j][k+1].p - 0.5*lb.mixture[i][j][k+2].p ) / delta;
            du = (-1.5*lb.mixture[i][j][k].u + 2.0*lb.mixture[i][j][k+1].u - 0.5*lb.mixture[i][j][k+2].u ) / delta;
            dv = (-1.5*lb.mixture[i][j][k].v + 2.0*lb.mixture[i][j][k+1].v - 0.5*lb.mixture[i][j][k+2].v ) / delta;
            dw = (-1.5*lb.mixture[i][j][k].w + 2.0*lb.mixture[i][j][k+1].w - 0.5*lb.mixture[i][j][k+2].w ) / delta;
            drho = (-1.5*lb.mixture[i][j][k].rho + 2.0*lb.mixture[i][j][k+1].rho - 0.5*lb.mixture[i][j][k+2].rho ) / delta;
            for (size_t a = 0; a < nSpecies; ++a)
                dY[a] = (-1.5*lb.species[a][i][j][k].rho/lb.mixture[i][j][k].rho + 2.0*lb.species[a][i][j][k+1].rho/lb.mixture[i][j][k+1].rho - 0.5*lb.species[a][i][j][k+2].rho/lb.mixture[i][j][k+2].rho ) / delta;
        }
        else if (isign == -1){
            dp = (1.5*lb.mixture[i][j][k].p - 2.0*lb.mixture[i][j][k-1].p + 0.5*lb.mixture[i][j][k-2].p ) / delta;
            du = (1.5*lb.mixture[i][j][k].u - 2.0*lb.mixture[i][j][k-1].u + 0.5*lb.mixture[i][j][k-2].u ) / delta;
            dv = (1.5*lb.mixture[i][j][k].v - 2.0*lb.mixture[i][j][k-1].v + 0.5*lb.mixture[i][j][k-2].v ) / delta;
            dw = (1.5*lb.mixture[i][j][k].w - 2.0*lb.mixture[i][j][k-1].w + 0.5*lb.mixture[i][j][k-2].w ) / delta;
            drho = (1.5*lb.mixture[i][j][k].rho - 2.0*lb.mixture[i][j][k-1].rho + 0.5*lb.mixture[i][j][k-2].rho ) / delta;
            for (size_t a = 0; a < nSpecies; ++a)
                dY[a] = (1.5*lb.species[a][i][j][k].rho/lb.mixture[i][j][k].rho - 2.0*lb.species[a][i][j][k-1].rho/lb.mixture[i][j][k-1].rho + 0.5*lb.species[a][i][j][k-2].rho/lb.mixture[i][j][k-2].rho ) / delta;
        }
    }

}

void tangential_derivative(LBM lb, int i, int j, int k, int idir, double delta, double &dp, double &du, double &dv, double &dw, double &drho, double *dY, size_t nSpecies)
{
    if (idir == 1){
        if (lb.mixture[i-1][j][k].type != TYPE_F){
            dp = (-1.5*lb.mixture[i][j][k].p + 2.0*lb.mixture[i+1][j][k].p - 0.5*lb.mixture[i+2][j][k].p ) / delta;
            du = (-1.5*lb.mixture[i][j][k].u + 2.0*lb.mixture[i+1][j][k].u - 0.5*lb.mixture[i+2][j][k].u ) / delta;
            dv = (-1.5*lb.mixture[i][j][k].v + 2.0*lb.mixture[i+1][j][k].v - 0.5*lb.mixture[i+2][j][k].v ) / delta;
            dw = (-1.5*lb.mixture[i][j][k].w + 2.0*lb.mixture[i+1][j][k].w - 0.5*lb.mixture[i+2][j][k].w ) / delta;
            drho = (-1.5*lb.mixture[i][j][k].rho + 2.0*lb.mixture[i+1][j][k].rho - 0.5*lb.mixture[i+2][j][k].rho ) / delta;
            for (size_t a = 0; a < nSpecies; ++a)
                dY[a] = (-1.5*lb.species[a][i][j][k].rho/lb.mixture[i][j][k].rho + 2.0*lb.species[a][i+1][j][k].rho/lb.mixture[i+1][j][k].rho - 0.5*lb.species[a][i+2][j][k].rho/lb.mixture[i+2][j][k].rho ) / delta;
        }
        else if (lb.mixture[i+1][j][k].type != TYPE_F){
            dp = (1.5*lb.mixture[i][j][k].p - 2.0*lb.mixture[i-1][j][k].p + 0.5*lb.mixture[i-2][j][k].p ) / delta;
            du = (1.5*lb.mixture[i][j][k].u - 2.0*lb.mixture[i-1][j][k].u + 0.5*lb.mixture[i-2][j][k].u ) / delta;
            dv = (1.5*lb.mixture[i][j][k].v - 2.0*lb.mixture[i-1][j][k].v + 0.5*lb.mixture[i-2][j][k].v ) / delta;
            dw = (1.5*lb.mixture[i][j][k].w - 2.0*lb.mixture[i-1][j][k].w + 0.5*lb.mixture[i-2][j][k].w ) / delta;
            drho = (1.5*lb.mixture[i][j][k].rho - 2.0*lb.mixture[i-1][j][k].rho + 0.5*lb.mixture[i-2][j][k].rho ) / delta;
            for (size_t a = 0; a < nSpecies; ++a)
                dY[a] = (1.5*lb.species[a][i][j][k].rho/lb.mixture[i][j][k].rho - 2.0*lb.species[a][i-1][j][k].rho/lb.mixture[i-1][j][k].rho + 0.5*lb.species[a][i-2][j][k].rho/lb.mixture[i-2][j][k].rho ) / delta;
        }
        else{
            dp = (lb.mixture[i+1][j][k].p - lb.mixture[i-1][j][k].p ) / (2.0 * delta);
            du = (lb.mixture[i+1][j][k].u - lb.mixture[i-1][j][k].u ) / (2.0 * delta);
            dv = (lb.mixture[i+1][j][k].v - lb.mixture[i-1][j][k].v ) / (2.0 * delta);
            dw = (lb.mixture[i+1][j][k].w - lb.mixture[i-1][j][k].w ) / (2.0 * delta);
            drho = (lb.mixture[i+1][j][k].rho - lb.mixture[i-1][j][k].rho ) / (2.0 * delta);
            for (size_t a = 0; a < nSpecies; ++a)
                dY[a] = (lb.species[a][i+1][j][k].rho/lb.mixture[i+1][j][k].rho - lb.species[a][i-1][j][k].rho/lb.mixture[i-1][j][k].rho ) / (2.0 * delta);
        }
    }
    else if (idir == 2){
        if (lb.mixture[i][j-1][k].type != TYPE_F){
            dp = (-1.5*lb.mixture[i][j][k].p + 2.0*lb.mixture[i][j+1][k].p - 0.5*lb.mixture[i][j+2][k].p ) / delta;
            du = (-1.5*lb.mixture[i][j][k].u + 2.0*lb.mixture[i][j+1][k].u - 0.5*lb.mixture[i][j+2][k].u ) / delta;
            dv = (-1.5*lb.mixture[i][j][k].v + 2.0*lb.mixture[i][j+1][k].v - 0.5*lb.mixture[i][j+2][k].v ) / delta;
            dw = (-1.5*lb.mixture[i][j][k].w + 2.0*lb.mixture[i][j+1][k].w - 0.5*lb.mixture[i][j+2][k].w ) / delta;
            drho = (-1.5*lb.mixture[i][j][k].rho + 2.0*lb.mixture[i][j+1][k].rho - 0.5*lb.mixture[i][j+2][k].rho ) / delta;
            for (size_t a = 0; a < nSpecies; ++a)
                dY[a] = (-1.5*lb.species[a][i][j][k].rho/lb.mixture[i][j][k].rho + 2.0*lb.species[a][i][j+1][k].rho/lb.mixture[i][j+1][k].rho - 0.5*lb.species[a][i][j+2][k].rho/lb.mixture[i][j+2][k].rho ) / delta;
        }
        else if (lb.mixture[i][j+1][k].type != TYPE_F){
            dp = (1.5*lb.mixture[i][j][k].p - 2.0*lb.mixture[i][j-1][k].p + 0.5*lb.mixture[i][j-2][k].p ) / delta;
            du = (1.5*lb.mixture[i][j][k].u - 2.0*lb.mixture[i][j-1][k].u + 0.5*lb.mixture[i][j-2][k].u ) / delta;
            dv = (1.5*lb.mixture[i][j][k].v - 2.0*lb.mixture[i][j-1][k].v + 0.5*lb.mixture[i][j-2][k].v ) / delta;
            dw = (1.5*lb.mixture[i][j][k].w - 2.0*lb.mixture[i][j-1][k].w + 0.5*lb.mixture[i][j-2][k].w ) / delta;
            drho = (1.5*lb.mixture[i][j][k].rho - 2.0*lb.mixture[i][j-1][k].rho + 0.5*lb.mixture[i][j-2][k].rho ) / delta;
            for (size_t a = 0; a < nSpecies; ++a)
                dY[a] = (1.5*lb.species[a][i][j][k].rho/lb.mixture[i][j][k].rho - 2.0*lb.species[a][i][j-1][k].rho/lb.mixture[i][j-1][k].rho + 0.5*lb.species[a][i][j-2][k].rho/lb.mixture[i][j-2][k].rho ) / delta;
        }
        else {
            dp = (lb.mixture[i][j+1][k].p - lb.mixture[i][j-1][k].p ) / (2.0 * delta);
            du = (lb.mixture[i][j+1][k].u - lb.mixture[i][j-1][k].u ) / (2.0 * delta);
            dv = (lb.mixture[i][j+1][k].v - lb.mixture[i][j-1][k].v ) / (2.0 * delta);
            dw = (lb.mixture[i][j+1][k].w - lb.mixture[i][j-1][k].w ) / (2.0 * delta);
            drho = (lb.mixture[i][j+1][k].rho - lb.mixture[i][j-1][k].rho ) / (2.0 * delta);
            for (size_t a = 0; a < nSpecies; ++a)
                dY[a] = (lb.species[a][i][j+1][k].rho/lb.mixture[i][j+1][k].rho - lb.species[a][i][j-1][k].rho/lb.mixture[i][j-1][k].rho ) / (2.0 * delta);
        }
    }
    else if (idir == 3){
        if (lb.mixture[i][j][k-1].type != TYPE_F){
            dp = (-1.5*lb.mixture[i][j][k].p + 2.0*lb.mixture[i][j][k+1].p - 0.5*lb.mixture[i][j][k+2].p ) / delta;
            du = (-1.5*lb.mixture[i][j][k].u + 2.0*lb.mixture[i][j][k+1].u - 0.5*lb.mixture[i][j][k+2].u ) / delta;
            dv = (-1.5*lb.mixture[i][j][k].v + 2.0*lb.mixture[i][j][k+1].v - 0.5*lb.mixture[i][j][k+2].v ) / delta;
            dw = (-1.5*lb.mixture[i][j][k].w + 2.0*lb.mixture[i][j][k+1].w - 0.5*lb.mixture[i][j][k+2].w ) / delta;
            drho = (-1.5*lb.mixture[i][j][k].rho + 2.0*lb.mixture[i][j][k+1].rho - 0.5*lb.mixture[i][j][k+2].rho ) / delta;
            for (size_t a = 0; a < nSpecies; ++a)
                dY[a] = (-1.5*lb.species[a][i][j][k].rho/lb.mixture[i][j][k].rho + 2.0*lb.species[a][i][j][k+1].rho/lb.mixture[i][j][k+1].rho - 0.5*lb.species[a][i][j][k+2].rho/lb.mixture[i][j][k+2].rho ) / delta;
        }
        else if (lb.mixture[i][j][k-1].type != TYPE_F){
            dp = (1.5*lb.mixture[i][j][k].p - 2.0*lb.mixture[i][j][k-1].p + 0.5*lb.mixture[i][j][k-2].p ) / delta;
            du = (1.5*lb.mixture[i][j][k].u - 2.0*lb.mixture[i][j][k-1].u + 0.5*lb.mixture[i][j][k-2].u ) / delta;
            dv = (1.5*lb.mixture[i][j][k].v - 2.0*lb.mixture[i][j][k-1].v + 0.5*lb.mixture[i][j][k-2].v ) / delta;
            dw = (1.5*lb.mixture[i][j][k].w - 2.0*lb.mixture[i][j][k-1].w + 0.5*lb.mixture[i][j][k-2].w ) / delta;
            drho = (1.5*lb.mixture[i][j][k].rho - 2.0*lb.mixture[i][j][k-1].rho + 0.5*lb.mixture[i][j][k-2].rho ) / delta;
            for (size_t a = 0; a < nSpecies; ++a)
                dY[a] = (1.5*lb.species[a][i][j][k].rho/lb.mixture[i][j][k].rho - 2.0*lb.species[a][i][j][k-1].rho/lb.mixture[i][j][k-1].rho + 0.5*lb.species[a][i][j][k-2].rho/lb.mixture[i][j][k-2].rho ) / delta;
        }
        else {
            dp = (lb.mixture[i][j][k+1].p - lb.mixture[i][j][k-1].p ) / (2.0 * delta);
            du = (lb.mixture[i][j][k+1].u - lb.mixture[i][j][k-1].u ) / (2.0 * delta);
            dv = (lb.mixture[i][j][k+1].v - lb.mixture[i][j][k-1].v ) / (2.0 * delta);
            dw = (lb.mixture[i][j][k+1].w - lb.mixture[i][j][k-1].w ) / (2.0 * delta);
            drho = (lb.mixture[i][j][k+1].rho - lb.mixture[i][j][k-1].rho ) / (2.0 * delta);
            for (size_t a = 0; a < nSpecies; ++a)
                dY[a] = (lb.species[a][i][j][k+1].rho/lb.mixture[i][j][k+1].rho - lb.species[a][i][j][k-1].rho/lb.mixture[i][j][k-1].rho ) / (2.0 * delta);
        }
    }

}

void compute_tranverse_terms(LBM lb, int i, int j, int k, int idir, double& T1, double& T2, double& T3, double& T4, double& T5,
double dpdx, double dudx, double dvdx, double dwdx, double drhodx,
    double dpdy, double dudy, double dvdy, double dwdy, double drhody,
    double dpdz, double dudz, double dvdz, double dwdz, double drhodz)
{
    size_t nSpecies = lb.get_nSpecies();
    std::vector<std::string> speciesName = lb.get_speciesName();

    int rank = omp_get_thread_num();
    auto gas = sols[rank]->thermo();   
    std::vector <double> Y (gas->nSpecies());
    for(size_t a = 0; a < nSpecies; ++a) Y[gas->speciesIndex(speciesName[a])] = lb.species[a][i][j][k].rho / lb.mixture[i][j][k].rho;
    gas->setMassFractions(&Y[0]);
    gas->setState_DP(units.si_rho(lb.mixture[i][j][k].rho), units.si_p(lb.mixture[i][j][k].p));

    double soundSpeed = units.u(gas->soundSpeed());
    double gamma = gas->cp_mass() / gas->cv_mass();

    double inv_rho = 1.0 / lb.mixture[i][j][k].rho;

    if (idir == 1) {
        T1 = (lb.mixture[i][j][k].v * (dpdy - lb.mixture[i][j][k].rho * soundSpeed * dudy)) +
             (lb.mixture[i][j][k].w * (dpdz - lb.mixture[i][j][k].rho * soundSpeed * dudz)) +
             (gamma * lb.mixture[i][j][k].p * (dvdy + dwdz));

        T2 = (lb.mixture[i][j][k].v * ((soundSpeed * soundSpeed * drhody) - dpdy)) +
             (lb.mixture[i][j][k].w * ((soundSpeed * soundSpeed * drhodz) - dpdz));

        T3 = lb.mixture[i][j][k].v * dvdy + lb.mixture[i][j][k].w * dvdz + dpdy * inv_rho;

        T4 = lb.mixture[i][j][k].v * dwdy + lb.mixture[i][j][k].w * dwdz + dpdz * inv_rho;

        T5 = (lb.mixture[i][j][k].v * (dpdy + lb.mixture[i][j][k].rho * soundSpeed * dudy)) +
             (lb.mixture[i][j][k].w * (dpdz + lb.mixture[i][j][k].rho * soundSpeed * dudz)) +
             (gamma * lb.mixture[i][j][k].p * (dvdy + dwdz));

    } else if (idir == 2) {
        T1 = (lb.mixture[i][j][k].u * (dpdx - lb.mixture[i][j][k].rho * soundSpeed * dvdx)) +
             (lb.mixture[i][j][k].w * (dpdz - lb.mixture[i][j][k].rho * soundSpeed * dvdz)) +
             (gamma * lb.mixture[i][j][k].p * (dudx + dwdz));

        T2 = lb.mixture[i][j][k].u * dudx + lb.mixture[i][j][k].w * dudz + dpdx * inv_rho;

        T3 = (lb.mixture[i][j][k].u * ((soundSpeed * soundSpeed * drhodx) - dpdx)) +
             (lb.mixture[i][j][k].w * ((soundSpeed * soundSpeed * drhodz) - dpdz));

        T4 = lb.mixture[i][j][k].u * dwdx + lb.mixture[i][j][k].w * dwdz + dpdz * inv_rho;

        T5 = (lb.mixture[i][j][k].u * (dpdx + lb.mixture[i][j][k].rho * soundSpeed * dvdx)) +
             (lb.mixture[i][j][k].w * (dpdz + lb.mixture[i][j][k].rho * soundSpeed * dvdz)) +
             (gamma * lb.mixture[i][j][k].p * (dudx + dwdz));

    } else if (idir == 3) {
        T1 = (lb.mixture[i][j][k].u * (dpdx - lb.mixture[i][j][k].rho * soundSpeed * dwdx)) +
             (lb.mixture[i][j][k].v * (dpdy - lb.mixture[i][j][k].rho * soundSpeed * dwdy)) +
             (gamma * lb.mixture[i][j][k].p * (dudx + dvdy));

        T2 = lb.mixture[i][j][k].u * dudx + lb.mixture[i][j][k].v * dudy + dpdx * inv_rho;

        T3 = lb.mixture[i][j][k].u * dvdx + lb.mixture[i][j][k].v * dvdy + dpdy * inv_rho;

        T4 = (lb.mixture[i][j][k].u * ((soundSpeed * soundSpeed * drhodx) - dpdx)) +
             (lb.mixture[i][j][k].v * ((soundSpeed * soundSpeed * drhody) - dpdy));

        T5 = (lb.mixture[i][j][k].u * (dpdx + lb.mixture[i][j][k].rho * soundSpeed * dwdx)) +
             (lb.mixture[i][j][k].v * (dpdy + lb.mixture[i][j][k].rho * soundSpeed * dwdy)) +
             (gamma * lb.mixture[i][j][k].p * (dudx + dvdy));

    } else {
        throw std::invalid_argument("Invalid idir in compute_transverse_terms");
    }
}

void compute_waves(LBM lb,
    int i, int j, int k, int idir, int isign,
    double T1, double T2, double T3, double T4, double T5,
    double& L1, double& L2, double& L3, double& L4, double& L5, double *L6,
    double dp, double du, double dv, double dw, double drho, double dYdx[], double dYdy[], double dYdz[])
{

    // Validate idir and isign
    if (idir < 1 || idir > 3) {
        throw std::invalid_argument("Problem with idir in compute_waves");
    }
    if (isign != 1 && isign != -1) {
        throw std::invalid_argument("Problem with isign in compute_waves");
    }

    size_t nSpecies = lb.get_nSpecies();
    std::vector<std::string> speciesName = lb.get_speciesName();

    int rank = omp_get_thread_num();
    auto gas = sols[rank]->thermo();   
    std::vector <double> Y (gas->nSpecies());
    for(size_t a = 0; a < nSpecies; ++a) Y[gas->speciesIndex(speciesName[a])] = lb.species[a][i][j][k].rho / lb.mixture[i][j][k].rho;
    gas->setMassFractions(&Y[0]);
    gas->setState_DP(units.si_rho(lb.mixture[i][j][k].rho), units.si_p(lb.mixture[i][j][k].p));

    double soundSpeed = units.u(gas->soundSpeed());
    double mach = v_mag(lb.mixture[i][j][k].u, lb.mixture[i][j][k].v, lb.mixture[i][j][k].w) / soundSpeed;

    // Recasting target values and numerical parameters
    double TARGET_VX = 0.0;
    double TARGET_VY = 0.0;
    double TARGET_VZ = 0.0;
    double TARGET_TEMPERATURE = 0.0;
    double TARGET_PRESSURE = 0.0;

    double relax_T = 0.0;
    double relax_U = 0.0;
    double relax_V = 0.0;
    double relax_W = 0.0;
    double sigma_out = 0.0;//0.25;

    int bc_type = 0.0;
    int len_char = 0.0;

    if (idir == 1) {
        TARGET_VX = lb.mixture[i-isign][j][k].u;
        TARGET_VY = lb.mixture[i-isign][j][k].v;
        TARGET_VZ = lb.mixture[i-isign][j][k].w;
        TARGET_TEMPERATURE = lb.mixture[i-isign][j][k].temp;
        TARGET_PRESSURE = lb.mixture[i-isign][j][k].p;
        bc_type = lb.mixture[i-isign][j][k].type;
        len_char = lb.get_Nx();
    } else if (idir == 2) {
        TARGET_VX = lb.mixture[i][j-isign][k].u;
        TARGET_VY = lb.mixture[i][j-isign][k].v;
        TARGET_VZ = lb.mixture[i][j-isign][k].w;
        TARGET_TEMPERATURE = lb.mixture[i][j-isign][k].temp;
        TARGET_PRESSURE = lb.mixture[i][j-isign][k].p;
        bc_type = lb.mixture[i][j-isign][k].type;
        len_char = lb.get_Ny();
    } else if (idir == 3) {
        TARGET_VX = lb.mixture[i][j][k-isign].u;
        TARGET_VY = lb.mixture[i][j][k-isign].v;
        TARGET_VZ = lb.mixture[i][j][k-isign].w;
        TARGET_TEMPERATURE = lb.mixture[i][j][k-isign].temp;
        TARGET_PRESSURE = lb.mixture[i][j][k-isign].p;
        bc_type = lb.mixture[i][j][k-isign].type;
        len_char = lb.get_Nz();
    }

    double beta = mach;

    // Computing known numerical LODI waves
    if (idir == 1) {
        L1 = (lb.mixture[i][j][k].u - soundSpeed) * (dp - (lb.mixture[i][j][k].rho * soundSpeed) * du);
        L2 = lb.mixture[i][j][k].u * (std::pow(soundSpeed, 2.0) * drho - dp);
        L3 = lb.mixture[i][j][k].u * dv;
        L4 = lb.mixture[i][j][k].u * dw;
        L5 = (lb.mixture[i][j][k].u + soundSpeed) * (dp + (lb.mixture[i][j][k].rho * soundSpeed) * du);
    } else if (idir == 2) {
        L1 = (lb.mixture[i][j][k].v - soundSpeed) * (dp - (lb.mixture[i][j][k].rho * soundSpeed) * dv);
        L2 = lb.mixture[i][j][k].v * du;
        L3 = lb.mixture[i][j][k].v * (std::pow(soundSpeed, 2.0) * drho - dp);
        L4 = lb.mixture[i][j][k].v * dw;
        L5 = (lb.mixture[i][j][k].v + soundSpeed) * (dp + (lb.mixture[i][j][k].rho * soundSpeed) * dv);
    } else if (idir == 3) {
        L1 = (lb.mixture[i][j][k].w - soundSpeed) * (dp - (lb.mixture[i][j][k].rho * soundSpeed) * dw);
        L2 = lb.mixture[i][j][k].w * du;
        L3 = lb.mixture[i][j][k].w * dv;
        L4 = lb.mixture[i][j][k].w * (std::pow(soundSpeed, 2.0) * drho - dp);
        L5 = (lb.mixture[i][j][k].w + soundSpeed) * (dp + (lb.mixture[i][j][k].rho * soundSpeed) * dw);
    }

    for (size_t a = 0; a < nSpecies; ++a)
        L6[a] = lb.mixture[i][j][k].u * dYdx[a];// + lb.mixture[i][j][k].v * dYdy[a] + lb.mixture[i][j][k].w * dYdz[a];

    // Additional computations for specific boundary conditions (Inflow, SlipWall, Outflow, etc.)
    if (bc_type == TYPE_I_C /* Inflow */) {
        // Add Inflow boundary-specific logic here
    } else if (bc_type == TYPE_O_C /* Outflow */) {
        double Kout = sigma_out * (1.0 - std::pow(mach, 2.0)) * (soundSpeed / len_char);
        if (isign == 1) {
            L5 = Kout * (lb.mixture[i][j][k].p - TARGET_PRESSURE) - (1.0 - beta) * T5;
        } else {
            L1 = Kout * (lb.mixture[i][j][k].p - TARGET_PRESSURE) - (1.0 - beta) * T1;
        }
    } else {
        throw std::runtime_error("Error: Unsupported boundary condition");
    }

    // // Shaping the waves
    // if (idir == 1) {
    //     L1 /= (lb.mixture[i][j][k].u - soundSpeed);
    //     L5 /= (lb.mixture[i][j][k].u + soundSpeed);
    //     if (lb.mixture[i][j][k].u == 0.0) {
    //         L2 = L3 = L4 = 0.0;
    //     } else {
    //         L2 /= lb.mixture[i][j][k].u;
    //         L3 /= lb.mixture[i][j][k].u;
    //         L4 /= lb.mixture[i][j][k].u;
    //     }
    // }
    // else if (idir == 2) {
    //     L1 /= (lb.mixture[i][j][k].v - soundSpeed);
    //     L5 /= (lb.mixture[i][j][k].v + soundSpeed);
    //     if (lb.mixture[i][j][k].v == 0.0) {
    //         L2 = L3 = L4 = 0.0;
    //     } else {
    //         L2 /= lb.mixture[i][j][k].v;
    //         L3 /= lb.mixture[i][j][k].v;
    //         L4 /= lb.mixture[i][j][k].v;
    //     }
    // }
    // else if (idir == 3) {
    //     L1 /= (lb.mixture[i][j][k].w - soundSpeed);
    //     L5 /= (lb.mixture[i][j][k].w + soundSpeed);
    //     if (lb.mixture[i][j][k].w == 0.0) {
    //         L2 = L3 = L4 = 0.0;
    //     } else {
    //         L2 /= lb.mixture[i][j][k].w;
    //         L3 /= lb.mixture[i][j][k].w;
    //         L4 /= lb.mixture[i][j][k].w;
    //     }
    // }
}

void update_bc_cells(LBM lb,
    int i, int j, int k, int idir, int isign,
    double L1, double L2, double L3, double L4, double L5, double L6[],
    double &rho_out, double rhoa_out[], double vel_out[], double &T_out)
{

    // Validate idir and isign
    if (idir < 1 || idir > 3) {
        throw std::invalid_argument("Problem with idir in compute_waves");
    }
    if (isign != 1 && isign != -1) {
        throw std::invalid_argument("Problem with isign in compute_waves");
    }

    size_t nSpecies = lb.get_nSpecies();
    std::vector<std::string> speciesName = lb.get_speciesName();

    int rank = omp_get_thread_num();
    auto gas = sols[rank]->thermo();   
    std::vector <double> Y (gas->nSpecies());
    for(size_t a = 0; a < nSpecies; ++a) Y[gas->speciesIndex(speciesName[a])] = lb.species[a][i][j][k].rho / lb.mixture[i][j][k].rho;
    gas->setMassFractions(&Y[0]);
    gas->setState_DP(units.si_rho(lb.mixture[i][j][k].rho), units.si_p(lb.mixture[i][j][k].p));

    double soundSpeed = units.u(gas->soundSpeed());

    double drho = 0.0, du = 0.0, dv = 0.0, dw = 0.0, dp = 0.0;

    double dt_sim = lb.get_dtsim();

    if (idir == 1) {
        drho = (L2 + 0.5 * (L1 + L5)) / (soundSpeed * soundSpeed);
        du = (L5 - L1) / (2.0 * soundSpeed * lb.mixture[i][j][k].rho);
        dv = L3;
        dw = L4;
        dp = 0.5 * (L1 + L5);
    } else if (idir == 2) {
        drho = (L3 + 0.5 * (L1 + L5)) / (soundSpeed * soundSpeed);
        du = L2;
        dv = (L5 - L1) / (2.0 * soundSpeed * lb.mixture[i][j][k].rho);
        dw = L4;
        dp = 0.5 * (L1 + L5);
    } else if (idir == 3) {
        drho = (L4 + 0.5 * (L1 + L5)) / (soundSpeed * soundSpeed);
        du = L2;
        dv = L3;
        dw = (L5 - L1) / (2.0 * soundSpeed * lb.mixture[i][j][k].rho);
        dp = 0.5 * (L1 + L5);
    }

    rho_out = lb.mixture[i][j][k].rho - dt_sim * drho;
    vel_out[0] = lb.mixture[i][j][k].u - dt_sim * du;
    vel_out[1] = lb.mixture[i][j][k].v - dt_sim * dv;
    vel_out[2] = lb.mixture[i][j][k].w - dt_sim * dw;
    double p_out = lb.mixture[i][j][k].p - dt_sim * dp;

    for(size_t a = 0; a < nSpecies; ++a)
        rhoa_out[a] = (lb.species[a][i][j][k].rho/lb.mixture[i][j][k].rho - dt_sim*L6[a]) * rho_out ;
    
    std::fill(Y.begin(), Y.end(), 0.0);
    for(size_t a = 0; a < nSpecies; ++a) Y[gas->speciesIndex(speciesName[a])] = rhoa_out[a] / rho_out;
    gas->setMassFractions(&Y[0]);
    gas->setState_DP(units.si_rho(rho_out), units.si_p(p_out));

    T_out = units.temp(gas->temperature());


}