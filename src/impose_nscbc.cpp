#include <iostream>
#include <omp.h>
#include "lbm.hpp"
#include "units.hpp"
#include "impose_nscbc.hpp"

// ============================================================================
// Navier-Stokes Characteristic Boundary Conditions (NSCBC) for multicomponent
// mixtures.
//
// Formulation:
//  - LODI wave framework: T. Poinsot & S. Lele, "Boundary conditions for
//    direct simulations of compressible viscous flows", JCP 101 (1992).
//  - Multicomponent extension (species advection waves, closure through the
//    mixture equation of state): M. Baum, T. Poinsot & D. Thevenin, "Accurate
//    boundary conditions for multicomponent reactive flows", JCP 116 (1995).
//  - Partially non-reflecting outflow with pressure relaxation
//    K = sigma (1-M^2) c / L  (Rudy & Strikwerda constant, sigma = 0.25 as
//    recommended by Poinsot & Lele) and transverse-term damping with
//    beta = Mach: C.S. Yoo & H.G. Im, Combust. Theory Modelling 11 (2007).
//
// Wave convention (idir = 1:x, 2:y, 3:z; un = boundary-normal velocity,
// c = frozen mixture sound speed from Cantera; d/dn = one-sided derivative
// along the boundary normal, pointing INTO the domain for isign=+1):
//   L1    = (un - c) (dp/dn - rho c dun/dn)   acoustic wave at speed un-c
//   Lt    = un dut/dn                         tangential velocity waves
//   Ls    = un (c^2 drho/dn - dp/dn)          entropy wave
//   L5    = (un + c) (dp/dn + rho c dun/dn)   acoustic wave at speed un+c
//   L6[a] = un dYa/dn                         species waves (Baum et al.)
// The boundary cell is advanced with the LODI relations
//   dp/dt   = -1/2 (L1 + L5)
//   dun/dt  = -(L5 - L1)/(2 rho c)
//   drho/dt = -[Ls + 1/2 (L1 + L5)]/c^2
//   dut/dt  = -Lt          dYa/dt = -L6[a]
// and the temperature is then recovered from (rho, p, Y) through the mixture
// equation of state (Cantera), which closes the multicomponent system in the
// sense of Baum et al. (1995).
//
// isign = +1: boundary on the minus side of the domain (outgoing wave: L1)
// isign = -1: boundary on the plus  side of the domain (outgoing wave: L5)
// ============================================================================

#ifdef MULTICOMP

void impose_NSCBC(LBM& lb, int i, int j, int k, int l_interface, double &rho_out, double rhoa_out[], double vel_out[], double &T_out )
{
    size_t nSpecies = lb.get_nSpecies();

    // derivatives along each coordinate direction (normal + transverse)
    double dpdx = 0.0, dudx = 0.0, dvdx = 0.0, dwdx = 0.0, drhodx = 0.0;
    double dpdy = 0.0, dudy = 0.0, dvdy = 0.0, dwdy = 0.0, drhody = 0.0;
    double dpdz = 0.0, dudz = 0.0, dvdz = 0.0, dwdz = 0.0, drhodz = 0.0;
    std::vector<double> dYdx(nSpecies, 0.0), dYdy(nSpecies, 0.0), dYdz(nSpecies, 0.0);
    int dx = lb.get_dx(); int dy = lb.get_dy(); int dz = lb.get_dz();

    double T1 = 0.0, T2 = 0.0, T3 = 0.0, T4 = 0.0, T5 = 0.0;
    double L1 = 0.0, L2 = 0.0, L3 = 0.0, L4 = 0.0, L5 = 0.0;
    std::vector<double> L6(nSpecies, 0.0);

    if (cx[l_interface] != 0){
        // x-normal boundary
        normal_derivative(lb, i, j, k, 1, cx[l_interface], dx, dpdx, dudx, dvdx, dwdx, drhodx, dYdx.data(), nSpecies);

        // transverse terms (Yoo & Im) only when the domain is resolved in the
        // tangential directions
        if (lb.get_NY() > 3)
            tangential_derivative(lb, i, j, k, 2, dy, dpdy, dudy, dvdy, dwdy, drhody, dYdy.data(), nSpecies);
        if (lb.get_NZ() > 3)
            tangential_derivative(lb, i, j, k, 3, dz, dpdz, dudz, dvdz, dwdz, drhodz, dYdz.data(), nSpecies);
        if (lb.get_NY() > 3 || lb.get_NZ() > 3)
            compute_tranverse_terms(lb, i, j, k, 1, T1, T2, T3, T4, T5, dpdx, dudx, dvdx, dwdx, drhodx, dpdy, dudy, dvdy, dwdy, drhody, dpdz, dudz, dvdz, dwdz, drhodz);

        compute_waves(lb, i, j, k, 1, cx[l_interface], T1, T2, T3, T4, T5, L1, L2, L3, L4, L5, L6.data(), dpdx, dudx, dvdx, dwdx, drhodx, dYdx.data());
        update_bc_cells(lb, i, j, k, 1, cx[l_interface], L1, L2, L3, L4, L5, L6.data(), rho_out, rhoa_out, vel_out, T_out);
    }
    else if (cy[l_interface] != 0){
        // y-normal boundary
        normal_derivative(lb, i, j, k, 2, cy[l_interface], dy, dpdy, dudy, dvdy, dwdy, drhody, dYdy.data(), nSpecies);

        if (lb.get_NX() > 3)
            tangential_derivative(lb, i, j, k, 1, dx, dpdx, dudx, dvdx, dwdx, drhodx, dYdx.data(), nSpecies);
        if (lb.get_NZ() > 3)
            tangential_derivative(lb, i, j, k, 3, dz, dpdz, dudz, dvdz, dwdz, drhodz, dYdz.data(), nSpecies);
        if (lb.get_NX() > 3 || lb.get_NZ() > 3)
            compute_tranverse_terms(lb, i, j, k, 2, T1, T2, T3, T4, T5, dpdx, dudx, dvdx, dwdx, drhodx, dpdy, dudy, dvdy, dwdy, drhody, dpdz, dudz, dvdz, dwdz, drhodz);

        compute_waves(lb, i, j, k, 2, cy[l_interface], T1, T2, T3, T4, T5, L1, L2, L3, L4, L5, L6.data(), dpdy, dudy, dvdy, dwdy, drhody, dYdy.data());
        update_bc_cells(lb, i, j, k, 2, cy[l_interface], L1, L2, L3, L4, L5, L6.data(), rho_out, rhoa_out, vel_out, T_out);
    }
    else if (cz[l_interface] != 0){
        // z-normal boundary
        normal_derivative(lb, i, j, k, 3, cz[l_interface], dz, dpdz, dudz, dvdz, dwdz, drhodz, dYdz.data(), nSpecies);

        if (lb.get_NX() > 3)
            tangential_derivative(lb, i, j, k, 1, dx, dpdx, dudx, dvdx, dwdx, drhodx, dYdx.data(), nSpecies);
        if (lb.get_NY() > 3)
            tangential_derivative(lb, i, j, k, 2, dy, dpdy, dudy, dvdy, dwdy, drhody, dYdy.data(), nSpecies);
        if (lb.get_NX() > 3 || lb.get_NY() > 3)
            compute_tranverse_terms(lb, i, j, k, 3, T1, T2, T3, T4, T5, dpdx, dudx, dvdx, dwdx, drhodx, dpdy, dudy, dvdy, dwdy, drhody, dpdz, dudz, dvdz, dwdz, drhodz);

        compute_waves(lb, i, j, k, 3, cz[l_interface], T1, T2, T3, T4, T5, L1, L2, L3, L4, L5, L6.data(), dpdz, dudz, dvdz, dwdz, drhodz, dYdz.data());
        update_bc_cells(lb, i, j, k, 3, cz[l_interface], L1, L2, L3, L4, L5, L6.data(), rho_out, rhoa_out, vel_out, T_out);
    }
}


void normal_derivative(LBM& lb, int i, int j, int k, int idir, int isign, double delta, double &dp, double &du, double &dv, double &dw, double &drho, double *dY, size_t nSpecies)
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

void tangential_derivative(LBM& lb, int i, int j, int k, int idir, double delta, double &dp, double &du, double &dv, double &dw, double &drho, double *dY, size_t nSpecies)
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

void compute_tranverse_terms(LBM& lb, int i, int j, int k, int idir, double& T1, double& T2, double& T3, double& T4, double& T5,
double dpdx, double dudx, double dvdx, double dwdx, double drhodx,
    double dpdy, double dudy, double dvdy, double dwdy, double drhody,
    double dpdz, double dudz, double dvdz, double dwdz, double drhodz)
{
    size_t nSpecies = lb.get_nSpecies();
    const std::vector<size_t>& speciesIdx = lb.get_speciesIdx();

    int rank = omp_get_thread_num();
    auto gas = sols[rank]->thermo();   
    std::vector <double> Y (lb.get_nSpeciesCantera());
    for(size_t a = 0; a < nSpecies; ++a) Y[speciesIdx[a]] = lb.species[a][i][j][k].rho / lb.mixture[i][j][k].rho;
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

void compute_waves(LBM& lb,
    int i, int j, int k, int idir, int isign,
    double T1, double T2, double T3, double T4, double T5,
    double& L1, double& L2, double& L3, double& L4, double& L5, double *L6,
    double dp, double du, double dv, double dw, double drho, double dYn[])
{
    // dp/du/dv/dw/drho/dYn are the one-sided derivatives along the boundary
    // normal (idir), as filled by normal_derivative().

    // Validate idir and isign
    if (idir < 1 || idir > 3) {
        throw std::invalid_argument("Problem with idir in compute_waves");
    }
    if (isign != 1 && isign != -1) {
        throw std::invalid_argument("Problem with isign in compute_waves");
    }

    size_t nSpecies = lb.get_nSpecies();
    const std::vector<size_t>& speciesIdx = lb.get_speciesIdx();
    const std::vector<double>& molarMass = lb.get_molarMass();

    int rank = omp_get_thread_num();
    auto gas = sols[rank]->thermo();
    std::vector <double> Y (lb.get_nSpeciesCantera());
    for(size_t a = 0; a < nSpecies; ++a) Y[speciesIdx[a]] = lb.species[a][i][j][k].rho / lb.mixture[i][j][k].rho;
    gas->setMassFractions(&Y[0]);
    gas->setState_DP(units.si_rho(lb.mixture[i][j][k].rho), units.si_p(lb.mixture[i][j][k].p));

    const double rho = lb.mixture[i][j][k].rho;
    const double soundSpeed = units.u(gas->soundSpeed());
    const double gamma = gas->cp_mass() / gas->cv_mass();
    const double gasconstant = units.cp(Cantera::GasConstant / gas->meanMolecularWeight());
    const double mach = v_mag(lb.mixture[i][j][k].u, lb.mixture[i][j][k].v, lb.mixture[i][j][k].w) / soundSpeed;

    // boundary (ghost) cell holding the target state, and geometry info
    int ib = i, jb = j, kb = k;
    double un = 0.0, un_target = 0.0, dun = 0.0;
    int len_char = 1;
    if (idir == 1)      { ib = i - isign; un = lb.mixture[i][j][k].u; dun = du; len_char = lb.get_Nx(); }
    else if (idir == 2) { jb = j - isign; un = lb.mixture[i][j][k].v; dun = dv; len_char = lb.get_Ny(); }
    else                { kb = k - isign; un = lb.mixture[i][j][k].w; dun = dw; len_char = lb.get_Nz(); }

    const int bc_type = lb.mixture[ib][jb][kb].type;
    const double TARGET_PRESSURE    = lb.mixture[ib][jb][kb].p;
    const double TARGET_TEMPERATURE = lb.mixture[ib][jb][kb].temp;
    if (idir == 1)      un_target = lb.mixture[ib][jb][kb].u;
    else if (idir == 2) un_target = lb.mixture[ib][jb][kb].v;
    else                un_target = lb.mixture[ib][jb][kb].w;

    // transverse damping factor (Yoo & Im 2007: beta = Mach number)
    const double beta = mach;

    // ------- interior (numerical) LODI waves along the boundary normal -------
    // acoustic waves
    L1 = (un - soundSpeed) * (dp - rho * soundSpeed * dun);
    L5 = (un + soundSpeed) * (dp + rho * soundSpeed * dun);
    // entropy wave (slot depends on direction, consistent with update_bc_cells)
    const double Ls_interior = un * (soundSpeed * soundSpeed * drho - dp);
    // tangential velocity waves
    if (idir == 1) {
        L2 = Ls_interior;
        L3 = un * dv;
        L4 = un * dw;
    } else if (idir == 2) {
        L2 = un * du;
        L3 = Ls_interior;
        L4 = un * dw;
    } else {
        L2 = un * du;
        L3 = un * dv;
        L4 = Ls_interior;
    }

    // species waves advect with the boundary-NORMAL velocity (Baum et al. 1995)
    for (size_t a = 0; a < nSpecies; ++a)
        L6[a] = un * dYn[a];

    // ---------------- boundary-specific incoming-wave modelling ----------------
    if (bc_type == TYPE_O_C /* partially non-reflecting subsonic outflow */) {
        if (mach < 1.0) {
            // Rudy & Strikwerda pressure relaxation, sigma = NSCBC_SIGMA_OUT
            // (0.25 recommended by Poinsot & Lele); the -(1-beta)T term is the
            // transverse correction of Yoo & Im (T is zero in 1D).
            double Kout = NSCBC_SIGMA_OUT * (1.0 - mach * mach) * (soundSpeed / len_char);
            if (isign == 1) {
                L5 = Kout * (lb.mixture[i][j][k].p - TARGET_PRESSURE) - (1.0 - beta) * T5;
            } else {
                L1 = Kout * (lb.mixture[i][j][k].p - TARGET_PRESSURE) - (1.0 - beta) * T1;
            }
        }
        // supersonic outflow: all characteristics leave the domain,
        // keep every wave from the interior
    } else if (bc_type == TYPE_I_C /* relaxed (soft) subsonic inflow */) {
        // Relax normal velocity, temperature and composition toward the target
        // state stored in the ghost cell. Signs follow from the LODI update
        // relations quoted at the top of this file:
        //   dun/dt = -(L5-L1)/(2 rho c)  and  drho/dt|entropy = -Ls/c^2.
        double Ku = NSCBC_RELAX_U * soundSpeed * (1.0 - mach * mach) / len_char;
        if (isign == 1) {
            L5 = L1 + 2.0 * rho * soundSpeed * Ku * (un - un_target) - (1.0 - beta) * T5;
        } else {
            L1 = L5 - 2.0 * rho * soundSpeed * Ku * (un - un_target) - (1.0 - beta) * T1;
        }

        // entropy wave: relaxes T toward the target at (approximately) constant
        // pressure; negative sign so that dT/dt = -(...)(T - T_target)
        const double Ls_relax = -NSCBC_RELAX_T * gamma * rho * gasconstant * soundSpeed / len_char
                                * (lb.mixture[i][j][k].temp - TARGET_TEMPERATURE);

        // tangential velocities relax toward the ghost-cell values
        const double Kt = NSCBC_RELAX_U * soundSpeed / len_char;
        if (idir == 1) {
            L2 = Ls_relax;
            L3 = Kt * (lb.mixture[i][j][k].v - lb.mixture[ib][jb][kb].v);
            L4 = Kt * (lb.mixture[i][j][k].w - lb.mixture[ib][jb][kb].w);
        } else if (idir == 2) {
            L2 = Kt * (lb.mixture[i][j][k].u - lb.mixture[ib][jb][kb].u);
            L3 = Ls_relax;
            L4 = Kt * (lb.mixture[i][j][k].w - lb.mixture[ib][jb][kb].w);
        } else {
            L2 = Kt * (lb.mixture[i][j][k].u - lb.mixture[ib][jb][kb].u);
            L3 = Kt * (lb.mixture[i][j][k].v - lb.mixture[ib][jb][kb].v);
            L4 = Ls_relax;
        }

        // species: relax mass fractions toward the target composition given by
        // the ghost-cell mole fractions (converted with the molar masses)
        double sumXW = 0.0;
        for (size_t a = 0; a < nSpecies; ++a)
            sumXW += lb.species[a][ib][jb][kb].X * molarMass[a];
        if (sumXW > 0.0) {
            const double Ky = NSCBC_RELAX_Y * soundSpeed / len_char;
            for (size_t a = 0; a < nSpecies; ++a) {
                double Y_target = lb.species[a][ib][jb][kb].X * molarMass[a] / sumXW;
                L6[a] = Ky * (lb.species[a][i][j][k].rho / lb.mixture[i][j][k].rho - Y_target);
            }
        }
    } else {
        throw std::runtime_error("Error: Unsupported boundary condition");
    }
}

void update_bc_cells(LBM& lb,
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
    const std::vector<size_t>& speciesIdx = lb.get_speciesIdx();

    int rank = omp_get_thread_num();
    auto gas = sols[rank]->thermo();   
    std::vector <double> Y (lb.get_nSpeciesCantera());
    for(size_t a = 0; a < nSpecies; ++a) Y[speciesIdx[a]] = lb.species[a][i][j][k].rho / lb.mixture[i][j][k].rho;
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
    for(size_t a = 0; a < nSpecies; ++a) Y[speciesIdx[a]] = rhoa_out[a] / rho_out;
    gas->setMassFractions(&Y[0]);
    gas->setState_DP(units.si_rho(rho_out), units.si_p(p_out));

    T_out = units.temp(gas->temperature());


}
#endif // MULTICOMP
