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

// Accumulate the LODI time-derivatives contributed by ONE boundary-normal
// direction (idir = 1,2,3; isign points from the boundary into the domain).
// Faces (1 normal), edges (2) and corners (3) are then handled uniformly by
// summing the locally-1D characteristic contributions of every incident face
// (Poinsot & Lele 1992; the summation corner treatment). The fully-coupled
// edge/corner formulation of Lodato, Domingo & Vervisch, JCP 227 (2008) would
// additionally modify the transverse terms at edges/corners; the summation used
// here is the simpler, robust variant.
static void nscbc_accumulate(LBM& lb, int i, int j, int k, int idir, int isign,
    size_t nSpecies, int bmask,
    double& drho_dt, double& du_dt, double& dv_dt, double& dw_dt, double& dp_dt, double* dYdt)
{
    const int dxyz[3] = { lb.get_dx(), lb.get_dy(), lb.get_dz() };
    const int Ntan[3] = { lb.get_NX(), lb.get_NY(), lb.get_NZ() };
    const int nrm = idir - 1;   // 0,1,2 normal-direction index

    // derivatives in all three directions; normal one is filled by
    // normal_derivative, the two tangential ones by tangential_derivative
    double dP[3] = {0.,0.,0.}, dU[3] = {0.,0.,0.}, dV[3] = {0.,0.,0.}, dW[3] = {0.,0.,0.}, dR[3] = {0.,0.,0.};
    std::vector<double> dYn(nSpecies, 0.0), dYt(nSpecies, 0.0);

    normal_derivative(lb, i, j, k, idir, isign, dxyz[nrm], dP[nrm], dU[nrm], dV[nrm], dW[nrm], dR[nrm], dYn.data(), nSpecies);

    bool have_transverse = false;
    #if NSCBC_TRANSVERSE
    for (int t = 0; t < 3; ++t){
        if (t == nrm) continue;
        if (Ntan[t] > 3){
            tangential_derivative(lb, i, j, k, t+1, dxyz[t], dP[t], dU[t], dV[t], dW[t], dR[t], dYt.data(), nSpecies);
            have_transverse = true;
            // Lodato edge/corner coupling: if this tangential direction is ITSELF a
            // boundary normal (edge/corner cell), its flux is already handled by its
            // own characteristic treatment, so weight its transverse contribution.
            if (bmask & (1 << t)){
                dP[t] *= NSCBC_EDGE_FACTOR; dU[t] *= NSCBC_EDGE_FACTOR;
                dV[t] *= NSCBC_EDGE_FACTOR; dW[t] *= NSCBC_EDGE_FACTOR; dR[t] *= NSCBC_EDGE_FACTOR;
            }
        }
    }
    #endif

    double T1 = 0.0, T2 = 0.0, T3 = 0.0, T4 = 0.0, T5 = 0.0;
    if (have_transverse)
        compute_tranverse_terms(lb, i, j, k, idir, T1, T2, T3, T4, T5,
            dP[0], dU[0], dV[0], dW[0], dR[0],
            dP[1], dU[1], dV[1], dW[1], dR[1],
            dP[2], dU[2], dV[2], dW[2], dR[2]);

    double L1 = 0.0, L2 = 0.0, L3 = 0.0, L4 = 0.0, L5 = 0.0;
    std::vector<double> L6(nSpecies, 0.0);
    compute_waves(lb, i, j, k, idir, isign, T1, T2, T3, T4, T5, L1, L2, L3, L4, L5, L6.data(),
                  dP[nrm], dU[nrm], dV[nrm], dW[nrm], dR[nrm], dYn.data());

    // frozen mixture sound speed of the boundary cell
    int rank = omp_get_thread_num();
    auto gas = sols[rank]->thermo();
    const std::vector<size_t>& speciesIdx = lb.get_speciesIdx();
    std::vector<double> Y (lb.get_nSpeciesCantera(), 0.0);
    for(size_t a = 0; a < nSpecies; ++a) Y[speciesIdx[a]] = lb.species[a][i][j][k].rho / lb.mixture[i][j][k].rho;
    gas->setMassFractions(&Y[0]);
    gas->setState_DP(units.si_rho(lb.mixture[i][j][k].rho), units.si_p(lb.mixture[i][j][k].p));
    const double c   = units.u(gas->soundSpeed());
    const double rho = lb.mixture[i][j][k].rho;

    // characteristic waves -> conservative time-derivatives (Poinsot & Lele)
    const double dun = (L5 - L1) / (2.0 * c * rho);
    if (idir == 1){
        drho_dt += (L2 + 0.5*(L1 + L5)) / (c*c);
        du_dt   += dun;   dv_dt += L3;   dw_dt += L4;
    } else if (idir == 2){
        drho_dt += (L3 + 0.5*(L1 + L5)) / (c*c);
        du_dt   += L2;    dv_dt += dun;  dw_dt += L4;
    } else {
        drho_dt += (L4 + 0.5*(L1 + L5)) / (c*c);
        du_dt   += L2;    dv_dt += L3;   dw_dt += dun;
    }
    dp_dt += 0.5*(L1 + L5);
    for (size_t a = 0; a < nSpecies; ++a) dYdt[a] += L6[a];
}

void impose_NSCBC(LBM& lb, int i, int j, int k, int l_interface, double &rho_out, double rhoa_out[], double vel_out[], double &T_out )
{
    (void) l_interface;   // direction(s) are detected here from the neighbourhood,
                          // so the same routine serves faces, edges and corners
    size_t nSpecies = lb.get_nSpecies();

    // accumulated LODI time-derivatives over every incident characteristic face
    double drho_dt = 0.0, du_dt = 0.0, dv_dt = 0.0, dw_dt = 0.0, dp_dt = 0.0;
    std::vector<double> dYdt(nSpecies, 0.0);

    // scan the six axis neighbours; each characteristic ghost (TYPE_O_C / TYPE_I_C)
    // is one boundary-normal contribution. isign points from the boundary into
    // the domain: ghost at cell + s*e  ->  isign = -s.
    const int off[6][4] = {  // {di,dj,dk, idir}
        {-1,0,0,1}, {+1,0,0,1},
        {0,-1,0,2}, {0,+1,0,2},
        {0,0,-1,3}, {0,0,+1,3} };
    // first pass: which axes carry a characteristic boundary (for edge/corner coupling)
    int bmask = 0;
    for (int n = 0; n < 6; ++n){
        int in = i+off[n][0], jn = j+off[n][1], kn = k+off[n][2];
        short t = lb.mixture[in][jn][kn].type;
        if (t == TYPE_O_C || t == TYPE_I_C) bmask |= (1 << (off[n][3]-1));
    }
    int ndir = 0;
    for (int n = 0; n < 6; ++n){
        int in = i+off[n][0], jn = j+off[n][1], kn = k+off[n][2];
        short t = lb.mixture[in][jn][kn].type;
        if (t == TYPE_O_C || t == TYPE_I_C){
            int s = off[n][0] + off[n][1] + off[n][2];   // -1 or +1 (the non-zero offset)
            nscbc_accumulate(lb, i, j, k, off[n][3], -s, nSpecies, bmask,
                             drho_dt, du_dt, dv_dt, dw_dt, dp_dt, dYdt.data());
            ++ndir;
        }
    }
    if (ndir == 0){   // should not happen, but stay safe: hold the current state
        rho_out = lb.mixture[i][j][k].rho;
        vel_out[0] = lb.mixture[i][j][k].u; vel_out[1] = lb.mixture[i][j][k].v; vel_out[2] = lb.mixture[i][j][k].w;
        T_out = lb.mixture[i][j][k].temp;
        for (size_t a = 0; a < nSpecies; ++a) rhoa_out[a] = lb.species[a][i][j][k].rho;
        return;
    }

    // single integration of the summed contributions
    const double dt_sim = lb.get_dtsim();
    rho_out    = lb.mixture[i][j][k].rho - dt_sim * drho_dt;
    vel_out[0] = lb.mixture[i][j][k].u   - dt_sim * du_dt;
    vel_out[1] = lb.mixture[i][j][k].v   - dt_sim * dv_dt;
    vel_out[2] = lb.mixture[i][j][k].w   - dt_sim * dw_dt;
    const double p_out = lb.mixture[i][j][k].p - dt_sim * dp_dt;
    for (size_t a = 0; a < nSpecies; ++a)
        rhoa_out[a] = (lb.species[a][i][j][k].rho/lb.mixture[i][j][k].rho - dt_sim*dYdt[a]) * rho_out;

    // temperature from the mixture equation of state (closure of Baum et al.)
    int rank = omp_get_thread_num();
    auto gas = sols[rank]->thermo();
    const std::vector<size_t>& speciesIdx = lb.get_speciesIdx();
    std::vector<double> Y (lb.get_nSpeciesCantera(), 0.0);
    for(size_t a = 0; a < nSpecies; ++a) Y[speciesIdx[a]] = rhoa_out[a] / rho_out;
    gas->setMassFractions(&Y[0]);
    gas->setState_DP(units.si_rho(rho_out), units.si_p(p_out));
    T_out = units.temp(gas->temperature());
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
    // transverse-term relaxation coefficient (Yoo & Im): the incoming wave keeps
    // (1-beta) of the transverse terms. beta = local Mach is the Yoo-Im default;
    // a fixed value can be forced through NSCBC_BETA (>= 0).
    const double beta = (NSCBC_BETA < 0.0) ? mach : NSCBC_BETA;

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

#endif // MULTICOMP
