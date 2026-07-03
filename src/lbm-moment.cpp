#include "lbm.hpp"
#include "math_util.hpp"
#include "units.hpp"
#include "omp.h"
#ifdef MULTICOMP
#include <Eigen/Dense>
#endif


#ifndef MULTICOMP
void LBM::calculate_moment()
{
    #ifdef PARALLEL
    #pragma omp parallel for schedule(dynamic)
    #endif
    for (int i = 0; i < Nx; ++i){
        for (int j = 0; j < Ny; ++j){
            for (int k = 0; k < Nz; ++k){
                if (mixture[i][j][k].type==TYPE_F ){ //
                    // defined variables
                    double rho = 0.0;
                    double rhou = 0.0;
                    double rhov = 0.0;
                    double rhow = 0.0;
                    double rhoe = 0.0;

                    double heat_flux_x=0.0;
                    double heat_flux_y=0.0;
                    double heat_flux_z=0.0;

                    for(int p=0; p < 3; ++p)
                        for(int q=0; q < 3; ++q)
                            mixture[i][j][k].p_tensor[p][q] = 0.0;
                    
                    for (int l = 0; l < npop; ++l){
                        rho+=mixture[i][j][k].f[l];
                        rhou+=mixture[i][j][k].f[l]*cx[l];
                        rhov+=mixture[i][j][k].f[l]*cy[l];
                        rhow+=mixture[i][j][k].f[l]*cz[l];
                        
                        double velocity_set[3] = {(double)cx[l], (double)cy[l], (double)cz[l]};

                        for(int p=0; p < 3; ++p)   // tensor is symmetric: accumulate the upper triangle only
                            for(int q=p; q < 3; ++q)
                                mixture[i][j][k].p_tensor[p][q] += mixture[i][j][k].f[l]*velocity_set[p]*velocity_set[q];

                        #ifndef ISOTHERM
                        rhoe += mixture[i][j][k].g[l];
                        heat_flux_x += mixture[i][j][k].g[l]*cx[l];
                        heat_flux_y += mixture[i][j][k].g[l]*cy[l];
                        heat_flux_z += mixture[i][j][k].g[l]*cz[l];
                        #endif
                    }

                    mixture[i][j][k].p_tensor[1][0] = mixture[i][j][k].p_tensor[0][1];
                    mixture[i][j][k].p_tensor[2][0] = mixture[i][j][k].p_tensor[0][2];
                    mixture[i][j][k].p_tensor[2][1] = mixture[i][j][k].p_tensor[1][2];

                    mixture[i][j][k].rho = rho;
                    mixture[i][j][k].u = rhou / rho;
                    mixture[i][j][k].v = rhov / rho;
                    mixture[i][j][k].w = rhow / rho;

                    #ifndef ISOTHERM
                    mixture[i][j][k].rhoe = rhoe;
                    mixture[i][j][k].energy_flux[0]=heat_flux_x;
                    mixture[i][j][k].energy_flux[1]=heat_flux_y;
                    mixture[i][j][k].energy_flux[2]=heat_flux_z;

                    double velocity[3] = {  mixture[i][j][k].u,
                                            mixture[i][j][k].v, 
                                            mixture[i][j][k].w};
                    double internalEnergy=mixture[i][j][k].rhoe/mixture[i][j][k].rho - 0.5*v_sqr(velocity[0], velocity[1], velocity[2]);
                    double cv = gas_const / (gamma - 1.0);
                    mixture[i][j][k].temp = internalEnergy / cv;
                    #endif
                    mixture[i][j][k].p = mixture[i][j][k].rho*gas_const*mixture[i][j][k].temp;

                    // const double G_lambda = Ra*(nu*nu/prtl) / (0.1*(Ny-2.0)*(Ny-2.0)*(Ny-2.0));
                    // mixture[i][j][k].v = mixture[i][j][k].v + 0.5*dt_sim*G_lambda*(mixture[i][j][k].temp-0.15);
                  
                       
                }

                #ifdef CONJUGATE
                else if (mixture[i][j][k].type==TYPE_S ){ //
                    // defined variables

                    double rhoe = 0.0;
                    double heat_flux_x=0.0;
                    double heat_flux_y=0.0;
                    double heat_flux_z=0.0;
                    
                    for (int l = 0; l < npop; ++l){
                        #ifndef ISOTHERM
                        rhoe += mixture[i][j][k].g[l];
                        heat_flux_x += mixture[i][j][k].g[l]*cx[l];
                        heat_flux_y += mixture[i][j][k].g[l]*cy[l];
                        heat_flux_z += mixture[i][j][k].g[l]*cz[l];
                        #endif
                    }

                    #ifndef ISOTHERM
                    mixture[i][j][k].rhoe = rhoe;
                    mixture[i][j][k].energy_flux[0]=heat_flux_x;
                    mixture[i][j][k].energy_flux[1]=heat_flux_y;
                    mixture[i][j][k].energy_flux[2]=heat_flux_z;

                    double internalEnergy=mixture[i][j][k].rhoe/mixture[i][j][k].rho;
                    double cv = gas_const / (gamma - 1.0);
                    mixture[i][j][k].temp = internalEnergy / cv;
                    #endif
                    mixture[i][j][k].p = mixture[i][j][k].rho*gas_const*mixture[i][j][k].temp;

                    // const double G_lambda = Ra*(nu*nu/prtl) / (0.1*(Ny-2.0)*(Ny-2.0)*(Ny-2.0));
                    // mixture[i][j][k].v = mixture[i][j][k].v + 0.5*dt_sim*G_lambda*(mixture[i][j][k].temp-0.15); 
                }
                #endif

            }
        }
    }

    // fill_BC();
}

void LBM::calculate_moment_point(int i, int j, int k)
{
    // defined variables
    double rho = 0.0;
    double rhou = 0.0;
    double rhov = 0.0;
    double rhow = 0.0;
    double rhoe = 0.0;

    double heat_flux_x=0.0;
    double heat_flux_y=0.0;
    double heat_flux_z=0.0;

    for(int p=0; p < 3; ++p)
        for(int q=0; q < 3; ++q)
            mixture[i][j][k].p_tensor[p][q] = 0.0;
    
    for (int l = 0; l < npop; ++l){
        rho+=mixture[i][j][k].f[l];
        rhou+=mixture[i][j][k].f[l]*cx[l];
        rhov+=mixture[i][j][k].f[l]*cy[l];
        rhow+=mixture[i][j][k].f[l]*cz[l];
        
        double velocity_set[3] = {(double)cx[l], (double)cy[l], (double)cz[l]};

        for(int p=0; p < 3; ++p)   // tensor is symmetric: accumulate the upper triangle only
            for(int q=p; q < 3; ++q)
                mixture[i][j][k].p_tensor[p][q] += mixture[i][j][k].f[l]*velocity_set[p]*velocity_set[q];

        #ifndef ISOTHERM
        rhoe += mixture[i][j][k].g[l];
        heat_flux_x += mixture[i][j][k].g[l]*cx[l];
        heat_flux_y += mixture[i][j][k].g[l]*cy[l];
        heat_flux_z += mixture[i][j][k].g[l]*cz[l];
        #endif
    }

    mixture[i][j][k].p_tensor[1][0] = mixture[i][j][k].p_tensor[0][1];
    mixture[i][j][k].p_tensor[2][0] = mixture[i][j][k].p_tensor[0][2];
    mixture[i][j][k].p_tensor[2][1] = mixture[i][j][k].p_tensor[1][2];

    mixture[i][j][k].rho = rho;
    mixture[i][j][k].u = rhou / rho;
    mixture[i][j][k].v = rhov / rho;
    mixture[i][j][k].w = rhow / rho;

    #ifndef ISOTHERM
    mixture[i][j][k].rhoe = rhoe;
    mixture[i][j][k].energy_flux[0]=heat_flux_x;
    mixture[i][j][k].energy_flux[1]=heat_flux_y;
    mixture[i][j][k].energy_flux[2]=heat_flux_z;

    double velocity[3] = {  mixture[i][j][k].u,
                            mixture[i][j][k].v, 
                            mixture[i][j][k].w};
    double internalEnergy=mixture[i][j][k].rhoe/mixture[i][j][k].rho - 0.5*v_sqr(velocity[0], velocity[1], velocity[2]);
    double cv = gas_const / (gamma - 1.0);
    mixture[i][j][k].temp = internalEnergy / cv;
    #endif
    mixture[i][j][k].p = mixture[i][j][k].rho*gas_const*mixture[i][j][k].temp;

    // const double G_lambda = Ra*(nu*nu/prtl) / (0.1*(Ny-2.0)*(Ny-2.0)*(Ny-2.0));
    // mixture[i][j][k].v = mixture[i][j][k].v + 0.5*dt_sim*G_lambda*(mixture[i][j][k].temp-0.15);
}

#elif defined MULTICOMP
void LBM::calculate_moment()
{
    #ifdef PARALLEL
    #pragma omp parallel
    #endif
    {
    // thread-local Cantera objects and work buffers, reused for every cell
    int rank = omp_get_thread_num();
    auto gas = sols[rank]->thermo();
    auto trans = sols[rank]->transport();
    const int ld = nSpeciesCantera;
    std::vector<double> Y (nSpeciesCantera);
    std::vector<double> X (nSpeciesCantera);
    std::vector<double> d (ld * ld);
    Eigen::MatrixXd matrix_A(nSpecies, nSpecies);
    Eigen::VectorXd vector_rx(nSpecies);
    Eigen::VectorXd vector_ry(nSpecies);
    Eigen::VectorXd vector_rz(nSpecies);

    #ifdef PARALLEL
    #pragma omp for collapse(2) schedule(dynamic)
    #endif
    for (int i = 0; i < Nx; ++i){
        for (int j = 0; j < Ny; ++j){
            for (int k = 0; k < Nz; ++k){
                if (mixture[i][j][k].type==TYPE_F){
                    std::fill(Y.begin(), Y.end(), 0.0);
                    for(size_t a = 0; a < nSpecies; ++a) Y[speciesIdx[a]] = species[a][i][j][k].rho / mixture[i][j][k].rho;
                    gas->setMassFractions(&Y[0]);
                    gas->setState_TD(units.si_temp(mixture[i][j][k].temp), units.si_rho(mixture[i][j][k].rho));

                    // defined variables
                    double rho = 0.0;
                    double rhou = 0.0;
                    double rhov = 0.0;
                    double rhow = 0.0;
                    for(int p=0; p < 3; ++p)
                        for(int q=0; q < 3; ++q)
                            mixture[i][j][k].p_tensor[p][q] = 0.0;
                    
                    #ifndef ISOTHERM
                    double rhoe = 0.0;
                    double heat_flux_x=0.0;
                    double heat_flux_y=0.0;
                    double heat_flux_z=0.0;

                    for (int l = 0; l < npop; ++l){
                        rhoe += mixture[i][j][k].g[l];
                        heat_flux_x += mixture[i][j][k].g[l]*cx[l];
                        heat_flux_y += mixture[i][j][k].g[l]*cy[l];
                        heat_flux_z += mixture[i][j][k].g[l]*cz[l];
                    }

                    mixture[i][j][k].rhoe = rhoe;
                    mixture[i][j][k].energy_flux[0]=heat_flux_x;
                    mixture[i][j][k].energy_flux[1]=heat_flux_y;
                    mixture[i][j][k].energy_flux[2]=heat_flux_z;
                    #endif

                    double rho_a[nSpecies] = {0.0};
                    double rhou_a[nSpecies] = {0.0};
                    double rhov_a[nSpecies] = {0.0};
                    double rhow_a[nSpecies] = {0.0};

                    for(size_t a = 0; a < nSpecies; ++a){
                        for (int l = 0; l < npop; ++l){
                            rho_a[a]  += species[a][i][j][k].f[l];
                            rhou_a[a] += species[a][i][j][k].f[l]*cx[l];
                            rhov_a[a] += species[a][i][j][k].f[l]*cy[l];
                            rhow_a[a] += species[a][i][j][k].f[l]*cz[l];

                            double velocity_set[3] = {(double)cx[l], (double)cy[l], (double)cz[l]};

                            for(int p=0; p < 3; ++p)   // tensor is symmetric: accumulate the upper triangle only
                                for(int q=p; q < 3; ++q)
                                    mixture[i][j][k].p_tensor[p][q] += species[a][i][j][k].f[l]*velocity_set[p]*velocity_set[q];
                        }
                        if (rho_a[a] > SPECIES_MIN)
                        {
                            species[a][i][j][k].rho = std::max(rho_a[a], SPECIES_MIN);
                            species[a][i][j][k].u = rhou_a[a] / rho_a[a];
                            species[a][i][j][k].v = rhov_a[a] / rho_a[a];
                            species[a][i][j][k].w = rhow_a[a] / rho_a[a];
                        }
                        else 
                        {
                            species[a][i][j][k].rho = SPECIES_MIN;
                            species[a][i][j][k].u = 0.0;
                            species[a][i][j][k].v = 0.0;
                            species[a][i][j][k].w = 0.0;
                        }
                        // std::cout << a << " | " << units.si_u(species[a][i][j][k].u) << std::endl;

                        rho += rho_a[a];
                        rhou += rhou_a[a];
                        rhov += rhov_a[a];
                        rhow += rhow_a[a];
                    }

                    mixture[i][j][k].p_tensor[1][0] = mixture[i][j][k].p_tensor[0][1];
                    mixture[i][j][k].p_tensor[2][0] = mixture[i][j][k].p_tensor[0][2];
                    mixture[i][j][k].p_tensor[2][1] = mixture[i][j][k].p_tensor[1][2];

                    mixture[i][j][k].rho = rho;
                    mixture[i][j][k].u = rhou / rho;
                    mixture[i][j][k].v = rhov / rho;
                    mixture[i][j][k].w = rhow / rho;
                                        
                    // -------------------------------------------------------------------------------------------------------------------------------
                    trans->getBinaryDiffCoeffs(ld, &d[0]);

                    double D_ab[nSpecies][nSpecies];
                    for(size_t a = 0; a < nSpecies; ++a)
                        for(size_t b = 0; b <= a; ++b){
                            D_ab[a][b] =  units.nu( d[ld*speciesIdx[b] + speciesIdx[a]] );
                            D_ab[b][a] = D_ab[a][b];
                        }

                    // dense symmetric system (the matrix is fully populated anyway,
                    // so a dense Cholesky factorization is much faster than a per-cell sparse solver)
                    for(size_t a = 0; a < nSpecies; ++a)
                        for(size_t b = 0; b < nSpecies; ++b)
                            if(b != a) matrix_A(a, b) = -1.0 * dt_sim/2.0 * mixture[i][j][k].p * species[a][i][j][k].X*species[b][i][j][k].X/D_ab[a][b];
                            else {
                                if(species[a][i][j][k].rho == 0.0) matrix_A(a, b) = 1.0;
                                else matrix_A(a, b) = species[a][i][j][k].rho;
                            }

                    for(size_t a = 0; a < nSpecies; ++a)
                        for(size_t b = 0; b < nSpecies; ++b)
                            if (b != a) matrix_A(a, a) -= matrix_A(a,b);

                    for(size_t a = 0; a < nSpecies; ++a){
                        vector_rx(a) = rhou_a[a];
                        vector_ry(a) = rhov_a[a];
                        vector_rz(a) = rhow_a[a];
                    }

                    Eigen::LLT<Eigen::MatrixXd> solver(matrix_A);

                    Eigen::VectorXd vector_u = solver.solve(vector_rx);
                    Eigen::VectorXd vector_v = solver.solve(vector_ry);
                    Eigen::VectorXd vector_w = solver.solve(vector_rz);

                    for(size_t a = 0; a < nSpecies; ++a){
                        species[a][i][j][k].u = vector_u(a);
                        species[a][i][j][k].v = vector_v(a);
                        species[a][i][j][k].w = vector_w(a);
                    }

                    // ------------------------------------------------------------------------------------------


                    std::fill(Y.begin(), Y.end(), 0.0);
                    for(size_t a = 0; a < nSpecies; ++a) Y[speciesIdx[a]] = species[a][i][j][k].rho / mixture[i][j][k].rho;
                    gas->setMassFractions(&Y[0]);

                    #ifndef ISOTHERM
                        double kineticEnergy =  0.0;
                        for (size_t a = 0; a < nSpecies; ++a)
                            kineticEnergy += 0.5 * species[a][i][j][k].rho * v_sqr(species[a][i][j][k].u, species[a][i][j][k].v, species[a][i][j][k].w);
                                                
                       double internal_energy = mixture[i][j][k].rhoe/mixture[i][j][k].rho; 
                        for (size_t a = 0; a < nSpecies; ++a){
                            internal_energy = internal_energy - 0.5 * species[a][i][j][k].rho/mixture[i][j][k].rho * v_sqr(species[a][i][j][k].u, species[a][i][j][k].v, species[a][i][j][k].w);
                        }
                        gas->setState_UV(units.si_energy_mass(internal_energy), 1.0/units.si_rho(mixture[i][j][k].rho),1.0e-15);
                    #elif defined ISOTHERM
                            gas->setState_TD(units.si_temp(mixture[i][j][k].temp), units.si_rho(mixture[i][j][k].rho));
                    #endif
                    
                    gas->getMoleFractions(&X[0]);
                    for(size_t a = 0; a < nSpecies; ++a) species[a][i][j][k].X = X[speciesIdx[a]];

                    mixture[i][j][k].temp = units.temp(gas->temperature());
                    mixture[i][j][k].p = units.p(gas->pressure());

                }
            }
        }
    }
    } // end parallel region

    // fill_BC();
}

void LBM::calculate_moment_point(int i, int j, int k)
{                   
    // defined variables
    double rho = 0.0;
    double rhou = 0.0;
    double rhov = 0.0;
    double rhow = 0.0;
    for(int p=0; p < 3; ++p)
        for(int q=0; q < 3; ++q)
            mixture[i][j][k].p_tensor[p][q] = 0.0;
    
    #ifndef ISOTHERM
    double rhoe = 0.0;
    double heat_flux_x=0.0;
    double heat_flux_y=0.0;
    double heat_flux_z=0.0;

    for (int l = 0; l < npop; ++l){
        rhoe += mixture[i][j][k].g[l];
        heat_flux_x += mixture[i][j][k].g[l]*cx[l];
        heat_flux_y += mixture[i][j][k].g[l]*cy[l];
        heat_flux_z += mixture[i][j][k].g[l]*cz[l];
    }

    mixture[i][j][k].rhoe = rhoe;
    mixture[i][j][k].energy_flux[0]=heat_flux_x;
    mixture[i][j][k].energy_flux[1]=heat_flux_y;
    mixture[i][j][k].energy_flux[2]=heat_flux_z;
    #endif

    double rho_a[nSpecies] = {0.0};
    double rhou_a[nSpecies] = {0.0};
    double rhov_a[nSpecies] = {0.0};
    double rhow_a[nSpecies] = {0.0};

    for(size_t a = 0; a < nSpecies; ++a){
        for (int l = 0; l < npop; ++l){
            rho_a[a]  += species[a][i][j][k].f[l];
            rhou_a[a] += species[a][i][j][k].f[l]*cx[l];
            rhov_a[a] += species[a][i][j][k].f[l]*cy[l];
            rhow_a[a] += species[a][i][j][k].f[l]*cz[l];

            double velocity_set[3] = {(double)cx[l], (double)cy[l], (double)cz[l]};

            for(int p=0; p < 3; ++p)   // tensor is symmetric: accumulate the upper triangle only
                for(int q=p; q < 3; ++q)
                    mixture[i][j][k].p_tensor[p][q] += species[a][i][j][k].f[l]*velocity_set[p]*velocity_set[q];
        }
        if (rho_a[a] != 0)
        {
            species[a][i][j][k].rho = rho_a[a];
            species[a][i][j][k].u = rhou_a[a] / rho_a[a];
            species[a][i][j][k].v = rhov_a[a] / rho_a[a];
            species[a][i][j][k].w = rhow_a[a] / rho_a[a];
        }
        else 
        {
            species[a][i][j][k].rho = 0.0;
            species[a][i][j][k].u = 0.0;
            species[a][i][j][k].v = 0.0;
            species[a][i][j][k].w = 0.0;
        }
        // std::cout << a << " | " << units.si_u(species[a][i][j][k].u) << std::endl;

        rho += rho_a[a];
        rhou += rhou_a[a];
        rhov += rhov_a[a];
        rhow += rhow_a[a];
    }

    mixture[i][j][k].p_tensor[1][0] = mixture[i][j][k].p_tensor[0][1];
    mixture[i][j][k].p_tensor[2][0] = mixture[i][j][k].p_tensor[0][2];
    mixture[i][j][k].p_tensor[2][1] = mixture[i][j][k].p_tensor[1][2];

    mixture[i][j][k].rho = rho;
    mixture[i][j][k].u = rhou / rho;
    mixture[i][j][k].v = rhov / rho;
    mixture[i][j][k].w = rhow / rho;

    // create Cantera Object
    int rank = omp_get_thread_num();
    auto gas = sols[rank]->thermo();
    std::vector<double> Y (nSpeciesCantera);
    std::vector<double> X (nSpeciesCantera);
    for(size_t a = 0; a < nSpecies; ++a) Y[speciesIdx[a]] = species[a][i][j][k].rho / mixture[i][j][k].rho;
    gas->setMassFractions(&Y[0]);
    gas->setState_TD(units.si_temp(mixture[i][j][k].temp), units.si_rho(mixture[i][j][k].rho));

    // -------------------------------------------------------------------------------------------------------------------------------
    const int ld = nSpeciesCantera;
    std::vector<double> d(ld * ld);
    auto trans=sols[rank]->transport();
    trans->getBinaryDiffCoeffs(ld, &d[0]);

    double D_ab[nSpecies][nSpecies];
    for(size_t a = 0; a < nSpecies; ++a)
        for(size_t b = 0; b < nSpecies; ++b)
            D_ab[a][b] =  units.nu( d[ld*speciesIdx[b] + speciesIdx[a]] );

    Eigen::MatrixXd matrix_A(nSpecies, nSpecies);
    Eigen::VectorXd vector_rx(nSpecies);
    Eigen::VectorXd vector_ry(nSpecies);
    Eigen::VectorXd vector_rz(nSpecies);

    for(size_t a = 0; a < nSpecies; ++a)
        for(size_t b = 0; b < nSpecies; ++b)
            if(b != a) matrix_A(a, b) = -1.0 * dt_sim/2.0 * mixture[i][j][k].p * species[a][i][j][k].X*species[b][i][j][k].X/D_ab[a][b];
            else {
                if(species[a][i][j][k].rho == 0.0) matrix_A(a, b) = 1.0;
                else matrix_A(a, b) = species[a][i][j][k].rho;
            }

    for(size_t a = 0; a < nSpecies; ++a)
        for(size_t b = 0; b < nSpecies; ++b)
            if (b != a) matrix_A(a, a) -= matrix_A(a,b);

    for(size_t a = 0; a < nSpecies; ++a){
        vector_rx(a) = rhou_a[a];
        vector_ry(a) = rhov_a[a];
        vector_rz(a) = rhow_a[a];
    }

    Eigen::PartialPivLU<Eigen::MatrixXd> solver(matrix_A);

    Eigen::VectorXd vector_u = solver.solve(vector_rx);
    Eigen::VectorXd vector_v = solver.solve(vector_ry);
    Eigen::VectorXd vector_w = solver.solve(vector_rz);

    for(size_t a = 0; a < nSpecies; ++a){
        species[a][i][j][k].u = vector_u(a);
        species[a][i][j][k].v = vector_v(a);
        species[a][i][j][k].w = vector_w(a);
    }

    // ------------------------------------------------------------------------------------------

    std::fill(Y.begin(), Y.end(), 0.0);
    for(size_t a = 0; a < nSpecies; ++a) Y[speciesIdx[a]] = species[a][i][j][k].rho / mixture[i][j][k].rho;
    gas->setMassFractions(&Y[0]);

    #ifndef ISOTHERM
        double kineticEnergy =  0.0;
        for (size_t a = 0; a < nSpecies; ++a)
            kineticEnergy += 0.5 * species[a][i][j][k].rho * v_sqr(species[a][i][j][k].u, species[a][i][j][k].v, species[a][i][j][k].w);
                                
        double internalEnergy = mixture[i][j][k].rhoe/mixture[i][j][k].rho - kineticEnergy/mixture[i][j][k].rho;
        gas->setState_UV(units.si_energy_mass(internalEnergy), 1.0/units.si_rho(mixture[i][j][k].rho),1.0e-15);
    #elif defined ISOTHERM
            gas->setState_TD(units.si_temp(mixture[i][j][k].temp), units.si_rho(mixture[i][j][k].rho));
    #endif
    
    gas->getMoleFractions(&X[0]);
    for(size_t a = 0; a < nSpecies; ++a) species[a][i][j][k].X = X[speciesIdx[a]];

    mixture[i][j][k].temp = units.temp(gas->temperature());
    mixture[i][j][k].p = units.p(gas->pressure());
    // gas_const = units.cp(Cantera::GasConstant/gas->meanMolecularWeight());

}


double LBM::calculate_temp(double U, double rho, double rho_a[])
{
    int rank = omp_get_thread_num();
    auto gas = sols[rank]->thermo();
    std::vector <double> Y (nSpeciesCantera);

    for(size_t a = 0; a < nSpecies; ++a){
        Y[speciesIdx[a]] = rho_a[a] / rho;
    }
    gas->setMassFractions(&Y[0]);
    gas->setState_UV(units.si_energy_mass(U), 1.0/units.si_rho(rho));

    return units.temp(gas->temperature());
}
#endif