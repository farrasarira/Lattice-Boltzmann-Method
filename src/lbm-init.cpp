#include "lbm.hpp"
#include "math_util.hpp"
#include "units.hpp"
#include "omp.h"


#ifndef MULTICOMP
void LBM::Init()
{
    #ifdef PARALLEL 
        #pragma omp parallel for schedule(static, 1) 
    #endif
    for(int i = 0; i < Nx ; ++i)
    {
        for(int j = 0; j < Ny; ++j)
        {
            for(int k = 0; k < Nz; ++k)
            {    
                if (mixture[i][j][k].type == TYPE_F || mixture[i][j][k].type == TYPE_I || mixture[i][j][k].type == TYPE_O)     
                {            
                    double velocity[3] = {  mixture[i][j][k].u,
                                            mixture[i][j][k].v, 
                                            mixture[i][j][k].w};
                    mixture[i][j][k].rho = mixture[i][j][k].p / (gas_const*mixture[i][j][k].temp);
                    double theta = gas_const*mixture[i][j][k].temp;   
                    double eq_p_tensor[3][3] = {{0., 0., 0.},    // pressure tensor
                                                {0., 0., 0.},
                                                {0., 0., 0.}};

                    for(int p=0; p < 3; ++p)
                    {
                        for(int q=0; q < 3; ++q)
                        {
                            eq_p_tensor[p][q] = (p==q) ? mixture[i][j][k].p+mixture[i][j][k].rho*velocity[p]*velocity[q] : mixture[i][j][k].rho*velocity[p]*velocity[q]; 
                            mixture[i][j][k].p_tensor[p][q] = eq_p_tensor[p][q]; 
                        }  
                    }

                    #ifndef ISOTHERM
                    double cv = gas_const / (gamma - 1.0);
                    double internal_energy = cv * mixture[i][j][k].temp;
                    mixture[i][j][k].rhoe = mixture[i][j][k].rho*(internal_energy + 0.5 * v_sqr(velocity[0], velocity[1], velocity[2]));
                    #endif


                    double corr[3] = {0, 0, 0}; 
                    for (int l = 0; l < npop; ++l)
                    {
                        // ------------- Mass and Momentum Distribution Function Initialization -----------------------------
                        mixture[i][j][k].f[l] = calculate_feq(l, mixture[i][j][k].rho, velocity, theta, corr);        
                        #ifndef ISOTHERM              
                        mixture[i][j][k].g[l] = calculate_geq(l, mixture[i][j][k].rho, internal_energy, theta, velocity);
                        #endif

                    }     

                    
                }
            }
        }
    }
    // fill_BC();
}

#elif defined MULTICOMP
void LBM::Init()
{
    #ifdef PARALLEL 
        #pragma omp parallel for schedule(static, 1) 
    #endif
    for(int i = 0; i < Nx ; ++i)
    {
        for(int j = 0; j < Ny; ++j)
        {
            for(int k = 0; k < Nz; ++k)
            {    
                if (mixture[i][j][k].type == TYPE_F || mixture[i][j][k].type == TYPE_I || mixture[i][j][k].type == TYPE_O || mixture[i][j][k].type == TYPE_O_C) // || mixture[i][j][k].type == TYPE_O     
                {            
                    // initiate Cantera object
                    int rank = omp_get_thread_num();
                    auto gas = sols[rank]->thermo();
                    std::vector <double> X (gas->nSpecies());
                    for(size_t a = 0; a < nSpecies; ++a) X[gas->speciesIndex(speciesName[a])] = species[a][i][j][k].X;
                    gas->setMoleFractions(&X[0]);
                    gas->setState_TP(units.si_temp(mixture[i][j][k].temp), units.si_p(mixture[i][j][k].p));

                    // calculate other macropscopic properties [pressure, internal energy, enthalpy, total energy]
                    mixture[i][j][k].rho = units.rho(gas->density());
                    double velocity[3] = {  mixture[i][j][k].u,
                                            mixture[i][j][k].v, 
                                            mixture[i][j][k].w};
                    #ifndef ISOTHERM
                    double internal_energy = units.energy_mass(gas->intEnergy_mass());
                    mixture[i][j][k].rhoe = mixture[i][j][k].rho*(internal_energy + 0.5 * v_sqr(velocity[0], velocity[1], velocity[2]));
                    #endif
                    double theta = units.energy_mass(gas->RT()/gas->meanMolecularWeight());

                    double eq_p_tensor[3][3] = {{0., 0., 0.},    // pressure tensor
                                                {0., 0., 0.},
                                                {0., 0., 0.}};

                    double corr[3] = {0, 0, 0}; 

                    // ------------- Energy Distribution Function Initialization -----------------------------
                    #ifndef ISOTHERM
                    for (int l = 0; l < npop; ++l)
                    {
                        mixture[i][j][k].g[l] = calculate_geq(l, mixture[i][j][k].rho, internal_energy, theta, velocity); 
                    }
                    #endif

                    // Species distribution function Initialization
                    for(size_t a = 0; a < nSpecies; ++a)
                    {
                        species[a][i][j][k].rho = gas->massFraction(gas->speciesIndex(speciesName[a])) * mixture[i][j][k].rho;
                        species[a][i][j][k].u = mixture[i][j][k].u;
                        species[a][i][j][k].v = mixture[i][j][k].v;
                        species[a][i][j][k].w = mixture[i][j][k].w;
                        theta = units.energy_mass(gas->RT() / gas->molecularWeight(gas->speciesIndex(speciesName[a])) );

                        for(int p=0; p < 3; ++p)
                        {
                            for(int q=0; q < 3; ++q)
                            {
                                eq_p_tensor[p][q] = (p==q) ? species[a][i][j][k].X*mixture[i][j][k].p+species[a][i][j][k].rho*velocity[p]*velocity[q] : species[a][i][j][k].rho*velocity[p]*velocity[q]; 
                                species[a][i][j][k].p_tensor[p][q] = eq_p_tensor[p][q]; 
                            }  
                        }

                        for (int l = 0; l < npop; ++l)
                            species[a][i][j][k].f[l]=calculate_feq(l, species[a][i][j][k].rho, velocity, theta, corr);

                    }
                }
            }
        }
    }

    // fill_BC();
}
#endif