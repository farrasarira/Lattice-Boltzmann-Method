#include "lbm.hpp"
#include "math_util.hpp"
#include "units.hpp"
#include "omp.h"


#ifndef MULTICOMP
void LBM::Collide()
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
                if (mixture[i][j][k].type == TYPE_F || mixture[i][j][k].type == TYPE_E)     
                {   
                    double cv = gas_const / (gamma - 1.0);
                    double cp = cv + gas_const;
                    double theta = gas_const*mixture[i][j][k].temp;
                    double mu = nu*mixture[i][j][k].rho;
                    double conduc_coeff = mu*cp/prtl;
                    double enthalpy = cp * mixture[i][j][k].temp; // H = Cp * T = (Cv + 1) * T

                    double velocity[3] = {  mixture[i][j][k].u,
                                            mixture[i][j][k].v, 
                                            mixture[i][j][k].w};
                
                    double omega1 = 2*mixture[i][j][k].p*cp*dt_sim / (mixture[i][j][k].p*cp*dt_sim + 2*conduc_coeff);
                    double omega = 2*mixture[i][j][k].p*dt_sim / (mixture[i][j][k].p*dt_sim + 2*mu);

                    double internal_energy = mixture[i][j][k].rhoe / mixture[i][j][k].rho - 0.5 * v_sqr(velocity[0], velocity[1], velocity[2]);
                    double total_enthalpy = enthalpy + 0.5 * v_sqr(velocity[0], velocity[1], velocity[2]);
                    double eq_heat_flux[3] = {  total_enthalpy*mixture[i][j][k].rho*mixture[i][j][k].u,
                                                total_enthalpy*mixture[i][j][k].rho*mixture[i][j][k].v,
                                                total_enthalpy*mixture[i][j][k].rho*mixture[i][j][k].w};
                    double eq_p_tensor[3][3] = {{0., 0., 0.},    // pressure tensor
                                                {0., 0., 0.},
                                                {0., 0., 0.}};

                    for(int p=0; p < 3; ++p)
                        for(int q=0; q < 3; ++q)
                            eq_p_tensor[p][q] = (p==q) ? mixture[i][j][k].p+mixture[i][j][k].rho*velocity[p]*velocity[q] : mixture[i][j][k].rho*velocity[p]*velocity[q]; 

                    double delQdevx, delQdevy, delQdevz;
                    delQdevx = FD_limiterVanleer( mixture[i+1][j][k].rho*mixture[i+1][j][k].u*(1-3*gas_const*mixture[i+1][j][k].temp)-mixture[i+1][j][k].rho*cb(mixture[i+1][j][k].u), mixture[i][j][k].rho*mixture[i][j][k].u*(1-3*gas_const*mixture[i][j][k].temp)-mixture[i][j][k].rho*cb(mixture[i][j][k].u), mixture[i-1][j][k].rho*mixture[i-1][j][k].u*(1-3*gas_const*mixture[i-1][j][k].temp)-mixture[i-1][j][k].rho*cb(mixture[i-1][j][k].u), dx) ;
                    delQdevy = FD_limiterVanleer( mixture[i][j+1][k].rho*mixture[i][j+1][k].v*(1-3*gas_const*mixture[i][j+1][k].temp)-mixture[i][j+1][k].rho*cb(mixture[i][j+1][k].v), mixture[i][j][k].rho*mixture[i][j][k].v*(1-3*gas_const*mixture[i][j][k].temp)-mixture[i][j][k].rho*cb(mixture[i][j][k].v), mixture[i][j-1][k].rho*mixture[i][j-1][k].v*(1-3*gas_const*mixture[i][j-1][k].temp)-mixture[i][j-1][k].rho*cb(mixture[i][j-1][k].v), dy) ;
                    delQdevz = FD_limiterVanleer( mixture[i][j][k+1].rho*mixture[i][j][k+1].w*(1-3*gas_const*mixture[i][j][k+1].temp)-mixture[i][j][k+1].rho*cb(mixture[i][j][k+1].w), mixture[i][j][k].rho*mixture[i][j][k].w*(1-3*gas_const*mixture[i][j][k].temp)-mixture[i][j][k].rho*cb(mixture[i][j][k].w), mixture[i][j][k-1].rho*mixture[i][j][k-1].w*(1-3*gas_const*mixture[i][j][k-1].temp)-mixture[i][j][k-1].rho*cb(mixture[i][j][k-1].w), dz) ;

                    double str_heat_flux[3] ={  eq_heat_flux[0] + velocity[0]*(mixture[i][j][k].p_tensor[0][0]-eq_p_tensor[0][0]) + velocity[1]*(mixture[i][j][k].p_tensor[1][0]-eq_p_tensor[1][0]) + velocity[2]*(mixture[i][j][k].p_tensor[2][0]-eq_p_tensor[2][0]) + 0.5*dt_sim*velocity[0]*delQdevx ,   // + 0.5*dt_sim*velocity[0]*mixture[i][j][k].dQdevx
                                                eq_heat_flux[1] + velocity[0]*(mixture[i][j][k].p_tensor[0][1]-eq_p_tensor[0][1]) + velocity[1]*(mixture[i][j][k].p_tensor[1][1]-eq_p_tensor[1][1]) + velocity[2]*(mixture[i][j][k].p_tensor[2][1]-eq_p_tensor[2][1]) + 0.5*dt_sim*velocity[1]*delQdevy ,   // + 0.5*dt_sim*velocity[1]*mixture[i][j][k].dQdevy
                                                eq_heat_flux[2] + velocity[0]*(mixture[i][j][k].p_tensor[0][2]-eq_p_tensor[0][2]) + velocity[1]*(mixture[i][j][k].p_tensor[1][2]-eq_p_tensor[1][2]) + velocity[2]*(mixture[i][j][k].p_tensor[2][2]-eq_p_tensor[2][2]) + 0.5*dt_sim*velocity[2]*delQdevz } ; // + 0.5*dt_sim*velocity[2]*mixture[i][j][k].dQdevz

                    double corr[3] = {  dt_sim*(2-omega)/(2*mixture[i][j][k].rho*omega)*delQdevx,
                                        dt_sim*(2-omega)/(2*mixture[i][j][k].rho*omega)*delQdevy,
                                        dt_sim*(2-omega)/(2*mixture[i][j][k].rho*omega)*delQdevz};

                    for(size_t l = 0; l < npop; ++l){
                        double feq = calculate_feq(l, mixture[i][j][k].rho, velocity, theta, corr);
                        mixture[i][j][k].fpc[l] = (1.0-omega)*mixture[i][j][k].f[l] + omega*feq;
                        
                        double geq = calculate_geq(l, mixture[i][j][k].rho, internal_energy, theta, velocity);
                        double gstr = calculate_gstr(l, geq, str_heat_flux, eq_heat_flux);
                        mixture[i][j][k].gpc[l] = mixture[i][j][k].g[l] + omega1*(geq-mixture[i][j][k].g[l]) + (omega-omega1)*(gstr-mixture[i][j][k].g[l]);
                    }
                        
                }
            }
        }
    }
}

#elif defined MULTICOMP
void LBM::Collide()
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
                if (mixture[i][j][k].type == TYPE_F || mixture[i][j][k].type == TYPE_E)     
                {   
                    // create Cantera object
                    int rank = omp_get_thread_num();
                    auto gas = sols[rank]->thermo();
                    std::vector <double> Y (gas->nSpecies());
                    for(size_t a = 0; a < nSpecies; ++a) Y [gas->speciesIndex(speciesName[a])] = species[a][i][j][k].rho / mixture[i][j][k].rho;
                    gas->setMassFractions(&Y[0]);
                    gas->setState_TD(units.si_temp(mixture[i][j][k].temp), units.si_rho(mixture[i][j][k].rho));

                    // get macropscopic properties
                    double cp = units.cp(gas->cp_mass());
                    double enthalpy = units.energy_mass(gas->enthalpy_mass());
                    std::vector <double> part_enthalpy(gas->nSpecies());
                    gas->getPartialMolarEnthalpies(&part_enthalpy[0]);

                    auto trans = sols[rank]->transport();
                    double mu = units.mu(trans->viscosity());
                    double conduc_coeff = units.thermalConductivity(trans->thermalConductivity());
                    
                    double theta = units.energy_mass(gas->RT()/gas->meanMolecularWeight());

                    double velocity[3] = {  mixture[i][j][k].u,
                                            mixture[i][j][k].v, 
                                            mixture[i][j][k].w};
                
                    double omega1 = 2*mixture[i][j][k].p*cp*dt_sim / (mixture[i][j][k].p*cp*dt_sim + 2*conduc_coeff);
                    double omega = 2*mixture[i][j][k].p*dt_sim / (mixture[i][j][k].p*dt_sim + 2*mu);

                    double internal_energy = mixture[i][j][k].rhoe / mixture[i][j][k].rho - 0.5 * v_sqr(velocity[0], velocity[1], velocity[2]);
                    double total_enthalpy = enthalpy + 0.5 * v_sqr(velocity[0], velocity[1], velocity[2]);
                    double eq_heat_flux[3] = {  total_enthalpy*mixture[i][j][k].rho*mixture[i][j][k].u,
                                                total_enthalpy*mixture[i][j][k].rho*mixture[i][j][k].v,
                                                total_enthalpy*mixture[i][j][k].rho*mixture[i][j][k].w};
                    double eq_p_tensor[3][3] = {{0., 0., 0.},    // pressure tensor
                                                {0., 0., 0.},
                                                {0., 0., 0.}};

                    for(int p=0; p < 3; ++p)
                    {
                        for(int q=0; q < 3; ++q)
                        {
                            eq_p_tensor[p][q] = (p==q) ? mixture[i][j][k].p+mixture[i][j][k].rho*velocity[p]*velocity[q] : mixture[i][j][k].rho*velocity[p]*velocity[q]; 
                        }  
                    }

                    double q_diff [3] = {0.0, 0.0, 0.0};
                    double q_corr [3] = {0.0, 0.0, 0.0};
                    for (size_t a = 0; a < nSpecies; ++a)
                    {
                        int speciesIdx = gas->speciesIndex(speciesName[a]);
                        double mmass = gas->molecularWeight(speciesIdx);
                        q_diff[0] += omega1/(omega-omega1)  * units.energy_mass(part_enthalpy[speciesIdx]/mmass) * mixture[i][j][k].rho * gas->massFraction(speciesIdx) * (species[a][i][j][k].Vdiff_x);
                        q_diff[1] += omega1/(omega-omega1)  * units.energy_mass(part_enthalpy[speciesIdx]/mmass) * mixture[i][j][k].rho * gas->massFraction(speciesIdx) * (species[a][i][j][k].Vdiff_y);
                        q_diff[2] += omega1/(omega-omega1)  * units.energy_mass(part_enthalpy[speciesIdx]/mmass) * mixture[i][j][k].rho * gas->massFraction(speciesIdx) * (species[a][i][j][k].Vdiff_z);
                    
                        q_corr[0] += 0.5 * (omega1-2.0)/(omega1-omega) * dt_sim * mixture[i][j][k].p * units.energy_mass(part_enthalpy[speciesIdx]/mmass) * (species[a][i+1][j][k].rho/mixture[i+1][j][k].rho - species[a][i-1][j][k].rho/mixture[i-1][j][k].rho) / (2*dx);//species[a][i][j][k].delYx;
                        q_corr[1] += 0.5 * (omega1-2.0)/(omega1-omega) * dt_sim * mixture[i][j][k].p * units.energy_mass(part_enthalpy[speciesIdx]/mmass) * (species[a][i][j+1][k].rho/mixture[i][j+1][k].rho - species[a][i][j-1][k].rho/mixture[i][j-1][k].rho) / (2*dy);//species[a][i][j][k].delYy;
                        q_corr[2] += 0.5 * (omega1-2.0)/(omega1-omega) * dt_sim * mixture[i][j][k].p * units.energy_mass(part_enthalpy[speciesIdx]/mmass) * (species[a][i][j][k+1].rho/mixture[i][j][k+1].rho - species[a][i][j][k-1].rho/mixture[i][j][k-1].rho) / (2*dz);//species[a][i][j][k].delYz;              
                    }

                    double delQdevx, delQdevy, delQdevz;
                    delQdevx = FD_limiterVanleer( mixture[i+1][j][k].rho*mixture[i+1][j][k].u*(1-3*gas_const*mixture[i+1][j][k].temp)-mixture[i+1][j][k].rho*cb(mixture[i+1][j][k].u), mixture[i][j][k].rho*mixture[i][j][k].u*(1-3*gas_const*mixture[i][j][k].temp)-mixture[i][j][k].rho*cb(mixture[i][j][k].u), mixture[i-1][j][k].rho*mixture[i-1][j][k].u*(1-3*gas_const*mixture[i-1][j][k].temp)-mixture[i-1][j][k].rho*cb(mixture[i-1][j][k].u), dx) ;
                    delQdevy = FD_limiterVanleer( mixture[i][j+1][k].rho*mixture[i][j+1][k].v*(1-3*gas_const*mixture[i][j+1][k].temp)-mixture[i][j+1][k].rho*cb(mixture[i][j+1][k].v), mixture[i][j][k].rho*mixture[i][j][k].v*(1-3*gas_const*mixture[i][j][k].temp)-mixture[i][j][k].rho*cb(mixture[i][j][k].v), mixture[i][j-1][k].rho*mixture[i][j-1][k].v*(1-3*gas_const*mixture[i][j-1][k].temp)-mixture[i][j-1][k].rho*cb(mixture[i][j-1][k].v), dy) ;
                    delQdevz = FD_limiterVanleer( mixture[i][j][k+1].rho*mixture[i][j][k+1].w*(1-3*gas_const*mixture[i][j][k+1].temp)-mixture[i][j][k+1].rho*cb(mixture[i][j][k+1].w), mixture[i][j][k].rho*mixture[i][j][k].w*(1-3*gas_const*mixture[i][j][k].temp)-mixture[i][j][k].rho*cb(mixture[i][j][k].w), mixture[i][j][k-1].rho*mixture[i][j][k-1].w*(1-3*gas_const*mixture[i][j][k-1].temp)-mixture[i][j][k-1].rho*cb(mixture[i][j][k-1].w), dz) ;

                    double str_heat_flux[3] ={  eq_heat_flux[0] + velocity[0]*(mixture[i][j][k].p_tensor[0][0]-eq_p_tensor[0][0]) + velocity[1]*(mixture[i][j][k].p_tensor[1][0]-eq_p_tensor[1][0]) + velocity[2]*(mixture[i][j][k].p_tensor[2][0]-eq_p_tensor[2][0]) + 0.5*dt_sim*velocity[0]*delQdevx ,   // + 0.5*dt_sim*velocity[0]*mixture[i][j][k].dQdevx
                                                eq_heat_flux[1] + velocity[0]*(mixture[i][j][k].p_tensor[0][1]-eq_p_tensor[0][1]) + velocity[1]*(mixture[i][j][k].p_tensor[1][1]-eq_p_tensor[1][1]) + velocity[2]*(mixture[i][j][k].p_tensor[2][1]-eq_p_tensor[2][1]) + 0.5*dt_sim*velocity[1]*delQdevy ,   // + 0.5*dt_sim*velocity[1]*mixture[i][j][k].dQdevy
                                                eq_heat_flux[2] + velocity[0]*(mixture[i][j][k].p_tensor[0][2]-eq_p_tensor[0][2]) + velocity[1]*(mixture[i][j][k].p_tensor[1][2]-eq_p_tensor[1][2]) + velocity[2]*(mixture[i][j][k].p_tensor[2][2]-eq_p_tensor[2][2]) + 0.5*dt_sim*velocity[2]*delQdevz } ; // + 0.5*dt_sim*velocity[2]*mixture[i][j][k].dQdevz

                    for (int l = 0; l < npop; ++l)
                    {
                        // ------------- Mass and Momentum collision -----------------------------
                        double geq = calculate_geq(l, mixture[i][j][k].rho, internal_energy, theta, velocity);
                        double gstr = calculate_gstr(l, geq, str_heat_flux, eq_heat_flux);
                        mixture[i][j][k].gpc[l] = mixture[i][j][k].g[l] + omega1*(geq-mixture[i][j][k].g[l]) + (omega-omega1)*(gstr-mixture[i][j][k].g[l]);

                    }     

                }
            }
        }
    }
}

void LBM::Collide_Species()
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
                if (mixture[i][j][k].type == TYPE_F || mixture[i][j][k].type == TYPE_E)     
                {   
                    int rank = omp_get_thread_num();
                    auto gas = sols[rank]->thermo();               
                    std::vector<double> Y (gas->nSpecies());
                    for(size_t a = 0; a < nSpecies; ++a) Y[gas->speciesIndex(speciesName[a])] = species[a][i][j][k].rho / mixture[i][j][k].rho;
                    gas->setMassFractions(&Y[0]);
                    gas->setState_TD(units.si_temp(mixture[i][j][k].temp), units.si_rho(mixture[i][j][k].rho));

                    // double w_dot[gas->nSpecies()];   // mole density rate [kmol/m3/s]
                    // // auto kinetics = sols[rank]->kinetics();
                    // // kinetics->getNetProductionRates(w_dot); 
                    // for (int a = 0; a < (int) gas->nSpecies(); ++a)
                    // {
                    //     if (w_dot[a] != 0)  // Check, rate != 0.
                    //     {
                    //         double idx_species;
                    //         bool new_species = true;
                    //         for (int b = 0; b < nSpecies; ++b)
                    //             if(gas->speciesName(a) == speciesName[b])
                    //             {
                    //                 new_species = false;
                    //                 idx_species = b;
                    //             }
                            
                    //         if (new_species) // adding product to the LBM variables
                    //         {
                    //             speciesName.push_back(gas->speciesName(a));
                    //             nSpecies++;
                    //             idx_species = nSpecies - 1;

                    //             // allocate memory for new species
                    //             this->species.resize(nSpecies);
                    //             this->species[idx_species] = new SPECIES **[this->Nx];
                    //             for (int i = 0; i < this->Nx; ++i)
                    //             {
                    //                 this->species[idx_species][i] = new SPECIES *[this->Ny];
                    //                 for (int j = 0; j < this->Ny; ++j)
                    //                 {
                    //                     this->species[idx_species][i][j] = new SPECIES [this->Nz];
                    //                 }
                    //             }                            
                    //         }
                    //         species[idx_species][i][j][k].rho_dot = units.rho_dot(w_dot[a] * gas->molecularWeight(a));                          
                    //     }
                    // }
                    
                    auto trans = sols[rank]->transport();
                    int ld = gas->nSpecies();
                    std::vector<double> d(ld * ld);
                    trans->getBinaryDiffCoeffs(ld, &d[0]);

                    std::vector<double> spec_visc(ld);
                    trans->getSpeciesViscosities(&spec_visc[0]);
                    double visc_a[nSpecies] = {0.0};
                    double tau_visc[nSpecies] = {0.0};

                    double mmass[nSpecies];
                    for(size_t a = 0; a < nSpecies; ++a) 
                        mmass[a] = gas->molecularWeight(gas->speciesIndex(speciesName[a]));

                    double phi[nSpecies][nSpecies] = {0.0};
                    for(size_t a = 0; a < nSpecies; ++a){                        
                        for(size_t b = a; b < nSpecies; ++b){
                            size_t idx_a = gas->speciesIndex(speciesName[a]);
                            size_t idx_b = gas->speciesIndex(speciesName[b]);

                            double factor1 = 1.0 + sqrt(spec_visc[idx_a]/spec_visc[idx_b]) * sqrt(sqrt(mmass[b]/mmass[a]));
                            phi[a][b] = factor1*factor1 / sqrt(8.0*(1.0 + mmass[a]/mmass[b]));
                            phi[b][a] = spec_visc[idx_b]/spec_visc[idx_a] * mmass[a]/mmass[b] * phi[a][b];
                        }
                    }        

                    double omega[nSpecies];  

                    for(size_t a = 0; a < nSpecies; ++a){
                        size_t idx_a = gas->speciesIndex(speciesName[a]);
                        
                        double denom = 0.0;
                        for(size_t b = 0; b < nSpecies; ++b)
                            denom = denom + species[b][i][j][k].X * phi[a][b];
                        
                        visc_a[a] = species[a][i][j][k].X * units.mu(spec_visc[idx_a]) / denom;
                        tau_visc[a] = visc_a[a]/(species[a][i][j][k].X * mixture[i][j][k].p) + dt_sim/2.0;
                        omega[a] = 2*mixture[i][j][k].p*dt_sim / (mixture[i][j][k].p*dt_sim + 2*visc_a[a]);
                    }               

                    // std::cout << visc_a[0]+visc_a[1]+visc_a[2] << " | " << units.mu(trans->viscosity()) << std::endl; 

                    double D_ab[nSpecies][nSpecies];
                    for(size_t a = 0; a < nSpecies; ++a)
                    {                       
                        for(size_t b = 0; b < nSpecies; ++b)
                        {
                            // get diffusion coefficient
                            D_ab[a][b] =  units.nu( d[ld*gas->speciesIndex(speciesName[b]) + gas->speciesIndex(speciesName[a])] );
                            // std::cout << "D_ab " << gas->speciesName(gas->speciesIndex(speciesName[a])) << "-" << gas->speciesName(gas->speciesIndex(speciesName[b])) << " : " << D_ab[a][b] << std::endl;
                        }
                    }

                    double ux[nSpecies];
                    double uy[nSpecies];
                    double uz[nSpecies];
                    double fx[nSpecies];
                    double fy[nSpecies];
                    double fz[nSpecies];
                    for(size_t a = 0; a < nSpecies; ++a){    
                        ux[a] = 0.0;
                        uy[a] = 0.0;
                        uz[a] = 0.0;   
                        for(size_t l = 0; l < npop; ++l){
                            ux[a] += species[a][i][j][k].f[l]*cx[l];
                            uy[a] += species[a][i][j][k].f[l]*cy[l];
                            uz[a] += species[a][i][j][k].f[l]*cz[l];
                        }     
                        ux[a] = ux[a] / species[a][i][j][k].rho;  
                        uy[a] = uy[a] / species[a][i][j][k].rho;
                        uz[a] = uz[a] / species[a][i][j][k].rho;

                        fx[a] = 0.0;
                        fy[a] = 0.0;
                        fz[a] = 0.0;
                        for(size_t b = 0; b < nSpecies; ++b){
                            fx[a] += -1.0 * mixture[i][j][k].p * species[a][i][j][k].X*species[b][i][j][k].X/D_ab[a][b] * (species[a][i][j][k].u - species[b][i][j][k].u);
                            fy[a] += -1.0 * mixture[i][j][k].p * species[a][i][j][k].X*species[b][i][j][k].X/D_ab[a][b] * (species[a][i][j][k].v - species[b][i][j][k].v);
                            fz[a] += -1.0 * mixture[i][j][k].p * species[a][i][j][k].X*species[b][i][j][k].X/D_ab[a][b] * (species[a][i][j][k].w - species[b][i][j][k].w);
                        }
                        fx[a] = fx[a]*dt_sim/species[a][i][j][k].rho;
                        fy[a] = fy[a]*dt_sim/species[a][i][j][k].rho;
                        fz[a] = fz[a]*dt_sim/species[a][i][j][k].rho;
                    }
                            
                    
                    for (int l = 0; l < npop; ++l)
                    {
                        double feq[nSpecies];
                        double feq_du[nSpecies];
                        // double freact[nSpecies];

                        for (size_t a = 0; a < nSpecies; ++a)
                        {
                            double velocity[3] = {  ux[a],
                                                    uy[a], 
                                                    uz[a]};
                            double velocity_du[3] = {   ux[a] + fx[a],
                                                        uy[a] + fy[a], 
                                                        uz[a] + fz[a]};
                            double theta = units.energy_mass( gas->RT()/gas->molecularWeight(gas->speciesIndex(speciesName[a])) );
                            // double gas_const = ;
                            
                            // double delQdevx, delQdevy, delQdevz;
                            // delQdevx = FD_limiterVanleer( species[a][i+1][j][k].rho*species[a][i+1][j][k].u*(1-3*gas-> )-species[a][i+1][j][k].rho*cb(species[a][i+1][j][k].u), species[a][i][j][k].rho*species[a][i][j][k].u*(1-3*theta)-species[a][i][j][k].rho*cb(species[a][i][j][k].u), species[a][i-1][j][k].rho*species[a][i-1][j][k].u*(1-3*gas_const*mixture[i-1][j][k].temp)-species[a][i-1][j][k].rho*cb(species[a][i-1][j][k].u), dx) ;
                            // delQdevy = FD_limiterVanleer( species[a][i][j+1][k].rho*species[a][i][j+1][k].v*(1-3*gas-> )-species[a][i][j+1][k].rho*cb(species[a][i][j+1][k].v), species[a][i][j][k].rho*species[a][i][j][k].v*(1-3*theta)-species[a][i][j][k].rho*cb(species[a][i][j][k].v), species[a][i][j-1][k].rho*species[a][i][j-1][k].v*(1-3*gas_const*mixture[i][j-1][k].temp)-species[a][i][j-1][k].rho*cb(species[a][i][j-1][k].v), dy) ;
                            // delQdevz = FD_limiterVanleer( species[a][i][j][k+1].rho*species[a][i][j][k+1].w*(1-3*gas-> )-species[a][i][j][k+1].rho*cb(species[a][i][j][k+1].w), species[a][i][j][k].rho*species[a][i][j][k].w*(1-3*theta)-species[a][i][j][k].rho*cb(species[a][i][j][k].w), species[a][i][j][k-1].rho*species[a][i][j][k-1].w*(1-3*gas_const*mixture[i][j][k-1].temp)-species[a][i][j][k-1].rho*cb(species[a][i][j][k-1].w), dz) ;
                            
                            // // double corr[3] = {  dt_sim*(2-omega)/(2*mixture[i][j][k].rho*omega)*delQdevx,
                            // //                     dt_sim*(2-omega)/(2*mixture[i][j][k].rho*omega)*delQdevy,
                            // //                     dt_sim*(2-omega)/(2*mixture[i][j][k].rho*omega)*delQdevz};

                            double corr[3] = {0, 0, 0};
                            
                            feq[a] = calculate_feq(l, species[a][i][j][k].rho, velocity, theta, corr);
                            feq_du[a] = calculate_feq(l, species[a][i][j][k].rho, velocity_du, theta, corr);
                            // freact[a] = feq[a] / species[a][i][j][k].rho * species[a][i][j][k].rho_dot;
                        }
                        
                        for(size_t a = 0; a < nSpecies; ++a)
                        {
                            if (tau_visc[a] != tau_visc[a]) continue;

                            species[a][i][j][k].fpc[l] = species[a][i][j][k].f[l] + dt_sim/tau_visc[a]*(feq[a]-species[a][i][j][k].f[l]) + (feq_du[a] - feq[a]);
                            // species[a][i][j][k].fpc[l] =(1.0-omega[a])*species[a][i][j][k].f[l] + omega[a]*feq[a];
                            // std::cout << a << " | " << (feq_du[a] - feq[a]) << " | " << species[a][i][j][k].fpc[l] << std::endl;
                            // std::cout << tau_visc[a] << std::endl;
                            // if (tau_visc[a] < 0.6) std::cout << "unstable !!!" << std::endl;
                        }
                    }  
                }
            }
        }
    }
}
#endif