#include "lbm.hpp"
#include "math_util.hpp"
#include "units.hpp"
#include "omp.h"
#include "FD.hpp"


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
                if (mixture[i][j][k].type == TYPE_F)     
                {   
                    double cv = gas_const / (gamma - 1.0);
                    double cp = cv + gas_const;
                    double theta = gas_const*mixture[i][j][k].temp;
                    double mu = nu*mixture[i][j][k].rho;
                    double conduc_coeff = mu*cp/prtl;

                    double omega1 = 2*mixture[i][j][k].p*cp*dt_sim / (mixture[i][j][k].p*cp*dt_sim + 2*conduc_coeff);
                    double omega = 2*mixture[i][j][k].p*dt_sim / (mixture[i][j][k].p*dt_sim + 2*mu);

                    // std::cout << 1.0/omega << " | " << 1.0/omega1 << std::endl;
                    if (1/omega <= 0.5 || 1/omega1 <=0.5)
                        std::cout << 1.0/omega << " | " << 1.0/omega1 << std::endl;

                    double velocity[3] = {  mixture[i][j][k].u,
                                            mixture[i][j][k].v, 
                                            mixture[i][j][k].w};
                    #ifndef ISOTHERM
                    double internal_energy = mixture[i][j][k].rhoe / mixture[i][j][k].rho - 0.5 * v_sqr(velocity[0], velocity[1], velocity[2]);
                    #endif
                    double eq_p_tensor[3][3] = {0.};
                    // double eq_R_tensor[3][3] = {{0., 0., 0.},    // second-order moment of g
                    //                             {0., 0., 0.},
                    //                             {0., 0., 0.}};
                    
                    // double enthalpy = cp * mixture[i][j][k].temp; // H = Cp * T = (Cv + 1) * T
                    // double total_enthalpy = enthalpy + 0.5 * v_sqr(velocity[0], velocity[1], velocity[2]);

                    for(int p=0; p < 3; ++p)
                        for(int q=0; q < 3; ++q){
                            eq_p_tensor[p][q] = (p==q) ? mixture[i][j][k].p+mixture[i][j][k].rho*velocity[p]*velocity[q] : mixture[i][j][k].rho*velocity[p]*velocity[q]; 
                            // eq_R_tensor[p][q] = total_enthalpy*eq_p_tensor[p][q] + mixture[i][j][k].p*velocity[p]*velocity[q];
                        }

                    // double eq_heat_flux[3] = {  total_enthalpy*mixture[i][j][k].rho*mixture[i][j][k].u,
                    //                             total_enthalpy*mixture[i][j][k].rho*mixture[i][j][k].v,
                    //                             total_enthalpy*mixture[i][j][k].rho*mixture[i][j][k].w};
                    
                    // double str_heat_flux[3] ={  eq_heat_flux[0] + velocity[0]*(mixture[i][j][k].p_tensor[0][0]-eq_p_tensor[0][0]) + velocity[1]*(mixture[i][j][k].p_tensor[1][0]-eq_p_tensor[1][0]) + velocity[2]*(mixture[i][j][k].p_tensor[2][0]-eq_p_tensor[2][0]) + 0.5*dt_sim*velocity[0]*mixture[i][j][k].dQdevx,   // - 0.5*dt_sim*velocity[0]*mixture[i][j][k].dQdevx
                    //                             eq_heat_flux[1] + velocity[0]*(mixture[i][j][k].p_tensor[0][1]-eq_p_tensor[0][1]) + velocity[1]*(mixture[i][j][k].p_tensor[1][1]-eq_p_tensor[1][1]) + velocity[2]*(mixture[i][j][k].p_tensor[2][1]-eq_p_tensor[2][1]) + 0.5*dt_sim*velocity[1]*mixture[i][j][k].dQdevy,   // - 0.5*dt_sim*velocity[1]*mixture[i][j][k].dQdevy
                    //                             eq_heat_flux[2] + velocity[0]*(mixture[i][j][k].p_tensor[0][2]-eq_p_tensor[0][2]) + velocity[1]*(mixture[i][j][k].p_tensor[1][2]-eq_p_tensor[1][2]) + velocity[2]*(mixture[i][j][k].p_tensor[2][2]-eq_p_tensor[2][2]) + 0.5*dt_sim*velocity[2]*mixture[i][j][k].dQdevz} ; // - 0.5*dt_sim*velocity[2]*mixture[i][j][k].dQdevz    

                    double delQdevx, delQdevy, delQdevz;
                    delQdevx = fd_central(mixture[i-1][j][k].rho*mixture[i-1][j][k].u*(1.0-3.0*gas_const*mixture[i-1][j][k].temp)-mixture[i-1][j][k].rho*cb(mixture[i-1][j][k].u), mixture[i][j][k].rho*mixture[i][j][k].u*(1.0-3.0*gas_const*mixture[i][j][k].temp)-mixture[i][j][k].rho*cb(mixture[i][j][k].u), mixture[i+1][j][k].rho*mixture[i+1][j][k].u*(1.0-3.0*gas_const*mixture[i+1][j][k].temp)-mixture[i+1][j][k].rho*cb(mixture[i+1][j][k].u), dx, mixture[i][j][k].u, mixture[i-1][j][k].type, mixture[i+1][j][k].type) ;
                    delQdevy = fd_central(mixture[i][j-1][k].rho*mixture[i][j-1][k].v*(1.0-3.0*gas_const*mixture[i][j-1][k].temp)-mixture[i][j-1][k].rho*cb(mixture[i][j-1][k].v), mixture[i][j][k].rho*mixture[i][j][k].v*(1.0-3.0*gas_const*mixture[i][j][k].temp)-mixture[i][j][k].rho*cb(mixture[i][j][k].v), mixture[i][j+1][k].rho*mixture[i][j+1][k].v*(1.0-3.0*gas_const*mixture[i][j+1][k].temp)-mixture[i][j+1][k].rho*cb(mixture[i][j+1][k].v), dy, mixture[i][j][k].v, mixture[i][j-1][k].type, mixture[i][j+1][k].type) ;
                    delQdevz = fd_central(mixture[i][j][k-1].rho*mixture[i][j][k-1].w*(1.0-3.0*gas_const*mixture[i][j][k-1].temp)-mixture[i][j][k-1].rho*cb(mixture[i][j][k-1].w), mixture[i][j][k].rho*mixture[i][j][k].w*(1.0-3.0*gas_const*mixture[i][j][k].temp)-mixture[i][j][k].rho*cb(mixture[i][j][k].w), mixture[i][j][k+1].rho*mixture[i][j][k+1].w*(1.0-3.0*gas_const*mixture[i][j][k+1].temp)-mixture[i][j][k+1].rho*cb(mixture[i][j][k+1].w), dz, mixture[i][j][k].w, mixture[i][j][k-1].type, mixture[i][j][k+1].type) ;

                    // std::cout << j << " | " << delQdevy << " | " << mixture[i][j][k].dQdevy << std::endl;
                    // std::cout << mixture[i][j][k].dQdevy << std::endl;

                    double d_str_heat_flux[3] ={  velocity[0]*(mixture[i][j][k].p_tensor[0][0]-eq_p_tensor[0][0]) + velocity[1]*(mixture[i][j][k].p_tensor[1][0]-eq_p_tensor[1][0]) + velocity[2]*(mixture[i][j][k].p_tensor[2][0]-eq_p_tensor[2][0]) + 0.5*dt_sim*velocity[0]*delQdevx,   //+ 0.5*dt_sim*velocity[0]*delQdevx + 0.5*dt_sim*velocity[0]*mixture[i][j][k].dQdevx
                                                  velocity[0]*(mixture[i][j][k].p_tensor[0][1]-eq_p_tensor[0][1]) + velocity[1]*(mixture[i][j][k].p_tensor[1][1]-eq_p_tensor[1][1]) + velocity[2]*(mixture[i][j][k].p_tensor[2][1]-eq_p_tensor[2][1]) + 0.5*dt_sim*velocity[1]*delQdevy,   //+ 0.5*dt_sim*velocity[1]*delQdevy + 0.5*dt_sim*velocity[1]*mixture[i][j][k].dQdevy
                                                  velocity[0]*(mixture[i][j][k].p_tensor[0][2]-eq_p_tensor[0][2]) + velocity[1]*(mixture[i][j][k].p_tensor[1][2]-eq_p_tensor[1][2]) + velocity[2]*(mixture[i][j][k].p_tensor[2][2]-eq_p_tensor[2][2]) + 0.5*dt_sim*velocity[2]*delQdevz} ; //+ 0.5*dt_sim*velocity[2]*delQdevz + 0.5*dt_sim*velocity[2]*mixture[i][j][k].dQdevz

                    double corr[3] = {  dt_sim*(2-omega)/(2*mixture[i][j][k].rho*omega)*delQdevx,
                                        dt_sim*(2-omega)/(2*mixture[i][j][k].rho*omega)*delQdevy,
                                        dt_sim*(2-omega)/(2*mixture[i][j][k].rho*omega)*delQdevz};
                    // double corr[3] = {0.0, 0.0, 0.0};

                    // const double G_lambda = Ra*(nu*nu/prtl) / (0.025*(Ny-2.0)*(Ny-2.0)*(Ny-2.0));
                    // double velocity_du[3] = {   mixture[i][j][k].u + 0.0,
                    //                             mixture[i][j][k].v + dt_sim*G_lambda*(mixture[i][j][k].temp-0.0375), 
                    //                             mixture[i][j][k].w + 0.0};

                    for(size_t l = 0; l < npop; ++l){
                        double feq = calculate_feq(l, mixture[i][j][k].rho, velocity, theta, corr);
                        // double feq_du = calculate_feq(l, mixture[i][j][k].rho, velocity_du, theta, corr);
                        mixture[i][j][k].fpc[l] = (1.0-omega)*mixture[i][j][k].f[l] + omega*feq ;//+ (feq_du-feq);
                        #ifndef ISOTHERM
                        double geq = calculate_geq(l, mixture[i][j][k].rho, internal_energy, theta, velocity);
                        double gstr = calculate_gstr(l, geq, d_str_heat_flux);
                        // double geq = calculate_geq(l, mixture[i][j][k].rhoe, eq_heat_flux, eq_R_tensor, theta);
                        // double gstr = calculate_geq(l, mixture[i][j][k].rhoe, str_heat_flux, eq_R_tensor, theta);
                        mixture[i][j][k].gpc[l] = mixture[i][j][k].g[l] + omega1*(geq-mixture[i][j][k].g[l]) + (omega-omega1)*(geq-gstr);
                        #endif
                    }
                        
                }
            }
        }
    }
}

#elif defined MULTICOMP
void LBM::Collide()
{   
    #ifndef ISOTHERM

    #ifdef PARALLEL 
        #pragma omp parallel for schedule(static, 1) 
    #endif
    for(int i = 0; i < Nx ; ++i)
    {
        for(int j = 0; j < Ny; ++j)
        {
            for(int k = 0; k < Nz; ++k)
            {    
                if (mixture[i][j][k].type == TYPE_F)     
                {   
                    // create Cantera object
                    int rank = omp_get_thread_num();
                    auto gas = sols[rank]->thermo();
                    std::vector <double> Y (gas->nSpecies());
                    for(size_t a = 0; a < nSpecies; ++a) Y[gas->speciesIndex(speciesName[a])] = species[a][i][j][k].rho / mixture[i][j][k].rho;
                    gas->setMassFractions(&Y[0]);
                    gas->setState_TD(units.si_temp(mixture[i][j][k].temp), units.si_rho(mixture[i][j][k].rho));

                    // get macropscopic properties
                    double cp = units.cp(gas->cp_mass());
                    double gasconstant = units.cp(Cantera::GasConstant/gas->meanMolecularWeight());
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

                    // std::cout << 1/omega << " | " << 1/omega1 << std::endl;

                    double internal_energy = mixture[i][j][k].rhoe / mixture[i][j][k].rho - 0.5 * v_sqr(velocity[0], velocity[1], velocity[2]);
                    double eq_p_tensor[3][3] = {{0., 0., 0.},    // pressure tensor
                                                {0., 0., 0.},
                                                {0., 0., 0.}};

                    for(int p=0; p < 3; ++p)
                        for(int q=0; q < 3; ++q){
                            eq_p_tensor[p][q] = (p==q) ? mixture[i][j][k].p+mixture[i][j][k].rho*velocity[p]*velocity[q] : mixture[i][j][k].rho*velocity[p]*velocity[q]; 
                            // eq_R_tensor[p][q] = total_enthalpy*eq_p_tensor[p][q] + mixture[i][j][k].p*velocity[p]*velocity[q];
                        }               
                   

                    double q_diff [3] = {0.0, 0.0, 0.0};
                    double q_corr [3] = {0.0, 0.0, 0.0};
                    for (size_t a = 0; a < nSpecies; ++a)
                    {
                        int speciesIdx = gas->speciesIndex(speciesName[a]);
                        double mmass = gas->molecularWeight(speciesIdx);
                        q_diff[0] += omega1/(omega-omega1)  * units.energy_mass(part_enthalpy[speciesIdx]/mmass) * mixture[i][j][k].rho * gas->massFraction(speciesIdx) * (species[a][i][j][k].u-mixture[i][j][k].u);
                        q_diff[1] += omega1/(omega-omega1)  * units.energy_mass(part_enthalpy[speciesIdx]/mmass) * mixture[i][j][k].rho * gas->massFraction(speciesIdx) * (species[a][i][j][k].v-mixture[i][j][k].v);
                        q_diff[2] += omega1/(omega-omega1)  * units.energy_mass(part_enthalpy[speciesIdx]/mmass) * mixture[i][j][k].rho * gas->massFraction(speciesIdx) * (species[a][i][j][k].w-mixture[i][j][k].w);

                        double dY_x = fd_fuw(species[a][i-1][j][k].rho/mixture[i-1][j][k].rho, species[a][i][j][k].rho/mixture[i][j][k].rho, species[a][i+1][j][k].rho/mixture[i+1][j][k].rho, dx, species[a][i][j][k].u-mixture[i][j][k].u, mixture[i-1][j][k].type, mixture[i+1][j][k].type);
                        double dY_y = fd_fuw(species[a][i][j-1][k].rho/mixture[i][j-1][k].rho, species[a][i][j][k].rho/mixture[i][j][k].rho, species[a][i][j+1][k].rho/mixture[i][j+1][k].rho, dy, species[a][i][j][k].v-mixture[i][j][k].v, mixture[i][j-1][k].type, mixture[i][j+1][k].type);
                        double dY_z = fd_fuw(species[a][i][j][k-1].rho/mixture[i][j][k-1].rho, species[a][i][j][k].rho/mixture[i][j][k].rho, species[a][i][j][k+1].rho/mixture[i][j][k+1].rho, dz, species[a][i][j][k].w-mixture[i][j][k].w, mixture[i][j][k-1].type, mixture[i][j][k+1].type);

                        q_corr[0] += 0.5 * (2.0-omega1)/(omega-omega1) * dt_sim * mixture[i][j][k].p * units.energy_mass(part_enthalpy[speciesIdx]/mmass) * dY_x;//species[a][i][j][k].delYx;
                        q_corr[1] += 0.5 * (2.0-omega1)/(omega-omega1) * dt_sim * mixture[i][j][k].p * units.energy_mass(part_enthalpy[speciesIdx]/mmass) * dY_y;//species[a][i][j][k].delYy;
                        q_corr[2] += 0.5 * (2.0-omega1)/(omega-omega1) * dt_sim * mixture[i][j][k].p * units.energy_mass(part_enthalpy[speciesIdx]/mmass) * dY_z;//species[a][i][j][k].delYz;  
                    }

                    double delQdevx, delQdevy, delQdevz;
                    delQdevx = fd_central(mixture[i-1][j][k].rho*mixture[i-1][j][k].u*(1-3*gasconstant*mixture[i-1][j][k].temp)-mixture[i-1][j][k].rho*cb(mixture[i-1][j][k].u), mixture[i][j][k].rho*mixture[i][j][k].u*(1-3*gasconstant*mixture[i][j][k].temp)-mixture[i][j][k].rho*cb(mixture[i][j][k].u), mixture[i+1][j][k].rho*mixture[i+1][j][k].u*(1-3*gasconstant*mixture[i+1][j][k].temp)-mixture[i+1][j][k].rho*cb(mixture[i+1][j][k].u), dx, mixture[i][j][k].u, mixture[i-1][j][k].type, mixture[i+1][j][k].type) ;
                    delQdevy = fd_central(mixture[i][j-1][k].rho*mixture[i][j-1][k].v*(1-3*gasconstant*mixture[i][j-1][k].temp)-mixture[i][j-1][k].rho*cb(mixture[i][j-1][k].v), mixture[i][j][k].rho*mixture[i][j][k].v*(1-3*gasconstant*mixture[i][j][k].temp)-mixture[i][j][k].rho*cb(mixture[i][j][k].v), mixture[i][j+1][k].rho*mixture[i][j+1][k].v*(1-3*gasconstant*mixture[i][j+1][k].temp)-mixture[i][j+1][k].rho*cb(mixture[i][j+1][k].v), dy, mixture[i][j][k].v, mixture[i][j-1][k].type, mixture[i][j+1][k].type) ;
                    delQdevz = fd_central(mixture[i][j][k-1].rho*mixture[i][j][k-1].w*(1-3*gasconstant*mixture[i][j][k-1].temp)-mixture[i][j][k-1].rho*cb(mixture[i][j][k-1].w), mixture[i][j][k].rho*mixture[i][j][k].w*(1-3*gasconstant*mixture[i][j][k].temp)-mixture[i][j][k].rho*cb(mixture[i][j][k].w), mixture[i][j][k+1].rho*mixture[i][j][k+1].w*(1-3*gasconstant*mixture[i][j][k+1].temp)-mixture[i][j][k+1].rho*cb(mixture[i][j][k+1].w), dz, mixture[i][j][k].w, mixture[i][j][k-1].type, mixture[i][j][k+1].type) ;

                    double d_str_heat_flux[3] ={    velocity[0]*(mixture[i][j][k].p_tensor[0][0]-eq_p_tensor[0][0]) + velocity[1]*(mixture[i][j][k].p_tensor[1][0]-eq_p_tensor[1][0]) + velocity[2]*(mixture[i][j][k].p_tensor[2][0]-eq_p_tensor[2][0]) - q_diff[0] - q_corr[0] + 0.5*dt_sim*velocity[0]*delQdevx,   //+ q_diff[0] + q_corr[0]  + 0.5*dt_sim*velocity[0]*delQdevx
                                                    velocity[0]*(mixture[i][j][k].p_tensor[0][1]-eq_p_tensor[0][1]) + velocity[1]*(mixture[i][j][k].p_tensor[1][1]-eq_p_tensor[1][1]) + velocity[2]*(mixture[i][j][k].p_tensor[2][1]-eq_p_tensor[2][1]) - q_diff[1] - q_corr[1] + 0.5*dt_sim*velocity[1]*delQdevy,   //+ q_diff[1] + q_corr[1]  + 0.5*dt_sim*velocity[1]*delQdevy
                                                    velocity[0]*(mixture[i][j][k].p_tensor[0][2]-eq_p_tensor[0][2]) + velocity[1]*(mixture[i][j][k].p_tensor[1][2]-eq_p_tensor[1][2]) + velocity[2]*(mixture[i][j][k].p_tensor[2][2]-eq_p_tensor[2][2]) - q_diff[2] - q_corr[2] + 0.5*dt_sim*velocity[2]*delQdevz} ; //+ q_diff[2] + q_corr[2]  + 0.5*dt_sim*velocity[2]*delQdevz
                
                    // double str_heat_flux[3] ={  eq_heat_flux[0] + velocity[0]*(mixture[i][j][k].p_tensor[0][0]-eq_p_tensor[0][0]) + velocity[1]*(mixture[i][j][k].p_tensor[1][0]-eq_p_tensor[1][0]) + velocity[2]*(mixture[i][j][k].p_tensor[2][0]-eq_p_tensor[2][0]) + q_diff[0] + q_corr[0],   // - 0.5*dt_sim*velocity[0]*mixture[i][j][k].dQdevx
                    //                             eq_heat_flux[1] + velocity[0]*(mixture[i][j][k].p_tensor[0][1]-eq_p_tensor[0][1]) + velocity[1]*(mixture[i][j][k].p_tensor[1][1]-eq_p_tensor[1][1]) + velocity[2]*(mixture[i][j][k].p_tensor[2][1]-eq_p_tensor[2][1]) + q_diff[1] + q_corr[1],   // - 0.5*dt_sim*velocity[1]*mixture[i][j][k].dQdevy
                    //                             eq_heat_flux[2] + velocity[0]*(mixture[i][j][k].p_tensor[0][2]-eq_p_tensor[0][2]) + velocity[1]*(mixture[i][j][k].p_tensor[1][2]-eq_p_tensor[1][2]) + velocity[2]*(mixture[i][j][k].p_tensor[2][2]-eq_p_tensor[2][2]) + q_diff[2] + q_corr[2]} ; // - 0.5*dt_sim*velocity[2]*mixture[i][j][k].dQdevz    
   

                    for (int l = 0; l < npop; ++l)
                    {
                        // ------------- Mass and Momentum collision -----------------------------
                       
                        double geq = calculate_geq(l, mixture[i][j][k].rho, internal_energy, theta, velocity);
                        double gstr = calculate_gstr(l, geq, d_str_heat_flux);
                        // double geq = calculate_geq(l, mixture[i][j][k].rhoe, eq_heat_flux, eq_R_tensor, theta);
                        // double gstr = calculate_geq(l, mixture[i][j][k].rhoe, str_heat_flux, eq_R_tensor, theta);
                        mixture[i][j][k].gpc[l] = mixture[i][j][k].g[l] + omega1*(geq-mixture[i][j][k].g[l]) + (omega-omega1)*(geq-gstr);
                    }     

                }
            }
        }
    }
    #endif
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
                if (mixture[i][j][k].type == TYPE_F)     
                {   
                    int rank = omp_get_thread_num();
                    auto gas = sols[rank]->thermo();   
                    std::vector <double> Y (gas->nSpecies());
                    for(size_t a = 0; a < nSpecies; ++a) Y[gas->speciesIndex(speciesName[a])] = (species[a][i][j][k].rho) / mixture[i][j][k].rho;
                    gas->setMassFractions(&Y[0]);   
                    gas->setState_TD(units.si_temp(mixture[i][j][k].temp), units.si_rho(mixture[i][j][k].rho));

                    double mmass[nSpecies];
                    size_t idx_species[nSpecies];
                    for(size_t a = 0; a < nSpecies; ++a){
                        idx_species[a] = gas->speciesIndex(speciesName[a]);
                        mmass[a] = gas->molecularWeight(idx_species[a]);
                    }
                                        
                    // Reaction ---------------------------------------------------------------------------------------------------------------------------
                    #ifdef REACTION
                        double velocity_mix[3] = {  mixture[i][j][k].u,
                                                    mixture[i][j][k].v, 
                                                    mixture[i][j][k].w};      
                        double internalEnergy = mixture[i][j][k].rhoe / mixture[i][j][k].rho - 0.5 * v_sqr(velocity_mix[0], velocity_mix[1], velocity_mix[2]);
                        std::vector <double> w_dot(gas->nSpecies());    // mole density rate [kmol/m3/s]
                        auto kinetics = sols[rank]->kinetics();

                        size_t p = 0;
                        double maxIter = 2.0;
                        double rho_a[gas->nSpecies()] = {0.0};
                        do{
                            kinetics->getNetProductionRates(&w_dot[0]); 
                            for (size_t a = 0; a < nSpecies; ++a){                                    
                                double rho_dot = units.rho_dot(w_dot[a] * gas->molecularWeight(a));
                                rho_a[a] += dt_sim/maxIter * rho_dot;
                            }
                    
                            std::fill(Y.begin(), Y.end(), 0.0);
                            for(size_t a = 0; a < nSpecies; ++a) Y[gas->speciesIndex(speciesName[a])] = (species[a][i][j][k].rho+rho_a[a]) / mixture[i][j][k].rho;
                            gas->setMassFractions(&Y[0]);
                            gas->setState_UV(units.si_energy_mass(internalEnergy), 1.0/units.si_rho(mixture[i][j][k].rho),1.0e-15);
                            p++;
                        }while(p < maxIter);
                    #endif

                    // Viscosity -----------------------------------------------------------------------------------------------------------------------------                    
                    std::fill(Y.begin(), Y.end(), 0.0);
                    for(size_t a = 0; a < nSpecies; ++a) Y[gas->speciesIndex(speciesName[a])] = (species[a][i][j][k].rho) / mixture[i][j][k].rho;
                    gas->setMassFractions(&Y[0]);
                    gas->setState_TD(units.si_temp(mixture[i][j][k].temp), units.si_rho(mixture[i][j][k].rho));
                    
                    auto trans = sols[rank]->transport();
                    int ld = gas->nSpecies();
                    std::vector<double> d(ld * ld);
                    trans->getBinaryDiffCoeffs(ld, &d[0]);

                    std::vector<double> spec_visc(ld);
                    trans->getSpeciesViscosities(&spec_visc[0]);
                    double visc_a[nSpecies] = {0.0};
                    double omega_a[nSpecies] = {0.0};

                    double phi[nSpecies][nSpecies] = {0.0};
                    for(size_t a = 0; a < nSpecies; ++a){                        
                        for(size_t b = a; b < nSpecies; ++b){
                            size_t idx_a = idx_species[a];
                            size_t idx_b = idx_species[b];

                            double factor1 = 1.0 + sqrt(spec_visc[idx_a]/spec_visc[idx_b]) * sqrt(sqrt(mmass[b]/mmass[a]));
                            phi[a][b] = factor1*factor1 / sqrt(8.0*(1.0 + mmass[a]/mmass[b]));
                            phi[b][a] = spec_visc[idx_b]/spec_visc[idx_a] * mmass[a]/mmass[b] * phi[a][b];
                        }
                    }        

                    for(size_t a = 0; a < nSpecies; ++a){
                        size_t idx_a = idx_species[a];
                        
                        double denom = 0.0;
                        for(size_t b = 0; b < nSpecies; ++b)
                            denom = denom + species[b][i][j][k].X * phi[a][b];
                        
                        visc_a[a] = species[a][i][j][k].X * units.mu(spec_visc[idx_a]) / denom;
                        omega_a[a] = 2*species[a][i][j][k].X*mixture[i][j][k].p*dt_sim / (species[a][i][j][k].X*mixture[i][j][k].p*dt_sim + 2*visc_a[a]);
                    }               

                    // std::cout << visc_a[1] << " | " << units.mu(trans->viscosity()) << std::endl; 

                    // Diffusion ----------------------------------------------------------------------------------------------------------------------------
                                        
                    double D_ab[nSpecies][nSpecies];
                    for(size_t a = 0; a < nSpecies; ++a)
                    {                       
                        for(size_t b = 0; b <= a; ++b)
                        {
                            // get diffusion coefficient
                            D_ab[a][b] =  units.nu( d[ld*gas->speciesIndex(speciesName[b]) + gas->speciesIndex(speciesName[a])] );
                            D_ab[b][a] =  D_ab[a][b];
                            // std::cout << "D_ab " << gas->speciesName(gas->speciesIndex(speciesName[a])) << "-" << gas->speciesName(gas->speciesIndex(speciesName[b])) << " : " << d[ld*gas->speciesIndex(speciesName[b]) + gas->speciesIndex(speciesName[a])] << std::endl;
                        }
                    }

                    double ux[nSpecies] = {};
                    double uy[nSpecies] = {};
                    double uz[nSpecies] = {};
                    double fx[nSpecies] = {};
                    double fy[nSpecies] = {};
                    double fz[nSpecies] = {};
                    for(size_t a = 0; a < nSpecies; ++a){     
                        for(size_t l = 0; l < npop; ++l){
                            ux[a] += species[a][i][j][k].f[l]*cx[l];
                            uy[a] += species[a][i][j][k].f[l]*cy[l];
                            uz[a] += species[a][i][j][k].f[l]*cz[l];
                        }     
                        ux[a] = ux[a] / species[a][i][j][k].rho;  
                        uy[a] = uy[a] / species[a][i][j][k].rho;
                        uz[a] = uz[a] / species[a][i][j][k].rho;

                        for(size_t b = 0; b < nSpecies; ++b){
                            fx[a] += -1.0 * mixture[i][j][k].p * species[a][i][j][k].X*species[b][i][j][k].X/D_ab[a][b] * (species[a][i][j][k].u - species[b][i][j][k].u);
                            fy[a] += -1.0 * mixture[i][j][k].p * species[a][i][j][k].X*species[b][i][j][k].X/D_ab[a][b] * (species[a][i][j][k].v - species[b][i][j][k].v);
                            fz[a] += -1.0 * mixture[i][j][k].p * species[a][i][j][k].X*species[b][i][j][k].X/D_ab[a][b] * (species[a][i][j][k].w - species[b][i][j][k].w);
                        }
                        fx[a] = fx[a]*dt_sim/species[a][i][j][k].rho;
                        fy[a] = fy[a]*dt_sim/species[a][i][j][k].rho;
                        fz[a] = fz[a]*dt_sim/species[a][i][j][k].rho;
                    }

                    // Bouyancy ------------------------------------------------------------------------------------------------------------------------------
                    // const double NU = units.nu(trans->viscosity()/gas->density());
                    // const double PR = trans->viscosity()*gas->cp_mass() / trans->thermalConductivity();       
                    // const double G_lambda = Ra*(NU*NU/PR) / (units.temp(50.0)*(Ny-2.0)*(Ny-2.0)*(Ny-2.0));

                    // Calculate Moment & Kinetic Equation ---------------------------------------------------------------------------------------------------
                    for (size_t a = 0; a < nSpecies; ++a)
                    {
                        double velocity[3] = {  ux[a],
                                                uy[a], 
                                                uz[a]   };
                        double velocity_du[3] = {   ux[a] + fx[a],
                                                    uy[a] + fy[a], //+ dt_sim*G_lambda*(mixture[i][j][k].temp-units.temp(325.0) 
                                                    uz[a] + fz[a]};

                        double gas_const_a = units.cp(Cantera::GasConstant/mmass[a]);
                        double theta = gas_const_a * mixture[i][j][k].temp;

                        double delQdevx, delQdevy, delQdevz;                            
                        delQdevx = fd_central(species[a][i-1][j][k].rho*species[a][i-1][j][k].u*(1-3*gas_const_a*mixture[i-1][j][k].temp)-species[a][i-1][j][k].rho*cb(species[a][i-1][j][k].u), species[a][i][j][k].rho*species[a][i][j][k].u*(1-3*gas_const_a*mixture[i][j][k].temp)-species[a][i][j][k].rho*cb(species[a][i][j][k].u), species[a][i+1][j][k].rho*species[a][i+1][j][k].u*(1-3*gas_const_a*mixture[i+1][j][k].temp)-species[a][i+1][j][k].rho*cb(species[a][i+1][j][k].u), dx, species[a][i][j][k].u, mixture[i-1][j][k].type, mixture[i+1][j][k].type);
                        delQdevy = fd_central(species[a][i][j-1][k].rho*species[a][i][j-1][k].v*(1-3*gas_const_a*mixture[i][j-1][k].temp)-species[a][i][j-1][k].rho*cb(species[a][i][j-1][k].v), species[a][i][j][k].rho*species[a][i][j][k].v*(1-3*gas_const_a*mixture[i][j][k].temp)-species[a][i][j][k].rho*cb(species[a][i][j][k].v), species[a][i][j+1][k].rho*species[a][i][j+1][k].v*(1-3*gas_const_a*mixture[i][j+1][k].temp)-species[a][i][j+1][k].rho*cb(species[a][i][j+1][k].v), dy, species[a][i][j][k].v, mixture[i][j-1][k].type, mixture[i][j+1][k].type);
                        delQdevz = fd_central(species[a][i][j][k-1].rho*species[a][i][j][k-1].w*(1-3*gas_const_a*mixture[i][j][k-1].temp)-species[a][i][j][k-1].rho*cb(species[a][i][j][k-1].w), species[a][i][j][k].rho*species[a][i][j][k].w*(1-3*gas_const_a*mixture[i][j][k].temp)-species[a][i][j][k].rho*cb(species[a][i][j][k].w), species[a][i][j][k+1].rho*species[a][i][j][k+1].w*(1-3*gas_const_a*mixture[i][j][k+1].temp)-species[a][i][j][k+1].rho*cb(species[a][i][j][k+1].w), dz, species[a][i][j][k].w, mixture[i][j][k-1].type, mixture[i][j][k+1].type);
                        double corr[3] = {  dt_sim*(2-omega_a[a])/(2*species[a][i][j][k].rho*omega_a[a])*delQdevx,
                                            dt_sim*(2-omega_a[a])/(2*species[a][i][j][k].rho*omega_a[a])*delQdevy,
                                            dt_sim*(2-omega_a[a])/(2*species[a][i][j][k].rho*omega_a[a])*delQdevz};

                        double corr_zero[3] = {0.0, 0.0, 0.0};

                        for (int l = 0; l < npop; ++l)
                        {
                            double feq_a = calculate_feq(l, species[a][i][j][k].rho, velocity, theta, corr);
                            double feq_du_a = calculate_feq(l, species[a][i][j][k].rho, velocity_du, theta, corr);
                            #ifdef REACTION
                            double freact_a = calculate_feq(l, rho_a[a], velocity_mix, theta, corr_zero);
                            #endif

                            if (species[a][i][j][k].X == 0.0 && freact_a == 0.0) 
                                continue;
                            else if (species[a][i][j][k].X == 0.0)
                                species[a][i][j][k].fpc[l] = freact_a;
                            else
                                species[a][i][j][k].fpc[l] = species[a][i][j][k].f[l] + omega_a[a]*(feq_a-species[a][i][j][k].f[l]) + (feq_du_a - feq_a) + freact_a;

                        }
                    }

                  

                }
            }
        }
    }
}
#endif