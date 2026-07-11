#include "lbm.hpp"
#include "math_util.hpp"
#include "units.hpp"
#include "omp.h"
#include "FD.hpp"


#ifndef MULTICOMP
void LBM::Collide()
{    
     #ifdef PARALLEL 
        #pragma omp parallel for schedule(dynamic) 
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

                #ifdef CONJUGATE
                if (mixture[i][j][k].type == TYPE_S)     
                {   
                    double cv = gas_const / (gamma - 1.0);
                    double cp = cv + gas_const;
                    double theta = gas_const*mixture[i][j][k].temp;
                    double mu = nu*mixture[i][j][k].rho;
                    double conduc_coeff = 10.0*mu*cp/prtl;

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
                    
                    for(int p=0; p < 3; ++p)
                        for(int q=0; q < 3; ++q)
                            eq_p_tensor[p][q] = (p==q) ? mixture[i][j][k].p+mixture[i][j][k].rho*velocity[p]*velocity[q] : mixture[i][j][k].rho*velocity[p]*velocity[q]; 
                            
                    double q_boundary [3] = {0.0, 0.0, 0.0};
                    // int i_nb, j_nb, k_nb;
                    // for (int l=0; l < npop; ++l){
                    //     i_nb = i - cx[l];
                    //     j_nb = j - cy[l];
                    //     k_nb = k - cz[l];

                    //     if((mixture[i_nb][j_nb][k_nb].type==TYPE_Q ) && mixture[(int)(i+cx[l])][(int)(j+cy[l])][(int)(k+cz[l])].type==TYPE_S){    //|| mixture[i_nb][j_nb][k_nb].type==TYPE_Q
                    //         if (cx[l] != 0)
                    //             q_boundary[0] = mixture[i_nb][j_nb][k_nb].energy_flux[0];
                    //         if (cy[l] != 0)
                    //             q_boundary[1] = mixture[i_nb][j_nb][k_nb].energy_flux[1];
                    //         if (cz[l] != 0)
                    //             q_boundary[2] = mixture[i_nb][j_nb][k_nb].energy_flux[2];
                    //         break;
                    //     }                            
                    // }
                    // q_boundary[0] *= omega1/(-omega+omega1) ;
                    // q_boundary[1] *= omega1/(-omega+omega1) ;
                    // q_boundary[2] *= omega1/(-omega+omega1) ;

                    double delQdevx, delQdevy, delQdevz;
                    delQdevx = fd_central(mixture[i-1][j][k].rho*mixture[i-1][j][k].u*(1.0-3.0*gas_const*mixture[i-1][j][k].temp)-mixture[i-1][j][k].rho*cb(mixture[i-1][j][k].u), mixture[i][j][k].rho*mixture[i][j][k].u*(1.0-3.0*gas_const*mixture[i][j][k].temp)-mixture[i][j][k].rho*cb(mixture[i][j][k].u), mixture[i+1][j][k].rho*mixture[i+1][j][k].u*(1.0-3.0*gas_const*mixture[i+1][j][k].temp)-mixture[i+1][j][k].rho*cb(mixture[i+1][j][k].u), dx, mixture[i][j][k].u, mixture[i-1][j][k].type, mixture[i+1][j][k].type) ;
                    delQdevy = fd_central(mixture[i][j-1][k].rho*mixture[i][j-1][k].v*(1.0-3.0*gas_const*mixture[i][j-1][k].temp)-mixture[i][j-1][k].rho*cb(mixture[i][j-1][k].v), mixture[i][j][k].rho*mixture[i][j][k].v*(1.0-3.0*gas_const*mixture[i][j][k].temp)-mixture[i][j][k].rho*cb(mixture[i][j][k].v), mixture[i][j+1][k].rho*mixture[i][j+1][k].v*(1.0-3.0*gas_const*mixture[i][j+1][k].temp)-mixture[i][j+1][k].rho*cb(mixture[i][j+1][k].v), dy, mixture[i][j][k].v, mixture[i][j-1][k].type, mixture[i][j+1][k].type) ;
                    delQdevz = fd_central(mixture[i][j][k-1].rho*mixture[i][j][k-1].w*(1.0-3.0*gas_const*mixture[i][j][k-1].temp)-mixture[i][j][k-1].rho*cb(mixture[i][j][k-1].w), mixture[i][j][k].rho*mixture[i][j][k].w*(1.0-3.0*gas_const*mixture[i][j][k].temp)-mixture[i][j][k].rho*cb(mixture[i][j][k].w), mixture[i][j][k+1].rho*mixture[i][j][k+1].w*(1.0-3.0*gas_const*mixture[i][j][k+1].temp)-mixture[i][j][k+1].rho*cb(mixture[i][j][k+1].w), dz, mixture[i][j][k].w, mixture[i][j][k-1].type, mixture[i][j][k+1].type) ;

                    double d_str_heat_flux[3] ={  velocity[0]*(mixture[i][j][k].p_tensor[0][0]-eq_p_tensor[0][0]) + velocity[1]*(mixture[i][j][k].p_tensor[1][0]-eq_p_tensor[1][0]) + velocity[2]*(mixture[i][j][k].p_tensor[2][0]-eq_p_tensor[2][0]) + 0.5*dt_sim*velocity[0]*delQdevx + q_boundary[0],   //+ 0.5*dt_sim*velocity[0]*delQdevx + 0.5*dt_sim*velocity[0]*mixture[i][j][k].dQdevx
                                                  velocity[0]*(mixture[i][j][k].p_tensor[0][1]-eq_p_tensor[0][1]) + velocity[1]*(mixture[i][j][k].p_tensor[1][1]-eq_p_tensor[1][1]) + velocity[2]*(mixture[i][j][k].p_tensor[2][1]-eq_p_tensor[2][1]) + 0.5*dt_sim*velocity[1]*delQdevy + q_boundary[1],   //+ 0.5*dt_sim*velocity[1]*delQdevy + 0.5*dt_sim*velocity[1]*mixture[i][j][k].dQdevy
                                                  velocity[0]*(mixture[i][j][k].p_tensor[0][2]-eq_p_tensor[0][2]) + velocity[1]*(mixture[i][j][k].p_tensor[1][2]-eq_p_tensor[1][2]) + velocity[2]*(mixture[i][j][k].p_tensor[2][2]-eq_p_tensor[2][2]) + 0.5*dt_sim*velocity[2]*delQdevz + q_boundary[2]} ; //+ 0.5*dt_sim*velocity[2]*delQdevz + 0.5*dt_sim*velocity[2]*mixture[i][j][k].dQdevz

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
                #endif


            }
        }
    }
}

#elif defined MULTICOMP
void LBM::Collide()
{
    #ifndef ISOTHERM

    #ifdef PARALLEL
        #pragma omp parallel
    #endif
    {
    // thread-local Cantera objects and work buffers, reused for every cell
    int rank = omp_get_thread_num();
    auto gas = sols[rank]->thermo();
    auto trans = sols[rank]->transport();
    std::vector <double> Y (nSpeciesCantera);
    std::vector <double> part_enthalpy(nSpeciesCantera);

    // specific gas constant R/W_mean from the local composition only;
    // equivalent to Cantera's meanMolecularWeight() without touching the thermo state
    auto mix_gasconst = [&](int ii, int jj, int kk) -> double {
        double sumYoverW = 0.0, sumY = 0.0;
        for(size_t a = 0; a < nSpecies; ++a){
            double Ya = species[a][ii][jj][kk].rho / mixture[ii][jj][kk].rho;
            sumYoverW += Ya / molarMass[a];
            sumY += Ya;
        }
        return units.cp(Cantera::GasConstant * sumYoverW / sumY);
    };

    auto is_fluidlike = [&](int ii, int jj, int kk) -> bool {
        short t = mixture[ii][jj][kk].type;
        return t == TYPE_F || t == TYPE_I || t == TYPE_O || t == TYPE_I_C || t == TYPE_O_C;
    };

    #ifdef PARALLEL
        #pragma omp for collapse(2) schedule(dynamic)
    #endif
    for(int i = 0; i < Nx ; ++i)
    {
        for(int j = 0; j < Ny; ++j)
        {
            for(int k = 0; k < Nz; ++k)
            {
                if (mixture[i][j][k].type == TYPE_F)
                {
                    std::fill(Y.begin(), Y.end(), 0.0);
                    for(size_t a = 0; a < nSpecies; ++a) Y[speciesIdx[a]] = species[a][i][j][k].rho / mixture[i][j][k].rho;
                    gas->setMassFractions(&Y[0]);
                    gas->setState_TD(units.si_temp(mixture[i][j][k].temp), units.si_rho(mixture[i][j][k].rho));

                    // get macropscopic properties
                    double cp = units.cp(gas->cp_mass());
                    double gasconstant = units.cp(Cantera::GasConstant/gas->meanMolecularWeight());
                    gas->getPartialMolarEnthalpies(&part_enthalpy[0]);

                    double mu = units.mu(trans->viscosity());
                    double conduc_coeff = units.thermalConductivity(trans->thermalConductivity());
                    
                    double theta = units.energy_mass(gas->RT()/gas->meanMolecularWeight());

                    double velocity[3] = {  mixture[i][j][k].u,
                                            mixture[i][j][k].v, 
                                            mixture[i][j][k].w};
                
                    double omega1 = 2*mixture[i][j][k].p*cp*dt_sim / (mixture[i][j][k].p*cp*dt_sim + 2*conduc_coeff);
                    double omega = 2*mixture[i][j][k].p*dt_sim / (mixture[i][j][k].p*dt_sim + 2*mu);

                    // std::cout << 1/omega << " | " << 1/omega1 << std::endl;

                    double internal_energy = mixture[i][j][k].rhoe/mixture[i][j][k].rho; 
                    for (size_t a = 0; a < nSpecies; ++a){
                        internal_energy = internal_energy - 0.5 * species[a][i][j][k].rho/mixture[i][j][k].rho * v_sqr(species[a][i][j][k].u, species[a][i][j][k].v, species[a][i][j][k].w);
                    }     
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
                        double h_a = units.energy_mass(part_enthalpy[speciesIdx[a]]/molarMass[a]);
                        q_diff[0] += omega1/(-omega+omega1)  * h_a * species[a][i][j][k].rho * (species[a][i][j][k].u-mixture[i][j][k].u);
                        q_diff[1] += omega1/(-omega+omega1)  * h_a * species[a][i][j][k].rho * (species[a][i][j][k].v-mixture[i][j][k].v);
                        q_diff[2] += omega1/(-omega+omega1)  * h_a * species[a][i][j][k].rho * (species[a][i][j][k].w-mixture[i][j][k].w);

                        double dY_x = fd_fuw(species[a][i-1][j][k].rho/mixture[i-1][j][k].rho, species[a][i][j][k].rho/mixture[i][j][k].rho, species[a][i+1][j][k].rho/mixture[i+1][j][k].rho, dx, species[a][i][j][k].u-mixture[i][j][k].u, mixture[i-1][j][k].type, mixture[i+1][j][k].type);
                        double dY_y = fd_fuw(species[a][i][j-1][k].rho/mixture[i][j-1][k].rho, species[a][i][j][k].rho/mixture[i][j][k].rho, species[a][i][j+1][k].rho/mixture[i][j+1][k].rho, dy, species[a][i][j][k].v-mixture[i][j][k].v, mixture[i][j-1][k].type, mixture[i][j+1][k].type);
                        double dY_z = fd_fuw(species[a][i][j][k-1].rho/mixture[i][j][k-1].rho, species[a][i][j][k].rho/mixture[i][j][k].rho, species[a][i][j][k+1].rho/mixture[i][j][k+1].rho, dz, species[a][i][j][k].w-mixture[i][j][k].w, mixture[i][j][k-1].type, mixture[i][j][k+1].type);

                        q_corr[0] += 0.5 * (2.0-omega1)/(-omega+omega1) * dt_sim * mixture[i][j][k].p * h_a * dY_x;//species[a][i][j][k].delYx;
                        q_corr[1] += 0.5 * (2.0-omega1)/(-omega+omega1) * dt_sim * mixture[i][j][k].p * h_a * dY_y;//species[a][i][j][k].delYy;
                        q_corr[2] += 0.5 * (2.0-omega1)/(-omega+omega1) * dt_sim * mixture[i][j][k].p * h_a * dY_z;//species[a][i][j][k].delYz;
                    }

                    double delQdevx, delQdevy, delQdevz;
                    double gasconst_in = 0.0, gasconst_ip = 0.0, gasconst_jn = 0.0, gasconst_jp = 0.0, gasconst_kn = 0.0, gasconst_kp = 0.0;

                    // gas constant of the neighbours, from composition only (no Cantera state updates)
                    if(is_fluidlike(i-1,j,k)) gasconst_in = mix_gasconst(i-1,j,k);
                    if(is_fluidlike(i+1,j,k)) gasconst_ip = mix_gasconst(i+1,j,k);
                    if(is_fluidlike(i,j-1,k)) gasconst_jn = mix_gasconst(i,j-1,k);
                    if(is_fluidlike(i,j+1,k)) gasconst_jp = mix_gasconst(i,j+1,k);
                    if(is_fluidlike(i,j,k-1)) gasconst_kn = mix_gasconst(i,j,k-1);
                    if(is_fluidlike(i,j,k+1)) gasconst_kp = mix_gasconst(i,j,k+1);

                    delQdevx = fd_central(mixture[i-1][j][k].rho*mixture[i-1][j][k].u*(1.0-3.0*gasconst_in*mixture[i-1][j][k].temp)-mixture[i-1][j][k].rho*cb(mixture[i-1][j][k].u), mixture[i][j][k].rho*mixture[i][j][k].u*(1.0-3.0*gasconstant*mixture[i][j][k].temp)-mixture[i][j][k].rho*cb(mixture[i][j][k].u), mixture[i+1][j][k].rho*mixture[i+1][j][k].u*(1.0-3.0*gasconst_ip*mixture[i+1][j][k].temp)-mixture[i+1][j][k].rho*cb(mixture[i+1][j][k].u), dx, mixture[i][j][k].u, mixture[i-1][j][k].type, mixture[i+1][j][k].type) ;
                    delQdevy = fd_central(mixture[i][j-1][k].rho*mixture[i][j-1][k].v*(1.0-3.0*gasconst_jn*mixture[i][j-1][k].temp)-mixture[i][j-1][k].rho*cb(mixture[i][j-1][k].v), mixture[i][j][k].rho*mixture[i][j][k].v*(1.0-3.0*gasconstant*mixture[i][j][k].temp)-mixture[i][j][k].rho*cb(mixture[i][j][k].v), mixture[i][j+1][k].rho*mixture[i][j+1][k].v*(1.0-3.0*gasconst_jp*mixture[i][j+1][k].temp)-mixture[i][j+1][k].rho*cb(mixture[i][j+1][k].v), dy, mixture[i][j][k].v, mixture[i][j-1][k].type, mixture[i][j+1][k].type) ;
                    delQdevz = fd_central(mixture[i][j][k-1].rho*mixture[i][j][k-1].w*(1.0-3.0*gasconst_kn*mixture[i][j][k-1].temp)-mixture[i][j][k-1].rho*cb(mixture[i][j][k-1].w), mixture[i][j][k].rho*mixture[i][j][k].w*(1.0-3.0*gasconstant*mixture[i][j][k].temp)-mixture[i][j][k].rho*cb(mixture[i][j][k].w), mixture[i][j][k+1].rho*mixture[i][j][k+1].w*(1.0-3.0*gasconst_kp*mixture[i][j][k+1].temp)-mixture[i][j][k+1].rho*cb(mixture[i][j][k+1].w), dz, mixture[i][j][k].w, mixture[i][j][k-1].type, mixture[i][j][k+1].type) ;

                    double d_str_heat_flux[3] ={    velocity[0]*(mixture[i][j][k].p_tensor[0][0]-eq_p_tensor[0][0]) + velocity[1]*(mixture[i][j][k].p_tensor[1][0]-eq_p_tensor[1][0]) + velocity[2]*(mixture[i][j][k].p_tensor[2][0]-eq_p_tensor[2][0]) + q_diff[0] + q_corr[0] + 0.5*dt_sim*velocity[0]*delQdevx,   //+ q_diff[0] + q_corr[0]  + 0.5*dt_sim*velocity[0]*delQdevx
                                                    velocity[0]*(mixture[i][j][k].p_tensor[0][1]-eq_p_tensor[0][1]) + velocity[1]*(mixture[i][j][k].p_tensor[1][1]-eq_p_tensor[1][1]) + velocity[2]*(mixture[i][j][k].p_tensor[2][1]-eq_p_tensor[2][1]) + q_diff[1] + q_corr[1] + 0.5*dt_sim*velocity[1]*delQdevy,   //+ q_diff[1] + q_corr[1]  + 0.5*dt_sim*velocity[1]*delQdevy
                                                    velocity[0]*(mixture[i][j][k].p_tensor[0][2]-eq_p_tensor[0][2]) + velocity[1]*(mixture[i][j][k].p_tensor[1][2]-eq_p_tensor[1][2]) + velocity[2]*(mixture[i][j][k].p_tensor[2][2]-eq_p_tensor[2][2]) + q_diff[2] + q_corr[2] + 0.5*dt_sim*velocity[2]*delQdevz} ; //+ q_diff[2] + q_corr[2]  + 0.5*dt_sim*velocity[2]*delQdevz
                
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
    } // end parallel region
    #endif
}

void LBM::Collide_Species()
{
    #ifdef PARALLEL
        #pragma omp parallel
    #endif
    {
    // thread-local Cantera objects and work buffers, reused for every cell
    int rank = omp_get_thread_num();
    auto gas = sols[rank]->thermo();
    auto trans = sols[rank]->transport();
    auto kinetics = sols[rank]->kinetics();
    const int ld = nSpeciesCantera;
    std::vector <double> Y (nSpeciesCantera);
    std::vector <double> w_dot(nSpeciesCantera);    // mole density rate [kmol/m3/s]
    std::vector <double> part_enthalpy(nSpeciesCantera);
    std::vector <double> spec_visc(nSpeciesCantera);
    std::vector <double> d(ld * ld);
    std::vector <double> rho_a(nSpecies);

    #ifdef PARALLEL
        #pragma omp for collapse(2) schedule(dynamic)
    #endif
    for(int i = 0; i < Nx ; ++i)
    {
        for(int j = 0; j < Ny; ++j)
        {
            for(int k = 0; k < Nz; ++k)
            {
                if (mixture[i][j][k].type == TYPE_F)
                {
                    std::fill(Y.begin(), Y.end(), 0.0);
                    for(size_t a = 0; a < nSpecies; ++a) Y[speciesIdx[a]] = species[a][i][j][k].rho / mixture[i][j][k].rho;
                    gas->setMassFractions(&Y[0]);
                    gas->setState_TD(units.si_temp(mixture[i][j][k].temp), units.si_rho(mixture[i][j][k].rho));
                                        
                    // Reaction ---------------------------------------------------------------------------------------------------------------------------
                    #ifdef REACTION
                        double velocity_mix[3] = {  mixture[i][j][k].u,
                                                    mixture[i][j][k].v, 
                                                    mixture[i][j][k].w};      
                        double internal_energy = mixture[i][j][k].rhoe/mixture[i][j][k].rho; 
                        for (size_t a = 0; a < nSpecies; ++a){
                            internal_energy = internal_energy - 0.5 * species[a][i][j][k].rho/mixture[i][j][k].rho * v_sqr(species[a][i][j][k].u, species[a][i][j][k].v, species[a][i][j][k].w);
                        }
                        size_t p = 0;
                        double maxIter = 2.0;
                        std::fill(rho_a.begin(), rho_a.end(), 0.0);
                        do{
                            kinetics->getNetProductionRates(&w_dot[0]);
                            for (size_t a = 0; a < nSpecies; ++a){
                                double rho_dot = units.rho_dot(w_dot[speciesIdx[a]] * molarMass[a]);
                                rho_a[a] += dt_sim/maxIter * rho_dot;
                            }

                            std::fill(Y.begin(), Y.end(), 0.0);
                            for(size_t a = 0; a < nSpecies; ++a) Y[speciesIdx[a]] = (species[a][i][j][k].rho+rho_a[a]) / mixture[i][j][k].rho;
                            gas->setMassFractions(&Y[0]);
                            gas->setState_UV(units.si_energy_mass(internal_energy), 1.0/units.si_rho(mixture[i][j][k].rho),1.0e-15);
                            p++;
                        }while(p < maxIter);

                        gas->getPartialMolarEnthalpies(&part_enthalpy[0]);
                        mixture[i][j][k].HRR = 0.0;
                        for (size_t a = 0; a < nSpecies; ++a)
                            mixture[i][j][k].HRR += -1.0 * units.energy_mass(part_enthalpy[speciesIdx[a]] / molarMass[a]) * rho_a[a];
                    #endif

                    // Viscosity -----------------------------------------------------------------------------------------------------------------------------
                    std::fill(Y.begin(), Y.end(), 0.0);
                    for(size_t a = 0; a < nSpecies; ++a) Y[speciesIdx[a]] = species[a][i][j][k].rho / mixture[i][j][k].rho;
                    gas->setMassFractions(&Y[0]);
                    gas->setState_TD(units.si_temp(mixture[i][j][k].temp), units.si_rho(mixture[i][j][k].rho));

                    trans->getSpeciesViscosities(&spec_visc[0]);
                    double visc_a[nSpecies]; std::fill_n(&visc_a[0], nSpecies, 0.0);   // (clang: VLAs cannot have initializers)
                    double omega_a[nSpecies]; std::fill_n(&omega_a[0], nSpecies, 0.0);   // (clang: VLAs cannot have initializers)

                    double phi[nSpecies][nSpecies]; std::fill_n(&phi[0][0], nSpecies*nSpecies, 0.0);   // (clang: VLAs cannot have initializers)
                    for(size_t a = 0; a < nSpecies; ++a){
                        for(size_t b = a; b < nSpecies; ++b){
                            size_t idx_a = speciesIdx[a];
                            size_t idx_b = speciesIdx[b];

                            double factor1 = 1.0 + sqrt(spec_visc[idx_a]/spec_visc[idx_b]) * sqrt(sqrt(molarMass[b]/molarMass[a]));
                            phi[a][b] = factor1*factor1 / sqrt(8.0*(1.0 + molarMass[a]/molarMass[b]));
                            phi[b][a] = spec_visc[idx_b]/spec_visc[idx_a] * molarMass[a]/molarMass[b] * phi[a][b];
                        }
                    }

                    for(size_t a = 0; a < nSpecies; ++a){
                        double denom = 0.0;
                        for(size_t b = 0; b < nSpecies; ++b)
                            denom = denom + species[b][i][j][k].X * phi[a][b];

                        visc_a[a] = species[a][i][j][k].X * units.mu(spec_visc[speciesIdx[a]]) / denom;
                        omega_a[a] = 2*species[a][i][j][k].X*mixture[i][j][k].p*dt_sim / (species[a][i][j][k].X*mixture[i][j][k].p*dt_sim + 2*visc_a[a]);
                    }

                    // Diffusion ----------------------------------------------------------------------------------------------------------------------------
                    trans->getBinaryDiffCoeffs(ld, &d[0]);

                    double D_ab[nSpecies][nSpecies];
                    for(size_t a = 0; a < nSpecies; ++a)
                    {
                        for(size_t b = 0; b <= a; ++b)
                        {
                            // get diffusion coefficient
                            D_ab[a][b] =  units.nu( d[ld*speciesIdx[b] + speciesIdx[a]] );
                            D_ab[b][a] =  D_ab[a][b];
                        }
                    }

                    double ux[nSpecies]; std::fill_n(&ux[0], nSpecies, 0.0);   // (clang: VLAs cannot have initializers)
                    double uy[nSpecies]; std::fill_n(&uy[0], nSpecies, 0.0);   // (clang: VLAs cannot have initializers)
                    double uz[nSpecies]; std::fill_n(&uz[0], nSpecies, 0.0);   // (clang: VLAs cannot have initializers)
                    double fx[nSpecies]; std::fill_n(&fx[0], nSpecies, 0.0);   // (clang: VLAs cannot have initializers)
                    double fy[nSpecies]; std::fill_n(&fy[0], nSpecies, 0.0);   // (clang: VLAs cannot have initializers)
                    double fz[nSpecies]; std::fill_n(&fz[0], nSpecies, 0.0);   // (clang: VLAs cannot have initializers)
                    // Cap for the species drift / diffusion velocity. The species
                    // equilibrium is centred on the species advection velocity; when a
                    // species vanishes (Y_a -> 0) the raw ratio (sum_l f_a c_l)/rho_a
                    // becomes an ill-posed 0/0 that can exceed the lattice range and
                    // drive the product-form equilibrium negative. Following Sawant &
                    // Karlin, a vanishing species is advected with the mixture velocity
                    // and every drift is capped near the sound speed, so the equilibrium
                    // never sees an unbounded velocity. Resolved species (X_a >= SPECIES_MIN)
                    // are unaffected: the cap is inactive for any stable state.
                    const double Vmax = 0.5/sqrt(3.0);   // ~ lattice sound speed
                    for(size_t a = 0; a < nSpecies; ++a){

                        if (species[a][i][j][k].X < SPECIES_MIN){
                            // trace species: advect with the mixture, no diffusion drift
                            ux[a] = mixture[i][j][k].u;
                            uy[a] = mixture[i][j][k].v;
                            uz[a] = mixture[i][j][k].w;
                            fx[a] = 0.0; fy[a] = 0.0; fz[a] = 0.0;
                            continue;
                        }

                        for(size_t l = 0; l < npop; ++l){
                            ux[a] += species[a][i][j][k].f[l]*cx[l];
                            uy[a] += species[a][i][j][k].f[l]*cy[l];
                            uz[a] += species[a][i][j][k].f[l]*cz[l];
                        }
                        ux[a] = ux[a] / species[a][i][j][k].rho;
                        uy[a] = uy[a] / species[a][i][j][k].rho;
                        uz[a] = uz[a] / species[a][i][j][k].rho;
                        // cap the drift relative to the mixture velocity (normally inactive)
                        ux[a] = mixture[i][j][k].u + fmax(-Vmax, fmin(Vmax, ux[a]-mixture[i][j][k].u));
                        uy[a] = mixture[i][j][k].v + fmax(-Vmax, fmin(Vmax, uy[a]-mixture[i][j][k].v));
                        uz[a] = mixture[i][j][k].w + fmax(-Vmax, fmin(Vmax, uz[a]-mixture[i][j][k].w));

                        for(size_t b = 0; b < nSpecies; ++b){
                            fx[a] += -1.0 * mixture[i][j][k].p * species[a][i][j][k].X*species[b][i][j][k].X/D_ab[a][b] * (species[a][i][j][k].u - species[b][i][j][k].u);
                            fy[a] += -1.0 * mixture[i][j][k].p * species[a][i][j][k].X*species[b][i][j][k].X/D_ab[a][b] * (species[a][i][j][k].v - species[b][i][j][k].v);
                            fz[a] += -1.0 * mixture[i][j][k].p * species[a][i][j][k].X*species[b][i][j][k].X/D_ab[a][b] * (species[a][i][j][k].w - species[b][i][j][k].w);
                        }
                        fx[a] = fmax(-Vmax, fmin(Vmax, fx[a]*dt_sim/species[a][i][j][k].rho));
                        fy[a] = fmax(-Vmax, fmin(Vmax, fy[a]*dt_sim/species[a][i][j][k].rho));
                        fz[a] = fmax(-Vmax, fmin(Vmax, fz[a]*dt_sim/species[a][i][j][k].rho));
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

                        double gas_const_a = units.cp(Cantera::GasConstant/molarMass[a]);
                        double theta = gas_const_a * mixture[i][j][k].temp;

                        double delQdevx, delQdevy, delQdevz;                            
                        delQdevx = fd_central(species[a][i-1][j][k].rho*species[a][i-1][j][k].u*(1.0-3.0*gas_const_a*mixture[i-1][j][k].temp)-species[a][i-1][j][k].rho*cb(species[a][i-1][j][k].u), species[a][i][j][k].rho*species[a][i][j][k].u*(1.0-3.0*gas_const_a*mixture[i][j][k].temp)-species[a][i][j][k].rho*cb(species[a][i][j][k].u), species[a][i+1][j][k].rho*species[a][i+1][j][k].u*(1.0-3.0*gas_const_a*mixture[i+1][j][k].temp)-species[a][i+1][j][k].rho*cb(species[a][i+1][j][k].u), dx, species[a][i][j][k].u, mixture[i-1][j][k].type, mixture[i+1][j][k].type);
                        delQdevy = fd_central(species[a][i][j-1][k].rho*species[a][i][j-1][k].v*(1.0-3.0*gas_const_a*mixture[i][j-1][k].temp)-species[a][i][j-1][k].rho*cb(species[a][i][j-1][k].v), species[a][i][j][k].rho*species[a][i][j][k].v*(1.0-3.0*gas_const_a*mixture[i][j][k].temp)-species[a][i][j][k].rho*cb(species[a][i][j][k].v), species[a][i][j+1][k].rho*species[a][i][j+1][k].v*(1.0-3.0*gas_const_a*mixture[i][j+1][k].temp)-species[a][i][j+1][k].rho*cb(species[a][i][j+1][k].v), dy, species[a][i][j][k].v, mixture[i][j-1][k].type, mixture[i][j+1][k].type);
                        delQdevz = fd_central(species[a][i][j][k-1].rho*species[a][i][j][k-1].w*(1.0-3.0*gas_const_a*mixture[i][j][k-1].temp)-species[a][i][j][k-1].rho*cb(species[a][i][j][k-1].w), species[a][i][j][k].rho*species[a][i][j][k].w*(1.0-3.0*gas_const_a*mixture[i][j][k].temp)-species[a][i][j][k].rho*cb(species[a][i][j][k].w), species[a][i][j][k+1].rho*species[a][i][j][k+1].w*(1.0-3.0*gas_const_a*mixture[i][j][k+1].temp)-species[a][i][j][k+1].rho*cb(species[a][i][j][k+1].w), dz, species[a][i][j][k].w, mixture[i][j][k-1].type, mixture[i][j][k+1].type);
                        // third-order-moment correction; division by rho_a makes it
                        // blow up for trace species, so it is disabled there (the
                        // species carries no dynamics worth correcting in that limit)
                        double corr[3] = {0.0, 0.0, 0.0};
                        if (species[a][i][j][k].X >= SPECIES_MIN){
                            corr[0] = dt_sim*(2.0-omega_a[a])/(2.0*species[a][i][j][k].rho*omega_a[a])*delQdevx;
                            corr[1] = dt_sim*(2.0-omega_a[a])/(2.0*species[a][i][j][k].rho*omega_a[a])*delQdevy;
                            corr[2] = dt_sim*(2.0-omega_a[a])/(2.0*species[a][i][j][k].rho*omega_a[a])*delQdevz;
                        }

                        #ifdef REACTION
                        double corr_zero[3] = {0.0, 0.0, 0.0};
                        #endif

                        for (int l = 0; l < npop; ++l)
                        {
                            double feq_a = calculate_feq(l, species[a][i][j][k].rho, velocity, theta, corr);
                            double feq_du_a = calculate_feq(l, species[a][i][j][k].rho, velocity_du, theta, corr);
                            double freact_a = 0.0;
                            #ifdef REACTION
                            freact_a = calculate_feq(l, rho_a[a], velocity_mix, theta, corr_zero);
                            #endif

                            if (species[a][i][j][k].rho <= SPECIES_MIN && freact_a <= 0.0) 
                                continue;
                            // else if (species[a][i][j][k].rho == SPECIES_MIN)
                            //     species[a][i][j][k].fpc[l] = freact_a;
                            else
                                species[a][i][j][k].fpc[l] = species[a][i][j][k].f[l] + omega_a[a]*(feq_a-species[a][i][j][k].f[l]) + (feq_du_a - feq_a) + freact_a;
                        }
                    }

                }
            }
        }
    }
    } // end parallel region
}
#endif
