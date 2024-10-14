#ifndef LBM_H
    #define LBM_H

    #include <math.h>
    #include <iostream>
    #include <vector>
    #include <string>
    #include <memory>
    #include "defines.hpp"
    // #include <eigen3/Eigen/Sparse>
    // #include <eigen3/Eigen/QR>
    #include "output.hpp"
    #include "cantera.hpp"
 

    #if defined D1Q3
        const int ndim = 1;
        const int npop = 3;
        const double cx[3] = {0.0,+1.0,-1.0};
        const double cy[3] = {0.0, 0.0, 0.0};
        const double cz[3] = {0.0, 0.0, 0.0};
        // D2Q9 weight factor | cs^2 = 1/3
        const double wi[3]  = {4./6,1./6,1./6};

    #elif defined D2Q9
        // D2Q9 velocity sets
        const int ndim = 2;
        const int npop = 9;
        const double cx[9] = {0.0,+1.0,-1.0, 0.0, 0.0,+1.0,-1.0,+1.0,-1.0};
        const double cy[9] = {0.0, 0.0, 0.0,+1.0,-1.0,+1.0,-1.0,-1.0,+1.0};
        const double cz[9] = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
        // D2Q9 weight factor | cs^2 = 1/3
        const double wi[9]  = {4./9,1./9,1./9,1./9,1./9,1./36,1./36,1./36,1./36};
                                                                                                                 
    #elif defined D3Q15
        const int ndim = 3;
        const int npop = 15;
        const double cx[15] = {0.0,+1.0,-1.0, 0.0, 0.0, 0.0, 0.0,+1.0,+1.0,+1.0,+1.0,-1.0,-1.0,-1.0,-1.0};
        const double cy[15] = {0.0, 0.0, 0.0,+1.0,-1.0, 0.0, 0.0,+1.0,+1.0,-1.0,-1.0,+1.0,+1.0,-1.0,-1.0};
        const double cz[15] = {0.0, 0.0, 0.0, 0.0, 0.0,+1.0,-1.0,+1.0,-1.0,+1.0,-1.0,+1.0,-1.0,+1.0,-1.0};
        const double wi[15] = {2/9, 1/9, 1/9, 1/9, 1/9, 1/9, 1/9,1/72,1/72,1/72,1/72,1/72,1/72,1/72,1/72};
    #elif defined D3Q27
        const int ndim = 3;
        const int npop = 27;
        const double cx[27] = {  0.0, +1.0, -1.0,  0.0,  0.0,  0.0,  0.0, +1.0, -1.0, +1.0, -1.0,  0.0,  0.0, +1.0, -1.0, +1.0, -1.0,  0.0,  0.0,  +1.0,  -1.0,  +1.0,  -1.0,  +1.0,  -1.0,  -1.0,  +1.0};
        const double cy[27] = {  0.0,  0.0,  0.0, +1.0, -1.0,  0.0,  0.0, +1.0, -1.0,  0.0,  0.0, +1.0, -1.0, -1.0, +1.0,  0.0,  0.0, +1.0, -1.0,  +1.0,  -1.0,  +1.0,  -1.0,  -1.0,  +1.0,  +1.0,  -1.0};
        const double cz[27] = {  0.0,  0.0,  0.0,  0.0,  0.0, +1.0, -1.0,  0.0,  0.0, +1.0, -1.0, +1.0, -1.0,  0.0,  0.0, -1.0, +1.0, -1.0, +1.0,  +1.0,  -1.0,  -1.0,  +1.0,  +1.0,  -1.0,  +1.0,  -1.0};
        const double wi[27] = {8./27,2./27,2./27,2./27,2./27,2./27,2./27,1./54,1./54,1./54,1./54,1./54,1./54,1./54,1./54,1./54,1./54,1./54,1./54,1./216,1./216,1./216,1./216,1./216,1./216,1./216,1./216};
        const int opposite[27] = {0,2,1,4,3,6,5,8,7,10,9,12,11,14,13,16,15,18,17,20,19,22,21,24,23,26,25};
    #endif

    class MIXTURE
    {
        public:
            short type = TYPE_F;    // type of lattice (FLUID, SOLID, ... see setup for more)
            
            // ###### Momentum Kinetic Equation Parameter ######
            #ifndef MULTICOMP
                double f[npop], fpc[npop];  // distribution function, distribution function post collistion
            #endif

            double p_tensor[3][3] = {{0., 0., 0.},    // pressure tensor
                                        {0., 0., 0.},
                                        {0., 0., 0.}}; 

            double rho;              // macroscopic quantity
            double u = 0.0;          // velocity in x-direction
            double v = 0.0;          // velocity in y-direction
            double w = 0.0;          // velocity in z-direction
            double p = 1./3.;        // pressure

            // ####### Energy Kinetic Equation Parameter #######
            #ifndef ISOTHERM
            double g[npop], gpc[npop];  // energy distribution function, energy distrbution function post collision
            double rhoe;                // density * total energy (E)
            double energy_flux[3] = {0., 0., 0.}; // heat flux
            #endif

            double temp = 1./3.;        // temperature

            #ifdef MULTICOMP
                double HRR = 0.0;
            #endif
    };

    class SPECIES
    {
        public:
            // ###### Momentum Kinetic Equation Parameter ######
            double X = 0.0;             // mole fraction
            double rho = 0.0;           // density

            double f[npop], fpc[npop];  // distribution function, distribution function post collistion  
            double u = 0.0;     // velocity in x-direction
            double v = 0.0;     // velocity in y-direction
            double w = 0.0;     // velocity in z-direction

            double p_tensor[3][3] = {{0., 0., 0.},    // pressure tensor
                                     {0., 0., 0.},
                                     {0., 0., 0.}};

    };

    class LBM
    {
        private:
            int dx = 1;
            int dy = 1;
            int dz = 1;
            int Nx = 1;
            int Ny = 1;
            int Nz = 1;
            int dt_sim = 1.0;
            int step = 0;

            #ifndef MULTICOMP
            double nu = 0.001;      // kinematic viscosity
            double gas_const = 1.0; // gas constant
            double gamma = 1.4;     // gamma (Cp/Cv)
            double Ra = 1.0;  // Reyleigh-Benard Constant
            double prtl = 0.5;      // prantdl number
            #else
            size_t nSpecies = 0;
            std::vector<std::string> speciesName;
            double permeability = 999999;  // Permeability
            #endif
        
        public:
            MIXTURE *** mixture;
            #ifdef MULTICOMP
                std::vector<SPECIES***> species;
            #endif

        public:
            // constructor
            LBM(int Nx, int Ny, int Nz, double nu);
            LBM(int Nx, int Ny, int Nz, std::vector<std::string> species);

            // calculate moment
            void calculate_moment();
            void calculate_moment_smoothing(){
                calculate_moment();
                Smoothing();
                calculate_moment();
            }
            double calculate_temp(double U, double rho, double Y[]);

            // calculate equlibrium density
            double calculate_feq(int l, double rho, double velocity[], double theta,  double corr[]);
            double calculate_geq(int l, double rho, double U, double theta, double v[]);
            double calculate_geq(int l, double rhoe, double eq_heat_flux[], double eq_R_tensor[][3], double theta);
            double calculate_gstr(int l, double geq, double d_str_heat_flux[]);
            void calculate_feq_geq(double f_tgt[], double g_tgt[], double rho_bb, double vel_tgt[], double temp_tgt);
            void calculate_feq_geq(double fa_tgt[][npop], double g_tgt[], double rho_bb, double rhoa_bb[], double vel_tgt[], double vela_tgt[][3], double temp_tgt);
            void calculate_feq_geq(double fa_tgt[][npop], double g_tgt[], double rho_bb, double rhoa_bb[], double vel_tgt[], double temp_tgt);

            // initialize
            void Init();    // initialize equilibrium  
            void Init_smooting(){
                Init();  
                Smoothing();
                calculate_moment();
            }

            // collision operator
            void Collide();
            void Collide_Species();
            void FD_species();
            
            // Boundary Conditions
            void fill_BC();
            void TMS_BC();
            void dirSlip(int l, int i, int j, int k, int &lp, int &ip, int &jp, int &kp);
            void FD_BC();

            // stream
            void Streaming();   

            // run simulation
            void run(int nstep, int tout);
            void loop(int nstep, int tout);

            // smoothing high gradient value
            void Smoothing();
            
            // get private data
            int get_Nx(){return Nx;};
            int get_Ny(){return Ny;};
            int get_Nz(){return Nz;};
            int get_NX(){return Nx-2;};
            int get_NY(){return Ny-2;};
            int get_NZ(){return Nz-2;};
            int get_dx(){return dx;};
            int get_dy(){return dy;};
            int get_dz(){return dz;};
            double get_dtsim(){return dt_sim;};
            double get_step(){return step;};
            
            #ifndef MULTICOMP
                double get_nu(){return nu;};
                double get_gasconst(){return gas_const;};
                double get_gamma(){return gamma;};
                double get_Ra(){return Ra;};
                double get_prtl(){return prtl;};
                size_t get_size(){
                    size_t size_scalar_int      = 9 * sizeof(int);
                    size_t size_scalar_double   = 5 * sizeof(double);
                    size_t size_mixture         = Nx*Ny*Nz * sizeof(MIXTURE);

                    size_t size_total = size_scalar_int + size_scalar_double + size_mixture;

                    return size_total;
                }

            #elif defined MULTICOMP
                int get_nSpecies(){return nSpecies;};
                std::vector<std::string> get_speciesName(){return speciesName;};
                double get_permeability(){return permeability;};
                size_t get_size(){
                    size_t size_scalar_int      = 9 * sizeof(int);
                    size_t size_scalar_species  = sizeof(size_t) + nSpecies*sizeof(size_t);
                    for(size_t a = 0; a < nSpecies; ++a)
                        size_scalar_species += speciesName[a].size();
                    size_scalar_species += sizeof(double); // permeability
                    size_t size_mixture         = Nx*Ny*Nz * sizeof(MIXTURE);
                    size_t size_species         = Nx*Ny*Nz*nSpecies * sizeof(SPECIES);

                    size_t size_total = size_scalar_int + size_scalar_species + size_mixture + size_species;

                    return size_total;
                }
            #endif


            // set private data
            void set_dx(int dx){this->dx = dx;};
            void set_dy(int dy){this->dy = dy;};
            void set_dz(int dz){this->dz = dz;};
            void set_Nx(int Nx){this->Nx = Nx;};
            void set_Ny(int Ny){this->Ny = Ny;};
            void set_Nz(int Nz){this->Nz = Nz;};
            void set_dtsim(int dt_sim){this->dt_sim = dt_sim;};
            void set_step(int step){this->step = step;};
            void set_permeability(double permeability){this->permeability = permeability;};

            #ifndef MULTICOMP
                void set_nu(double nu){this->nu = nu;};
                void set_gasconst(double gas_const){this->gas_const = gas_const;};
                void set_prtl(double prtl){this->prtl = prtl;};
                void set_gamma(double gamma){this->gamma = gamma;};
                void set_Ra(double Ra){this->Ra = Ra;};
                double get_soundspeed(double temp){return sqrt(this->gamma*this->gas_const*temp);};
                double get_conduc_coeff(){return nu*1.0/prtl;};
            #endif

    };

#endif
