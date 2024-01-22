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

    class LATTICE
    {
        public:
            short type = TYPE_F;    // type of lattice (FLUID, SOLID, ... see setup for more)
            
            // ###### Momentum Kinetic Equation Parameter ######
            double f[npop], fpc[npop];  // distribution function, distribution function post collistion  
            double rho;                 // macroscopic quantity
            double u = 0.0;          // velocity in x-direction
            double v = 0.0;          // velocity in y-direction
            double w = 0.0;          // velocity in z-direction
            double p = 1./3.;           // pressure
            double p_tensor[3][3] = {{0., 0., 0.},    // pressure tensor
                                     {0., 0., 0.},
                                     {0., 0., 0.}};

            // ####### Energy Kinetic Equation Parameter #######
            double g[npop], gpc[npop];  // energy distribution function, energy distrbution function post collision
            double temp = 1./3.;        // temperature
            double rhoe;                // density * total energy (E)
            double energy_flux[3] = {0., 0., 0.}; // heat flux

            // Third-order Moment Deviation
            double dQdevx = 0.0;
            double dQdevy = 0.0;
            double dQdevz = 0.0;

    };

    class SPECIES
    {
        public:
            // ###### Momentum Kinetic Equation Parameter ######
            double X = 0.0;             // mole fraction
            double rho = 0.0;           // density
            // double rho_dot = 0.0;   // rate of formation/destruction during chemical reaction.
            double Vdiff_x = 0.0;       // diffusion velocity in x direction               
            double Vdiff_y = 0.0;       // diffusion velocity in y direction       
            double Vdiff_z = 0.0;       // diffusion velocity in z direction
            double delYx = 0.0;         // gradient of mass fraction in x direction
            double delYy = 0.0;         // gradient of mass fraction in y direction
            double delYz = 0.0;         // gradient of mass fraction in z direction

            #ifndef FD
                double f[npop], fpc[npop];  // distribution function, distribution function post collistion  
                double u = 0.0;     // velocity in x-direction
                double v = 0.0;     // velocity in y-direction
                double w = 0.0;     // velocity in z-direction
            #elif defined FD
                double rho_n = 0.0;        // mixture density * mass fraction
            #endif

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
            int nstep = 1000;
            int tout = 100;

            double nu = 0.001;      // kinematic viscosity
            double gas_const = 1.0; // gas constant
            double prtl = 0.5;      // prantdl number
            double gamma = 1.4;     // gamma (Cp/Cv)

            size_t nSpecies = 0;
            std::vector<std::string> speciesName;
        
        public:
            LATTICE *** mixture;
            #ifdef MULTICOMP
                std::vector<SPECIES***> species;
            #endif

        public:
            // constructor
            LBM(int Nx, int Ny, int Nz, double nu);
            LBM(int Nx, int Ny, int Nz, std::vector<std::string> species);

            // calculate moment
            void calculate_moment();
            void calculate_moment(SPECIES ***species);

            // calculate equlibrium density
            double calculate_feq(int l, double rho, double velocity[], double theta,  double corr[]);
            double calculate_geq(int l, double rhoe, double eq_heat_flux[], double eq_R_tensor[][3], double theta);

            // initialize
            void Init();    // initialize equilibrium  
            void Init(SPECIES*** species, int idxSpecies); 

            // collision operator
            void Collide(); // BGK collision
            void Collide_Species();
            void FD_species();
            void FD_BC();
            
            // stream
            void Streaming();   

            // run simulation
            void run(int nstep, int tout);
            
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
            double get_nu(){return nu;};
            int get_nSpecies(){return nSpecies;};
            std::vector<std::string> get_speciesName(){return speciesName;};

            // set private data
            void set_nu(double nu){this->nu = nu;};
            void set_gasconst(double gas_const){this->gas_const = gas_const;};
            void set_prtl(double prtl){this->prtl = prtl;};
            void set_gamma(double gamma){this->gamma = gamma;};
    };

#endif