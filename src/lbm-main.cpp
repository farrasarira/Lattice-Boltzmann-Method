
#include "cantera.hpp"
#include "lbm.hpp"
#include "math_util.hpp"
#include "units.hpp"
#include "restart_file.hpp"
#include <omp.h>
#include <numeric>
#include <vector>

#ifdef MULTICOMP
    std::vector<std::shared_ptr<Cantera::Solution>> sols;
#endif

#ifndef MULTICOMP
LBM::LBM(int Nx, int Ny, int Nz, double nu)
{
    // Number lattice used in simulation
    this->Nx = Nx + 2;  // + 2 for ghost lattice in the boundary
    this->Ny = Ny + 2;
    this->Nz = Nz + 2;
    this->nu = nu;

    // allocate memory for lattice
    mixture = new MIXTURE **[this->Nx];
    for (int i = 0; i < this->Nx; ++i)
    {
        mixture[i] = new MIXTURE *[this->Ny];
        for (int j = 0; j < this->Ny; ++j)
        {
            mixture[i][j] = new MIXTURE [this->Nz];
        }
    }
}
#endif

#ifdef MULTICOMP
LBM::LBM(int Nx, int Ny, int Nz, std::vector<std::string> species)
{
    // Number lattice used in simulation
    this->Nx = Nx + 2;  // + 2 for ghost lattice in the boundary
    this->Ny = Ny + 2;
    this->Nz = Nz + 2;
    this->speciesName = species;
    this->nSpecies = species.size();

    // allocate memory for mixture
    mixture = new MIXTURE **[this->Nx];
    for (int i = 0; i < this->Nx; ++i)
    {
        mixture[i] = new MIXTURE *[this->Ny];
        for (int j = 0; j < this->Ny; ++j)
        {
            mixture[i][j] = new MIXTURE [this->Nz];
        }
    }

    // allocate memory for species
    this->species.resize(nSpecies);
    for(size_t p = 0; p < nSpecies; ++p)
    {
        this->species[p] = new SPECIES **[this->Nx];
        for (int i = 0; i < this->Nx; ++i)
        {
            this->species[p][i] = new SPECIES *[this->Ny];
            for (int j = 0; j < this->Ny; ++j)
            {
                this->species[p][i][j] = new SPECIES [this->Nz];
            }
        }
    }

    // Create Cantera's Solution object
    int nThreads = omp_get_max_threads();
    std::cout << "nThreads : " << nThreads << std::endl;
    for(int i = 0; i < nThreads; ++i)
    {
        // auto sol = Cantera::newSolution("gri30.yaml", "gri30","mixture-averaged");
        // auto sol = Cantera::newSolution("h2o2.yaml", "ohmech");
        // auto sol = Cantera::newSolution("gri30.yaml", "gri30", "multicomponent");
        // auto sol = Cantera::newSolution("./src/reaction-mech/one-step.yaml", "FakeGas");
        auto sol = Cantera::newSolution("./src/reaction-mech/CH4_2S.yaml", "CH4_BFER_multi");
        // auto sol = Cantera::newSolution("./src/reaction-mech/propane_mech.yaml");
        sols.push_back(sol);
    }

}
#endif

LBM::~LBM(){

}


void LBM::run(int nstep, int tout)
{ 
    std::cout << "  Setup Done" << std::endl;

    // initialize the distribution function 
    std::cout << "  Initialization ..." << std::endl;
    step = 0;
    #ifdef SMOOTHING
        Init_smooting();  
    #else
        Init(); 
    #endif
    std::cout << "  Initialization Done" << std::endl;
    
    // Save the macroscopic at t=0
    OutputVTK(step, this);
    OutputKeEns(step, this);

    // Simulation loop
    loop(nstep, tout);
    
}

void LBM::loop(int nstep, int tout)
{ 

    // Simulation loop
    for (int step = this->step+1; step <= nstep; ++step)
    {
        // double start = omp_get_wtime();
        #ifdef MULTICOMP
        Collide_Species();  // collide species distribution function
        std::cout << "  Species Collision Done" << std::endl;
        #endif
        // std::cout << "Collision species     : " << double(omp_get_wtime()-start) << " seconds" << std::endl;

        // start = omp_get_wtime();
        Collide();   // collision step
        std::cout << "  Mixture Collision Done" << std::endl;
        // std::cout << "Collision mixture     : " << double(omp_get_wtime()-start) << " seconds" << std::endl;

        // start = omp_get_wtime();
        Streaming();        // streaming step & BC
        std::cout << "  Streaming Done" << std::endl;
        // std::cout << "Streaming process     : " << double(omp_get_wtime()-start) << " seconds" << std::endl;

        // start = omp_get_wtime();
        TMS_BC();
        std::cout << "  Apply BC Done" << std::endl;  
        // std::cout << "TMS BC                : " << double(omp_get_wtime()-start) << " seconds" << std::endl;
                
        // start = omp_get_wtime();
        #ifdef SMOOTHING
            calculate_moment_smoothing(); // calculate moment
        #else
            calculate_moment(); // calculate moment
        #endif
        std::cout << "  Calculate Moment Done" << std::endl;
        // std::cout << "Moment                : " << double(omp_get_wtime()-start) << " seconds" << std::endl;
        

        if (step % tout == 0){
            OutputVTK(step, this);      // Save the macroscopic quantity
            OutputKeEns(step, this);
        }
        if (step % (20*tout) == 0)
            write_restart(step, this);

    }
    
}


// // --------------------------------------------------------------------------------------------------------------------------------------------------------------
