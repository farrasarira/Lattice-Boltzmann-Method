
#include "cantera.hpp"
#include "lbm.hpp"
#include "math_util.hpp"
#include "units.hpp"
#include <omp.h>
#include <numeric>
#include <vector>

#ifdef MULTICOMP
    std::vector<std::shared_ptr<Cantera::Solution>> sols;
#endif

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
        auto sol = Cantera::newSolution("gri30.yaml", "gri30","mixture-averaged");
        // auto sol = Cantera::newSolution("gri30.yaml", "gri30", "multicomponent");
        sols.push_back(sol);
    }

}
#endif


void LBM::run(int nstep, int tout)
{ 
    this->nstep = nstep;
    this->tout = tout;

    std::cout << "  Setup Done" << std::endl;

    // initialize the distribution function 
    std::cout << "  Initialization ..." << std::endl;
    Init();  
    std::cout << "  Initialization Done" << std::endl;
    
    // initialize time step variable
    int step = 0;

    // Save the macroscopic at t=0
    OutputVTK(step, this);
    OutputKeEns(step, this);

    // Simulation loop
    for (step = 1; step <= nstep; ++step)
    {
        #ifdef MULTICOMP
        Collide_Species();  // collide species distribution function
        // std::cout << "  Species Collision Done" << std::endl;
        #endif

        Collide();   // collision step
        // std::cout << "  Mixture Collision Done" << std::endl;
        Streaming();        // streaming step & BC
        // std::cout << "  Streaming Done" << std::endl;
        calculate_moment(); // calculate moment
        // std::cout << "  Calculate Moment Done" << std::endl;

        if (step % tout == 0)
        {
            OutputVTK(step, this); // Save the macroscopic quantity
            OutputKeEns(step, this);
        }

    }
}


// --------------------------------------------------------------------------------------------------------------------------------------------------------------