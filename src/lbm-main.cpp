
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

    // allocate memory for lattice (one contiguous block + pointer tables for [i][j][k] indexing)
    mixture_data = new MIXTURE [(size_t)this->Nx * this->Ny * this->Nz];
    mixture = new MIXTURE **[this->Nx];
    for (int i = 0; i < this->Nx; ++i)
    {
        mixture[i] = new MIXTURE *[this->Ny];
        for (int j = 0; j < this->Ny; ++j)
        {
            mixture[i][j] = mixture_data + ((size_t)i * this->Ny + j) * this->Nz;
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

    // allocate memory for mixture (one contiguous block + pointer tables for [i][j][k] indexing)
    const size_t ncell = (size_t)this->Nx * this->Ny * this->Nz;
    mixture_data = new MIXTURE [ncell];
    mixture = new MIXTURE **[this->Nx];
    for (int i = 0; i < this->Nx; ++i)
    {
        mixture[i] = new MIXTURE *[this->Ny];
        for (int j = 0; j < this->Ny; ++j)
        {
            mixture[i][j] = mixture_data + ((size_t)i * this->Ny + j) * this->Nz;
        }
    }

    // allocate memory for species (one contiguous block per species)
    this->species.resize(nSpecies);
    this->species_data.resize(nSpecies);
    for(size_t p = 0; p < nSpecies; ++p)
    {
        this->species_data[p] = new SPECIES [ncell];
        this->species[p] = new SPECIES **[this->Nx];
        for (int i = 0; i < this->Nx; ++i)
        {
            this->species[p][i] = new SPECIES *[this->Ny];
            for (int j = 0; j < this->Ny; ++j)
            {
                this->species[p][i][j] = this->species_data[p] + ((size_t)i * this->Ny + j) * this->Nz;
            }
        }
    }

    // Create Cantera's Solution object
    int nThreads = omp_get_max_threads();
    std::cout << "nThreads : " << nThreads << std::endl;
    for(int i = 0; i < nThreads; ++i)
    {
        // auto sol = Cantera::newSolution("gri30.yaml", "gri30","mixture-averaged");
        auto sol = Cantera::newSolution("h2o2.yaml", "ohmech");
        // auto sol = Cantera::newSolution("gri30.yaml", "gri30", "multicomponent");
        // auto sol = Cantera::newSolution("./src/reaction-mech/one-step.yaml", "FakeGas");
        // auto sol = Cantera::newSolution("./src/reaction-mech/CH4_2S.yaml", "CH4_BFER_multi");
        // auto sol = Cantera::newSolution("./src/reaction-mech/propane_mech.yaml");
        sols.push_back(sol);
    }

    // Precompute the Cantera species indices and molecular weights once,
    // so the hot loops don't have to do string lookups for every lattice cell
    auto gas = sols[0]->thermo();
    nSpeciesCantera = gas->nSpecies();
    speciesIdx.resize(nSpecies);
    molarMass.resize(nSpecies);
    for(size_t a = 0; a < nSpecies; ++a)
    {
        speciesIdx[a] = gas->speciesIndex(speciesName[a]);
        // speciesIndex() returns Cantera::npos when the name is not in the mechanism.
        // This happens when the reaction mechanism loaded by the LBM constructor (above)
        // does not match the one used by main_setup(). Fail loudly with an actionable
        // message instead of letting molecularWeight(npos) abort deep inside Cantera.
        if (speciesIdx[a] == Cantera::npos)
        {
            std::cerr << "\nERROR (LBM constructor): species \"" << speciesName[a]
                      << "\" is not present in the reaction mechanism loaded here.\n"
                      << "The mechanism in src/lbm-main.cpp (currently loading "
                      << gas->nSpecies() << " species) must match the one used in your "
                      << "main_setup() case in src/setup.cpp.\n"
                      << "Update the Cantera::newSolution(...) call in the LBM constructor "
                      << "to the same mechanism file/phase.\n" << std::endl;
            std::abort();
        }
        molarMass[a] = gas->molecularWeight(speciesIdx[a]);
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
        // std::cout << "  Species Collision Done" << std::endl;
        #endif
        // std::cout << "Collision species     : " << double(omp_get_wtime()-start) << " seconds" << std::endl;

        // start = omp_get_wtime();
        Collide();   // collision step
        // std::cout << "  Mixture Collision Done" << std::endl;
        // std::cout << "Collision mixture     : " << double(omp_get_wtime()-start) << " seconds" << std::endl;

        // start = omp_get_wtime();
        Streaming();        // streaming step & BC
        // std::cout << "  Streaming Done" << std::endl;
        // std::cout << "Streaming process     : " << double(omp_get_wtime()-start) << " seconds" << std::endl;

        // start = omp_get_wtime();
        TMS_BC();
        // std::cout << "  Apply BC Done" << std::endl;  
        // std::cout << "TMS BC                : " << double(omp_get_wtime()-start) << " seconds" << std::endl;
                
        // start = omp_get_wtime();
        #ifdef SMOOTHING
            calculate_moment_smoothing(); // calculate moment
        #else
            calculate_moment(); // calculate moment
        #endif
        // std::cout << "  Calculate Moment Done" << std::endl;
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
