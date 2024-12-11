#include "lbm.hpp"
#include "restart_file.hpp"
#include <vector>
#include <fcntl.h>     // For open, O_CREAT, O_WRONLY, O_DIRECT
#include <unistd.h>    // For write, close
#include <iostream>    // For cout
#include <cstring>     // For memset
#include <cstdlib>     // For posix_memalign
#include <sys/stat.h>  // For fstat

#include <fstream>
#include <string>

#ifdef RESTART_O_DIRECT
// Write restart file subroutine
void write_restart(int &step, LBM *lbm)
{
    LBM lb = *lbm;

    // Create name of file
	char	filename[128];
    sprintf(filename,"./restart%06d.dat",step);

    // Define block_size
    const size_t blockSize = 512;  // Example block size, this may vary
    const size_t lbmSize = lb.get_size();
    const size_t dataSize = (lbmSize/blockSize + 1)*blockSize;  // Write 4096 bytes (8 blocks of 512 bytes)

    // Open the file with O_DIRECT flag
    int fd = open(filename, O_WRONLY | O_CREAT | O_TRUNC | O_DIRECT, 0644);
    if (fd == -1) {
        std::cerr << "Failed to open file: " << strerror(errno) << std::endl;
    }

    // Allocate aligned memory for O_DIRECT
    char *buffer;
    if (posix_memalign(reinterpret_cast<void**> (&buffer), blockSize, dataSize) != 0) {
        std::cerr << "Failed to allocate aligned memory" << std::endl;
        close(fd);
    }

    // Fill the buffer with data (for example, set everything to 'A')
    memset(buffer, 0, dataSize);
     
    // Write scalar data first to buffer
    int dx = lb.get_dx();                   memcpy(buffer, &dx, sizeof(int));
    int dy = lb.get_dy();                   memcpy(buffer + sizeof(int), &dy, sizeof(int));
    int dz = lb.get_dz();                   memcpy(buffer + 2 * sizeof(int), &dz, sizeof(int));
    int Nx = lb.get_Nx();                   memcpy(buffer + 3 * sizeof(int), &Nx, sizeof(int));
    int Ny = lb.get_Ny();                   memcpy(buffer + 4 * sizeof(int), &Ny, sizeof(int));
    int Nz = lb.get_Nz();                   memcpy(buffer + 5 * sizeof(int), &Nz, sizeof(int));
    int dt_sim = lb.get_dtsim();            memcpy(buffer + 6 * sizeof(int), &dt_sim, sizeof(int));
                                            memcpy(buffer + 7 * sizeof(int), &step, sizeof(int));
    unsigned long size_scalar_int = 8 * sizeof(int);

    #ifndef MULTICOMP
    double nu = lb.get_nu();                memcpy(buffer + size_scalar_int                     , &nu, sizeof(double));
    double gas_const = lb.get_gasconst();   memcpy(buffer + size_scalar_int + 1 * sizeof(double), &gas_const, sizeof(double));
    double gamma = lb.get_gamma();          memcpy(buffer + size_scalar_int + 2 * sizeof(double), &gamma, sizeof(double));
    double Ra = lb.get_Ra();                memcpy(buffer + size_scalar_int + 3 * sizeof(double), &Ra, sizeof(double));
    double prtl = lb.get_prtl();            memcpy(buffer + size_scalar_int + 4 * sizeof(double), &prtl, sizeof(double));
    unsigned long size_scalar_double = 5 * sizeof(double);
    unsigned long size_scalar = size_scalar_int + size_scalar_double;

    for (int i = 0; i < Nx; ++i){
        for (int j = 0; j < Ny; ++j){
            for (int k = 0; k < Nz; ++k){                
                // Write mixture data
                memcpy(buffer + size_scalar + (i*Ny*Nz+j*Nz+k) * (sizeof(MIXTURE)), &lb.mixture[i][j][k], sizeof(MIXTURE));
            }
        }
    }

    #else
    size_t nSpecies = lb.get_nSpecies();                            memcpy(buffer + size_scalar_int, &nSpecies, sizeof(size_t));
    std::vector<std::string> speciesName = lb.get_speciesName();  

    unsigned long size_scalar_species = sizeof(size_t);
    for (size_t a = 0; a < nSpecies; ++a){
        size_t size_species_a = speciesName[a].size();
        memcpy(buffer + size_scalar_int + size_scalar_species , &size_species_a, sizeof(size_t));
        memcpy(buffer + size_scalar_int + size_scalar_species + sizeof(size_t), speciesName[a].c_str(), size_species_a);
        size_scalar_species += sizeof(size_t) + size_species_a;
    }

    double peremeability = lb.get_permeability();                   memcpy(buffer + size_scalar_int + size_scalar_species, &peremeability, sizeof(double));
    size_scalar_species += sizeof(double);
    
    unsigned long size_scalar = size_scalar_int + size_scalar_species;

    unsigned long size_field = 0;
    for (int i = 0; i < Nx; ++i){
        for (int j = 0; j < Ny; ++j){
            for (int k = 0; k < Nz; ++k){                
                // Write mixture data
                memcpy(buffer + size_scalar + size_field, &lb.mixture[i][j][k], sizeof(MIXTURE));
                size_field += sizeof(MIXTURE);

                // Write species data
                for (size_t a = 0; a < nSpecies; ++a){
                    memcpy(buffer + size_scalar + size_field, &lb.species[a][i][j][k], sizeof(SPECIES));
                    size_field += sizeof(SPECIES);
                }
            }
        }
    }

    #endif

    // Write the data to the file
    ssize_t bytesWritten = write(fd, buffer, dataSize);
    if (bytesWritten == -1) {
        std::cerr << "Write failed: " << strerror(errno) << std::endl;
        free(buffer);
        close(fd);
    }

    std::cout << "Wrote " << bytesWritten << std::endl;

    // Cleanup
    free(buffer); // Free the aligned memory
    close(fd);    // Close the file

}


// Read restart file subroutine
LBM read_restart(const std::string& filename)
{

    // Define block_size
    const size_t blockSize = 512;  // Example block size, this may vary
    const size_t chunkSize = 512 * 1024 * 1024;
    const size_t dataSize = (chunkSize/blockSize + 1)*blockSize;  // Write 4096 bytes (8 blocks of 512 bytes)

    // Open the file with O_DIRECT for direct I/O
    int fd = open(filename.c_str(), O_RDONLY | O_DIRECT);
    if (fd < 0) {
        std::cerr << "Failed to open file: " << strerror(errno) << std::endl;
    }

    // Allocate aligned memory for O_DIRECT
    char *buffer;
    if (posix_memalign(reinterpret_cast<void**> (&buffer), blockSize, dataSize) != 0) {
        std::cerr << "Failed to allocate aligned memory" << std::endl;
        close(fd);
    }

    // Read the scalar data first
    memset(buffer, 0, dataSize);
    ssize_t bytes_read = read(fd, buffer, dataSize);
    if (bytes_read < 0) {
        std::cerr << "Failed to read scalar data from file: " << strerror(errno) << std::endl;
        free(buffer);
        close(fd);
    }

    // Start reading data from the buffer
    size_t offset = 0;
     
    // Get scalar data first from buffer
    int dx;     memcpy(&dx,     buffer + offset, sizeof(int)); offset += sizeof(int);
    int dy;     memcpy(&dy,     buffer + offset, sizeof(int)); offset += sizeof(int);
    int dz;     memcpy(&dz,     buffer + offset, sizeof(int)); offset += sizeof(int);
    int Nx;     memcpy(&Nx,     buffer + offset, sizeof(int)); offset += sizeof(int);
    int Ny;     memcpy(&Ny,     buffer + offset, sizeof(int)); offset += sizeof(int);
    int Nz;     memcpy(&Nz,     buffer + offset, sizeof(int)); offset += sizeof(int);
    int dt_sim; memcpy(&dt_sim, buffer + offset, sizeof(int)); offset += sizeof(int);
    int step;   memcpy(&step,   buffer + offset, sizeof(int)); offset += sizeof(int);
    
    // create LBM object
    #ifndef MULTICOMP
        double nu;          memcpy(&nu, buffer + size_scalar_int, sizeof(double));
        double gas_const;   memcpy(&gas_const  , buffer + size_scalar_int + 1 * sizeof(double), sizeof(double)); 
        double gamma;       memcpy(&gamma      , buffer + size_scalar_int + 2 * sizeof(double), sizeof(double)); 
        double Ra;          memcpy(&Ra         , buffer + size_scalar_int + 3 * sizeof(double), sizeof(double)); 
        double prtl;        memcpy(&prtl       , buffer + size_scalar_int + 4 * sizeof(double), sizeof(double)); 
        unsigned long size_scalar_double = 5 * sizeof(double);

        
        LBM lb(Nx-2, Ny-2, Nz-2, nu);
        lb.set_dx(dx);
        lb.set_dy(dy);
        lb.set_dz(dz);
        lb.set_dtsim(dt_sim);
        lb.set_nstep(nstep);
        lb.set_tout(tout);

        lb.set_gasconst(gas_const);
        lb.set_gamma(gamma);
        lb.set_Ra(Ra);
        lb.set_prtl(prtl);

        unsigned long size_scalar = size_scalar_int + size_scalar_double;

        for (int i = 0; i < Nx; ++i){
            for (int j = 0; j < Ny; ++j){
                for (int k = 0; k < Nz; ++k){
                    // Read mixture data
                    memcpy(&lb.mixture[i][j][k], buffer + size_scalar + (i*Ny*Nz+j*Nz+k) * sizeof(MIXTURE), sizeof(MIXTURE));
                }
            }
        }



    #else

        size_t nSpecies;
        memcpy(&nSpecies, buffer + offset, sizeof(size_t)); offset += sizeof(size_t);

        std::vector<std::string> speciesName(nSpecies);
        for (size_t a = 0; a < nSpecies; ++a) {
            size_t size_species_a;
            memcpy(&size_species_a, buffer + offset, sizeof(size_t)); offset += sizeof(size_t);

            char *speciesBuffer = new char[size_species_a + 1];
            memcpy(speciesBuffer, buffer + offset, size_species_a); offset += size_species_a;
            speciesBuffer[size_species_a] = '\0'; // Null-terminate the string

            speciesName[a] = std::string(speciesBuffer);
            delete[] speciesBuffer;
        }

        double permeability;                    memcpy(&permeability, buffer + offset, sizeof(double)); offset += sizeof(double);



        LBM lb(Nx-2, Ny-2, Nz-2, speciesName);
        lb.set_dx(dx);
        lb.set_dy(dy);
        lb.set_dz(dz);
        lb.set_dtsim(dt_sim);
        lb.set_step(step);
        lb.set_permeability(permeability);

        for (int i = 0; i < Nx; ++i){
            for (int j = 0; j < Ny; ++j){
                for (int k = 0; k < Nz; ++k){
                    // Read mixture data
                    memcpy( &lb.mixture[i][j][k], buffer + offset, sizeof(MIXTURE)); offset += sizeof(MIXTURE);

                    // Read species data
                    for (size_t a = 0; a < nSpecies; ++a){
                        memcpy( &lb.species[a][i][j][k], buffer + offset, sizeof(SPECIES)); offset += sizeof(SPECIES);
                    }

                }
            }
        }

    #endif

    // Clean up and close the file
    free(buffer);
    close(fd);
    std::cout << "Restart file read successfully!" << std::endl;

    return lb;
}

size_t getFileSize(int fd) {
    struct stat stat_buf;
    if (fstat(fd, &stat_buf) == 0) {
        return stat_buf.st_size;  // Return the file size in bytes
    } else {
        std::cerr << "Failed to get file size: " << strerror(errno) << std::endl;
        return -1;  // Error case
    }
}

#else

void write_restart(int &step, LBM *lbm) {
    LBM lb = *lbm;

    // Create name of file
    char filename[128];
    sprintf(filename, "./restart%06d.dat", step);

    // Open file in binary mode
    std::ofstream outFile(filename, std::ios::binary);
    if (!outFile) {
        throw std::runtime_error("Failed to open file: " + std::string(filename));
    }

    // Write scalar data
    int dx = lb.get_dx(), dy = lb.get_dy(), dz = lb.get_dz();
    int Nx = lb.get_Nx(), Ny = lb.get_Ny(), Nz = lb.get_Nz();
    int dt_sim = lb.get_dtsim();

    outFile.write(reinterpret_cast<char*>(&dx), sizeof(dx));
    outFile.write(reinterpret_cast<char*>(&dy), sizeof(dy));
    outFile.write(reinterpret_cast<char*>(&dz), sizeof(dz));
    outFile.write(reinterpret_cast<char*>(&Nx), sizeof(Nx));
    outFile.write(reinterpret_cast<char*>(&Ny), sizeof(Ny));
    outFile.write(reinterpret_cast<char*>(&Nz), sizeof(Nz));
    outFile.write(reinterpret_cast<char*>(&dt_sim), sizeof(dt_sim));
    outFile.write(reinterpret_cast<char*>(&step), sizeof(step));

#ifndef MULTICOMP
    // Write double scalar data
    double nu = lb.get_nu(), gas_const = lb.get_gasconst(), gamma = lb.get_gamma();
    double Ra = lb.get_Ra(), prtl = lb.get_prtl();

    outFile.write(reinterpret_cast<char*>(&nu), sizeof(nu));
    outFile.write(reinterpret_cast<char*>(&gas_const), sizeof(gas_const));
    outFile.write(reinterpret_cast<char*>(&gamma), sizeof(gamma));
    outFile.write(reinterpret_cast<char*>(&Ra), sizeof(Ra));
    outFile.write(reinterpret_cast<char*>(&prtl), sizeof(prtl));

    // Write mixture data
    for (int i = 0; i < Nx; ++i) {
        for (int j = 0; j < Ny; ++j) {
            for (int k = 0; k < Nz; ++k) {
                outFile.write(reinterpret_cast<char*>(&lb.mixture[i][j][k]), sizeof(MIXTURE));
            }
        }
    }
#else
    size_t nSpecies = lb.get_nSpecies();
    outFile.write(reinterpret_cast<char*>(&nSpecies), sizeof(nSpecies));

    // Write species names
    std::vector<std::string> speciesName = lb.get_speciesName();
    for (const auto &name : speciesName) {
        size_t size_species_a = name.size();
        outFile.write(reinterpret_cast<char*>(&size_species_a), sizeof(size_species_a));
        outFile.write(name.c_str(), size_species_a);
    }

    double permeability = lb.get_permeability();
    outFile.write(reinterpret_cast<char*>(&permeability), sizeof(permeability));

    // Write mixture and species data
    for (int i = 0; i < Nx; ++i) {
        for (int j = 0; j < Ny; ++j) {
            for (int k = 0; k < Nz; ++k) {
                outFile.write(reinterpret_cast<char*>(&lb.mixture[i][j][k]), sizeof(MIXTURE));
                for (size_t a = 0; a < nSpecies; ++a) {
                    outFile.write(reinterpret_cast<char*>(&lb.species[a][i][j][k]), sizeof(SPECIES));
                }
            }
        }
    }
#endif

    outFile.close();
    std::cout << "Restart file written successfully!" << std::endl;
}

LBM read_restart(const std::string& filename) {
    std::ifstream inFile(filename, std::ios::binary);
    if (!inFile) {
        throw std::runtime_error("Failed to open file: " + filename);
    }

    // Read scalar data
    int dx, dy, dz, Nx, Ny, Nz, dt_sim, step;
    inFile.read(reinterpret_cast<char*>(&dx), sizeof(dx));
    inFile.read(reinterpret_cast<char*>(&dy), sizeof(dy));
    inFile.read(reinterpret_cast<char*>(&dz), sizeof(dz));
    inFile.read(reinterpret_cast<char*>(&Nx), sizeof(Nx));
    inFile.read(reinterpret_cast<char*>(&Ny), sizeof(Ny));
    inFile.read(reinterpret_cast<char*>(&Nz), sizeof(Nz));
    inFile.read(reinterpret_cast<char*>(&dt_sim), sizeof(dt_sim));
    inFile.read(reinterpret_cast<char*>(&step), sizeof(step));

#ifndef MULTICOMP
    // Read double scalar data
    double nu, gas_const, gamma, Ra, prtl;
    inFile.read(reinterpret_cast<char*>(&nu), sizeof(nu));
    inFile.read(reinterpret_cast<char*>(&gas_const), sizeof(gas_const));
    inFile.read(reinterpret_cast<char*>(&gamma), sizeof(gamma));
    inFile.read(reinterpret_cast<char*>(&Ra), sizeof(Ra));
    inFile.read(reinterpret_cast<char*>(&prtl), sizeof(prtl));

    // Initialize LBM object
    LBM lb(Nx - 2, Ny - 2, Nz - 2, nu);
    lb.set_dx(dx);
    lb.set_dy(dy);
    lb.set_dz(dz);
    lb.set_dtsim(dt_sim);
    lb.set_step(step);
    lb.set_gasconst(gas_const);
    lb.set_gamma(gamma);
    lb.set_Ra(Ra);
    lb.set_prtl(prtl);

    // Read mixture data
    for (int i = 0; i < Nx; ++i) {
        for (int j = 0; j < Ny; ++j) {
            for (int k = 0; k < Nz; ++k) {
                inFile.read(reinterpret_cast<char*>(&lb.mixture[i][j][k]), sizeof(MIXTURE));
            }
        }
    }
#else
    size_t nSpecies;
    inFile.read(reinterpret_cast<char*>(&nSpecies), sizeof(nSpecies));

    std::vector<std::string> speciesName(nSpecies);
    for (size_t a = 0; a < nSpecies; ++a) {
        size_t size_species_a;
        inFile.read(reinterpret_cast<char*>(&size_species_a), sizeof(size_species_a));

        std::string speciesBuffer(size_species_a, '\0');
        inFile.read(&speciesBuffer[0], size_species_a);
        speciesName[a] = speciesBuffer;
    }

    double permeability;
    inFile.read(reinterpret_cast<char*>(&permeability), sizeof(permeability));

    // Initialize LBM object
    LBM lb(Nx - 2, Ny - 2, Nz - 2, speciesName);
    lb.set_dx(dx);
    lb.set_dy(dy);
    lb.set_dz(dz);
    lb.set_dtsim(dt_sim);
    lb.set_step(step);
    lb.set_permeability(permeability);

    // Read mixture and species data
    for (int i = 0; i < Nx; ++i) {
        for (int j = 0; j < Ny; ++j) {
            for (int k = 0; k < Nz; ++k) {
                inFile.read(reinterpret_cast<char*>(&lb.mixture[i][j][k]), sizeof(MIXTURE));
                for (size_t a = 0; a < nSpecies; ++a) {
                    inFile.read(reinterpret_cast<char*>(&lb.species[a][i][j][k]), sizeof(SPECIES));
                }
            }
        }
    }
#endif

    inFile.close();
    std::cout << "Restart file read successfully!" << std::endl;
    return lb;
}

#endif