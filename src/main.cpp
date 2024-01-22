#include<iostream>

#include "./headers/lbm.hpp"
#include "./headers/setup.hpp"
#include "./headers/output.hpp"
#include <iostream>
#include <omp.h>

int main(int argc, char** argv)
{
    #ifdef PARALLEL
        omp_set_num_threads(NUM_THREADS);
    #endif

    printLogo();
    double start;
    start = omp_get_wtime();

    // LBM
    main_setup();  

    std::cout << "Comp Time  : " << double(omp_get_wtime()-start) << " seconds" << std::endl;
    std::cout << '\a' ;

    return 0;
}