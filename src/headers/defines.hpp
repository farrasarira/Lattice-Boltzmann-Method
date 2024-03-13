#ifndef DEFINES_H
    #define DEFINES_H

    #include<math.h>

    // !! ALWAYS DO "MAKE CLEAN" AND THEN "MAKE" AFTER MODIFY THIS CODE !!

    //_______________________________________________________________________________________________________
    // Choose the simulation condition

    // ################ Flow Cases ####################
    // Choose one of the following Flow Cases
    // #define CYLINDER_2D
    // #define TAYLOR_GREEN_2D
    // #define TAsYLOR_GREEN_3D
    // #define CHANNEL_FLOW_3D
    // #define CYLINDER_3D
    // #define VISCOSITY_TEST
    // #define SOD_SHOCK_1D
    // #define SOD_SHOCK
    // #define SOD_SHOCK_SIUNIT
    #define TERNARY_DIFFUSION
    // #define SHEAR_LAYER_MULTICOMP
    // #define PERFECTLY_STIRRED_REACTOR_3D
    // #define CONDUCTION_1D
    // #define COUETTE_FLOW
    // #define COUETTE_FLOW_MULTICOMP
    // #define TAYLOR_GREEN_3D_MULTICOMP

    
    // ####### USE FD for species conservation #######
    // Uncomment for using LBM instead
    // #define FD 

    // ########## Parallel Computation ###############
    // Uncomment 2 lines of code below to utilize parallel computation using OpenMP
    #define PARALLEL
    #define NUM_THREADS 7

    // ############### OUTPUT UNIT ###################
    #define OUTPUT_SI               // uncomment for SI UNIT, comment for LATTICE UNIT

    // ############### Limiter Type ##################
    #define LIMITER_TYPE limiterVanleer
    // #define LIMITER_TYPE limiterMinmod
    // #define LIMITER_TYPE limiterMaxmod
    // #define LIMITER_TYPE limiterSuperbee
    // #define LIMITER_TYPE limiterMC


    // _____________________________________________________________________________________________________
    // Definition For code, don not change!

    // ############# Flow Parameters ##################
    // Reference Moment Value [in Lattice Unit]
    #define VEL0 0.1    // velocity    
    #define RHO0 1.0    // Density
    #define TEMP0 0.025   // Temperature

    // ################# Velocity Set ##################
    #if defined SOD_SHOCK_1D
        #define D1Q3
        #define NDIM 1    
    #elif defined TAYLOR_GREEN_2D 
        #define D2Q9
        #define NDIM 2
    #elif defined TAYLOR_GREEN_3D || defined CHANNEL_FLOW_3D || defined CYLINDER_3D || defined VISCOSITY_TEST || defined SOD_SHOCK || defined TERNARY_DIFFUSION || defined SOD_SHOCK_SIUNIT || defined SHEAR_LAYER_MULTICOMP || defined PERFECTLY_STIRRED_REACTOR_3D || defined CYLINDER_2D || defined CONDUCTION_1D || defined BPVT_1 || defined COUETTE_FLOW || defined COUETTE_FLOW_MULTICOMP || defined TAYLOR_GREEN_3D_MULTICOMP
        #define D3Q27
        #define NDIM 3
    #endif

    #if defined TERNARY_DIFFUSION || defined SOD_SHOCK_SIUNIT || defined SHEAR_LAYER_MULTICOMP || defined PERFECTLY_STIRRED_REACTOR_3D || defined CONDUCTION_1D || defined COUETTE_FLOW_MULTICOMP || defined TAYLOR_GREEN_3D_MULTICOMP
        #define MULTICOMP
    #endif

    // ############## Boundary Condition ###############
    #define TYPE_F 0 // fluid domain
    #define TYPE_S 1 // (stationary or moving) solid boundary
    #define TYPE_E 2 // equilibrium boundary (inflow/outflow)
    #define TYPE_P 3 // periodic boundary
    #define TYPE_A 4 // Adiabatic Wall

    

#endif