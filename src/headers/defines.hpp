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
    // #define TAYLOR_GREEN_3D
    // #define CHANNEL_FLOW_3D
    // #define CYLINDER_3D
    // #define VISCOSITY_TEST
    // #define VISCOSITY_TEST_ROTATED
    // #define CONDUCTIVITY_TEST
    // #define SOUNDSPEED_TEST
    // #define SOD_SHOCK_1D
    // #define SOD_SHOCK
    // #define SOD_SHOCK_SIUNIT
    // #define TERNARY_DIFFUSION
    // #define SHEAR_LAYER_MULTICOMP
    // #define PERFECTLY_STIRRED_REACTOR_3D
    // #define CONDUCTION_1D
    // #define COUETTE_FLOW
    // #define COUETTE_FLOW_MULTICOMP
    // #define TAYLOR_GREEN_3D_MULTICOMP
    // #define SHEAR_LAYER
    // #define SHOCK_VORTEX_INTERAC
    // #define RB_INSTABILITY
    // #define RB_INSTABILITY_MULTICOMP
    // #define CONDUCTION_BLACK
    // #define VISCOSITY_TEST_MULTICOMP
    // #define OPPOSED_JET
    // #define OPPOSED_JET_MULTICOMP
    // #define POINT_COMBUSTION_2D
    // #define FLAME_SPEED
    // #define PREMIXED_LAMINAR_FLAME_2D
    // #define CONVECTED_VORTEX
    // #define COFLOW_FLAME
    // #define OPPOSED_JET_FLAME
    // #define PREMIXED_JET_FLAME
    // #define MICROCHANNEL_FLAME
    // #define CIRCULAR_FLAME
    // #define HEATED_OBSTACLE
    // #define HEAT_LID_DRIVEN
    // #define HEAT_LID_DRIVEN_MULTICOMP
    // #define DARIUS_LANDAW
    #define FLAME_VORTEX
    
    // ################ ISOTHERMAL ###################
    // #define ISOTHERM

    // ################# REACTION ####################
    #define REACTION

    // ##### CONJUGATE HEAT TRANSFER FLUID-SOLID #####
    // #define CONJUGATE
    
    // ################ SMOOTHING ####################
    // #define SMOOTHING

    // ####### USE FD for species conservation #######
    // Uncomment for using LBM instead
    // #define FD 

    // ########## Parallel Computation ###############
    // Uncomment 2 lines of code below to utilize parallel computation using OpenMP
    #define PARALLEL
    #define NUM_THREADS 10

    // ############### OUTPUT UNIT ###################
    // #define OUTPUT_SI               // uncomment for SI UNIT, comment for LATTICE UNIT

    // ########### RESTART FILE METHOD ###############
    // #define RESTART_O_DIRECT         // For fast computation but cannot be used for large file

    // ############### Limiter Type ##################
    #define LIMITER_TYPE limiterVanleer
    // #define LIMITER_TYPE limiterMinmod
    // #define LIMITER_TYPE limiterMaxmod
    // #define LIMITER_TYPE limiterSuperbee
    // #define LIMITER_TYPE limiterMC


    // _____________________________________________________________________________________________________
    // Definition For code, don not change!

    #define SPECIES_MIN 1E-8

    // ############# Flow Parameters ##################
    // Reference Moment Value [in Lattice Unit]
    #define VEL0 0.1    // velocity    
    #define RHO0 1.0    // Density
    #define TEMP0 0.025   // Temperature

    // ################# Velocity Set ##################
    #if defined SOD_SHOCK_1D
        #define D1Q3
        #define NDIM 1    
    #elif defined AA
        #define D2Q9
        #define NDIM 2
    #elif defined TAYLOR_GREEN_2D || defined TAYLOR_GREEN_3D || defined CHANNEL_FLOW_3D || defined CYLINDER_3D || defined VISCOSITY_TEST || defined CONDUCTIVITY_TEST || defined SOUNDSPEED_TEST || defined SOD_SHOCK || defined TERNARY_DIFFUSION || defined SOD_SHOCK_SIUNIT || defined SHEAR_LAYER_MULTICOMP || defined PERFECTLY_STIRRED_REACTOR_3D || defined CYLINDER_2D || defined CONDUCTION_1D || defined BPVT_1 || defined COUETTE_FLOW || defined COUETTE_FLOW_MULTICOMP || defined TAYLOR_GREEN_3D_MULTICOMP || defined SHEAR_LAYER || defined SHOCK_VORTEX_INTERAC || defined RB_INSTABILITY || defined CONDUCTION_BLACK || defined RB_INSTABILITY_MULTICOMP || defined VISCOSITY_TEST_MULTICOMP || defined VISCOSITY_TEST_ROTATED || defined OPPOSED_JET || defined POINT_COMBUSTION_2D || defined FLAME_SPEED || defined PREMIXED_LAMINAR_FLAME_2D || defined OPPOSED_JET_MULTICOMP || defined CONVECTED_VORTEX || defined COFLOW_FLAME || defined OPPOSED_JET_FLAME || defined PREMIXED_JET_FLAME || defined MICROCHANNEL_FLAME || defined CIRCULAR_FLAME || defined HEATED_OBSTACLE || defined HEAT_LID_DRIVEN || defined HEAT_LID_DRIVEN_MULTICOMP || defined DARIUS_LANDAW || defined FLAME_VORTEX
        #define D3Q27
        #define NDIM 3
    #endif

    #if defined TERNARY_DIFFUSION || defined SOD_SHOCK_SIUNIT || defined SHEAR_LAYER_MULTICOMP || defined PERFECTLY_STIRRED_REACTOR_3D || defined CONDUCTION_1D || defined COUETTE_FLOW_MULTICOMP || defined TAYLOR_GREEN_3D_MULTICOMP || defined RB_INSTABILITY_MULTICOMP || defined VISCOSITY_TEST_MULTICOMP || defined POINT_COMBUSTION_2D || defined FLAME_SPEED || defined PREMIXED_LAMINAR_FLAME_2D || defined OPPOSED_JET_MULTICOMP || defined CONVECTED_VORTEX || defined COFLOW_FLAME || defined OPPOSED_JET_FLAME || defined PREMIXED_JET_FLAME || defined MICROCHANNEL_FLAME || defined CIRCULAR_FLAME || defined HEAT_LID_DRIVEN_MULTICOMP || defined DARIUS_LANDAW || defined FLAME_VORTEX
        #define MULTICOMP
    #endif

    // ############## Boundary Condition ###############
    #define TYPE_F 0    // fluid domain
    #define TYPE_S 1    // dirichlet boundary condition
    #define TYPE_P 2    // periodic boundary
    #define TYPE_A 3    // Adiabatic No-Slip Wall
    #define TYPE_FS 4   // Adiabatic Free-Slip Wall
    #define TYPE_I 5    // Inflow boundary condition
    #define TYPE_O 6    // Outflow boundary condition (Neumann Boundary Condition | Zero gradient for all variables)
    #define TYPE_I_C 7  // Inflow boundary condition
    #define TYPE_O_C 8  // Outflow boundary condition (Neumann Boundary Condition | Zero gradient for all variables)
    #define TYPE_Q 9    // Heat boundary condition



    

#endif