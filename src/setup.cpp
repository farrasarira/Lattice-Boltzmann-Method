#include "lbm.hpp"
#include "geometry.hpp"
#include "units.hpp"
#include "cantera.hpp"
#include <iostream>
#include <math.h>
#include <omp.h>
#include "restart_file.hpp"

Units_Conv units;

#if defined CYLINDER_2D
void main_setup() // 2D Flow over cylinder --------------------------------------------------------
{
    double RE = 200;
    LBM lb(1000, 500, 1, 0.01);
    int Nx = lb.get_Nx(); int Ny = lb.get_Ny(); int Nz = lb.get_Nz();
    double nu = lb.get_nu();

    int D = Ny/10; 

    double u_max = RE*nu/D;    // [Velocity in Lattice Unit]
    double rho0 = 1.0;      // [Density in Lattice Unit]
    double T0 = 0.1;       // [Temperature in Lattice Unit]

    // double si_u_max = 0.05; // [m/s]
    // double si_rho = 1.225; // [kg/m^3]
    // double si_temp = 273.15; // [K]
    // double si_cylinder_temp = 500.0; // [K]

    std::cout << "umax : " << u_max << std::endl;

    // auto sol = Cantera::newSolution("gri30.yaml","gri30");
    // auto gas = sol->thermo();
    // auto trans = sol->transport();
    // gas->setState_TR(si_temp, si_rho);
    
    // double si_D = RE*trans->viscosity()/si_rho/si_u_max; // [m]

    // units.set_m_kg_s(D, u_max, rho0, T0, si_D, si_u_max, si_rho, si_temp); // setting the conversion factor 

    // std::cout << gas->report() << std::endl;
    // std::cout << units.nu(trans->viscosity()/si_rho) << std::endl;

    #pragma omp parallel for schedule(dynamic)
    for(int i = 0; i < Nx ; ++i)
    {
        for(int j = 0; j < Ny; ++j)
        {
            for(int k = 0; k < Nz; ++k)
            {
                if (k == 0 || k == Nz-1)
                {
                    lb.mixture[i][j][k].type = TYPE_P;
                }
                if (j == 0 || j == Ny-1) 
                {
                    lb.mixture[i][j][k].type = TYPE_P;  // Solid Boundary in upper and lower side
                }
                if (i == 0 )
                {
                    lb.mixture[i][j][k].type = TYPE_I;  // Equilibrium inlet and outlet
                    lb.mixture[i][j][k].u = u_max;
                    lb.mixture[i][j][k].p = T0;
                    lb.mixture[i][j][k].temp = T0;
                }
                if (i == Nx-1)
                {
                    lb.mixture[i][j][k].type = TYPE_A;  // Equilibrium inlet and outlet
                    lb.mixture[i][j][k].u = u_max;
                    lb.mixture[i][j][k].p = T0;
                    lb.mixture[i][j][k].temp = T0;
                }
                if (lb.mixture[i][j][k].type == TYPE_F)  // Fluid Domain
                {
                    lb.mixture[i][j][k].rho = rho0;
                    lb.mixture[i][j][k].u = 0.0;
                    lb.mixture[i][j][k].p = T0;
                    lb.mixture[i][j][k].temp = T0;
                }
                if(lb.mixture[i][j][k].type == TYPE_A)
                {
                    lb.mixture[i][j][k].temp = T0;
                }
            }
        }
    }

    cylinder_generator(lb, D);

    lb.run(10000, 100);
}

#elif defined TAYLOR_GREEN_3D
void main_setup() // 3D Taylor-Green Vortex
{
    double RE = 200;
    LBM lb(100, 100, 100, 0.008);

    const double gasconst = 1.0;
    const double gamma = 1.4;
    const double rho0 = 1.0;
    const double p0 = 1.0;
    const double temp0 = 1.0/3.0;//p0/(gasconst * rho0); 
    lb.set_gasconst(gasconst);
    lb.set_gamma(gamma);
    lb.set_prtl(0.71);
    
    
    int Nx = lb.get_Nx(); int Ny = lb.get_Ny(); int Nz = lb.get_Nz();
    int NX = lb.get_NX(); int NY = lb.get_NY(); int NZ = lb.get_NZ();
    double NU = lb.get_nu();

    double periodicity = 1.0;
    const double a = (double)NX/periodicity;
    const double b = (double)NY/periodicity;
    const double c = (double)NZ/periodicity;

    double u_max = RE * NU / (NX/(2.*M_PI));
    std::cout << "umax : " << u_max << std::endl;
    std::cout << "sound speed : " << lb.get_soundspeed(temp0) << std::endl; 
    std::cout << "Mach : " << u_max / lb.get_soundspeed(temp0) << std::endl; 

    #pragma omp parallel for schedule(dynamic)
    for(int i = 0; i < Nx ; ++i)
    {
        for(int j = 0; j < Ny; ++j)
        {
            for(int k = 0; k < Nz; ++k)
            {
                if (i==0 || i==Nx-1 ||  j==0 || j==Ny-1 || k==0 || k==Nz-1) // set periodic boundary condition
                {
                    lb.mixture[i][j][k].type = TYPE_P;
                }
                if (lb.mixture[i][j][k].type == TYPE_F)
                {
                    lb.mixture[i][j][k].u =  u_max * sin(2.0*M_PI/a*(double)i) * cos(2.0*M_PI/b*(double)j) * cos(2.0*M_PI/c*(double)k);
                    lb.mixture[i][j][k].v = -u_max * cos(2.0*M_PI/a*(double)i) * sin(2.0*M_PI/b*(double)j) * cos(2.0*M_PI/c*(double)k);
                    lb.mixture[i][j][k].w = 0.0;
                    lb.mixture[i][j][k].p = p0 + rho0*u_max*u_max/16.0 * (cos(4.0*M_PI/a*(double)i) + cos(4.0*M_PI/b*(double)j)) * (cos(4.0*M_PI/c*(double)k) + 2.0);
                    lb.mixture[i][j][k].temp = temp0;
                }
            }
        }
    }
    
    lb.run(10000000,10);
}
#elif defined TAYLOR_GREEN_2D
void main_setup() // 2D Taylor-Green Vortex
{
    double RE = 100;

    LBM lb(30, 30, 1, 0.01);
    int Nx = lb.get_Nx(); int Ny = lb.get_Ny(); int Nz = lb.get_Nz();
    double NX = lb.get_NX();
    double NU = lb.get_nu();
    lb.set_Ra(2E+5);

    const double u_max = RE * NU / NX;    // maximum initial velocity
    std::cout << "umax      : " << u_max << std::endl;
    const unsigned int periodicity = 1;
    const double a = (double)NX/periodicity;
    const double b = (double)NX/periodicity;
    
    #pragma omp parallel for schedule(dynamic) 
    for(int i = 0; i < Nx ; ++i)
    {
        for(int j = 0; j < Ny; ++j)
        {
            for(int k = 0; k < Nz; ++k)
            {
                if (i==0 || i==Nx-1 || k==0 || k==Nz-1) // set periodic boundary condition
                {
                    lb.mixture[i][j][k].type = TYPE_P;
                }
                if ( j==0 )
                {
                    lb.mixture[i][j][k].type = TYPE_S;
                    lb.mixture[i][j][k].temp = 0.025;
                }
                if ( j==Ny-1 )
                {
                    lb.mixture[i][j][k].type = TYPE_S;
                    lb.mixture[i][j][k].temp = 0.05;
                }

                if (lb.mixture[i][j][k].type == TYPE_F)
                {
                    double X = i ;
                    double Y = j ;
                    lb.mixture[i][j][k].u = -u_max * cos(2.0*M_PI/a*X) * sin(2.0*M_PI/b*Y) ;
                    lb.mixture[i][j][k].v =  u_max * sin(2.0*M_PI/a*X) * cos(2.0*M_PI/b*Y) ;
                    lb.mixture[i][j][k].rho = 1.0 - u_max*u_max * 3.0/4.0 * (cos(4.0*M_PI/a*X) + cos(4.0*M_PI/b*Y));
                    lb.mixture[i][j][k].temp = 0.03;
                    lb.mixture[i][j][k].p = lb.mixture[i][j][k].temp * (1.0 - u_max*u_max * 3.0/4.0 * (cos(4.0*M_PI/a*X) + cos(4.0*M_PI/b*Y)));
                }
            }
        }
    }
    
    lb.run(1000,10);
}

#elif defined CHANNEL_FLOW_3D
void main_setup() // 3D Channel Flow
{
    LBM lb(NX,NY,NY,NU);
    int Nx = lb.getNx(); int Ny = lb.getNy(); int Nz = lb.getNz();

    const double u_max = RE * NU / NY;    // maximum initial velocity
    std::cout << "umax      : " << u_max << std::endl;
    
    #pragma omp parallel for schedule(dynamic)
    for(int i = 0; i < Nx ; ++i)
    {
        for(int j = 0; j < Ny; ++j)
        {
            for(int k = 0; k < Nz; ++k)
            {
                if (j == 0 || j == Ny-1) 
                {
                    lb.mixture[i][j][k].type = TYPE_S;  // Solid Boundary in upper and lower side
                }
                if (i == 0 || i == Nx-1)
                {
                    lb.mixture[i][j][k].type = TYPE_E;  // Equilibrium inlet and outlet
                }
                if (k == 0 || k == Nz-1)
                {
                    lb.mixture[i][j][k].type = TYPE_P;
                }
                if (lb.mixture[i][j][k].type == TYPE_F || lb.mixture[i][j][k].type == TYPE_E)  // Fluid Domain
                {
                    lb.mixture[i][j][k].rho = 1.0;
                    lb.mixture[i][j][k].u = u_max;
                    lb.mixture[i][j][k].v = 0;
                    lb.mixture[i][j][k].w = 0;
                }
            }
        }
    }
    return lb;
}

#elif defined CYLINDER_3D
void main_setup() // 3D Flow over cylinder --------------------------------------------------------
{
    LBM lb(NX,NY,NZ,NU);
    int Nx = lb.getNx(); int Ny = lb.getNy(); int Nz = lb.getNz();

    int D = NY/10; 
    double nu = lb.getNu();
    double u_max = RE * nu / D;

    cylinder_generator(lb, D);

    #pragma omp parallel for schedule(dynamic)
    for(int i = 0; i < Nx ; ++i)
    {
        for(int j = 0; j < Ny; ++j)
        {
            for(int k = 0; k < Nz; ++k)
            {
                if (j == 0 || j == Ny-1 || k == 0 || k == Nz-1 ) 
                {
                    lb.mixture[i][j][k].type = TYPE_S;  // Solid Boundary in upper and lower side
                }
                if (i == 0 || i == Nx-1)
                {
                    lb.mixture[i][j][k].type = TYPE_E;  // Equilibrium inlet and outlet
                }
                if (lb.mixture[i][j][k].type == TYPE_F || lb.mixture[i][j][k].type == TYPE_E)  // Fluid Domain
                {
                    lb.mixture[i][j][k].rho = 1.0;
                    lb.mixture[i][j][k].u = u_max;
                    lb.mixture[i][j][k].v = 0;
                    lb.mixture[i][j][k].w = 0;
                }
            }
        }
    }
    std::cout << "setup selesai" << std::endl;
    return lb;
}

#elif defined VISCOSITY_TEST
void main_setup() // 2D Viscos Test --------------------------------------------------------
{
    double NX = 1;
    double NY = 200;
    double NZ = 1;
    double nu = 0.001;//0.1;

    LBM lb(NX,NY,NZ, nu);
    int Nx = lb.get_Nx(); int Ny = lb.get_Ny(); int Nz = lb.get_Nz();

    const double mach = 0.2;
    const double gamma = 1.162;//1.4;
    const double gas_const = 0.0283563;
    const double prtl = 0.7;
    const double T0 = 0.35;//0.1;
    const double a0 = 0.001*sqrt(gamma*gas_const*T0); // amplitude
    const double u0 = mach*sqrt(gamma*gas_const*T0); 

    std::cout << "mach  : " << mach << std::endl;
    std::cout << "u0    : " <<  u0 << std::endl;

    lb.set_gamma(gamma);
    lb.set_gasconst(gas_const);
    lb.set_prtl(prtl);


    #pragma omp parallel for schedule(dynamic)
    for(int i = 0; i < Nx ; ++i)
    {
        for(int j = 0; j < Ny; ++j)
        {
            for(int k = 0; k < Nz; ++k)
            {
                if (i==0 || i==Nx-1 ||  j==0 || j==Ny-1 || k==0 || k==Nz-1) // set periodic boundary condition
                {
                    lb.mixture[i][j][k].type = TYPE_P;
                }
                if (lb.mixture[i][j][k].type == TYPE_F)
                {
                    //double X = i ;
                    double Y = j ;
                    lb.mixture[i][j][k].u = a0*sin(Y/NY*2*M_PI);
                    lb.mixture[i][j][k].v = u0;// +
                    lb.mixture[i][j][k].w = 0.0;
                    lb.mixture[i][j][k].temp = T0;
                    lb.mixture[i][j][k].p = 1.0*gas_const*T0;
                    
                }
            }
        }
    }
    lb.run(10000,1000);
}

#elif defined VISCOSITY_TEST_ROTATED
void main_setup() // 2D Viscos Test --------------------------------------------------------
{
    double NX = 100;
    double NY = 100;
    double NZ = 1;
    double nu = 0.1;//0.1;

    LBM lb(NX,NY,NZ, nu);
    int Nx = lb.get_Nx(); int Ny = lb.get_Ny(); int Nz = lb.get_Nz();

    const double mach = 0.0;
    const double gamma = 1.4;//1.4;
    const double gas_const = 1.0;
    const double a0 = 0.001; // amplitude
    const double T0 = 0.05;//0.1;
    const double u0 = mach*sqrt(gamma*gas_const*T0); 

    std::cout << "mach  : " << mach << std::endl;
    std::cout << "u0    : " <<  u0 << std::endl;

    lb.set_gamma(gamma);
    lb.set_gasconst(gas_const);


    #pragma omp parallel for schedule(dynamic)
    for(int i = 0; i < Nx ; ++i)
    {
        for(int j = 0; j < Ny; ++j)
        {
            for(int k = 0; k < Nz; ++k)
            {
                if (i==0 || i==Nx-1 ||  j==0 || j==Ny-1 || k==0 || k==Nz-1) // set periodic boundary condition
                {
                    lb.mixture[i][j][k].type = TYPE_P;
                }
                if (lb.mixture[i][j][k].type == TYPE_F)
                {
                    //double X = i ;
                    double theta = 2.0*M_PI/360 * 45;
                    lb.mixture[i][j][k].u = (sin((400179318894815*j*cos(theta))/9007199254740992 - (400179318894815*i*sin(theta))/9007199254740992)*cos(theta))/1000;
                    lb.mixture[i][j][k].v = u0 + (sin((400179318894815*j*cos(theta))/9007199254740992 - (400179318894815*i*sin(theta))/9007199254740992)*sin(theta))/1000;// +
                    lb.mixture[i][j][k].w = 0.0;
                    lb.mixture[i][j][k].temp = T0;
                    lb.mixture[i][j][k].p = 1.0*gas_const*T0;
                    
                }
            }
        }
    }
    lb.run(10000,1000);
}


#elif defined CONDUCTIVITY_TEST
void main_setup() // 2D Viscos Test --------------------------------------------------------
{
    double NX = 1;
    double NY = 200;
    double NZ = 1;
    double nu = 0.1;

    LBM lb(NX,NY,NZ, nu);
    int Nx = lb.get_Nx(); int Ny = lb.get_Ny(); int Nz = lb.get_Nz();

    const double mach = 0.9;
    const double gamma = 1.4;
    const double a0 = 0.001; // amplitude
    const double T0 = 0.1;
    const double u0 = mach*sqrt(gamma*T0); 
    const double rho0 = 1.0;

    std::cout << "mach  : " << mach << std::endl;
    std::cout << "u0    : " << u0 << std::endl;

    lb.set_gamma(gamma);

    // double si_u_max = 1.0; // [m/s]
    // double si_rho = 1.225; // [kg/m^3]
    // double si_temp = 288.15; // [K]

    // auto gas = sols[0]->thermo();
    // auto trans = sols[0]->transport();
    // gas->setState_TR(si_temp, si_rho);

    // units.set_m_kg_s(NY, a0, rho0, T0, 1.0, si_u_max, si_rho, si_temp); // setting the conversion factor 

    // double si_pressure = gas->pressure();
    
    // double si_gas_constant = si_pressure/(si_rho*si_temp);
    // double gas_const = units.R(si_gas_constant);
    
    #pragma omp parallel for schedule(dynamic)
    for(int i = 0; i < Nx ; ++i)
    {
        for(int j = 0; j < Ny; ++j)
        {
            for(int k = 0; k < Nz; ++k)
            {
                if (i==0 || i==Nx-1 ||  j==0 || j==Ny-1 || k==0 || k==Nz-1) // set periodic boundary condition
                {
                    lb.mixture[i][j][k].type = TYPE_P;
                }
                if (lb.mixture[i][j][k].type == TYPE_F)
                {
                    //double X = i ;
                    double Y = j ;
                    lb.mixture[i][j][k].u = u0;
                    lb.mixture[i][j][k].v = 0.0;
                    lb.mixture[i][j][k].w = 0.0;
                    lb.mixture[i][j][k].rho = rho0 + a0*sin(Y/NY*2*M_PI);
                    lb.mixture[i][j][k].temp = rho0*0.1/lb.mixture[i][j][k].rho;
                    lb.mixture[i][j][k].p = lb.mixture[i][j][k].rho*lb.mixture[i][j][k].temp;
                    
                }
            }
        }
    }
    lb.run(10000,1000);
}

#elif defined SOUNDSPEED_TEST
void main_setup() // 2D Viscos Test --------------------------------------------------------
{
    double NX = 10000;
    double NY = 1;
    double NZ = 1;
    double nu = 0.0001;

    LBM lb(NX,NY,NZ, nu);
    int Nx = lb.get_Nx(); int Ny = lb.get_Ny(); int Nz = lb.get_Nz();

    const double gamma = 1.50601;
    const double gascont = 0.046897;
    
    const double temp = 0.025;
    const double dp = 1e-4;

    lb.set_gamma(gamma);
    lb.set_gasconst(gascont);

    std::cout << "Speed of Sound : " << lb.get_soundspeed(temp) << std::endl;

    // double si_u_max = 1.0; // [m/s]
    // double si_rho = 1.225; // [kg/m^3]
    // double si_temp = 288.15; // [K]

    // auto gas = sols[0]->thermo();
    // auto trans = sols[0]->transport();
    // gas->setState_TR(si_temp, si_rho);

    // units.set_m_kg_s(NY, a0, rho0, T0, 1.0, si_u_max, si_rho, si_temp); // setting the conversion factor 

    // double si_pressure = gas->pressure();
    
    // double si_gas_constant = si_pressure/(si_rho*si_temp);
    // double gas_const = units.R(si_gas_constant);
    
    #pragma omp parallel for schedule(dynamic)
    for(int i = 0; i < Nx ; ++i)
    {
        for(int j = 0; j < Ny; ++j)
        {
            for(int k = 0; k < Nz; ++k)
            {
                if (j==0 || j==Ny-1 || k==0 || k==Nz-1) // set periodic boundary condition
                {
                    lb.mixture[i][j][k].type = TYPE_P;
                }
                if ( i==0 || i==Nx-1 )
                {
                    lb.mixture[i][j][k].type = TYPE_A;
                }
                if (lb.mixture[i][j][k].type == TYPE_F)
                {
                    lb.mixture[i][j][k].u = 0.0;
                    lb.mixture[i][j][k].v = 0.0;
                    lb.mixture[i][j][k].w = 0.0;
                    lb.mixture[i][j][k].temp = temp;

                    lb.mixture[i][j][k].p = smooth(temp + dp, temp, i, 0.5*Nx, 0.07);   
                    
                }
            }
        }
    }
    lb.run(10000,1000);
}

#elif defined SOD_SHOCK
void main_setup() // Sos shock tube test case --------------------------------------------------------
{
    int NX = 3000; 
    int NY = 1; 
    int NZ = 1;
    
    double NU = 0.00227006;

    const double gamma = 1.38558;
    const double gasconst = 2.19691;
    const double T0 = 0.025;
    const double prantdl = 0.762829;

    LBM lb(NX, NY, NZ, NU);
    int Nx = lb.get_Nx(); int Ny = lb.get_Ny(); int Nz = lb.get_Nz();
    lb.set_gamma(gamma);
    lb.set_gasconst(gasconst);
    lb.set_prtl(prantdl);

    
    #pragma omp parallel for schedule(dynamic)
    for(int i = 0; i < Nx ; ++i)
    {
        for(int j = 0; j < Ny; ++j)
        {
            for(int k = 0; k < Nz; ++k)
            {
                if ( j==0 || j==Ny-1 || k==0 || k==Nz-1) // set periodic boundary condition
                {
                    lb.mixture[i][j][k].type = TYPE_P;
                }
                
                if (i==0 || i==Nx-1 )
                {
                    lb.mixture[i][j][k].type = TYPE_A;
                }

                if (lb.mixture[i][j][k].type == TYPE_F)
                {   
                    lb.mixture[i][j][k].p = smooth(0.00395, 0.00375, i, 0.5*Nx, 0.3);   
                    // std::cout << " p : " << units.si_p(2.0) << std::endl;
                    lb.mixture[i][j][k].temp = smooth(T0, T0, i, 0.5*Nx, 0.3);   
                    if ((float)i/(float)Nx <= 0.5 )
                    {
                        lb.mixture[i][j][k].u = 0.0;
                        lb.mixture[i][j][k].v = 0.0;
                        lb.mixture[i][j][k].w = 0.0;
                        // lb.mixture[i][j][k].p = 0.01;
                        // lb.mixture[i][j][k].temp = 0.2;
                    }
                    else
                    {
                        lb.mixture[i][j][k].u = 0.0;
                        lb.mixture[i][j][k].v = 0.0;
                        lb.mixture[i][j][k].w = 0.0;
                        // lb.mixture[i][j][k].p = 0.001;
                        // lb.mixture[i][j][k].temp = 0.025;
                    }                       
                }
            }
        }
    }

    lb.run(5000,100);
}

#elif defined SOD_SHOCK_1D
void main_setup() // 1D Sod shock tube --------------------------------------------------------
{
    int NX = 3000; 
    int NY = 1; 
    int NZ = 1;
    double NU = 0.025;
    
    double si_len = 0.1;    // [m]
    double si_u_max = 1.0;  // [m/s]
    double si_rho = 1.225;  // [kg/m^3]
    double si_temp = 288.15;// [K]

    units.set_m_kg_s(NX, VEL0, RHO0, TEMP0, si_len, si_u_max, si_rho, si_temp); // setting the conversion factor 

    LBM lb(NX, NY, NZ, NU);
    int Nx = lb.get_Nx(); int Ny = lb.get_Ny(); int Nz = lb.get_Nz();

    
    #pragma omp parallel for schedule(dynamic)
    for(int i = 0; i < Nx ; ++i)
    {
        for(int j = 0; j < Ny; ++j)
        {
            for(int k = 0; k < Nz; ++k)
            {
                if ( j==0 || j==Ny-1 || k==0 || k==Nz-1) // set periodic boundary condition
                {
                    lb.mixture[i][j][k].type = TYPE_P;
                }
                
                if (i==0 || i==Nx-1 )
                {
                    lb.mixture[i][j][k].type = TYPE_S;
                }

                if (lb.mixture[i][j][k].type == TYPE_F)
                {
                    if ((float)i/(float)Nx <= 0.5 )
                    {
                        lb.mixture[i][j][k].u = 0.0;
                        lb.mixture[i][j][k].v = 0.0;
                        lb.mixture[i][j][k].w = 0.0;
                        lb.mixture[i][j][k].rho = 0.5;
                        lb.mixture[i][j][k].temp = 0.2;
                    }
                    else
                    {
                        lb.mixture[i][j][k].u = 0.0;
                        lb.mixture[i][j][k].v = 0.0;
                        lb.mixture[i][j][k].w = 0.0;
                        lb.mixture[i][j][k].rho = 2.0;
                        lb.mixture[i][j][k].temp = 0.025;
                    }                       
                }
            }
        }
    }

    lb.run(1000,100);
}

#elif defined SOD_SHOCK_SIUNIT
void main_setup() // SOD SHOCK TUBE WITH SI UNIT --------------------------------------------------------
{
    std::cout << "-------------------- SOD SCHOCK SIUNIT ---------------------" << std::endl;
    int NX = 200; 
    int NY = 1; 
    int NZ = 1;
       
    //  std::vector<std::string> species = { "Ar", "H2" };
    // std::vector<std::string> species = {"H2"};
     std::vector<std::string> species = { "H2", "Ar", "CH4" };
    LBM lb(NX, NY, NZ, species);
    int Nx = lb.get_Nx(); int Ny = lb.get_Ny(); int Nz = lb.get_Nz();

    auto sol = Cantera::newSolution("gri30.yaml", "gri30");
    auto gas = sol->thermo();   
    std::vector <double> X (gas->nSpecies());
    X[gas->speciesIndex("H2")] = 1;
    X[gas->speciesIndex("Ar")] = 1;
    X[gas->speciesIndex("CH4")] = 1;

    // X[gas->speciesIndex("H2")] = 1.0;
    gas->setMoleFractions(&X[0]);
    gas->setState_TP(300.0, 1*Cantera::OneAtm);
    auto trans = sol->transport();

    double dx = 1e-6;
    units.set_m_kg_s(dx, 1e-4*dx);
    // units.set_m_kg_s(0.00001, gas->soundSpeed());
    // units.set_m_kg_s(0.001, 7*gas->soundSpeed());
    // units.set_m_kg_s(NX, VEL0, RHO0, 0.025, si_len, si_u_max, si_rho, si_temp); // setting the conversion factor 


    double NU = units.nu(trans->viscosity() / gas->density());
    double prandtl = trans->viscosity()*gas->cp_mass() / trans->thermalConductivity();
    
    std::cout << "nu : " << NU << std::endl;
    std::cout << "nu (unit) : " << trans->viscosity() / gas->density() << std::endl;
    std::cout << "prandtl number : " << prandtl << std::endl;
    std::cout << "gas const : " << units.cp(Cantera::GasConstant/gas->meanMolecularWeight()) << std::endl;
    std::cout << "gas const (unit) : " << Cantera::GasConstant/gas->meanMolecularWeight() << std::endl;
    std::cout << "temperature : " << units.temp(gas->temperature()) << std::endl;
    std::cout << "gamma : " << gas->cp_mass()/gas->cv_mass() << std::endl;
    std::cout << "sound speed : " << units.u(gas->soundSpeed()) << std::endl;
    // std::cout << "ux mean : " <<  u0 << std::endl;
    std::cout << "RT : " << units.cp(Cantera::GasConstant/gas->meanMolecularWeight())*units.temp(gas->temperature()) << std::endl;


    #pragma omp parallel for schedule(dynamic)
    for(int i = 0; i < Nx ; ++i)
    {
        for(int j = 0; j < Ny; ++j)
        {
            for(int k = 0; k < Nz; ++k)
            {
                if ( j==0 || j==Ny-1 || k==0 || k==Nz-1) // set periodic boundary condition
                {
                    lb.mixture[i][j][k].type = TYPE_P;
                }
                
                if (i==0 || i==Nx-1 )
                {
                    lb.mixture[i][j][k].type = TYPE_A;
                }

                if (lb.mixture[i][j][k].type == TYPE_F)
                {
                    lb.species[0][i][j][k].X = smooth(0.491, 0.000, i, 0.5*Nx, 1.0);
                    lb.species[1][i][j][k].X = smooth(0.509, 0.485, i, 0.5*Nx, 1.0);
                    lb.species[2][i][j][k].X = smooth(0.000, 0.515, i, 0.5*Nx, 1.0);

                    // lb.species[0][i][j][k].X = smooth(0.45, 0.55, i, 0.5*Nx, 0.07);
                    // lb.species[1][i][j][k].X = smooth(0.55, 0.45, i, 0.5*Nx, 0.07);
                    // lb.species[2][i][j][k].X = smooth(0.58, 0.43, i, 0.5*Nx, 0.07);

                    // lb.species[0][i][j][k].X = 1.0;
                    // lb.species[1][i][j][k].X = 1.0;


                    // lb.species[0][i][j][k].X = smooth(0.52, 0.48, i, 0.5*Nx, 0.3);
                    // lb.species[1][i][j][k].X = smooth(0.48, 0.52, i, 0.5*Nx, 0.3);

                    // lb.mixture[i][j][k].p = smooth(2*units.p(Cantera::OneAtm), units.p(Cantera::OneAtm), i, 0.5*Nx, 0.3);      
                    // lb.mixture[i][j][k].temp = smooth(units.temp(400.0), units.temp(400.0), i, 0.5*Nx, 0.3);                  

                    lb.mixture[i][j][k].p = smooth(1.00*units.p(Cantera::OneAtm), units.p(Cantera::OneAtm), i, 0.5*Nx, 0.1);      
                    lb.mixture[i][j][k].temp = smooth(units.temp(300.0), units.temp(300.0), i, 0.5*Nx, 0.1); 

                    if ((float)i/(float)Nx <= 0.5 )
                    {
                        lb.mixture[i][j][k].u = 0.0;
                        lb.mixture[i][j][k].v = 0.0;
                        lb.mixture[i][j][k].w = 0.0; 
                    }
                    else 
                    {
                        lb.mixture[i][j][k].u = 0.0;
                        lb.mixture[i][j][k].v = 0.0;
                        lb.mixture[i][j][k].w = 0.0;
                    }                       
                }
            }
        }
    }

    lb.run(100000000,1000);
}

#elif defined TERNARY_DIFFUSION

void main_setup() // Ternary Gas Diffusion --------------------------------------------------------
{
    int NX = 864; 
    int NY = 1; 
    int NZ = 1;
       
    //  std::vector<std::string> species = { "Ar", "N2" };
    //  std::vector<std::string> species = { "Ar" };
    // std::vector<std::string> species = {"H2"};
     std::vector<std::string> species = { "H2", "Ar", "CH4" };
    LBM lb(NX, NY, NZ, species);
    int Nx = lb.get_Nx(); int Ny = lb.get_Ny(); int Nz = lb.get_Nz();

    auto sol = Cantera::newSolution("gri30.yaml", "gri30");
    auto gas = sol->thermo();   
    std::vector <double> X (gas->nSpecies());
    X[gas->speciesIndex("H2")] = 1;
    X[gas->speciesIndex("Ar")] = 1;
    X[gas->speciesIndex("CH4")] = 1;

    gas->setMoleFractions(&X[0]);
    gas->setState_TP(300.0, 1*Cantera::OneAtm);
    auto trans = sol->transport();

    double dx = 1e-6;
    units.set_m_kg_s(dx, 1E-4*dx);
    // units.set_m_kg_s(0.000001, gas->soundSpeed());
    // units.set_m_kg_s(0.00001, gas->soundSpeed());
    // units.set_m_kg_s(0.001, 7*gas->soundSpeed());
    // units.set_m_kg_s(NX, VEL0, RHO0, 0.025, si_len, si_u_max, si_rho, si_temp); // setting the conversion factor 


    double NU = units.nu(trans->viscosity() / gas->density());
    double prandtl = trans->viscosity()*gas->cp_mass() / trans->thermalConductivity();
    
    std::cout << "nu : " << NU << std::endl;
    std::cout << "nu (unit) : " << trans->viscosity() / gas->density() << std::endl;
    std::cout << "prandtl number : " << prandtl << std::endl;
    std::cout << "gas const : " << units.cp(Cantera::GasConstant/gas->meanMolecularWeight()) << std::endl;
    std::cout << "gas const (unit) : " << Cantera::GasConstant/gas->meanMolecularWeight() << std::endl;
    std::cout << "temperature : " << units.temp(gas->temperature()) << std::endl;
    std::cout << "gamma : " << gas->cp_mass()/gas->cv_mass() << std::endl;
    std::cout << "sound speed : " << units.u(gas->soundSpeed()) << std::endl;
    std::cout << "RT : " << units.cp(Cantera::GasConstant/gas->meanMolecularWeight())*units.temp(gas->temperature()) << std::endl;

    #pragma omp parallel for schedule(dynamic)
    for(int i = 0; i < Nx ; ++i)
    {
        for(int j = 0; j < Ny; ++j)
        {
            for(int k = 0; k < Nz; ++k)
            {
                if ( j==0 || j==Ny-1 || k==0 || k==Nz-1) // set periodic boundary condition
                {
                    lb.mixture[i][j][k].type = TYPE_P;
                }
                
                if (i==0 || i==Nx-1 )
                {
                    lb.mixture[i][j][k].type = TYPE_A;
                }

                if (lb.mixture[i][j][k].type == TYPE_F)
                {
                    lb.species[0][i][j][k].X = smooth(0.491, 1e-6, i, 0.5*Nx, 0.3);
                    lb.species[1][i][j][k].X = smooth(0.509, 0.485, i, 0.5*Nx, 0.3);
                    lb.species[2][i][j][k].X = smooth(1e-6, 0.515, i, 0.5*Nx, 0.3);

                    // lb.species[0][i][j][k].X = smooth(0.491, 0.300, i, 0.5*Nx, 0.03);
                    // lb.species[1][i][j][k].X = smooth(0.209, 0.185, i, 0.5*Nx, 0.03);
                    // lb.species[2][i][j][k].X = smooth(0.300, 0.515, i, 0.5*Nx, 0.03);

                    // lb.species[0][i][j][k].X = smooth(1.0, 0.5, i, 0.5*Nx, 0.3);
                    // lb.species[1][i][j][k].X = smooth(0.5, 1.0, i, 0.5*Nx, 0.3);

                    // lb.species[0][i][j][k].X = smooth(0.9, 0.1, i, 0.5*Nx, 0.3);
                    // lb.species[1][i][j][k].X = smooth(0.1, 0.9, i, 0.5*Nx, 0.3);

                    // lb.species[0][i][j][k].X = smooth(1.0, 1.0, i, 0.5*Nx, 0.3);
                    
                    lb.mixture[i][j][k].p = smooth(units.p(Cantera::OneAtm), units.p(Cantera::OneAtm), i, 0.5*Nx, 0.03);      
                    lb.mixture[i][j][k].temp = smooth(units.temp(300.0), units.temp(300.0), i, 0.5*Nx, 0.01);                  
                  
                    lb.mixture[i][j][k].u = 0.0;
                    lb.mixture[i][j][k].v = 0.0;
                    lb.mixture[i][j][k].w = 0.0;
                     
                }
            }
        }
    }

    lb.run(10000000,1000);
}


#elif defined SHEAR_LAYER_MULTICOMP
void main_setup() // 3D Shear layer multicomponent
{
    int NX = 100; 
    int NY = 100; 
    int NZ = 1;
    
    // units.set_m_kg_s(NX, VEL0, RHO0, 0.025, si_len, si_u_max, si_rho, si_temp); // setting the conversion factor 

    double dx = 1e-5;

    units.set_m_kg_s(dx, 1e-4*dx);

    double Re = 50.0;
    double H = NX*dx;

    auto sol = Cantera::newSolution("gri30.yaml", "gri30");
    auto gas = sol->thermo();
    gas->setState_TP(300.0, Cantera::OneAtm);
    std::vector <double> X(gas->nSpecies());
    X[gas->speciesIndex("N2")] = 0.9;
    X[gas->speciesIndex("Ar")] = 0.1;
    gas->setMoleFractions(&X[0]);
    auto trans = sol->transport();
    double nu = trans->viscosity()/gas->density();
    std::cout << "nu = " << units.nu(nu) << std::endl;
    std::cout << "gamma = " << gas->cp_mass() / gas->cv_mass() << std::endl;
    std::cout << "gas_const = " << units.cp(gas->RT()/gas->meanMolecularWeight()/gas->temperature()) << std::endl;
    std::cout << "prandtl = " << trans->viscosity() * gas->cp_mass() / trans->thermalConductivity() << std::endl;


    // double u_max = 0.0;
    double u_max = units.u(Re*nu/H);
    std::cout << "u_max : " << u_max << " || " << Re*nu/H << std::endl;
    std::cout << "Speed of Sound : " << units.u(gas->soundSpeed()) << std::endl;
    std::cout << "Mach number : " << u_max / units.u(gas->soundSpeed()) << std::endl;

    // std::vector<std::string> species = { "N2" };
    std::vector<std::string> species = { "N2" , "Ar"};
    
    LBM lb(NX, NY, NZ, species);
    int Nx = lb.get_Nx(); int Ny = lb.get_Ny(); int Nz = lb.get_Nz();


    #pragma omp parallel for schedule(dynamic)
    for(int i = 0; i < Nx ; ++i)
    {
        for(int j = 0; j < Ny; ++j)
        {
            for(int k = 0; k < Nz; ++k)
            {
                if (i==0 || i==Nx-1 || j==0 || j==Ny-1 || k==0 || k==Nz-1) // set periodic boundary condition
                {
                    lb.mixture[i][j][k].type = TYPE_P;
                }
                
                if (lb.mixture[i][j][k].type == TYPE_F)
                {                    
                    lb.mixture[i][j][k].p = smooth(1*units.p(Cantera::OneAtm), units.p(Cantera::OneAtm), i, 0.5*Nx, 0.3);                        
                    lb.mixture[i][j][k].v = 0.01*u_max*sin(2*M_PI*i/Nx);
                    lb.mixture[i][j][k].w = 0.0;

                    if ((float)j/(float)Ny <= 0.5){
                        lb.mixture[i][j][k].u = smooth(u_max, -u_max, j, 0.25*Ny, 0.3);
                        lb.species[0][i][j][k].X = smooth(0.1, 0.9, j, 0.25*Ny, 0.3);
                        lb.species[1][i][j][k].X = smooth(0.9, 0.1, j, 0.25*Ny, 0.3);
                        lb.mixture[i][j][k].temp = smooth(units.temp(1400.0), units.temp(300.0), j, 0.25*Ny, 0.3);
                    } 
                    else{
                        lb.mixture[i][j][k].u = smooth(-u_max, u_max, j, 0.75*Ny, 0.3);
                        lb.species[0][i][j][k].X = smooth(0.9, 0.1, j, 0.75*Ny, 0.3);
                        lb.species[1][i][j][k].X = smooth(0.1, 0.9, j, 0.75*Ny, 0.3);
                        lb.mixture[i][j][k].temp = smooth(units.temp(300.0), units.temp(1400.0), j, 0.75*Ny, 0.3);
                    }
                     
                }
            }
        }
    }

    lb.run(100000000,1000);
}


#elif defined PERFECTLY_STIRRED_REACTOR_3D
void main_setup() // Perfectly stirred reactor ------------------------------------------------------------------------------------------
{
    int NX = 4; 
    int NY = 4; 
    int NZ = 4;
    
    double si_temp = 1400.0;// [K]

    units.set_m_kg_s(1e-5, 1e-9);

    auto sol = Cantera::newSolution("h2o2.yaml", "ohmech");
    auto gas = sol->thermo();

    std::vector<std::string> species = gas->speciesNames();
    

    LBM lb(NX, NY, NZ, species);
    int Nx = lb.get_Nx(); int Ny = lb.get_Ny(); int Nz = lb.get_Nz();

    #pragma omp parallel for schedule(dynamic)
    for(int i = 0; i < Nx ; ++i)
    {
        for(int j = 0; j < Ny; ++j)
        {
            for(int k = 0; k < Nz; ++k)
            {
                if ( i==0 || i==Nx-1 || j==0 || j==Ny-1 || k==0 || k==Nz-1) // set periodic boundary condition
                {
                    lb.mixture[i][j][k].type = TYPE_P;
                }

                if (lb.mixture[i][j][k].type == TYPE_F)
                {
                    lb.mixture[i][j][k].temp = units.temp(si_temp);
                    lb.mixture[i][j][k].p = 1.0*units.p(Cantera::OneAtm); 
                
                    lb.mixture[i][j][k].u = 0.0;
                    lb.mixture[i][j][k].v = 0.0;
                    lb.mixture[i][j][k].w = 0.0; 
                    lb.species[gas->speciesIndex("H2")][i][j][k].X = 0.29586;     // H2 Mole Fraction
                    lb.species[gas->speciesIndex("O2")][i][j][k].X = 0.14793;     // O2 Mole Fraction
                    lb.species[gas->speciesIndex("N2")][i][j][k].X = 0.55621;     // N2 Mole Fraction
                }
            }
        }
    }

    lb.run(100000,1000);
}

#elif defined CONDUCTION_1D

void main_setup() // 1D Heat conduction test case --------------------------------------------------------
{
    int NX = 3000; 
    int NY = 5; 
    int NZ = 5;
    
    double si_len = 1E-4;    // [m]
    double si_u_max = 1.0E+3;  // [m/s]
    double si_rho = 1.225;  // [kg/m^3]
    double si_temp = 300.0;// [K]
    

    units.set_m_kg_s(NX, VEL0, RHO0, TEMP0, si_len, si_u_max, si_rho, si_temp); // setting the conversion factor 

    std::vector<std::string> species = { "AR" };
    
    LBM lb(NX, NY, NZ, species);

    int Nx = lb.get_Nx(); int Ny = lb.get_Ny(); int Nz = lb.get_Nz(); int nSpecies = lb.get_nSpecies();

    auto sol = Cantera::newSolution("gri30.yaml", "gri30");
    auto gas = sol->thermo();
    
    // double gamma = gas->cp_mass() / gas->cv_mass();
    // double a_sound = sqrt(gamma*Cantera::GasConstant*);

    std::vector<double> XL (gas->nSpecies());
    std::vector<double> XR (gas->nSpecies()); 

    XL[gas->speciesIndex("AR")]  = 1.0;

    XR[gas->speciesIndex("AR")]  = 1.0;


    #pragma omp parallel for schedule(dynamic)
    for(int i = 0; i < Nx ; ++i)
    {
        for(int j = 0; j < Ny; ++j)
        {
            for(int k = 0; k < Nz; ++k)
            {
                if ( j==0 || j==Ny-1 || k==0 || k==Nz-1) // set periodic boundary condition
                {
                    lb.mixture[i][j][k].type = TYPE_P;
                }
                
                if (i==0 || i==Nx-1 )
                {
                    lb.mixture[i][j][k].type = TYPE_S;
                }

                if (lb.mixture[i][j][k].type == TYPE_F)
                {
                    lb.mixture[i][j][k].p = units.p(Cantera::OneAtm);
                    if ((float)i/(float)Nx <= 0.5 )
                    {
                        for (int a = 0; a < nSpecies; ++a) lb.species[a][i][j][k].X = XL[gas->speciesIndex(species[a])];
                        lb.mixture[i][j][k].temp = units.temp(300.0);
                    }
                    else
                    {
                        for (int a = 0; a < nSpecies; ++a) lb.species[a][i][j][k].X = XR[gas->speciesIndex(species[a])];
                        lb.mixture[i][j][k].temp = units.temp(350.0);
                    }        
                }
            }
        }
    }

    lb.run(10000,100);
}

#elif defined COUETTE_FLOW

void main_setup() // 2D Couette Flow --------------------------------------------------------
{
    int NX = 1; 
    int NY = 50; 
    int NZ = 1;
    
    const double Ec = 40.0;
    const double prandtl = 0.7;
    const double gasconst = 1.0;
    const double gamma = 1.4;
    const double Ma = 0.01;
    const double Re = 10.0;
    const double temp_c = 0.025;

    const double U = Ma*sqrt(gamma*gasconst*temp_c);
    const double nu = U*NY/Re;
    const double temp_h = temp_c + U*U / (gamma*gasconst/(gamma-1.0)*Ec);
    
    LBM lb(NX, NY, NZ, nu);
    lb.set_prtl(prandtl);
    lb.set_gasconst(gasconst);
    lb.set_gamma(gamma);

    int Nx = lb.get_Nx(); int Ny = lb.get_Ny(); int Nz = lb.get_Nz();

    std::cout << "Ec number : " << Ec << std::endl;
    std::cout << "prandtl number : " << prandtl << std::endl;
    std::cout << "Re number : " << Re << std::endl;
    std::cout << "nu : " << nu << std::endl;
    std::cout << "sound speed : " << lb.get_soundspeed(temp_c) << std::endl;
    std::cout << "hot temperature : " << temp_h << std::endl;
    std::cout << "U (lattice unit) : " << U << std::endl;
    std::cout << "gamma : " << gamma << std::endl;
    std::cout << "gas const : " << gasconst << std::endl;


    #pragma omp parallel for schedule(dynamic)
    for(int i = 0; i < Nx ; ++i)
    {
        for(int j = 0; j < Ny; ++j)
        {
            for(int k = 0; k < Nz; ++k)
            {                
                if ( i==0 || i==Nx-1 || k==0 || k==Nz-1) // set periodic boundary condition
                {
                    lb.mixture[i][j][k].type = TYPE_P;
                }

                if (j==0 )
                {
                    lb.mixture[i][j][k].type = TYPE_S;
                    lb.mixture[i][j][k].u = 0.0;
                    lb.mixture[i][j][k].temp = temp_c;
                }
                if (j==Ny-1 )
                {
                    lb.mixture[i][j][k].type = TYPE_S;
                    lb.mixture[i][j][k].u = U;
                    lb.mixture[i][j][k].temp = temp_h;
                }

                if (lb.mixture[i][j][k].type == TYPE_F)
                {
                    lb.mixture[i][j][k].u = 0.0;
                    lb.mixture[i][j][k].temp = temp_c;
                    lb.mixture[i][j][k].p = gasconst*temp_c;
                }
            }
        }
    }

    lb.run(150000,1000);
    for (int j = 0; j < Ny; ++j){
        std::cout << lb.mixture[1][j][1].temp << std::endl; 
    }
}

#elif defined COUETTE_FLOW_MULTICOMP
void main_setup() // 2D Couette flow (multicomponent)  --------------------------------------------------------
{  
    int NX = 1; 
    int NY = 49; 
    int NZ = 1;
    
    const double Ec = 10.0;
    const double Ma = 0.2;
    // const double Re = 100;

    const double temp0 = 300.0;
    const double p0 = 101325.0;

    std::vector<std::string> species = { "AR" };
    
    LBM lb(NX, NY, NZ, species);
    int Nx = lb.get_Nx(); int Ny = lb.get_Ny(); int Nz = lb.get_Nz();
    size_t nSpecies = lb.get_nSpecies();

    auto sol = Cantera::newSolution("gri30.yaml", "gri30");
    auto gas = sol->thermo();   
    std::vector <double> X (gas->nSpecies());
    X[gas->speciesIndex("AR")] = 1.0;
    gas->setMoleFractions(&X[0]);
    gas->setState_TP(temp0, p0);
    auto trans = sol->transport();
    const double NU = trans->viscosity() / gas->density();
    const double U = Ma*gas->soundSpeed();
    const double temp_h = temp0 + U*U / (gas->cp_mass()*Ec);

    // units.set_m_kg_s(NY, 1.0, 1.0, 0.025, 1E-3, 1.0, 1.0, 300.0); // setting the conversion factor 
    double dx = 1e-6;
    units.set_m_kg_s(dx, 1e-3*dx);

    std::cout << "Ec number : " << Ec << std::endl;
    std::cout << "prandtl number : " << trans->viscosity()*gas->cp_mass()/trans->thermalConductivity() << std::endl;
    std::cout << "Re number : " << U * units.si_x(NY) / NU << std::endl;
    std::cout << "nu : " << units.nu(NU) << std::endl;
    std::cout << "sound speed : " << units.u(gas->soundSpeed()) << std::endl;
    std::cout << "hot temperature : " << units.temp(temp_h) << std::endl;
    std::cout << "U (lattice unit) : " << units.u(U) << std::endl;
    std::cout << "gamma : " << gas->cp_mass()/gas->cv_mass() << std::endl;
    std::cout << "gas const : " << units.cp(gas->cp_mass() - gas->cv_mass()) << std::endl;

    #pragma omp parallel for schedule(dynamic)
    for(int i = 0; i < Nx ; ++i)
    {
        for(int j = 0; j < Ny; ++j)
        {
            for(int k = 0; k < Nz; ++k)
            {                
                if ( i==0 || i==Nx-1 || k==0 || k==Nz-1) // set periodic boundary condition
                {
                    lb.mixture[i][j][k].type = TYPE_P;
                }

                if (j==0 )
                {
                    lb.mixture[i][j][k].type = TYPE_S;
                    lb.mixture[i][j][k].u = 0.0;
                    lb.mixture[i][j][k].temp = units.temp(temp0);
                }
                if (j==Ny-1 )
                {
                    lb.mixture[i][j][k].type = TYPE_S;
                    lb.mixture[i][j][k].u = units.u(U);
                    lb.mixture[i][j][k].temp = units.temp(temp_h);
                }

                if (lb.mixture[i][j][k].type == TYPE_F)
                {
                    for(size_t a = 0; a < nSpecies; ++a)
                        lb.species[a][i][j][k].X = 1.0;

                    lb.mixture[i][j][k].u = 0.0;
                    lb.mixture[i][j][k].temp = units.temp(temp0);
                    lb.mixture[i][j][k].p = units.p(p0);
                }
            }
        }
    }

    lb.run(150000,1000);
    for (int j = 0; j < Ny; ++j){
        std::cout << lb.mixture[1][j][1].temp << std::endl; 
    }
}

#elif defined TAYLOR_GREEN_3D_MULTICOMP
void main_setup() // 3D Taylor-Green Vortex (multicomponent) ------------------------------------------------------------------------------------
{
    double RE = 200;
    int NX = 100; 
    int NY = 100; 
    int NZ = 100;
    
    units.set_m_kg_s(1E-5, 1e-9);

    // std::vector<std::string> species = { "CO2" , "AR" };
    // std::vector<std::string> species = { "AR" };
    std::vector<std::string> species = { "H2O", "N2" };
    double pressure = Cantera::OneAtm;
    double temp = 300.0;

    LBM lb(NX, NY, NZ, species);

    int Nx = lb.get_Nx(); int Ny = lb.get_Ny(); int Nz = lb.get_Nz();
    int nSpecies = lb.get_nSpecies();

    double periodicity = 1.0;
    const double a = (double)NX/periodicity;
    const double b = (double)NY/periodicity;
    const double c = (double)NZ/periodicity;

    auto sol = Cantera::newSolution("gri30.yaml", "gri30");
    auto gas = sol->thermo();   
    std::vector <double> X (gas->nSpecies());
    // X[gas->speciesIndex("CO2")] = 1.0;
    // X[gas->speciesIndex("AR")] = 1.0;
    X[gas->speciesIndex("H2O")] = 1.0;
    X[gas->speciesIndex("N2")] = 1.0;   

    gas->setMoleFractions(&X[0]);
    gas->setState_TP(temp, pressure);
    auto trans = sol->transport();
    double NU = units.nu(trans->viscosity() / gas->density());
    std::cout << "NU : " << NU << std::endl;

    double u_max = RE * NU / (NX/(2.*M_PI));
    std::cout << "umax : " << u_max << std::endl; 

    double spd_sound = units.u(gas->soundSpeed());
    std::cout << "Mach : " << u_max/spd_sound << std::endl;

    double gas_const = units.cp(Cantera::GasConstant/gas->meanMolecularWeight());
    std::cout << "gas const : " << gas_const << std::endl;

    const double rho0 = units.rho(gas->density());
    const double p0 = units.p(gas->pressure());
    const double temp0 = units.temp(gas->temperature()); 
    

    
    #pragma omp parallel for schedule(dynamic)
    for(int i = 0; i < Nx ; ++i)
    {
        for(int j = 0; j < Ny; ++j)
        {
            for(int k = 0; k < Nz; ++k)
            {
                if (i==0 || i==Nx-1 ||  j==0 || j==Ny-1 || k==0 || k==Nz-1) // set periodic boundary condition
                {
                    lb.mixture[i][j][k].type = TYPE_P;
                }
                if (lb.mixture[i][j][k].type == TYPE_F)
                {
                    lb.species[0][i][j][k].X = 1.0;
                    lb.species[1][i][j][k].X = 1.0;

                    lb.mixture[i][j][k].u =  u_max * sin(2.0*M_PI/a*(double)i) * cos(2.0*M_PI/b*(double)j) * cos(2.0*M_PI/c*(double)k);
                    lb.mixture[i][j][k].v = -u_max * cos(2.0*M_PI/a*(double)i) * sin(2.0*M_PI/b*(double)j) * cos(2.0*M_PI/c*(double)k);
                    lb.mixture[i][j][k].w = 0.0;
                    lb.mixture[i][j][k].p = p0 + rho0*u_max*u_max/16.0 * (cos(4.0*M_PI/a*(double)i) + cos(4.0*M_PI/b*(double)j)) * (cos(4.0*M_PI/c*(double)k) + 2.0);
                    lb.mixture[i][j][k].temp = temp0; 
                }
            }
        }
    }
    
    lb.run(1000000,1);
}

#elif defined SHEAR_LAYER
void main_setup() // 3D Shear layer
{
    int NX = 400; 
    int NY = 400; 
    int NZ = 1;
    
    double Re = 1000.0;
    double u_max = 0.1;
    double nu = u_max*NX/Re;
    std::cout << nu << std::endl;

    LBM lb(NX, NY, NZ, nu);
    lb.set_prtl(1.0);

    int Nx = lb.get_Nx(); int Ny = lb.get_Ny(); int Nz = lb.get_Nz();

    double dx = 3e-5;
    units.set_m_kg_s(dx, 1e-9);


    #pragma omp parallel for schedule(dynamic)
    for(int i = 0; i < Nx ; ++i)
    {
        for(int j = 0; j < Ny; ++j)
        {
            for(int k = 0; k < Nz; ++k)
            {
                if (i==0 || i==Nx-1 || j==0 || j==Ny-1 || k==0 || k==Nz-1) // set periodic boundary condition
                {
                    lb.mixture[i][j][k].type = TYPE_P;
                }
                
                if (lb.mixture[i][j][k].type == TYPE_F)
                {                    
                    lb.mixture[i][j][k].p = 0.025; //units.p(101325);      
                    lb.mixture[i][j][k].temp = 0.025; // units.temp(450);                  
                  
                    lb.mixture[i][j][k].v = 0.05*u_max*sin(2*M_PI*i/Nx);
                    lb.mixture[i][j][k].w = 0.0;

                    if ((float)j/(float)Ny <= 0.5){
                        lb.mixture[i][j][k].u = smooth(u_max, -u_max, j, 0.25*Ny, 0.7);
                    } 
                    else{
                        lb.mixture[i][j][k].u = smooth(-u_max, u_max, j, 0.75*Ny, 0.7);
                    }
                     
                }
            }
        }
    }

    lb.run(100000,100);
}


#elif defined SHOCK_VORTEX_INTERAC
void main_setup() // 3D Taylor-Green Vortex
{
    const int NX = 1680;
    const int NY = 1440;
    const int NZ = 1;

    const double gamma = 1.4;
    const double gas_const = 1.0;
    const double prantdl = 0.75;
    const double Re = 800;

    const double Ma_v = 0.25;
    const double r_v = NX/28;
    const double x_v = 6*r_v;
    const double y_v = NY/2;

    const double Ma_s = 1.2;
    const double x_s = 8.0*r_v;

    const double rho_l = 1.0;
    const double temp_l = 0.05;
    const double ux_l = Ma_s*sqrt(gamma*temp_l);
    const double uy = 0.0;

    const double rho_r = (gamma+1.0)*pow(Ma_s,2.0) / ((gamma-1.0)*pow(Ma_s,2.0)+2.0) * rho_l;
    const double temp_r = (1.0+(gamma-1.0)/2.0*pow(Ma_s,2.0)) * (2.0*gamma/(gamma-1.0)*pow(Ma_s,2.0)-1.0) / (pow(Ma_s,2.0)*((2.0*gamma)/(gamma-1.0)+(gamma-1.0)/2.0)) * temp_l;
    const double Ma_s_r = ((gamma-1.0)*pow(Ma_s,2.0) + 2.0) / (2.0*gamma*pow(Ma_s,2.0)-(gamma-1.0));
    
    const double nu = sqrt(gamma*temp_l) * r_v / Re;

    LBM lb(NX, NY, NZ, nu);
    lb.set_gamma(gamma);
    lb.set_gasconst(gas_const);
    lb.set_prtl(prantdl);
    int Nx = lb.get_Nx(); int Ny = lb.get_Ny(); int Nz = lb.get_Nz();
   
    std::cout << "umax : " << ux_l << std::endl;

    #pragma omp parallel for schedule(dynamic)
    for(int i = 0; i < Nx ; ++i)
    {
        for(int j = 0; j < Ny; ++j)
        {
            for(int k = 0; k < Nz; ++k)
            {
                if (i==0 || i==Nx-1 ||  j==0 || j==Ny-1 || k==0 || k==Nz-1) // set periodic boundary condition
                {
                    lb.mixture[i][j][k].type = TYPE_P;
                }
                if (lb.mixture[i][j][k].type == TYPE_F)
                {
                    lb.mixture[i][j][k].u =  u_max * sin(2.0*M_PI/a*(double)i) * cos(2.0*M_PI/b*(double)j) * cos(2.0*M_PI/c*(double)k);
                    lb.mixture[i][j][k].v = -u_max * cos(2.0*M_PI/a*(double)i) * sin(2.0*M_PI/b*(double)j) * cos(2.0*M_PI/c*(double)k);
                    lb.mixture[i][j][k].w = 0.0;
                    lb.mixture[i][j][k].p = p0 + rho0*u_max*u_max/16.0 * (cos(4.0*M_PI/a*(double)i) + cos(4.0*M_PI/b*(double)j)) * (cos(4.0*M_PI/c*(double)k) + 2.0);
                    lb.mixture[i][j][k].temp = temp0;
                }
            }
        }
    }
    
    lb.run(10000000,10);
}

#elif defined RB_INSTABILITY
void main_setup() // 2D Rayleigh-Benard Instability ------------------------------------------------------
{
    int NX = 80; 
    int NY = 40; 
    int NZ = 1;
    
    const double Ra = 2E4;
    const double prandtl = 0.7;
    const double gasconst = 1.0;
    const double gamma = 1.4;
    const double nu = 0.02;

    const double rho0 = 1.0;
    const double Temp_h = 0.050;
    const double Temp_c = 0.025;   
    const double delta_Temp = Temp_h-Temp_c;
    
    LBM lb(NX, NY, NZ, nu);
    lb.set_prtl(prandtl);
    lb.set_gasconst(gasconst);
    lb.set_gamma(gamma);
    lb.set_Ra(Ra);
    
    const double G_lambda = Ra*(nu*nu/prandtl) / (delta_Temp*NY*NY*NY);

    int Nx = lb.get_Nx(); int Ny = lb.get_Ny(); int Nz = lb.get_Nz();

    std::cout << "Ra number : " << Ra << std::endl;
    std::cout << "prandtl number : " << prandtl << std::endl;
    std::cout << "nu : " << nu << std::endl;
    std::cout << "sound speed : " << lb.get_soundspeed(Temp_c) << std::endl;
    std::cout << "gamma : " << gamma << std::endl;
    std::cout << "gas const : " << gasconst << std::endl;


    #pragma omp parallel for schedule(dynamic)
    for(int i = 0; i < Nx ; ++i)
    {
        for(int j = 0; j < Ny; ++j)
        {
            for(int k = 0; k < Nz; ++k)
            {
                if (i==0 || i==Nx-1 || k==0 || k==Nz-1) // set periodic boundary condition
                {
                    lb.mixture[i][j][k].type = TYPE_P;
                }
                if ( j==0 )
                {
                    lb.mixture[i][j][k].type = TYPE_S;
                    lb.mixture[i][j][k].temp = Temp_h;
                }
                if ( j==Ny-1 )
                {
                    lb.mixture[i][j][k].type = TYPE_S;
                    lb.mixture[i][j][k].temp = Temp_c;
                }

                if (lb.mixture[i][j][k].type == TYPE_F)
                {
                    lb.mixture[i][j][k].temp = Temp_h - delta_Temp*j/Ny;
                    lb.mixture[i][j][k].p = (1.0 + rho0*G_lambda*delta_Temp*j/2.0*(1.0-j/Ny))*(1.0+0.001*cos(2.0*M_PI*i/(Nx-2.0)));
                }
            }
        }
    }

    lb.run(50000,1000);
}

#elif defined RB_INSTABILITY_MULTICOMP
void main_setup() // 2D Multicomp Rayleigh-Benard Instability ------------------------------------------------------
{
    int NX = 80; 
    int NY = 40; 
    int NZ = 1;

    const double Ra = 2E4;
    
    // std::vector<std::string> species = {"Ar", "H2", "CH4"};
    std::vector<std::string> species = {"Ar", "CH4"};
    LBM lb(NX, NY, NZ, species);
    int Nx = lb.get_Nx(); int Ny = lb.get_Ny(); int Nz = lb.get_Nz();
    lb.set_Ra(Ra);

    auto sol = Cantera::newSolution("gri30.yaml", "gri30");
    auto gas = sol->thermo();   
    std::vector <double> X (gas->nSpecies());
    X[gas->speciesIndex("Ar")] = 1.0;
    // X[gas->speciesIndex("H2")] = 1.0;
    X[gas->speciesIndex("CH4")] = 1.0;
    gas->setMoleFractions(&X[0]);
    gas->setState_TP(300.0, Cantera::OneAtm);
    auto trans = sol->transport();

    units.set_m_kg_s(NY, 1E-4, gas->soundSpeed());
    // units.set_m_kg_s(1e-6,1e-9);
    
    double NU = units.nu(trans->viscosity() / gas->density());
    double prandtl = trans->viscosity()*gas->cp_mass() / trans->thermalConductivity();
    double rho0 = units.rho(gas->density());

    const double Temp_h = units.temp(350.0); // [K]
    const double Temp_c = units.temp(300.0); // [K]
    const double delta_Temp = Temp_h-Temp_c;
    const double pressure = units.p(Cantera::OneAtm);
    const double G_lambda = Ra*(NU*NU/prandtl) / (delta_Temp*NY*NY*NY);

    std::cout << "Ra number : " << Ra << std::endl;
    std::cout << "nu : " << NU << std::endl;
    std::cout << "prandtl number : " << prandtl << std::endl;
    std::cout << "gas const : " << units.cp(gas->RT()/gas->temperature()) << std::endl;
    std::cout << "gamma : " << gas->cp_mass()/gas->cv_mass() << std::endl;
    std::cout << "sound speed : " << units.u(gas->soundSpeed()) << std::endl;


    #pragma omp parallel for schedule(dynamic)
    for(int i = 0; i < Nx ; ++i)
    {
        for(int j = 0; j < Ny; ++j)
        {
            for(int k = 0; k < Nz; ++k)
            {
                if (i==0 || i==Nx-1 || k==0 || k==Nz-1) // set periodic boundary condition
                {
                    lb.mixture[i][j][k].type = TYPE_P;
                }
                if ( j==0 )
                {
                    lb.mixture[i][j][k].type = TYPE_A;
                    // lb.mixture[i][j][k].temp = Temp_h;
                }
                if ( j==Ny-1 )
                {
                    lb.mixture[i][j][k].type = TYPE_A;
                    // lb.mixture[i][j][k].temp = Temp_c;
                }

                if (lb.mixture[i][j][k].type == TYPE_F)
                {
                    lb.species[0][i][j][k].X = smooth(1.0, 0.0, j, 0.5*Ny, 0.7);
                    lb.species[1][i][j][k].X = 1.0 - lb.species[0][i][j][k].X;
                    // lb.species[2][i][j][k].X = 1.0 - lb.species[0][i][j][k].X;
                    lb.mixture[i][j][k].temp =Temp_h - delta_Temp*j/Ny;
                    lb.mixture[i][j][k].p = pressure;//(pressure+rho0*0.001*delta_Temp*j/2.0*(1.0-j/Ny))*(1.0+0.001*cos(2.0*M_PI*i/(Nx-2.0)));//(1.0 + rho0*G_lambda*delta_Temp*j/2.0*(1.0-j/Ny))*(1.0+0.001*cos(2.0*M_PI*i/(Nx-2.0)));
                }
            }
        }
    }

    lb.run(10000000,1000);
}

#elif defined CONDUCTION_BLACK
void main_setup() // 2D Rayleigh-Benard Instability ------------------------------------------------------
{
    int NX = 200; 
    int NY = 1; 
    int NZ = 1;
    
    const double prandtl = 0.7;
    const double gasconst = 1.0;
    const double gamma = 1.4;
    const double nu = 0.025;

    const double rho0 = 1.0;
    const double Temp_h = 0.2;
    const double Temp_c = 0.1;
    
    LBM lb(NX, NY, NZ, nu);
    lb.set_prtl(prandtl);
    lb.set_gasconst(gasconst);
    lb.set_gamma(gamma);
    
    int Nx = lb.get_Nx(); int Ny = lb.get_Ny(); int Nz = lb.get_Nz();

    std::cout << "prandtl number : " << prandtl << std::endl;
    std::cout << "nu : " << nu << std::endl;
    std::cout << "sound speed : " << lb.get_soundspeed(Temp_c) << std::endl;
    std::cout << "gamma : " << gamma << std::endl;
    std::cout << "gas const : " << gasconst << std::endl;


    #pragma omp parallel for schedule(dynamic)
    for(int i = 0; i < Nx ; ++i)
    {
        for(int j = 0; j < Ny; ++j)
        {
            for(int k = 0; k < Nz; ++k)
            {                
                if ( j==0 || j==Ny-1 || k==0 || k==Nz-1) // set periodic boundary condition
                {
                    lb.mixture[i][j][k].type = TYPE_P;
                }
                if ( i==0 || i==Nx-1 ) // set periodic boundary condition
                {
                    lb.mixture[i][j][k].type = TYPE_S;
                    if (i == 0)
                        lb.mixture[i][j][k].temp = Temp_h;
                    else
                        lb.mixture[i][j][k].temp = Temp_c;
                }

                if (lb.mixture[i][j][k].type == TYPE_F)
                {
                    lb.mixture[i][j][k].u = 0.0;
                    lb.mixture[i][j][k].rho = rho0;
                    lb.mixture[i][j][k].p = 0.1;

                    if (i < Nx/2)
                        lb.mixture[i][j][k].temp = Temp_h;
                    else
                        lb.mixture[i][j][k].temp = Temp_c;
                }
            }
        }
    }

    lb.run(1000000,1000);
}

#elif defined VISCOSITY_TEST_MULTICOMP
void main_setup() // 2D Viscos Test --------------------------------------------------------
{
    double NX = 1;
    double NY = 200;
    double NZ = 1;


    units.set_m_kg_s(1e-5, 1e-8);

    const double a0 = 0.001;
    const double mach = 0.1;
   
    // std::vector<std::string> species = { "H" };
    

    auto sol = Cantera::newSolution("h2o2.yaml", "ohmech");
    auto gas = sol->thermo();   
    std::vector <double> X (gas->nSpecies());
    X[gas->speciesIndex("H")] = 1.0;
    gas->setMoleFractions(&X[0]);
    gas->setState_TP(300.0, 1*Cantera::OneAtm);
    auto trans = sol->transport();

    std::vector<std::string> species = gas->speciesNames();
    LBM lb(NX, NY, NZ, species);
    int Nx = lb.get_Nx(); int Ny = lb.get_Ny(); int Nz = lb.get_Nz();

    double NU = units.nu(trans->viscosity() / gas->density());
    double prandtl = trans->viscosity()*gas->cp_mass() / trans->thermalConductivity();
    double u0 = mach*units.u(gas->soundSpeed());
    
    std::cout << "nu : " << NU << std::endl;
    std::cout << "nu (unit) : " << trans->viscosity() / gas->density() << std::endl;
    std::cout << "prandtl number : " << prandtl << std::endl;
    std::cout << "gas const : " << units.cp(Cantera::GasConstant/gas->meanMolecularWeight()) << std::endl;
    std::cout << "temperature : " << units.temp(gas->temperature()) << std::endl;
    std::cout << "gamma : " << gas->cp_mass()/gas->cv_mass() << std::endl;
    std::cout << "sound speed : " << units.u(gas->soundSpeed()) << std::endl;
    std::cout << "ux mean : " <<  u0 << std::endl;
    std::cout << "RT : " << units.cp(Cantera::GasConstant/gas->meanMolecularWeight())*units.temp(gas->temperature()) << std::endl;
    
    #pragma omp parallel for schedule(dynamic)
    for(int i = 0; i < Nx ; ++i)
    {
        for(int j = 0; j < Ny; ++j)
        {
            for(int k = 0; k < Nz; ++k)
            {
                if (i==0 || i==Nx-1 ||  j==0 || j==Ny-1 || k==0 || k==Nz-1) // set periodic boundary condition
                {
                    lb.mixture[i][j][k].type = TYPE_P;
                }
                if (lb.mixture[i][j][k].type == TYPE_F)
                {
                    //double X = i ;
                    double Y = j ;
                    lb.species[0][i][j][k].X = 1.0;
                    // lb.species[1][i][j][k].X = 1.0;
                    lb.mixture[i][j][k].u = a0*u0*sin(Y/NY*2*M_PI);// + u0*cos(M_PI/6);
                    lb.mixture[i][j][k].v = u0;//*sin(M_PI/6);
                    lb.mixture[i][j][k].w = 0.0;
                    lb.mixture[i][j][k].temp = units.temp(300.0);
                    lb.mixture[i][j][k].p = units.p(Cantera::OneAtm);
                    
                }
            }
        }
    }
    lb.run(10000,1);
}

#elif defined OPPOSED_JET
void main_setup() // 2D Opposed Jet --------------------------------------------------------
{
    double NX = 200;
    double NY = 400;
    double NZ = 1;

    double jet_width = 0.4*NY; 
   
    // std::vector<std::string> species = {"Ar"};
    // LBM lb(NX, NY, NZ, species);
    LBM lb(NX, NY, NZ, 0.1);
    int Nx = lb.get_Nx(); int Ny = lb.get_Ny(); int Nz = lb.get_Nz();

    // auto sol = Cantera::newSolution("gri30.yaml", "gri30");
    // auto gas = sol->thermo();   
    // std::vector <double> X (gas->nSpecies());
    // X[gas->speciesIndex("Ar")] = 1.0;
    // gas->setMoleFractions(&X[0]);
    // gas->setState_TP(300.0, 1*Cantera::OneAtm);
    // auto trans = sol->transport();

    // units.set_m_kg_s(1e-5, 1e-9);
 
    // double NU = units.nu(trans->viscosity() / gas->density());
    // double prandtl = trans->viscosity()*gas->cp_mass() / trans->thermalConductivity();
    double u0 = 0.02;
    
    // std::cout << "nu : " << NU << std::endl;
    // std::cout << "nu (unit) : " << trans->viscosity() / gas->density() << std::endl;
    // std::cout << "prandtl number : " << prandtl << std::endl;
    // std::cout << "gas const : " << units.cp(Cantera::GasConstant/gas->meanMolecularWeight()) << std::endl;
    // std::cout << "temperature : " << units.temp(gas->temperature()) << std::endl;
    // std::cout << "gamma : " << gas->cp_mass()/gas->cv_mass() << std::endl;
    // std::cout << "sound speed : " << units.u(gas->soundSpeed()) << std::endl;
    // std::cout << "ux mean : " <<  u0 << std::endl;
    // std::cout << "RT : " << units.cp(Cantera::GasConstant/gas->meanMolecularWeight())*units.temp(gas->temperature()) << std::endl;
    
    #pragma omp parallel for schedule(dynamic)
    for(int i = 0; i < Nx ; ++i){
        for(int j = 0; j < Ny; ++j){
            for(int k = 0; k < Nz; ++k){
                // set periodic boundary condition
                if ( k==0 || k==Nz-1 ) 
                    lb.mixture[i][j][k].type = TYPE_P;
                
                // Left and Right Boundary Condition
                if ( i==0 || i==Nx-1)
                    lb.mixture[i][j][k].type = TYPE_A;

                // Upper and Lower Boundary Condition
                if ( j==0 || j==Ny-1 )
                {
                    lb.mixture[i][j][k].type = TYPE_O;
                    lb.mixture[i][j][k].temp = 0.1;
                    lb.mixture[i][j][k].p = 0.1; 
                }

                // Left Inlet Boundary Condition
                if ( i==0 && j > NY/2 - jet_width && j < NY/2 + jet_width ) 
                {
                    lb.mixture[i][j][k].type = TYPE_I;
                    lb.mixture[i][j][k].u = u0;
                    lb.mixture[i][j][k].temp = 0.1;
                    lb.mixture[i][j][k].p = 0.1; 
                }
                // Right Inlet Boundary Condition
                if ( i==Nx-1 && j > NY/2 - jet_width && j < NY/2 + jet_width)
                {
                    lb.mixture[i][j][k].type = TYPE_I;
                    lb.mixture[i][j][k].u = -u0;
                    lb.mixture[i][j][k].temp = 0.1;
                    lb.mixture[i][j][k].p = 0.1; 
                }
                
                if (lb.mixture[i][j][k].type == TYPE_F)
                {
                    // lb.species[0][i][j][k].X = 1.0;
                    lb.mixture[i][j][k].u = 0.0;
                    lb.mixture[i][j][k].v = 0;
                    lb.mixture[i][j][k].w = 0.0;
                    lb.mixture[i][j][k].temp = 0.1;
                    lb.mixture[i][j][k].p = 0.1;                    
                }
            }
        }
    }
    lb.run(15000,100);
}

#elif defined OPPOSED_JET_MULTICOMP
void main_setup() // 2D Opposed Jet --------------------------------------------------------
{
    double NX = 20;
    double NY = 40;
    double NZ = 1;
   
    double jet_width = 0.2*NY; 
    double Mach = 0.1;
    double temp0 = 300.0;
    double press0 = Cantera::OneAtm;
    std::vector<std::string> species = {"H2", "N2" , "O2", "H2O"};
    
    LBM lb(NX, NY, NZ, species);
    int Nx = lb.get_Nx(); int Ny = lb.get_Ny(); int Nz = lb.get_Nz();

    auto sol = Cantera::newSolution("gri30.yaml", "gri30");
    auto gas = sol->thermo();   
    std::vector <double> X (gas->nSpecies());
    X[gas->speciesIndex("H2")] = 0.10;
    X[gas->speciesIndex("N2")] = 0.85;
    X[gas->speciesIndex("O2")] = 0.0;
    X[gas->speciesIndex("H2O")] = 0.05;
    gas->setMoleFractions(&X[0]);
    gas->setState_TP(temp0, press0);
    auto trans = sol->transport();

    units.set_m_kg_s(1e-6, 4e-10);
 
    double NU = units.nu(trans->viscosity() / gas->density());
    double prandtl = trans->viscosity()*gas->cp_mass() / trans->thermalConductivity();
    double sound_spd = units.u(gas->soundSpeed());

    double u0 = Mach * sound_spd;;
    
    std::cout << "nu : " << NU << std::endl;
    std::cout << "nu (unit) : " << trans->viscosity() / gas->density() << std::endl;
    std::cout << "prandtl number : " << prandtl << std::endl;
    std::cout << "gas const : " << units.cp(Cantera::GasConstant/gas->meanMolecularWeight()) << std::endl;
    std::cout << "temperature : " << units.temp(gas->temperature()) << std::endl;
    std::cout << "gamma : " << gas->cp_mass()/gas->cv_mass() << std::endl;
    std::cout << "sound speed : " << units.u(gas->soundSpeed()) << std::endl;
    std::cout << "ux mean : " <<  u0 << std::endl;
    std::cout << "RT : " << units.cp(Cantera::GasConstant/gas->meanMolecularWeight())*units.temp(gas->temperature()) << std::endl;
    
    #pragma omp parallel for schedule(dynamic)
    for(int i = 0; i < Nx ; ++i){
        for(int j = 0; j < Ny; ++j){
            for(int k = 0; k < Nz; ++k){
                // set periodic boundary condition
                if ( k==0 || k==Nz-1 ) 
                    lb.mixture[i][j][k].type = TYPE_P;
               
                // Upper and Lower Boundary Condition
                if ( j==0 || j==Ny-1 )
                {
                    lb.mixture[i][j][k].type = TYPE_O;
                    lb.mixture[i][j][k].temp = units.temp(temp0);
                    lb.mixture[i][j][k].p = units.p(press0); 

                    lb.species[1][i][j][k].X = 0.90;
                }
               
                if ( i==0 || i==Nx-1)
                {
                    lb.mixture[i][j][k].type = TYPE_I;
                    lb.mixture[i][j][k].u = 0.0;
                    lb.mixture[i][j][k].temp = units.temp(temp0);
                    lb.mixture[i][j][k].p = units.p(press0); 

                    lb.species[0][i][j][k].X = 0.05;
                    lb.species[1][i][j][k].X = 0.90;
                    lb.species[2][i][j][k].X = 0.05;
                    lb.species[3][i][j][k].X = 0.05;
                }

                // Left Inlet Boundary Condition
                if ( i==0 && j > NY/2 - jet_width && j < NY/2 + jet_width ) 
                {
                    lb.mixture[i][j][k].type = TYPE_I;
                    lb.mixture[i][j][k].u = u0;
                    lb.mixture[i][j][k].temp = units.temp(temp0);
                    lb.mixture[i][j][k].p = units.p(press0); 

                    lb.species[0][i][j][k].X = 0.10;
                    lb.species[1][i][j][k].X = 0.85;
                    lb.species[2][i][j][k].X = 0.001;
                    lb.species[3][i][j][k].X = 0.05;
                }
                // Right Inlet Boundary Condition
                if ( i==Nx-1 && j > NY/2 - jet_width && j < NY/2 + jet_width)
                {
                    lb.mixture[i][j][k].type = TYPE_I;
                    lb.mixture[i][j][k].u = -u0;
                    lb.mixture[i][j][k].temp = units.temp(temp0);
                    lb.mixture[i][j][k].p = units.p(press0); 

                    lb.species[0][i][j][k].X = 0.001;
                    lb.species[1][i][j][k].X = 0.90;
                    lb.species[2][i][j][k].X = 0.10;
                    lb.species[3][i][j][k].X = 0.001;
                }
                
                if (lb.mixture[i][j][k].type == TYPE_F)
                {
                    lb.mixture[i][j][k].u = 0.0;
                    lb.mixture[i][j][k].v = 0;
                    lb.mixture[i][j][k].w = 0.0;
                    lb.mixture[i][j][k].temp = units.temp(temp0);
                    lb.mixture[i][j][k].p = units.p(press0); 

                    lb.species[0][i][j][k].X = 0.05;
                    lb.species[1][i][j][k].X = 0.90;
                    lb.species[2][i][j][k].X = 0.05;
                    lb.species[3][i][j][k].X = 0.05;
                }
            }
        }
    }

    lb.run(1000000000,100);
}

#elif defined POINT_COMBUSTION_2D
void main_setup() // Perfectly stirred reactor ------------------------------------------------------------------------------------------
{
    int NX = 200; 
    int NY = 200; 
    int NZ = 1;
    
    // units.set_m_kg_s(NX, VEL0, RHO0, 0.1, si_len, si_u_max, si_rho, si_temp); // setting the conversion factor 
    units.set_m_kg_s(1e-5, 1e-9);

    double temp_amb = units.temp(300.0); // [K]
    double temp_ign = units.temp(3000.0); // [K]
    double sigma = 100;

    auto sol = Cantera::newSolution("h2o2.yaml", "ohmech");
    auto gas = sol->thermo();

    std::vector<std::string> species = gas->speciesNames();

    LBM lb(NX, NY, NZ, species);
    int Nx = lb.get_Nx(); int Ny = lb.get_Ny(); int Nz = lb.get_Nz();

    #pragma omp parallel for schedule(dynamic)
    for(int i = 0; i < Nx ; ++i)
    {
        for(int j = 0; j < Ny; ++j)
        {
            for(int k = 0; k < Nz; ++k)
            {
                if ( k==0 || k==Nz-1) // set periodic boundary condition
                {
                    lb.mixture[i][j][k].type = TYPE_P;
                }
                if (i==0 || i==Nx-1 || j==0 || j==Ny-1)
                {
                    lb.mixture[i][j][k].type = TYPE_A;
                }

                if (lb.mixture[i][j][k].type == TYPE_F)
                {
                    double X = i - Nx/2;
                    double Y = j - Ny/2;
                    lb.mixture[i][j][k].temp = temp_amb + (temp_ign - temp_amb) * exp(-((X*X + Y*Y) / (2.0 * sigma*sigma)));
                    lb.mixture[i][j][k].p = 1.0*units.p(Cantera::OneAtm); 
                
                    lb.mixture[i][j][k].u = 0.0;
                    lb.mixture[i][j][k].v = 0.0;
                    lb.mixture[i][j][k].w = 0.0; 
                    lb.species[gas->speciesIndex("H2")][i][j][k].X = 0.29586;     // H2 Mole Fraction
                    lb.species[gas->speciesIndex("O2")][i][j][k].X = 0.14793;     // O2 Mole Fraction
                    lb.species[gas->speciesIndex("N2")][i][j][k].X = 0.55621;     // N2 Mole Fraction
                }
            }
        }
    }

    lb.run(1000000,1000);
}

#elif defined FLAME_SPEED
void main_setup() // Perfectly stirred reactor ------------------------------------------------------------------------------------------
{
    units.set_m_kg_s(5.0e-6, 5.0e-10 ); //6.4e-6, 6.4e-10 (0.75-1.0)

    int NX = 4000; 
    int NY = 1; 
    int NZ = 1;
    
    auto sol = Cantera::newSolution("h2o2.yaml", "ohmech");
    auto gas = sol->thermo();

    std::vector<std::string> species = gas->speciesNames();

    LBM lb(NX, NY, NZ, species);
    int Nx = lb.get_Nx(); int Ny = lb.get_Ny(); int Nz = lb.get_Nz();

    double temp_u = units.temp(300.0); // [K]
    double smoothn = 0.02;
    double minim = 1e-13;
    double midpoint = 0.55*Nx;

    #pragma omp parallel for schedule(dynamic)
    for(int i = 0; i < Nx ; ++i)
    {
        for(int j = 0; j < Ny; ++j)
        {
            for(int k = 0; k < Nz; ++k)
            {
                if (j==0 || j==Ny-1 || k==0 || k==Nz-1) // set periodic boundary condition
                {
                    lb.mixture[i][j][k].type = TYPE_P;
                }
                if (i==0)
                {
                    lb.mixture[i][j][k].type = TYPE_A;
                }
                if (i==Nx-1)
                {
                    lb.mixture[i][j][k].type = TYPE_A;
                    lb.mixture[i][j][k].p = 1.0*units.p(Cantera::OneAtm);                 
                }

                if (lb.mixture[i][j][k].type == TYPE_F)
                {
                    lb.mixture[i][j][k].u = 0.0;
                    lb.mixture[i][j][k].v = 0.0;
                    lb.mixture[i][j][k].w = 0.0; 

                    // phi = 0.5
                    // double temp_b = units.temp(1646.2); // [K]
                    // lb.species[gas->speciesIndex("H2")][i][j][k].X = smooth(0.17361,  6.9363e-06, i, 0.8*Nx, smoothn);
                    // lb.species[gas->speciesIndex("H")][i][j][k].X = smooth(minim, 2.3043e-07, i, 0.8*Nx, smoothn);
                    // lb.species[gas->speciesIndex("O")][i][j][k].X = smooth(minim, 7.6421e-06, i, 0.8*Nx, smoothn);
                    // lb.species[gas->speciesIndex("O2")][i][j][k].X = smooth(0.17361, 0.094979 , i, 0.8*Nx, smoothn);
                    // lb.species[gas->speciesIndex("OH")][i][j][k].X = smooth(minim, 0.00028188, i, 0.8*Nx, smoothn);
                    // lb.species[gas->speciesIndex("H2O")][i][j][k].X = smooth(minim, 0.18995, i, 0.8*Nx, smoothn);
                    // lb.species[gas->speciesIndex("HO2")][i][j][k].X = smooth(minim, 3.8302e-07, i, 0.8*Nx, smoothn);
                    // lb.species[gas->speciesIndex("H2O2")][i][j][k].X = smooth(minim, 2.6553e-08, i, 0.8*Nx, smoothn);
                    // lb.species[gas->speciesIndex("N2")][i][j][k].X = smooth(0.65278, 0.71477, i, 0.8*Nx, smoothn);

                    // phi = 0.75
                    // double temp_b = units.temp(2101.8); // [K]
                    // lb.species[gas->speciesIndex("H2")][i][j][k].X      = smooth(0.23962, 0.00076445    , i, 0.8*Nx, smoothn);
                    // lb.species[gas->speciesIndex("H")][i][j][k].X       = smooth(minim  , 8.7175e-05    , i, 0.8*Nx, smoothn);
                    // lb.species[gas->speciesIndex("O")][i][j][k].X       = smooth(minim  , 0.00029721    , i, 0.8*Nx, smoothn);
                    // lb.species[gas->speciesIndex("O2")][i][j][k].X      = smooth(0.15974, 0.044642      , i, 0.8*Nx, smoothn);
                    // lb.species[gas->speciesIndex("OH")][i][j][k].X      = smooth(minim  , 0.0036599     , i, 0.8*Nx, smoothn);
                    // lb.species[gas->speciesIndex("H2O")][i][j][k].X     = smooth(minim  , 0.26918       , i, 0.8*Nx, smoothn);
                    // lb.species[gas->speciesIndex("HO2")][i][j][k].X     = smooth(minim  , 2.159e-06     , i, 0.8*Nx, smoothn);
                    // lb.species[gas->speciesIndex("H2O2")][i][j][k].X    = smooth(minim  , 1.4662e-07    , i, 0.8*Nx, smoothn);
                    // lb.species[gas->speciesIndex("N2")][i][j][k].X      = smooth(0.60064, 0.68136       , i, 0.8*Nx, smoothn);

                    // phi = 1.0
                    // double temp_b = units.temp(2387.637); // [K]
                    // lb.species[gas->speciesIndex("H2")][i][j][k].X = smooth(0.29586, 0.014565, i, 0.7*Nx, smoothn);
                    // lb.species[gas->speciesIndex("H")][i][j][k].X = smooth(minim, 0.00181, i, 0.7*Nx, smoothn);
                    // lb.species[gas->speciesIndex("O")][i][j][k].X = smooth(minim, 0.00060738, i, 0.7*Nx, smoothn);
                    // lb.species[gas->speciesIndex("O2")][i][j][k].X = smooth(0.14793, 0.0056094, i, 0.7*Nx, smoothn);
                    // lb.species[gas->speciesIndex("OH")][i][j][k].X = smooth(minim, 0.0072848, i, 0.7*Nx, smoothn);
                    // lb.species[gas->speciesIndex("H2O")][i][j][k].X = smooth(minim, 0.32437, i, 0.7*Nx, smoothn);
                    // lb.species[gas->speciesIndex("HO2")][i][j][k].X = smooth(minim, 1.2571e-06, i, 0.7*Nx, smoothn);
                    // lb.species[gas->speciesIndex("H2O2")][i][j][k].X = smooth(minim, 1.3391e-07, i, 0.7*Nx, smoothn);
                    // lb.species[gas->speciesIndex("N2")][i][j][k].X = smooth(0.55621, 0.64575, i, 0.7*Nx, smoothn);
                    
                    // phi = 1.25
                    // double temp_b = units.temp(2348.7); // [K]
                    // lb.species[gas->speciesIndex("H2")][i][j][k].X = smooth(0.34435, 0.079375, i, 0.6*Nx, smoothn);
                    // lb.species[gas->speciesIndex("H")][i][j][k].X = smooth(minim, 0.0034937, i, 0.6*Nx, smoothn);
                    // lb.species[gas->speciesIndex("O")][i][j][k].X = smooth(minim, 7.099e-05, i, 0.6*Nx, smoothn);
                    // lb.species[gas->speciesIndex("O2")][i][j][k].X = smooth(0.13774, 0.00011754, i, 0.6*Nx, smoothn);
                    // lb.species[gas->speciesIndex("OH")][i][j][k].X = smooth(minim, 0.0023875, i, 0.6*Nx, smoothn);
                    // lb.species[gas->speciesIndex("H2O")][i][j][k].X = smooth(minim,  0.3158, i, 0.6*Nx, smoothn);
                    // lb.species[gas->speciesIndex("HO2")][i][j][k].X = smooth(minim, 6.1036e-08, i, 0.6*Nx, smoothn);
                    // lb.species[gas->speciesIndex("H2O2")][i][j][k].X = smooth(minim, 1.7195e-08, i, 0.6*Nx, smoothn);
                    // lb.species[gas->speciesIndex("N2")][i][j][k].X = smooth(0.51791, 0.59876, i, 0.6*Nx, smoothn);

                    // phi = 1.5
                    // double temp_b = units.temp(2246.4); // [K]
                    // lb.species[gas->speciesIndex("H2")][i][j][k].X = smooth(0.3866, 0.14671 , i, 0.6*Nx, smoothn);
                    // lb.species[gas->speciesIndex("H")][i][j][k].X = smooth(minim, 0.0027918, i, 0.6*Nx, smoothn);
                    // lb.species[gas->speciesIndex("O")][i][j][k].X = smooth(minim, 1.0969e-05, i, 0.6*Nx, smoothn);
                    // lb.species[gas->speciesIndex("O2")][i][j][k].X = smooth(0.12887, 9.2513e-06, i, 0.6*Nx, smoothn);
                    // lb.species[gas->speciesIndex("OH")][i][j][k].X = smooth(minim, 0.00083597, i, 0.6*Nx, smoothn);
                    // lb.species[gas->speciesIndex("H2O")][i][j][k].X = smooth(minim,  0.29445, i, 0.6*Nx, smoothn);
                    // lb.species[gas->speciesIndex("HO2")][i][j][k].X = smooth(minim, 6.3987e-09, i, 0.6*Nx, smoothn);
                    // lb.species[gas->speciesIndex("H2O2")][i][j][k].X = smooth(minim, 3.4716e-09, i, 0.6*Nx, smoothn);
                    // lb.species[gas->speciesIndex("N2")][i][j][k].X = smooth(0.48454, 0.55519, i, 0.6*Nx, smoothn);
                    
                    // phi = 1.75
                    // double temp_b = units.temp(2149.6); // [K]
                    // lb.species[gas->speciesIndex("H2")][i][j][k].X = smooth(0.42373, 0.20559 , i, 0.8*Nx, smoothn);
                    // lb.species[gas->speciesIndex("H")][i][j][k].X = smooth(minim, 0.0019092, i, 0.8*Nx, smoothn);
                    // lb.species[gas->speciesIndex("O")][i][j][k].X = smooth(minim, 2.1507e-06, i, 0.8*Nx, smoothn);
                    // lb.species[gas->speciesIndex("O2")][i][j][k].X = smooth(0.12107, 1.2201e-06, i, 0.8*Nx, smoothn);
                    // lb.species[gas->speciesIndex("OH")][i][j][k].X = smooth(minim, 0.00032887, i, 0.8*Nx, smoothn);
                    // lb.species[gas->speciesIndex("H2O")][i][j][k].X = smooth(minim,  0.27484, i, 0.8*Nx, smoothn);
                    // lb.species[gas->speciesIndex("HO2")][i][j][k].X = smooth(minim, 9.7833e-10, i, 0.8*Nx, smoothn);
                    // lb.species[gas->speciesIndex("H2O2")][i][j][k].X = smooth(minim, 9.0087e-10, i, 0.8*Nx, smoothn);
                    // lb.species[gas->speciesIndex("N2")][i][j][k].X = smooth(0.45521, 0.51733, i, 0.8*Nx, smoothn);
                    
                    // phi = 2.0
                    // double temp_b = units.temp(2060.4); // [K]
                    // lb.species[gas->speciesIndex("H2")][i][j][k].X = smooth(0.45662, 0.25701 , i, 0.8*Nx, smoothn);
                    // lb.species[gas->speciesIndex("H")][i][j][k].X = smooth(minim, 0.0012315, i, 0.8*Nx, smoothn);
                    // lb.species[gas->speciesIndex("O")][i][j][k].X = smooth(minim, 4.7242e-07, i, 0.8*Nx, smoothn);
                    // lb.species[gas->speciesIndex("O2")][i][j][k].X = smooth(0.11416, 2.0277e-07, i, 0.8*Nx, smoothn);
                    // lb.species[gas->speciesIndex("OH")][i][j][k].X = smooth(minim, 0.00013708, i, 0.8*Nx, smoothn);
                    // lb.species[gas->speciesIndex("H2O")][i][j][k].X = smooth(minim,  0.25742, i, 0.8*Nx, smoothn);
                    // lb.species[gas->speciesIndex("HO2")][i][j][k].X = smooth(minim, 1.7805e-10, i, 0.8*Nx, smoothn);
                    // lb.species[gas->speciesIndex("H2O2")][i][j][k].X = smooth(minim, 2.6322e-10, i, 0.8*Nx, smoothn);
                    // lb.species[gas->speciesIndex("N2")][i][j][k].X = smooth(0.42922, 0.4842, i, 0.8*Nx, smoothn);

                    // phi = 2.25
                    double temp_b = units.temp(1978.6); // [K]
                    lb.species[gas->speciesIndex("H2")][i][j][k].X = smooth(0.48596, 0.30218, i, midpoint, smoothn);
                    lb.species[gas->speciesIndex("H")][i][j][k].X = smooth(minim, 0.00077213, i, midpoint, smoothn);
                    lb.species[gas->speciesIndex("O")][i][j][k].X = smooth(minim, 1.112e-07, i, midpoint, smoothn);
                    lb.species[gas->speciesIndex("O2")][i][j][k].X = smooth(0.10799, 3.8523e-08, i, midpoint, smoothn);
                    lb.species[gas->speciesIndex("OH")][i][j][k].X = smooth(minim, 5.9246e-05, i, midpoint, smoothn);
                    lb.species[gas->speciesIndex("H2O")][i][j][k].X = smooth(minim,  0.24197, i, midpoint, smoothn);
                    lb.species[gas->speciesIndex("HO2")][i][j][k].X = smooth(minim, 3.5935e-11, i, midpoint, smoothn);
                    lb.species[gas->speciesIndex("H2O2")][i][j][k].X = smooth(minim, 8.265e-11, i, midpoint, smoothn);
                    lb.species[gas->speciesIndex("N2")][i][j][k].X = smooth(0.40605, 0.45502, i, midpoint, smoothn);
                    
                    
                    lb.mixture[i][j][k].p = smooth(1.0*units.p(Cantera::OneAtm), 1.2*units.p(Cantera::OneAtm), i, midpoint, smoothn); ;//1.0*units.p(Cantera::OneAtm);                 
                    lb.mixture[i][j][k].temp = smooth(temp_u, temp_b, i, midpoint, smoothn); 



                }
            }
        }
    }

    lb.run(10000000000,1000);

    // LBM lb = read_restart("restart004000.dat");
    // lb.loop(10000000000, 1000);
}


#elif defined PREMIXED_LAMINAR_FLAME_2D
void main_setup() // Perfectly stirred reactor ------------------------------------------------------------------------------------------
{
    int NX = 200; 
    int NY = 200; 
    int NZ = 1;
    
    // units.set_m_kg_s(NX, VEL0, RHO0, 0.1, si_len, si_u_max, si_rho, si_temp); // setting the conversion factor 
    units.set_m_kg_s(1e-5, 1e-9);

    double temp_amb = units.temp(300.0); // [K]
    double temp_ign = units.temp(3000.0); // [K]

    auto sol = Cantera::newSolution("h2o2.yaml", "ohmech");
    auto gas = sol->thermo();

    std::vector<std::string> species = gas->speciesNames();

    LBM lb(NX, NY, NZ, species);
    int Nx = lb.get_Nx(); int Ny = lb.get_Ny(); int Nz = lb.get_Nz();

    #pragma omp parallel for schedule(dynamic)
    for(int i = 0; i < Nx ; ++i)
    {
        for(int j = 0; j < Ny; ++j)
        {
            for(int k = 0; k < Nz; ++k)
            {
                if (k==0 || k==Nz-1) // set periodic boundary condition
                {
                    lb.mixture[i][j][k].type = TYPE_P;
                }
                if (i==0 || i==Nx-1 || j==0 || j==Ny-1)
                {
                    lb.mixture[i][j][k].type = TYPE_A;
                }
                if (i==0 && j > 0.4*Ny && j < 0.6*Ny)
                {
                    lb.mixture[i][j][k].type = TYPE_I;
                    lb.mixture[i][j][k].temp = temp_amb;
                    lb.mixture[i][j][k].p = 1.0*units.p(Cantera::OneAtm); 
                    
                    lb.mixture[i][j][k].u = 0.01;
                    lb.species[gas->speciesIndex("Ar")][i][j][k].X = 1.0;     // H2 Mole Fraction
                }
                if (i==Nx-1)
                {
                    lb.mixture[i][j][k].type = TYPE_O;
                }

                if (lb.mixture[i][j][k].type == TYPE_F)
                {
                    lb.mixture[i][j][k].temp = temp_amb;
                    lb.mixture[i][j][k].p = 1.0*units.p(Cantera::OneAtm); 
                
                    lb.mixture[i][j][k].u = 0.0;
                    lb.mixture[i][j][k].v = 0.0;
                    lb.mixture[i][j][k].w = 0.0; 
                    lb.species[gas->speciesIndex("Ar")][i][j][k].X = 1.0;     // H2 Mole Fraction
                }
            }
        }
    }

    lb.run(10000000,100);
}

#elif defined CONVECTED_VORTEX
void main_setup()
{
    units.set_m_kg_s(1e-5, 1e-8);

    double NX = 100;
    double NY = 100;
    double NZ = 1;

    auto sol = Cantera::newSolution("h2o2.yaml", "ohmech");
    // auto sol = Cantera::newSolution("./src/reaction-mech/propane_mech.yaml");
    auto gas = sol->thermo();
    auto trans = sol->transport();
    std::vector <double> X (gas->nSpecies());
    X[gas->speciesIndex("N2")] = 0.9;
    X[gas->speciesIndex("H2")] = 0.1;
    gas->setState_TP(300, Cantera::OneAtm);
    gas->setMoleFractions(&X[0]);

    std::vector<std::string> species = { "N2", "H2" };
    LBM lb(NX, NY, NZ, species);
    int Nx = lb.get_Nx(); int Ny = lb.get_Ny(); int Nz = lb.get_Nz();


    // double Mach_inf = 0.1;
    // double u_inf = Mach_inf * units.u(gas->soundSpeed());
    // double b = 0.005;
    double core_rad = 20;
    double x_center = NX/2.0;
    double y_center = NY/2.0;
    
    
    std::cout << "nu (lu) : " << units.nu(trans->viscosity() / gas->density()) << std::endl;
    std::cout << "gas const : " << units.cp(Cantera::GasConstant/gas->meanMolecularWeight()) << std::endl;
    std::cout << "temperature : " << units.temp(gas->temperature()) << std::endl;
    std::cout << "gamma : " << gas->cp_mass()/gas->cv_mass() << std::endl;
    std::cout << "sound speed : " << units.u(gas->soundSpeed()) << " " << gas->soundSpeed() <<  std::endl;
    std::cout << "RT : " << units.cp(Cantera::GasConstant/gas->meanMolecularWeight())*units.temp(gas->temperature()) << std::endl;
    
    #pragma omp parallel for schedule(dynamic)
    for(int i = 0; i < Nx ; ++i){
        for(int j = 0; j < Ny; ++j){
            for(int k = 0; k < Nz; ++k){
                // set periodic boundary condition
                if ( k==0 || k==Nz-1 ) 
                    lb.mixture[i][j][k].type = TYPE_P;
                
                // if ( i==0 || i ==Nx-1 )
                //     lb.mixture[i][j][k].type = TYPE_O_C;

                // if ( j==0 || j==Ny-1)
                //     lb.mixture[i][j][k].type = TYPE_O_C;    

                if ( j==Ny-1 )
                    lb.mixture[i][j][k].type = TYPE_O;  

                if ( j==0 )
                    lb.mixture[i][j][k].type = TYPE_O;

                if ( i==Nx-1 )
                    lb.mixture[i][j][k].type = TYPE_O;

                if ( i==0 ) {
                    lb.mixture[i][j][k].type = TYPE_O;
                    lb.mixture[i][j][k].u       = units.u(10e-2);
                    lb.mixture[i][j][k].v       = 0.0;
                    lb.mixture[i][j][k].w       = 0.0;
                    lb.mixture[i][j][k].temp = units.temp(300);
                    lb.mixture[i][j][k].p = units.p(Cantera::OneAtm);

                    lb.species[0][i][j][k].X = 0.1;
                    lb.species[1][i][j][k].X = 0.9; 

                }

                                
                                
                // if ( i==0 && (j < 0.3*Nx || j > 0.7*Nx)){
                //     lb.mixture[i][j][k].type = TYPE_I_C;
                    
                //     lb.mixture[i][j][k].u = units.u(0);
                //     // lb.mixture[i][j][k].v = 0.0;
                //     // lb.mixture[i][j][k].w = 0.0;
                //     lb.mixture[i][j][k].temp = units.temp(300);
                //     lb.mixture[i][j][k].p = units.p(Cantera::OneAtm);  //  *
                    
                //     lb.species[0][i][j][k].X = 0.5;   
                //     // lb.species[1][i][j][k].X = 0.5;
                // }

                                
                if (lb.mixture[i][j][k].type == TYPE_F ||  lb.mixture[i][j][k].type == TYPE_O  ||  lb.mixture[i][j][k].type == TYPE_O_C  )
                {

                    // double r_dist_2 = pow(i-x_center, 2) + pow(j-y_center, 2);
                    // double core_rad_2 = pow(core_rad, 2);

                    // lb.mixture[i][j][k].u = u_inf - b/(2*M_PI)*exp(0.5*(1.0-r_dist_2/core_rad_2))*(j-y_center);
                    // lb.mixture[i][j][k].v = b/(2*M_PI)*exp(0.5*(1.0-r_dist_2/core_rad_2))*(i-x_center);
                    // lb.mixture[i][j][k].w = 0.0;
                    // lb.mixture[i][j][k].rho = pow( 1.0-(gamma-1.0)*b*b/(8.0*gamma*M_PI*M_PI) * exp(1.0-r_dist_2/core_rad_2) , 1.0/(gamma-1.0));
                    // lb.mixture[i][j][k].temp = temp0;
                    // lb.mixture[i][j][k].p = lb.mixture[i][j][k].rho*gas_const*lb.mixture[i][j][k].temp; 

                    // -----------------------------------------------------------------------------------------------------------------------------

                    double r_dist_2 = pow(i-x_center, 2) + pow(j-y_center, 2);
                    double core_rad_2 = pow(core_rad, 2);

                    lb.mixture[i][j][k].u = 0.0;//0.0018;
                    lb.mixture[i][j][k].v = 0.0;
                    lb.mixture[i][j][k].w = 0.0;
                    lb.mixture[i][j][k].temp = units.temp(300);
                    lb.mixture[i][j][k].p = units.p(Cantera::OneAtm) * (1.0 + 1.0e-1 * exp(-r_dist_2/core_rad_2));//  ;  //    
                    
                    if (i < 0.1*Nx){
                        lb.species[0][i][j][k].X = 0.9;
                        lb.species[1][i][j][k].X = 0.1;   
                    }
                    else{
                        lb.species[0][i][j][k].X = 0.9;
                        lb.species[1][i][j][k].X = 0.1;   
                    }
                }
            }
        }
    }
    lb.run(5000000000,1);
}

#elif defined COFLOW_FLAME
void main_setup()
{
    // units.set_m_kg_s(3e-5, 2e-8);
    // units.set_m_kg_s(3e-5, 2e-8); //3e-5, 2e-8
    units.set_m_kg_s(3e-5, 8e-9);

    double NX = 187;
    double NY = 334;
    double NZ = 1;

    auto sol = Cantera::newSolution("./src/reaction-mech/CH4_2S.yaml", "CH4_BFER_multi");
    auto gas = sol->thermo();
    auto trans = sol->transport();
    std::vector <double> Y (gas->nSpecies());
    Y[gas->speciesIndex("CH4")] = 0.1624900877705;
    Y[gas->speciesIndex("N2")] = 1.0 - 0.1624900877705;
    gas->setMoleFractions(&Y[0]);
    gas->setState_TP(950, Cantera::OneAtm);
    std::vector <double> X (gas->nSpecies());
    gas->getMassFractions(&X[0]);
    // std::cout << X[gas->speciesIndex("CH4")] << " | " << X[gas->speciesIndex("N2")] << std::endl;

    std::vector<std::string> species = gas->speciesNames(); //{"CH4"};//
    LBM lb(NX, NY, NZ, species);
    int Nx = lb.get_Nx(); int Ny = lb.get_Ny(); int Nz = lb.get_Nz();

    // double Mach_inf = 0.1;
    // double u_inf = Mach_inf * units.u(gas->soundSpeed());
    // double b = 0.005;
    
    std::cout << "nu (lu) : " << units.nu(trans->viscosity() / gas->density()) << std::endl;
    std::cout << "gas const : " << units.cp(Cantera::GasConstant/gas->meanMolecularWeight()) << std::endl;
    std::cout << "temperature : " << units.temp(gas->temperature()) << std::endl;
    std::cout << "gamma : " << gas->cp_mass()/gas->cv_mass() << std::endl;
    std::cout << "sound speed : " << units.u(gas->soundSpeed()) << " " << gas->soundSpeed() <<  std::endl;
    std::cout << "RT : " << units.cp(Cantera::GasConstant/gas->meanMolecularWeight())*units.temp(gas->temperature()) << std::endl;
   
    double eta = 1e-1;
    
    #pragma omp parallel for schedule(dynamic)
    for(int i = 0; i < Nx ; ++i){
        for(int j = 0; j < Ny; ++j){
            for(int k = 0; k < Nz; ++k){
                // set periodic boundary condition
                if ( k==0 || k==Nz-1 ) 
                    lb.mixture[i][j][k].type = TYPE_P;

                
                if ( i==0 )
                    lb.mixture[i][j][k].type = TYPE_P;

                if ( i==Nx-1 )
                    lb.mixture[i][j][k].type = TYPE_P;

                if ( j==Ny-1 )
                    lb.mixture[i][j][k].type = TYPE_O;  

                if ( j==0 ) {
                    lb.mixture[i][j][k].type = TYPE_I;
                    lb.mixture[i][j][k].u       = 0.0;
                    lb.mixture[i][j][k].w       = 0.0;
                    lb.mixture[i][j][k].temp = units.temp(950);
                    lb.mixture[i][j][k].p = units.p(Cantera::OneAtm);

                    if(i > 0.366*Nx && i < 0.634*Nx){ // fuel in central inlet
                        lb.mixture[i][j][k].v       = units.u(0.8);
                        lb.species[gas->speciesIndex("CH4")][i][j][k].X = 0.1624900877705;
                        lb.species[gas->speciesIndex("N2")][i][j][k].X = 1.0 - 0.1624900877705;

                        // lb.species[0][i][j][k].X = 0.0598242;
                    }else{
                        lb.mixture[i][j][k].v       = units.u(0.5);
                        // lb.species[gas->speciesIndex("O2")][i][j][k].X = 0.2017291647356;
                        lb.species[gas->speciesIndex("N2")][i][j][k].X = 1.0 - 0.2017291647356;

                        // lb.species[0][i][j][k].X = 0.0598242;

                    }




                }

                     
                if (lb.mixture[i][j][k].type == TYPE_F ||  lb.mixture[i][j][k].type == TYPE_O  ||  lb.mixture[i][j][k].type == TYPE_O_C )
                {

                    // double r_dist_2 = pow(i-x_center, 2) + pow(j-y_center, 2);
                    // double core_rad_2 = pow(core_rad, 2);

                    // lb.mixture[i][j][k].u = u_inf - b/(2*M_PI)*exp(0.5*(1.0-r_dist_2/core_rad_2))*(j-y_center);
                    // lb.mixture[i][j][k].v = b/(2*M_PI)*exp(0.5*(1.0-r_dist_2/core_rad_2))*(i-x_center);
                    // lb.mixture[i][j][k].w = 0.0;
                    // lb.mixture[i][j][k].rho = pow( 1.0-(gamma-1.0)*b*b/(8.0*gamma*M_PI*M_PI) * exp(1.0-r_dist_2/core_rad_2) , 1.0/(gamma-1.0));
                    // lb.mixture[i][j][k].temp = temp0;
                    // lb.mixture[i][j][k].p = lb.mixture[i][j][k].rho*gas_const*lb.mixture[i][j][k].temp; 

                    // -----------------------------------------------------------------------------------------------------------------------------

                    lb.mixture[i][j][k].u = 0.0;
                    lb.mixture[i][j][k].v = units.u(0.5);
                    lb.mixture[i][j][k].w = 0.0;
                    lb.mixture[i][j][k].temp = units.temp(950);
                    lb.mixture[i][j][k].p = units.p(Cantera::OneAtm)  ;  //    * (1.0 + 1.0e-1 * exp(-r_dist_2/core_rad_2))

                    // lb.species[gas->speciesIndex("O2")][i][j][k].X = 0.2017291647356 - eta/2.0;
                    lb.species[gas->speciesIndex("N2")][i][j][k].X = 1.0 - 0.2017291647356 - eta/2.0;
                    lb.species[gas->speciesIndex("CH4")][i][j][k].X = eta;

                    // lb.species[0][i][j][k].X = 0.0598242;

                                       
                }
            }
        }
    }
    lb.run(6000,1);

    // LBM lb =read_restart("restart020000.dat");
    // lb.loop(5000000000,1);
}

#elif defined OPPOSED_JET_FLAME
void main_setup()
{
    units.set_m_kg_s(1e-5, 0.8e-9);

    double NX = 300;
    double NY = 180;
    double NZ = 1;

    auto sol = Cantera::newSolution("h2o2.yaml", "ohmech");
    auto gas = sol->thermo();
    auto trans = sol->transport();
    std::vector <double> X (gas->nSpecies());
    X[gas->speciesIndex("H2")] = 0.070312;
    X[gas->speciesIndex("O2")] = 0.19531;
    // X[gas->speciesIndex("H2O")] = 0.061;
    // X[gas->speciesIndex("CO2")] = 0.111;
    X[gas->speciesIndex("N2")] = 0.73438;
    gas->setMoleFractions(&X[0]);
    gas->setState_TP(300, Cantera::OneAtm);
    gas->getMoleFractions(&X[0]);
    std::cout << X[gas->speciesIndex("H2")] << " | " << X[gas->speciesIndex("O2")] << " | " << X[gas->speciesIndex("H2O")] << " | "  << X[gas->speciesIndex("CO2")] << " | "  << X[gas->speciesIndex("N2")] << std::endl;

    std::vector<std::string> species = gas->speciesNames(); //{"CH4"};//
    LBM lb(NX, NY, NZ, species);
    int Nx = lb.get_Nx(); int Ny = lb.get_Ny(); int Nz = lb.get_Nz();

    // double Mach_inf = 0.1;
    // double u_inf = Mach_inf * units.u(gas->soundSpeed());
    // double b = 0.005;
    
    std::cout << "nu (lu) : " << units.nu(trans->viscosity() / gas->density()) << std::endl;
    std::cout << "gas const : " << units.cp(Cantera::GasConstant/gas->meanMolecularWeight()) << std::endl;
    std::cout << "temperature : " << units.temp(gas->temperature()) << std::endl;
    std::cout << "gamma : " << gas->cp_mass()/gas->cv_mass() << std::endl;
    std::cout << "sound speed : " << units.u(gas->soundSpeed()) << " " << gas->soundSpeed() <<  std::endl;
    std::cout << "RT : " << units.cp(Cantera::GasConstant/gas->meanMolecularWeight())*units.temp(gas->temperature()) << std::endl;

    double smoothn = 0.1;
    double minim = 1e-40;

    #pragma omp parallel for schedule(dynamic)
    for(int i = 0; i < Nx ; ++i){
        for(int j = 0; j < Ny; ++j){
            for(int k = 0; k < Nz; ++k){
                // set periodic boundary condition
                if ( k==0 || k==Nz-1 ) 
                    lb.mixture[i][j][k].type = TYPE_P;

                if ( i==0 )
                    lb.mixture[i][j][k].type = TYPE_A;
                if ( j==0 )
                    lb.mixture[i][j][k].type = TYPE_A;

                if ( i==Nx-1 )
                    lb.mixture[i][j][k].type = TYPE_O;
                if ( j==Ny-1 ){
                    lb.mixture[i][j][k].type = TYPE_I;
                    lb.mixture[i][j][k].v = -units.u(0.01);
                    lb.mixture[i][j][k].p = units.p(Cantera::OneAtm)  ;
                    lb.mixture[i][j][k].temp = units.temp(300);
                    // lb.species[0][i][j][k].X = 1.0;

                    // lb.mixture[i][j][k].temp = units.temp(smooth(870.15, 300.0, j, 0.4*Ny, smoothn));
                    // lb.species[gas->speciesIndex("H2")][i][j][k].X      = smooth(1.8428e-13, 0.070312, j, 0.25*Ny, smoothn);
                    // lb.species[gas->speciesIndex("H")][i][j][k].X       = smooth(1.0000e-40, minim  , j, 0.25*Ny, smoothn);
                    // lb.species[gas->speciesIndex("O")][i][j][k].X       = smooth(6.8999e-13, minim  , j, 0.25*Ny, smoothn);
                    // lb.species[gas->speciesIndex("O2")][i][j][k].X      = smooth(0.16599   , 0.19531, j, 0.25*Ny, smoothn);
                    // lb.species[gas->speciesIndex("OH")][i][j][k].X      = smooth(5.016e-09 , minim  , j, 0.25*Ny, smoothn);
                    // lb.species[gas->speciesIndex("H2O")][i][j][k].X     = smooth(0.072874  , minim  , j, 0.25*Ny, smoothn);
                    // lb.species[gas->speciesIndex("HO2")][i][j][k].X     = smooth(6.1495e-11, minim  , j, 0.25*Ny, smoothn);
                    // lb.species[gas->speciesIndex("H2O2")][i][j][k].X    = smooth(1.2509e-11, minim  , j, 0.25*Ny, smoothn);
                    lb.species[gas->speciesIndex("N2")][i][j][k].X      = smooth(0.76113   , 0.73438, j, 0.25*Ny, smoothn);

                }

                // if (i > 0.7*Nx && j > 0.8*Ny)
                //     lb.mixture[i][j][k].type = TYPE_A;
                
                     
                if (lb.mixture[i][j][k].type == TYPE_F ||  lb.mixture[i][j][k].type == TYPE_O )
                {

                    
                    lb.mixture[i][j][k].u = 0.0;//smooth(units.u(0.01), 0.0, j, 0.85*Ny, smoothn);
                    lb.mixture[i][j][k].v = 0.0;//smooth(0.0, -units.u(0.01)    , j, 0.85*Ny, smoothn);
                    lb.mixture[i][j][k].w = 0.0;
                    lb.mixture[i][j][k].p = units.p(Cantera::OneAtm)  ;  //    * (1.0 + 1.0e-1 * exp(-r_dist_2/core_rad_2))
                    // lb.mixture[i][j][k].temp = units.temp(300);
                    
                    // lb.species[gas->speciesIndex("C3H8")][i][j][k].X    = smooth(1e-5      , 0.0553053  , j, 0.25*Ny, smoothn);       
                    // lb.species[gas->speciesIndex("O2")][i][j][k].X      = smooth(0.0964249 , 0.242956   , j, 0.25*Ny, smoothn);
                    // lb.species[gas->speciesIndex("H2O")][i][j][k].X     = smooth(0.0372084 , 1E-5       , j, 0.25*Ny, smoothn);
                    // lb.species[gas->speciesIndex("CO2")][i][j][k].X     = smooth(0.165402  , 1E-5       , j, 0.25*Ny, smoothn);
                    // lb.species[gas->speciesIndex("N2")][i][j][k].X      = smooth(0.700965  ,  0.701739  , j, 0.25*Ny, smoothn);

                    lb.mixture[i][j][k].temp = units.temp(smooth(870.15, 300.0, j, 0.4*Ny, smoothn));
                    // lb.species[gas->speciesIndex("H2")][i][j][k].X      = smooth(1.8428e-13, 0.070312, j, 0.25*Ny, smoothn);
                    // lb.species[gas->speciesIndex("H")][i][j][k].X       = smooth(1.0000e-40, minim  , j, 0.25*Ny, smoothn);
                    // lb.species[gas->speciesIndex("O")][i][j][k].X       = smooth(6.8999e-13, minim  , j, 0.25*Ny, smoothn);
                    // lb.species[gas->speciesIndex("O2")][i][j][k].X      = smooth(0.16599   , 0.19531, j, 0.25*Ny, smoothn);
                    // lb.species[gas->speciesIndex("OH")][i][j][k].X      = smooth(5.016e-09 , minim  , j, 0.25*Ny, smoothn);
                    // lb.species[gas->speciesIndex("H2O")][i][j][k].X     = smooth(0.072874  , minim  , j, 0.25*Ny, smoothn);
                    // lb.species[gas->speciesIndex("HO2")][i][j][k].X     = smooth(6.1495e-11, minim  , j, 0.25*Ny, smoothn);
                    // lb.species[gas->speciesIndex("H2O2")][i][j][k].X    = smooth(1.2509e-11, minim  , j, 0.25*Ny, smoothn);
                    lb.species[gas->speciesIndex("N2")][i][j][k].X      = smooth(0.76113   , 0.73438, j, 0.25*Ny, smoothn);

                                       
                }
            }
        }
    }
    lb.run(5000000000,100);

    // LBM lb =read_restart("restart020000.dat");
    // lb.loop(5000000000,1);
}

#elif defined PREMIXED_JET_FLAME
void main_setup()
{
    units.set_m_kg_s(1e-5, 1e-8);

    double NX = 300;
    double NY = 180;
    double NZ = 1;

    auto sol = Cantera::newSolution("./src/reaction-mech/CH4_2S.yaml");
    auto gas = sol->thermo();
    auto trans = sol->transport();
    std::vector <double> Y (gas->nSpecies());
    Y[gas->speciesIndex("CH4")] = 0.037;
    // Y[gas->speciesIndex("O2")] = 0.089;
    // Y[gas->speciesIndex("H2O")] = 0.061;
    // Y[gas->speciesIndex("CO2")] = 0.111;
    Y[gas->speciesIndex("N2")] = 0.739;
    gas->setMoleFractions(&Y[0]);
    gas->setState_TP(300, Cantera::OneAtm);
    std::vector <double> X (gas->nSpecies());
    gas->getMassFractions(&X[0]);
    std::cout << X[gas->speciesIndex("CH4")] << " | " << X[gas->speciesIndex("O2")] << " | " << X[gas->speciesIndex("H2O")] << " | "  << X[gas->speciesIndex("CO2")] << " | "  << X[gas->speciesIndex("N2")] << std::endl;

    std::vector<std::string> species = gas->speciesNames(); //{"CH4"};//
    LBM lb(NX, NY, NZ, species);
    int Nx = lb.get_Nx(); int Ny = lb.get_Ny(); int Nz = lb.get_Nz();

    // double Mach_inf = 0.1;
    // double u_inf = Mach_inf * units.u(gas->soundSpeed());
    // double b = 0.005;
    
    std::cout << "nu (lu) : " << units.nu(trans->viscosity() / gas->density()) << std::endl;
    std::cout << "gas const : " << units.cp(Cantera::GasConstant/gas->meanMolecularWeight()) << std::endl;
    std::cout << "temperature : " << units.temp(gas->temperature()) << std::endl;
    std::cout << "gamma : " << gas->cp_mass()/gas->cv_mass() << std::endl;
    std::cout << "sound speed : " << units.u(gas->soundSpeed()) << " " << gas->soundSpeed() <<  std::endl;
    std::cout << "RT : " << units.cp(Cantera::GasConstant/gas->meanMolecularWeight())*units.temp(gas->temperature()) << std::endl;

    double smoothn = 0.07;

    #pragma omp parallel for schedule(dynamic)
    for(int i = 0; i < Nx ; ++i){
        for(int j = 0; j < Ny; ++j){
            for(int k = 0; k < Nz; ++k){
                // set periodic boundary condition
                if ( k==0 || k==Nz-1 ) 
                    lb.mixture[i][j][k].type = TYPE_P;

                if ( i==0 )
                    lb.mixture[i][j][k].type = TYPE_A;
                if ( j==0 )
                    lb.mixture[i][j][k].type = TYPE_A;

                if ( i==Nx-1 )
                    lb.mixture[i][j][k].type = TYPE_O;
                if ( j==Ny-1 ){
                    lb.mixture[i][j][k].type = TYPE_I;
                    lb.mixture[i][j][k].v = -units.u(0.2);
                    lb.mixture[i][j][k].p = units.p(Cantera::OneAtm)  ;
                    lb.mixture[i][j][k].temp = units.temp(300);
                    // lb.species[0][i][j][k].X = 1.0;

                }

                if (i > 0.7*Nx && j > 0.8*Ny)
                    lb.mixture[i][j][k].type = TYPE_A;
                
                     
                if (lb.mixture[i][j][k].type == TYPE_F ||  lb.mixture[i][j][k].type == TYPE_O || lb.mixture[i][j][k].type == TYPE_I )
                {

                    
                    lb.mixture[i][j][k].u = smooth(units.u(0.2), 0.0, j, 0.85*Ny, smoothn);
                    lb.mixture[i][j][k].v = smooth(0.0, -units.u(0.2)    , j, 0.85*Ny, smoothn);
                    lb.mixture[i][j][k].w = 0.0;
                    lb.mixture[i][j][k].p = units.p(Cantera::OneAtm)  ;  //    * (1.0 + 1.0e-1 * exp(-r_dist_2/core_rad_2))

                    lb.mixture[i][j][k].temp = units.temp(smooth(300.0, 300.0, j, 0.4*Ny, smoothn)); // 1970
                    
                    lb.species[gas->speciesIndex("C3H8")][i][j][k].X    = smooth(1e-5      , 0.0553053  , j, 0.25*Ny, smoothn);       
                    lb.species[gas->speciesIndex("O2")][i][j][k].X      = smooth(0.0964249 , 0.242956   , j, 0.25*Ny, smoothn);
                    lb.species[gas->speciesIndex("H2O")][i][j][k].X     = smooth(0.0372084 , 1E-5       , j, 0.25*Ny, smoothn);
                    lb.species[gas->speciesIndex("CO2")][i][j][k].X     = smooth(0.165402  , 1E-5       , j, 0.25*Ny, smoothn);
                    lb.species[gas->speciesIndex("N2")][i][j][k].X      = smooth(0.700965  ,  0.701739  , j, 0.25*Ny, smoothn);

                                       
                }
            }
        }
    }
    lb.run(5000000000,100);

    // LBM lb =read_restart("restart020000.dat");
    // lb.loop(5000000000,1);
}

#elif defined MICROCHANNEL_FLAME
void main_setup() // Perfectly stirred reactor ------------------------------------------------------------------------------------------
{
    // units.set_m_kg_s(2.778e-5, 2.778e-8/5 ); //6.4e-6, 6.4e-10 (0.75-1.0)
    // units.set_m_kg_s(2.5e-5, 0.9e-8 );
    units.set_m_kg_s(20e-6, 1e-9 );

    int NX = 300; 
    int NY = 30; 
    int NZ = 1;
    
    auto sol = Cantera::newSolution("h2o2.yaml", "ohmech");
    auto gas = sol->thermo();

    std::vector<std::string> species = gas->speciesNames(); //{"N2"};//

    LBM lb(NX, NY, NZ, species);
    int Nx = lb.get_Nx(); int Ny = lb.get_Ny(); int Nz = lb.get_Nz();

    double temp_u = units.temp(300.0); // [K]
    double temp_b = units.temp(960.0); // [K]960
    double smoothn = 0.4;
    double minim = 1E-13;
    double midpoint = 0.55*Nx/20;

    auto trans = sol->transport();
    std::vector <double> X (gas->nSpecies());
    // X[0] = 0.17361;
    // X[gas->speciesIndex("O2")] = 0.17361;
    X[gas->speciesIndex("H2O")] = 0.061;
    // X[gas->speciesIndex("CO2")] = 0.111;
    // X[gas->speciesIndex("N2")] = 0.65278;
    gas->setMoleFractions(&X[0]);
    gas->setState_TP(960, Cantera::OneAtm);

    std::cout << "nu (lu) : " << units.nu(trans->viscosity() / gas->density()) << std::endl;
    std::cout << "gas const : " << units.cp(Cantera::GasConstant/gas->meanMolecularWeight()) << std::endl;
    std::cout << "temperature : " << units.temp(gas->temperature()) << std::endl;
    std::cout << "gamma : " << gas->cp_mass()/gas->cv_mass() << std::endl;
    std::cout << "sound speed : " << units.u(gas->soundSpeed()) << " " << gas->soundSpeed() <<  std::endl;
    std::cout << "RT : " << units.cp(Cantera::GasConstant/gas->meanMolecularWeight())*units.temp(gas->temperature()) << std::endl;

    #pragma omp parallel for schedule(dynamic)
    for(int i = 0; i < Nx ; ++i)
    {
        for(int j = 0; j < Ny; ++j)
        {
            for(int k = 0; k < Nz; ++k)
            {
                if ( k==0 || k==Nz-1) // set periodic boundary condition
                {
                    lb.mixture[i][j][k].type = TYPE_P;
                }

                if (i==0)
                {
                    lb.mixture[i][j][k].type = TYPE_I;
                    lb.mixture[i][j][k].p = 1.0*units.p(Cantera::OneAtm);  
                }
                if (i==Nx-1)
                {
                    lb.mixture[i][j][k].type = TYPE_O;
                    lb.mixture[i][j][k].p = 1.0*units.p(Cantera::OneAtm);                 
                }
                if (j==0 || j==Ny-1 )
                {
                    lb.mixture[i][j][k].type = TYPE_S;
                    lb.mixture[i][j][k].temp = smooth(temp_u, temp_b, i, midpoint, 0.7);
                }

                if (lb.mixture[i][j][k].type == TYPE_F || lb.mixture[i][j][k].type == TYPE_O  || lb.mixture[i][j][k].type == TYPE_O_C || lb.mixture[i][j][k].type == TYPE_I )
                {
                    lb.mixture[i][j][k].u = units.u(50e-2);
                    lb.mixture[i][j][k].v = 0.0;
                    lb.mixture[i][j][k].w = 0.0; 

                    // phi = 0.5
                    lb.species[gas->speciesIndex("H2")][i][j][k].X  = 0.17361;   
                    lb.species[gas->speciesIndex("H")][i][j][k].X   = minim  ;  
                    lb.species[gas->speciesIndex("O")][i][j][k].X   = minim  ;  
                    lb.species[gas->speciesIndex("O2")][i][j][k].X  = 0.17361;   
                    lb.species[gas->speciesIndex("OH")][i][j][k].X  = minim  ;  
                    lb.species[gas->speciesIndex("H2O")][i][j][k].X = minim  ;   
                    lb.species[gas->speciesIndex("HO2")][i][j][k].X = minim  ;  
                    lb.species[gas->speciesIndex("H2O2")][i][j][k].X= minim  ;  
                    lb.species[gas->speciesIndex("N2")][i][j][k].X  = 0.65278;   

                    // lb.species[0][i][j][k].X  = 0.65278;   

                    
                    lb.mixture[i][j][k].p = 1.0*units.p(Cantera::OneAtm);                 
                    // lb.mixture[i][j][k].temp = temp_b; 
                    lb.mixture[i][j][k].temp = smooth(temp_u, temp_b, i, midpoint, 0.7); 



                }
            }
        }
    }

    lb.run(1000000000,100);

    // LBM lb = read_restart("restart840000.dat");
    // lb.loop(840001, 1);
}

#elif defined CIRCULAR_FLAME
void main_setup() // Perfectly stirred reactor ------------------------------------------------------------------------------------------
{
    units.set_m_kg_s(2e-6, 1e-10 );

    int NX = 800; 
    int NY = 800; 
    int NZ = 1;
    
    auto sol = Cantera::newSolution("h2o2.yaml", "ohmech");
    auto gas = sol->thermo();

    std::vector<std::string> species = gas->speciesNames();

    LBM lb(NX, NY, NZ, species);
    int Nx = lb.get_Nx(); int Ny = lb.get_Ny(); int Nz = lb.get_Nz();

    auto trans = sol->transport();
    std::vector <double> X (gas->nSpecies());
    X[gas->speciesIndex("H2")] = 0.17361;
    // X[gas->speciesIndex("O2")] = 0.17361;
    // X[gas->speciesIndex("H2O")] = 0.061;
    // X[gas->speciesIndex("CO2")] = 0.111;
    // X[gas->speciesIndex("N2")] = 0.65278;
    gas->setMoleFractions(&X[0]);
    gas->setState_TP(1844.3, 5.0*Cantera::OneAtm);

    std::cout << "nu (lu) : " << units.nu(trans->viscosity() / gas->density()) << std::endl;
    std::cout << "gas const : " << units.cp(Cantera::GasConstant/gas->meanMolecularWeight()) << std::endl;
    std::cout << "temperature : " << units.temp(gas->temperature()) << std::endl;
    std::cout << "gamma : " << gas->cp_mass()/gas->cv_mass() << std::endl;
    std::cout << "sound speed : " << units.u(gas->soundSpeed()) << " " << gas->soundSpeed() <<  std::endl;
    std::cout << "RT : " << units.cp(Cantera::GasConstant/gas->meanMolecularWeight())*units.temp(gas->temperature()) << std::endl;
    
    double x_center = 0.0;
    double y_center = 0.0;
    double core_rad = 0.0;
    double smoothness = 0.04;

    double minim = 1E-13;

    #pragma omp parallel for schedule(dynamic)
    for(int i = 0; i < Nx ; ++i)
    {
        for(int j = 0; j < Ny; ++j)
        {
            for(int k = 0; k < Nz; ++k)
            {
                if ( k==0 || k==Nz-1) // set periodic boundary condition
                {
                    lb.mixture[i][j][k].type = TYPE_P;
                }

                if (i==0 || j==0 )
                {
                    lb.mixture[i][j][k].type = TYPE_FS;
                }
                if (i==Nx-1 || j==Ny-1 )
                {
                    lb.mixture[i][j][k].type = TYPE_O;
                }

                if (lb.mixture[i][j][k].type == TYPE_F || lb.mixture[i][j][k].type == TYPE_O  || lb.mixture[i][j][k].type == TYPE_O_C || lb.mixture[i][j][k].type == TYPE_I || lb.mixture[i][j][k].type == TYPE_S)
                {
                    lb.mixture[i][j][k].u = 0.0;
                    lb.mixture[i][j][k].v = 0.0;
                    lb.mixture[i][j][k].w = 0.0; 

                    lb.mixture[i][j][k].p = 5.0*units.p(Cantera::OneAtm);  

                    double theta = atan(j/(double)i);
                    double core_rad = 340.0 * (1.0 + 0.05*cos(4.0*4.0*theta) );

                    double temp_b = units.temp(1844.3); // [K]
                    double temp_u = units.temp(298.0); // [K]
                    lb.mixture[i][j][k].temp = smooth2D(temp_b, temp_u, i, j, core_rad, smoothness);
                    lb.species[gas->speciesIndex("H2")][i][j][k].X      = smooth2D(2.9513e-05   , 0.20134, i, j, core_rad, smoothness);
                    lb.species[gas->speciesIndex("H")][i][j][k].X       = smooth2D(1.2511e-06   , minim,   i, j, core_rad, smoothness);   
                    lb.species[gas->speciesIndex("O")][i][j][k].X       = smooth2D(2.2334e-05   , minim,   i, j, core_rad, smoothness);   
                    lb.species[gas->speciesIndex("O2")][i][j][k].X      = smooth2D(0.074442     , 0.16779, i, j, core_rad, smoothness);   
                    lb.species[gas->speciesIndex("OH")][i][j][k].X      = smooth2D(0.00069034   , minim,   i, j, core_rad, smoothness);   
                    lb.species[gas->speciesIndex("H2O")][i][j][k].X     = smooth2D(0.22346      , minim,   i, j, core_rad, smoothness);   
                    lb.species[gas->speciesIndex("HO2")][i][j][k].X     = smooth2D(1.4784e-06   , minim,   i, j, core_rad, smoothness);   
                    lb.species[gas->speciesIndex("H2O2")][i][j][k].X    = smooth2D(1.4581e-07   , minim,   i, j, core_rad, smoothness);   
                    lb.species[gas->speciesIndex("N2")][i][j][k].X      = smooth2D(0.70135      , 0.63087, i, j, core_rad, smoothness);                   


                }
            }
        }
    }

    lb.run(1000000000,1000);

    // LBM lb = read_restart("restart100000.dat");
    // lb.loop(100000000000, 1);
}

#elif defined HEATED_OBSTACLE
void main_setup() // Perfectly stirred reactor ------------------------------------------------------------------------------------------
{
    units.set_m_kg_s(1.25e-5, 1e-8 );

    int NX = 2100; 
    int NY = 80; 
    int NZ = 1;

    int h = 20;

    double nu = units.nu(1.568e-5) ;
    
    LBM lb(NX, NY, NZ, nu);
    int Nx = lb.get_Nx(); int Ny = lb.get_Ny(); int Nz = lb.get_Nz();

    lb.set_prtl(0.71);
    lb.set_gamma(1.4);
    lb.set_gasconst(units.cp(287));

    std::cout << "u_ave (lu)    : " << units.u(3.136) << std::endl;
    std::cout << "nu (lu)       : " << nu << std::endl;
    std::cout << "gas const (lu): " << lb.get_gasconst() << std::endl;
    std::cout << "RT (lu)       : " << lb.get_gasconst()*units.temp(300.0) << std::endl;




    #pragma omp parallel for schedule(dynamic)
    for(int i = 0; i < Nx ; ++i)
    {
        for(int j = 0; j < Ny; ++j)
        {
            for(int k = 0; k < Nz; ++k)
            {

                if ( k==0 || k==Nz-1) // set periodic boundary condition
                {
                    lb.mixture[i][j][k].type = TYPE_P;
                }

                if ( i==0 )
                {
                    lb.mixture[i][j][k].type = TYPE_I;
                    lb.mixture[i][j][k].u = units.u(3.136);
                    lb.mixture[i][j][k].v = 0.0;
                    lb.mixture[i][j][k].w = 0.0; 

                    lb.mixture[i][j][k].p = 1.0*units.p(Cantera::OneAtm);  
                    lb.mixture[i][j][k].temp = units.temp(300.0);                 
                }
                if ( i==Nx-1 )
                {
                    lb.mixture[i][j][k].type = TYPE_O;
                    lb.mixture[i][j][k].p = 1.0*units.p(Cantera::OneAtm);                 
                }

                if ( j==0 || j==Ny-1)
                {
                    lb.mixture[i][j][k].type = TYPE_A;
                }

                if ( i > 5*Ny && i < 5*Ny+h && j == 0){
                    lb.mixture[i][j][k].type = TYPE_Q;
                    lb.mixture[i][j][k].temp = units.temp(500.0);
                    lb.mixture[i][j][k].p = 1.0*units.p(Cantera::OneAtm); 
                    lb.mixture[i][j][k].energy_flux[1] = units.energy_flux(436.7); 
                }
                if ( i > 5*Ny && i < 5*Ny+h && j > 0 && j < h && lb.mixture[i][j][k].type != TYPE_P)
                {
                    lb.mixture[i][j][k].type = TYPE_S;
                    lb.mixture[i][j][k].temp = units.temp(300.0);
                    lb.mixture[i][j][k].p = 1.0*units.p(Cantera::OneAtm);  
                }   
                if ( i > 5.5*Ny && i < 5.5*Ny+h && j == Ny-1){
                    lb.mixture[i][j][k].type = TYPE_Q;
                    lb.mixture[i][j][k].temp = units.temp(500.0);
                    lb.mixture[i][j][k].p = 1.0*units.p(Cantera::OneAtm); 
                    lb.mixture[i][j][k].energy_flux[1] = units.energy_flux(436.7); 
                }             
                if ( i > 5.5*Ny && i < 5.5*Ny+h && j < Ny-1 && j > Ny-1-h && lb.mixture[i][j][k].type != TYPE_P)
                {
                    lb.mixture[i][j][k].type = TYPE_S;
                    lb.mixture[i][j][k].temp = units.temp(300.0);
                    lb.mixture[i][j][k].p = 1.0*units.p(Cantera::OneAtm);  
                }
                if ( i > 6*Ny && i < 6*Ny+h && j == 0){
                    lb.mixture[i][j][k].type = TYPE_Q;
                    lb.mixture[i][j][k].temp = units.temp(500.0);
                    lb.mixture[i][j][k].p = 1.0*units.p(Cantera::OneAtm); 
                    lb.mixture[i][j][k].energy_flux[1] = units.energy_flux(436.7); 
                }  
                if ( i > 6*Ny && i < 6*Ny+h && j > 0 && j < h && lb.mixture[i][j][k].type != TYPE_P)
                {
                    lb.mixture[i][j][k].type = TYPE_S;
                    lb.mixture[i][j][k].temp = units.temp(300.0);
                    lb.mixture[i][j][k].p = 1.0*units.p(Cantera::OneAtm);  
                }


                if (lb.mixture[i][j][k].type == TYPE_F || lb.mixture[i][j][k].type == TYPE_O )
                {
                    lb.mixture[i][j][k].u = units.u(3.136);
                    lb.mixture[i][j][k].v = 0.0;
                    lb.mixture[i][j][k].w = 0.0; 

                    lb.mixture[i][j][k].p = 1.0*units.p(Cantera::OneAtm);  
                    lb.mixture[i][j][k].temp = units.temp(300.0);
                }
               
            }
        }
    }

    lb.run(1000000000,1000);

    // LBM lb = read_restart("restart000001.dat");
    // lb.loop(840001, 1);
}

#elif defined HEAT_LID_DRIVEN
void main_setup() // Perfectly stirred reactor ------------------------------------------------------------------------------------------
{
    units.set_m_kg_s(1.0e-5, 0.5e-8 );

    int NX = 200; 
    int NY = 200; 
    int NZ = 1;

    double nu = units.nu(1.0e-4) ;
    
    LBM lb(NX, NY, NZ, nu);
    int Nx = lb.get_Nx(); int Ny = lb.get_Ny(); int Nz = lb.get_Nz();

    lb.set_prtl(0.71);
    lb.set_gamma(1.4);
    lb.set_gasconst(units.cp(287));

    double Re = 800;
    double u_wall = Re*nu/NX;

    double gas_const = lb.get_gasconst();
    double gamma = lb.get_gamma();
    double prtl = lb.get_prtl();
    double cv = gas_const / (gamma - 1.0);
    double cp = cv + gas_const;
    double mu = nu*(1.0);
    double conduc_coeff = mu*cp/prtl;

    std::cout << "u_ave (lu)    : " << u_wall << " " << units.si_u(u_wall) << std::endl;
    std::cout << "nu (lu)       : " << nu << std::endl;
    std::cout << "dynamic visc  : " << mu << " " << units.si_mu(mu) << std::endl;
    std::cout << "gas const (lu): " << lb.get_gasconst() << std::endl;
    std::cout << "cp (lu)       : " << cp << " " << units.si_cp(cp) << std::endl;
    std::cout << "thermal cond  : " << conduc_coeff << " " << units.si_thermalConductivity(conduc_coeff) << std::endl;
    std::cout << "RT (lu)       : " << lb.get_gasconst()*units.temp(300.0) << std::endl;

    #pragma omp parallel for schedule(dynamic)
    for(int i = 0; i < Nx ; ++i)
    {
        for(int j = 0; j < Ny; ++j)
        {
            for(int k = 0; k < Nz; ++k)
            {

                if ( k==0 || k==Nz-1) // set periodic boundary condition
                {
                    lb.mixture[i][j][k].type = TYPE_P;
                }

                if ( i==0 || i==Nx-1 || j==0 )
                {
                    lb.mixture[i][j][k].type = TYPE_S;
                    lb.mixture[i][j][k].u = 0.0;
                    lb.mixture[i][j][k].v = 0.0;
                    lb.mixture[i][j][k].w = 0.0; 

                    lb.mixture[i][j][k].p = 1.0*units.p(86100);  
                    lb.mixture[i][j][k].temp = units.temp(300.0);                 
                }
                if ( j==Ny-1 )
                {
                    lb.mixture[i][j][k].type = TYPE_S;
                    lb.mixture[i][j][k].u = u_wall;
                    lb.mixture[i][j][k].v = 0.0;
                    lb.mixture[i][j][k].w = 0.0; 

                    lb.mixture[i][j][k].p = 1.0*units.p(86100);  
                    lb.mixture[i][j][k].temp = units.temp(1500.0);                 
                }
                

                if (lb.mixture[i][j][k].type == TYPE_F || lb.mixture[i][j][k].type == TYPE_O )
                {
                    lb.mixture[i][j][k].u = 0.0;
                    lb.mixture[i][j][k].v = 0.0;
                    lb.mixture[i][j][k].w = 0.0; 

                    lb.mixture[i][j][k].p = 1.0*units.p(86100);  
                    lb.mixture[i][j][k].temp = units.temp(300.0);
                }
               
            }
        }
    }

    lb.run(1000000000,1000);

    // LBM lb = read_restart("restart000001.dat");
    // lb.loop(840001, 1);
}

#elif defined HEAT_LID_DRIVEN_MULTICOMP
void main_setup() // Perfectly stirred reactor ------------------------------------------------------------------------------------------
{
    units.set_m_kg_s(1.0e-6, 5e-10 );

    int NX = 200; 
    int NY = 200; 
    int NZ = 1;
    
    auto sol = Cantera::newSolution("h2o2.yaml", "ohmech");
    auto gas = sol->thermo();

    std::vector<std::string> species = { "N2" };

    LBM lb(NX, NY, NZ, species);
    int Nx = lb.get_Nx(); int Ny = lb.get_Ny(); int Nz = lb.get_Nz();

    auto trans = sol->transport();
    std::vector <double> X (gas->nSpecies());
    X[gas->speciesIndex("N2")] = 1.0;
    gas->setMoleFractions(&X[0]);
    gas->setState_TP(300.0, 1.0*Cantera::OneAtm);

    std::cout << "nu (lu) : " << units.nu(trans->viscosity() / gas->density()) << std::endl;
    std::cout << "gas const : " << units.cp(Cantera::GasConstant/gas->meanMolecularWeight()) << std::endl;
    std::cout << "temperature : " << units.temp(gas->temperature()) << std::endl;
    std::cout << "gamma : " << gas->cp_mass()/gas->cv_mass() << std::endl;
    std::cout << "sound speed : " << units.u(gas->soundSpeed()) << " " << gas->soundSpeed() <<  std::endl;
    std::cout << "RT : " << units.cp(Cantera::GasConstant/gas->meanMolecularWeight())*units.temp(gas->temperature()) << std::endl;

    double Re = 400;
    double u_wall = Re*units.nu(trans->viscosity() / gas->density()) / NX;

  

    #pragma omp parallel for schedule(dynamic)
    for(int i = 0; i < Nx ; ++i)
    {
        for(int j = 0; j < Ny; ++j)
        {
            for(int k = 0; k < Nz; ++k)
            {

                if ( k==0 || k==Nz-1) // set periodic boundary condition
                {
                    lb.mixture[i][j][k].type = TYPE_P;
                }

                if ( i==0 || i==Nx-1 || j==0 )
                {
                    lb.mixture[i][j][k].type = TYPE_S;
                    lb.mixture[i][j][k].u = 0.0;
                    lb.mixture[i][j][k].v = 0.0;
                    lb.mixture[i][j][k].w = 0.0; 

                    lb.mixture[i][j][k].p = 1.0*units.p(Cantera::OneAtm);  
                    lb.mixture[i][j][k].temp = units.temp(300.0);     
                    
                    lb.species[0][i][j][k].X = 1.0;                 
                }
                if ( j==Ny-1 )
                {
                    lb.mixture[i][j][k].type = TYPE_S;
                    lb.mixture[i][j][k].u = u_wall;
                    lb.mixture[i][j][k].v = 0.0;
                    lb.mixture[i][j][k].w = 0.0; 

                    lb.mixture[i][j][k].p = 1.0*units.p(Cantera::OneAtm);  
                    lb.mixture[i][j][k].temp = units.temp(1500.0);       
                    lb.species[0][i][j][k].X = 1.0;                 
                }
                

                if (lb.mixture[i][j][k].type == TYPE_F || lb.mixture[i][j][k].type == TYPE_O )
                {
                    lb.mixture[i][j][k].u = 0.0;
                    lb.mixture[i][j][k].v = 0.0;
                    lb.mixture[i][j][k].w = 0.0; 

                    lb.mixture[i][j][k].p = 1.0*units.p(Cantera::OneAtm);  
                    lb.mixture[i][j][k].temp = units.temp(300.0);

                    lb.species[0][i][j][k].X = 1.0;                 

                }
               
            }
        }
    }

    lb.run(1000000000,1000);

    // LBM lb = read_restart("restart000001.dat");
    // lb.loop(840001, 1);
}


#endif


