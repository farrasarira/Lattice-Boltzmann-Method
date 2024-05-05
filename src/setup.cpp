#include "lbm.hpp"
#include "geometry.hpp"
#include "units.hpp"
#include "cantera.hpp"
#include <iostream>
#include <math.h>
#include <omp.h>

Units_Conv units;

#if defined CYLINDER_2D
void main_setup() // 2D Flow over cylinder --------------------------------------------------------
{
    double RE = 1000;
    LBM lb(1000, 500, 5, 0.00001);
    int Nx = lb.get_Nx(); int Ny = lb.get_Ny(); int Nz = lb.get_Nz();
    double nu = lb.get_nu();

    int D = Ny/10; 
    cylinder_generator(lb, D);

    double u_max = RE*nu/D;    // [Velocity in Lattice Unit]
    double rho0 = 1.0;      // [Density in Lattice Unit]
    double T0 = 0.03;       // [Temperature in Lattice Unit]

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

    #pragma omp parallel for
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
                    lb.mixture[i][j][k].rho = rho0;
                    lb.mixture[i][j][k].u = u_max;
                    lb.mixture[i][j][k].v = 0;
                    lb.mixture[i][j][k].w = 0;
                    lb.mixture[i][j][k].temp = T0;
                }
                if(lb.mixture[i][j][k].type == TYPE_S)
                {
                    lb.mixture[i][j][k].temp = T0;
                }
            }
        }
    }

    lb.run(1000, 10);
}

#elif defined TAYLOR_GREEN_3D
void main_setup() // 3D Taylor-Green Vortex
{
    double RE = 200;
    LBM lb(100, 100, 100, 0.00297753);
    // LBM lb(100, 100, 100, 0.025);


    const double gasconst = 0.033059;
    const double gamma = 1.4;
    const double rho0 = 1.0;
    const double p0 = 0.1;
    const double temp0 = p0/(gasconst * rho0); 
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

    #pragma omp parallel for
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
    
    lb.run(10000000,100);
}
#elif defined TAYLOR_GREEN_2D
void main_setup() // 2D Taylor-Green Vortex
{
    double RE = 10;

    LBM lb(30, 30, 1, 0.01);
    int Nx = lb.get_Nx(); int Ny = lb.get_Ny(); int Nz = lb.get_Nz();
    double NX = lb.get_NX();
    double NU = lb.get_nu();

    const double u_max = RE * NU / NX;    // maximum initial velocity
    std::cout << "umax      : " << u_max << std::endl;
    const unsigned int periodicity = 1;
    const double a = (double)NX/periodicity;
    const double b = (double)NX/periodicity;
    
    #pragma omp parallel for
    for(int i = 0; i < Nx ; ++i)
    {
        for(int j = 0; j < Ny; ++j)
        {
            for(int k = 0; k < Nz; ++k)
            {
                if (i==0 || i==Nx-1 ||  j==0 || j==Ny-1) // set periodic boundary condition
                {
                    lb.mixture[i][j][k].type = TYPE_P;
                }
                if (lb.mixture[i][j][k].type == TYPE_F)
                {
                    double X = i ;
                    double Y = j ;
                    lb.mixture[i][j][k].u = -u_max * cos(2.0*M_PI/a*X) * sin(2.0*M_PI/b*Y) ;
                    lb.mixture[i][j][k].v =  u_max * sin(2.0*M_PI/a*X) * cos(2.0*M_PI/b*Y) ;
                    lb.mixture[i][j][k].rho = 1.0 - u_max*u_max * 3.0/4.0 * (cos(4.0*M_PI/a*X) + cos(4.0*M_PI/b*Y));
                    lb.mixture[i][j][k].p = 1.0/3.0;
                    lb.mixture[i][j][k].temp = lb.mixture[i][j][k].p / (1.0 - u_max*u_max * 3.0/4.0 * (cos(4.0*M_PI/a*X) + cos(4.0*M_PI/b*Y)));
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
    
    #pragma omp parallel for
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

    #pragma omp parallel for
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
    double nu = 0.1;

    LBM lb(NX,NY,NZ, nu);
    int Nx = lb.get_Nx(); int Ny = lb.get_Ny(); int Nz = lb.get_Nz();

    const double mach = 0.9;
    const double gamma = 1.4;
    const double a0 = 0.001; // amplitude
    const double T0 = 0.1;
    const double u0 = mach*sqrt(gamma*T0); 

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
    
    #pragma omp parallel for
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
                    lb.mixture[i][j][k].u = a0*sin(Y/NY*2*M_PI) + u0*cos(M_PI/6);
                    lb.mixture[i][j][k].v = u0*sin(M_PI/6);
                    lb.mixture[i][j][k].w = 0.0;
                    lb.mixture[i][j][k].temp = 0.1;
                    lb.mixture[i][j][k].p = 0.1;
                    
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
    
    #pragma omp parallel for
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
    double NX = 3000;
    double NY = 1;
    double NZ = 1;
    double nu = 0.1;

    LBM lb(NX,NY,NZ, nu);
    int Nx = lb.get_Nx(); int Ny = lb.get_Ny(); int Nz = lb.get_Nz();

    const double gamma = 1.4;
    const double gascont = 0.02;
    
    const double temp = 1.0;
    const double dp = 1e-4;

    lb.set_gamma(gamma);
    lb.set_gasconst(gascont);

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
    
    #pragma omp parallel for
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
    lb.run(1000,100);
}

#elif defined SOD_SHOCK
void main_setup() // Sos shock tube test case --------------------------------------------------------
{
    int NX = 3000; 
    int NY = 1; 
    int NZ = 1;
    
    double NU = 0.025;

    double si_len = 1E-0;    // [m]
    double si_u_max = 10.0;  // [m/s]
    double si_rho = 1.225;  // [kg/m^3]
    double si_temp = 300.0;// [K]

    units.set_m_kg_s(NX, VEL0, RHO0, TEMP0, si_len, si_u_max, si_rho, si_temp); // setting the conversion factor 

    LBM lb(NX, NY, NZ, NU);
    int Nx = lb.get_Nx(); int Ny = lb.get_Ny(); int Nz = lb.get_Nz();

    
    #pragma omp parallel for
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
                    lb.mixture[i][j][k].p = smooth(2.0, 1.0, i, 0.5*Nx, 0.03);   
                    // std::cout << " p : " << units.si_p(2.0) << std::endl;
                    lb.mixture[i][j][k].temp = smooth(0.2, 0.1, i, 0.5*Nx, 0.03);   
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

    lb.run(1000,100);
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

    
    #pragma omp parallel for
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
    int NX = 900; 
    int NY = 1; 
    int NZ = 1;
    
    double si_len = 1E-3;    // [m]
    double si_u_max = 1E+3;  // [m/s]
    double si_rho = 1.225;  // [kg/m^3]
    double si_temp = 300.0;// [K]

    // int NX = 6000; 
    // int NY = 1; 
    // int NZ = 1;
    
    // double si_len = 1E-3;    // [m]
    // double si_u_max = 1.0E+3;  // [m/s]
    // double si_rho = 1.225;  // [kg/m^3]
    // double si_temp = 400.0;// [K]

    units.set_m_kg_s(NX, VEL0, RHO0, 0.1, si_len, si_u_max, si_rho, si_temp); // setting the conversion factor 

    std::vector<std::string> species = { "Ar" };
    // std::vector<std::string> species = { "ch4", "o2"};
    // std::vector<std::string> species = {"Ar", "CH4", "NH3"};
    // std::vector<std::string> species = {"H2", "N2", "C2H6"};
    
    LBM lb(NX, NY, NZ, species);
    int Nx = lb.get_Nx(); int Ny = lb.get_Ny(); int Nz = lb.get_Nz();


    #pragma omp parallel for
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
                    // lb.species[0][i][j][k].X = smooth(0.491, 0.000, i, 0.5*Nx, 0.1);
                    // lb.species[1][i][j][k].X = smooth(0.509, 0.485, i, 0.5*Nx, 0.1);
                    // lb.species[2][i][j][k].X = smooth(0.000, 0.515, i, 0.5*Nx, 0.1);

                    // lb.species[0][i][j][k].X = smooth(0.45, 0.55, i, 0.5*Nx, 0.07);
                    // lb.species[1][i][j][k].X = smooth(0.55, 0.45, i, 0.5*Nx, 0.07);
                    // lb.species[2][i][j][k].X = smooth(0.58, 0.43, i, 0.5*Nx, 0.07);

                    lb.species[0][i][j][k].X = 1.0;

                    // lb.species[0][i][j][k].X = smooth(0.52, 0.48, i, 0.5*Nx, 0.3);
                    // lb.species[1][i][j][k].X = smooth(0.48, 0.52, i, 0.5*Nx, 0.3);

                    // lb.mixture[i][j][k].p = smooth(2*units.p(Cantera::OneAtm), units.p(Cantera::OneAtm), i, 0.5*Nx, 0.3);      
                    // lb.mixture[i][j][k].temp = smooth(units.temp(400.0), units.temp(400.0), i, 0.5*Nx, 0.3);                  

                    lb.mixture[i][j][k].p = smooth(2*units.p(Cantera::OneAtm), units.p(Cantera::OneAtm), i, 0.5*Nx, 0.3);      
                    lb.mixture[i][j][k].temp = smooth(units.temp(300.0), units.temp(300.0), i, 0.5*Nx, 0.3); 

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

    lb.run(1,1);
}

#elif defined TERNARY_DIFFUSION

void main_setup() // Ternary Gas Diffusion --------------------------------------------------------
{
    int NX = 200; 
    int NY = 1; 
    int NZ = 1;
    
    double si_len = 1E-4;   // [m]
    double si_u_max = 1E+3; // [m/s]
    double si_rho = 1.225;  // [kg/m^3]
    double si_temp = 300.0; // [K]

    units.set_m_kg_s(NX, VEL0, RHO0, 0.1, si_len, si_u_max, si_rho, si_temp); // setting the conversion factor 

    std::vector<std::string> species = { "H2", "Ar" , "CH4"};
    // std::vector<std::string> species = { "AR" , "N2"};
    
    LBM lb(NX, NY, NZ, species);
    int Nx = lb.get_Nx(); int Ny = lb.get_Ny(); int Nz = lb.get_Nz();


    #pragma omp parallel for
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
                    lb.species[0][i][j][k].X = smooth(0.491, 0.000, i, 0.5*Nx, 0.7);
                    lb.species[1][i][j][k].X = smooth(0.509, 0.485, i, 0.5*Nx, 0.7);
                    lb.species[2][i][j][k].X = smooth(0.000, 0.515, i, 0.5*Nx, 0.7);

                    // lb.species[0][i][j][k].X = smooth(0.491, 0.300, i, 0.5*Nx, 0.03);
                    // lb.species[1][i][j][k].X = smooth(0.209, 0.185, i, 0.5*Nx, 0.03);
                    // lb.species[2][i][j][k].X = smooth(0.300, 0.515, i, 0.5*Nx, 0.03);

                    // lb.species[0][i][j][k].X = smooth(1.0, 0.5, i, 0.5*Nx, 0.3);
                    // lb.species[1][i][j][k].X = smooth(0.5, 1.0, i, 0.5*Nx, 0.3);

                    
                    lb.mixture[i][j][k].p = smooth(1*units.p(Cantera::OneAtm), units.p(Cantera::OneAtm), i, 0.5*Nx, 0.3);      
                    lb.mixture[i][j][k].temp = smooth(units.temp(300.0), units.temp(300.0), i, 0.5*Nx, 0.3);                  
                  
                    lb.mixture[i][j][k].u = 0.0;
                    lb.mixture[i][j][k].v = 0.0;
                    lb.mixture[i][j][k].w = 0.0;
                     
                }
            }
        }
    }

    lb.run(1000000000,1000);
}


#elif defined SHEAR_LAYER_MULTICOMP
void main_setup() // 3D Shear layer multicomponent
{
    int NX = 100; 
    int NY = 100; 
    int NZ = 1;
    
    // units.set_m_kg_s(NX, VEL0, RHO0, 0.025, si_len, si_u_max, si_rho, si_temp); // setting the conversion factor 

    double dx = 3e-5;

    units.set_m_kg_s(dx, 1e-9);

    double Re = 500.0;
    double H = NX*dx;

    auto sol = Cantera::newSolution("gri30.yaml", "gri30","mixture-averaged");
    auto gas = sol->thermo();
    gas->setState_TP(300.0, Cantera::OneAtm);
    std::vector <double> X(gas->nSpecies());
    X[gas->speciesIndex("N2")] = 0.4;
    X[gas->speciesIndex("H2O")] = 0.6;
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
    std::vector<std::string> species = { "N2" , "H2O"};
    
    LBM lb(NX, NY, NZ, species);
    int Nx = lb.get_Nx(); int Ny = lb.get_Ny(); int Nz = lb.get_Nz();


    #pragma omp parallel for
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
                    lb.mixture[i][j][k].temp = smooth(units.temp(450.0), units.temp(450.0), i, 0.5*Nx, 0.3);                  
                  
                    lb.mixture[i][j][k].v = 0.05*u_max*sin(2*M_PI*i/Nx);
                    lb.mixture[i][j][k].w = 0.0;

                    if ((float)j/(float)Ny <= 0.5){
                        lb.mixture[i][j][k].u = smooth(u_max, -u_max, j, 0.25*Ny, 0.3);
                        lb.species[0][i][j][k].X = smooth(0.4, 0.4, j, 0.25*Ny, 0.3);
                        lb.species[1][i][j][k].X = smooth(0.6, 0.6, j, 0.25*Ny, 0.3);
                    } 
                    else{
                        lb.mixture[i][j][k].u = smooth(-u_max, u_max, j, 0.75*Ny, 0.3);
                        lb.species[0][i][j][k].X = smooth(0.4, 0.4, j, 0.75*Ny, 0.3);
                        lb.species[1][i][j][k].X = smooth(0.6, 0.6, j, 0.75*Ny, 0.3);
                    }
                     
                }
            }
        }
    }

    lb.run(100000000,1);
}


#elif defined PERFECTLY_STIRRED_REACTOR_3D
void main_setup() // Perfectly stirred reactor ------------------------------------------------------------------------------------------
{
    int NX = 4; 
    int NY = 4; 
    int NZ = 4;
    
    double si_len = 0.01;    // [m]
    double si_u_max = 1.0E+3;  // [m/s]
    double si_rho = 1.225;  // [kg/m^3]
    double si_temp = 1400.0;// [K]

    units.set_m_kg_s(NX, VEL0, RHO0, 0.3, si_len, si_u_max, si_rho, si_temp); // setting the conversion factor 

    std::vector<std::string> species = { "H2", "H2O" };

    LBM lb(NX, NY, NZ, species);
    int Nx = lb.get_Nx(); int Ny = lb.get_Ny(); int Nz = lb.get_Nz();

    #pragma omp parallel for
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
                    lb.species[0][i][j][k].X = 2./3.;     // H2 Mole Fraction
                    lb.species[1][i][j][k].X = 1./3.;     // H2O Mole Fraction
        
                }
            }
        }
    }

    lb.run(1000,10);
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


    #pragma omp parallel for
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
    int NX = 2; 
    int NY = 100; 
    int NZ = 1;
    
    double si_len = 1E-3;    // [m]
    double si_u_max = 1.0E+2;  // [m/s]
    double si_rho = 1.225;  // [kg/m^3]
    double si_temp = 300.0;// [K]
    double nu = 0.025;

    units.set_m_kg_s(NX, VEL0, RHO0, TEMP0, si_len, si_u_max, si_rho, si_temp); // setting the conversion factor 

    
    LBM lb(NX, NY, NZ, nu);

    int Nx = lb.get_Nx(); int Ny = lb.get_Ny(); int Nz = lb.get_Nz();


    #pragma omp parallel for
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
                    lb.mixture[i][j][k].temp = 0.1;
                }
                if (j==Ny-1 )
                {
                    lb.mixture[i][j][k].type = TYPE_S;
                    lb.mixture[i][j][k].u = 0.05;
                    lb.mixture[i][j][k].temp = 0.1003125;
                }

                if (lb.mixture[i][j][k].type == TYPE_F)
                {
                    lb.mixture[i][j][k].temp = 0.1;
                    lb.mixture[i][j][k].p = 1.0;
                    lb.mixture[i][j][k].u = 0.0;
                }
            }
        }
    }

    lb.set_prtl(0.25);
    lb.set_gasconst(0.4/1.4);
    lb.run(200000,10000);
}

#elif defined COUETTE_FLOW_MULTICOMP
void main_setup() // 2D Couette flow (multicomponent)  --------------------------------------------------------
{
    int NX = 2; 
    int NY = 50; 
    int NZ = 5;
    
    double si_len = 1E-6;    // [m]
    double si_u_max = 1.0E+3;  // [m/s]
    double si_rho = 0.5;  // [kg/m^3]
    double si_temp = 300.0;// [K]
    

    units.set_m_kg_s(NX, VEL0, RHO0, TEMP0, si_len, si_u_max, si_rho, si_temp); // setting the conversion factor 

    std::vector<std::string> species = { "AR" };
    
    LBM lb(NX, NY, NZ, species);

    int Nx = lb.get_Nx(); int Ny = lb.get_Ny(); int Nz = lb.get_Nz(); int nSpecies = lb.get_nSpecies();


    #pragma omp parallel for
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
                    lb.mixture[i][j][k].temp = units.si_temp(300);
                }
                if (j==Ny-1 )
                {
                    lb.mixture[i][j][k].type = TYPE_S;
                    lb.mixture[i][j][k].u = units.u(10);
                    lb.mixture[i][j][k].temp = units.temp(350);
                }

                if (lb.mixture[i][j][k].type == TYPE_F)
                {
                    lb.mixture[i][j][k].temp = units.temp(300);
                    lb.mixture[i][j][k].p = units.p(Cantera::OneAtm);
                    
                    for(int a=0; a < nSpecies; ++a)
                        lb.species[a][i][j][k].X = 1.0;
                }
            }
        }
    }
    lb.run(20000,1000);
}

#elif defined TAYLOR_GREEN_3D_MULTICOMP
void main_setup() // 3D Taylor-Green Vortex (multicomponent) ------------------------------------------------------------------------------------
{
    double RE = 200;
    int NX = 150; 
    int NY = 150; 
    int NZ = 150;
    
    units.set_m_kg_s(1e-7, 1e-11);

    std::vector<std::string> species = { "CO2" , "AR" };
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
    X[gas->speciesIndex("CO2")] = 1.0;
    X[gas->speciesIndex("AR")] = 2.0;
    // X[gas->speciesIndex("H2")] = 1.0;
    gas->setMoleFractions(&X[0]);
    gas->setState_TP(temp, pressure);
    auto trans = sol->transport();
    double NU = units.nu(trans->viscosity() / gas->density());
    std::cout << "NU : " << NU << std::endl;

    double u_max = RE * NU / (NX/(2.*M_PI));
    std::cout << "umax : " << u_max << std::endl; 
    

    
    #pragma omp parallel for
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
                    for(int spec_idx = 0; spec_idx < nSpecies; ++spec_idx)
                        lb.species[spec_idx][i][j][k].X = 1.0 + spec_idx;
                    lb.mixture[i][j][k].u =  u_max * sin(2.0*M_PI/a*(double)i) * cos(2.0*M_PI/b*(double)j) * cos(2.0*M_PI/c*(double)k);
                    lb.mixture[i][j][k].v = -u_max * cos(2.0*M_PI/a*(double)i) * sin(2.0*M_PI/b*(double)j) * cos(2.0*M_PI/c*(double)k);
                    lb.mixture[i][j][k].w = 0.0;
                    // lb.mixture[i][j][k].rho = units.rho(gas->density()) - u_max*u_max * 3.0/16.0 * (cos(4.0*M_PI/a*(double)i) + cos(4.0*M_PI/b*(double)j)) * (cos(4.0*M_PI/c*(double)k) + 2);
                    lb.mixture[i][j][k].temp = units.temp(temp);
                    lb.mixture[i][j][k].p = units.p(pressure);
                }
            }
        }
    }
    
    lb.run(100,10);
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


    #pragma omp parallel for
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

#endif


