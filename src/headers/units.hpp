#ifndef UNIT

    #define UNIT
    
    #include"math_util.hpp"
    #include<math.h>
    
    class Units_Conv {  // contains the 3 base units m, kg, s for unit conversions
    private:
        float m=1.0f, kg=1.0f, s=1.0f, K=1.0f;  // 1 lattice unit times m/kg/s is meter/kilogram/seconds
    public:
        void set_m_kg_s(const float x, const float u, const float rho, const float temp, const float si_x, const float si_u, const float si_rho, const float si_temp) { // length x, velocity u, density rho in both simulation and SI units
            m = si_x/x;                     // length si_x = x*[m]
            kg = si_rho/rho*cb((double)m);  // density si_rho = rho*[kg/m^3]
            s = m * u/si_u;                 // velocity si_u = u*[m/s]
            K = si_temp/temp;               // temperature

            std::cout << "=== UNIT CONVERSION ===" << std::endl;
            std::cout << "dt* : " << s << std::endl;
            std::cout << "dx* : " << m << std::endl;
            std::cout << "dT* : " << K << std::endl;
            std::cout << "dkg* : " << kg << std::endl;
            std::cout << "=======================" << std::endl;

        }
        // void set_m_kg_s(const float m, const float kg, const float s, const float K) { // do unit conversion manually
        //     this->m = m;
        //     this->kg = kg;
        //     this->s = s;
        //     this->K = K;
        // }
        void set_m_kg_s(const float dx, const float dt) { // do unit conversion manually
            this->m = dx/1.0;
            this->kg = 1.0/1.0*cb((double)m);
            this->s = dt/1.0;
            this->K = 12000.0;
            
            std::cout << "=== UNIT CONVERSION ===" << std::endl;
            std::cout << "dt* : " << s << std::endl;
            std::cout << "dx* : " << m << std::endl;
            std::cout << "dT* : " << K << std::endl;
            std::cout << "dkg* : " << kg << std::endl;
            std::cout << "=======================" << std::endl;
        }

        void set_m_kg_s(const float lbm_len, const float ref_len, const float sound_spd) { // do unit conversion manually
            this->m = ref_len/lbm_len;
            this->kg = 1.0/1.0*cb((double)m);
            this->s = m/sound_spd/10.0;
            this->K = 12000.0;
            
            std::cout << "=== UNIT CONVERSION ===" << std::endl;
            std::cout << "dt* : " << s << std::endl;
            std::cout << "dx* : " << m << std::endl;
            std::cout << "dT* : " << K << std::endl;
            std::cout << "dkg* : " << kg << std::endl;
            std::cout << "=======================" << std::endl;
        }

        // void set_m_kg_s(const float ref_len, const float sound_spd) { // do unit conversion manually
        //     this->m = ref_len;
        //     this->kg = 1.0/1.0*cb((double)m);
        //     this->s = ref_len/(1.35*sound_spd);
        //     this->K = 12000.0;
            
        //     std::cout << "=== UNIT CONVERSION ===" << std::endl;
        //     std::cout << "dt* : " << s << std::endl;
        //     std::cout << "dx* : " << m << std::endl;
        //     std::cout << "dT* : " << K << std::endl;
        //     std::cout << "dkg* : " << kg << std::endl;
        //     std::cout << "=======================" << std::endl;
        // }

        // the following methods convert SI units into simulation units (have to be called after set_m_kg_s(...);)
        inline float x(const float si_x) const { return si_x/m; } // length si_x = x*[m]
        // float M(const float si_M) const { return si_M/kg; } // mass si_M = M*[kg]
        // float t(const float si_t) const { return (si_t/s); } // time si_t = t*[s]
        inline float temp(const float si_temp) const {return si_temp/K;} // temperature
        inline float energy_mass(const float si_energy_mass) const {return si_energy_mass * (s*s)/(m*m);}
        // float energy(const float si_energy) const {return si_energy * (s*s)/(kg*m*m);}
        // float frequency(const float si_frequency) const { return si_frequency*s; } // frequency si_frequency = frequency*[1/s]
        // float omega(const float si_omega) const { return si_omega*s; } // frequency si_omega = omega/[s]
        inline float u(const float si_u) const { return si_u*s/m; } // velocity si_u = u*[m/s]
        inline float rho(const float si_rho) const { return si_rho*cb(m)/kg; } // density si_rho = rho*[kg/m^3]
        // float rho_u(const float si_rho_u) const { return si_rho_u* cb(m)/kg * s/m; } // mass flux si_rho_u = rho_u
        inline float rho_dot(const float si_rho) const { return si_rho*cb(m)*s/kg; } // density si_rho_dot = rho*[kg/m^3/s]
        // float Q(const float si_Q) const { return si_Q*s/cb(m); } // flow rate si_Q = Q*[m^3/s]
        inline float nu(const float si_nu) const { return si_nu*s/sq(m); } // kinematic shear viscosity si_nu = nu*[m^2/s]
        inline float mu(const float si_mu) const { return si_mu*s*m/kg; } // dynamic shear viscosity si_mu = mu*[kg/(m*s)]
        // float g(const float si_g) const { return si_g/m*sq(s); } // gravitational acceleration si_g = g*[m/s^2]
        // float f(const float si_f) const { return si_f*sq(m*s)/kg; } // force per volume si_f = f*[kg/(m*s)^2]
        // float f(const float si_rho, const float si_g) const { return si_rho*si_g*sq(m*s)/kg; } // force per volume f = rho*g = si_rho/[kg/m^3]*si_g/[m/s^2] = si_rho*si_g*[(m*s)^2/kg]
        // float F(const float si_F) const { return si_F*sq(s)/(kg*m); } // force si_F = F*[kg*m/s^2]
        // float T(const float si_T) const { return si_T*sq(s)/(kg*sq(m)); } // torque si_T = T*[kg*m^2/s^2]
        // float sigma(const float si_sigma) const { return si_sigma*sq(s)/kg; } // surface tension si_sigma = sigma*[kg/s^2]
        inline float thermalConductivity(const float si_lamda) const { return si_lamda*s*s*s*K / (kg*m) ; } // thermal conductivity
        inline float cp(const float si_cp) const { return si_cp*s*s*K / (m*m); }
        // float R(const float si_R) const{ return si_R*s*s*K / (m*m);}  
        inline float p(const float si_p) const { return si_p*(m*sq(s))/kg; } // pressure 
        // float bin_diff(const float si_bin_diff) const { return si_bin_diff * s / (m*m); } // pressure 
        // float diff_coeff(const float si_dif_coeff) const { return si_dif_coeff*s/sq(m); } // diffusion coefficient si_nu = nu*[m^2/s]
        inline float permeability(const float si_K) const { return si_K/sq(m); } // permeability K
        inline float HRR(const float si_HRR) const { return si_HRR / (kg/(m*cb(s))); } // Heat Relese Rate 


        
        // // the following methods convert simulation units into SI units (have to be called after set_m_kg_s(...);)
        inline float si_x(const int x) const { return (float)x*m; } // length si_x = x*[m]
        // float si_x(const float x) const { return x*m; } // length si_x = x*[m]
        // float si_M(const float M) const { return M*kg; } // mass si_M = M*[kg]
        // float si_t(const unsigned long t) const { return (float)t*s; } // time si_t = t*[s]
        inline float si_temp(const float temp) const {return temp*K;} // temperature
        inline float si_energy_mass(const float energy_mass) const {return energy_mass * (m*m)/(s*s);}
        inline float si_energy(const float energy_mass) const {return energy_mass * kg*(m*m)/(s*s);}
        // float si_mu(const float mu) const {return mu * kg / (m * s) ;} // dynamic viscosity
        // float si_frequency(const float frequency) const { return frequency/s; } // frequency si_frequency = frequency*[1/s]
        // float si_V(const float V) const { return V*cb(m); } // volume si_V = V*[m^3]
        inline float si_u(const float u) const { return u*m/s; } // velocity si_u = u*[m/s]
        inline float si_rho(const float rho) const { return rho*kg/cb(m); } // density si_rho = rho*[kg/m^3]
        inline float si_p(const float p) const { return p*kg/(m*sq(s)); } // pressure si_p = p*[kg/(m*s^2)]
        // float si_Q(const float Q) const { return Q*cb(m)/s; } // flow rate si_Q = Q*[m^3/s]
        // float si_nu(const float nu) const { return nu*sq(m)/s; } // kinematic shear viscosity si_nu = nu*[m^2/s]
        // float si_g(const float g) const { return g*m/sq(s); } // gravitational acceleration si_g = g*[m/s^2]
        // float si_f(const float f) const { return f*kg/sq(m*s); } // force per volume si_f = f*[kg/(m*s)^2]
        // float si_F(const float F) const { return F*kg*m/sq(s); } // force si_F = F*[kg*m/s^2]
        // float si_T(const float T) const { return T*kg*sq(m)/sq(s); } // torque si_T = T*[kg*m^2/s^2]
        // float si_sigma(const float sigma) const { return sigma*kg/sq(s); } // surface tension si_sigma = sigma*[kg/s^2]
        inline float si_HRR(const float HRR) const { return HRR * (kg/(m*cb(s))); } // Heat Relese Rate  


        // // other conversions in simulation units (can be called before set_m_kg_s(...);)
        // float Re(const float si_Re) const { return si_Re; } // Reynolds number Re = x*u/nu = [1] no unit
        // float Re(const float x, const float u, const float nu) const { return x*u/nu; } // Reynolds number Re = x*u/nu = [1] no unit
        // float Re(const float x, const float u, const float mu, const float rho) const { return x*u*rho/mu; } // Reynolds number Re = x*u*rho/mu = [1] no unit
        // float We(const float x, const float u, const float rho, const float sigma) const { return x*sq(u)*rho/sigma; } // Weber number We = x*u^2*rho/sigma = [1] no unit
        // float Fr(const float x, const float u, const float g) const { return u/sqrt(x*g); } // Froude number Fr = u/sqrt(x*g) = [1] no unit
        // float Ca(const float u, const float mu, const float sigma) const { return u*mu/sigma; } // Capillary number Ca = u*mu/sigma = [1] no unit
        // float Ca(const float u, const float rho, const float nu, const float sigma) const { return u*rho*nu/sigma; } // Capillary number Ca = u*rho*nu/sigma = [1] no unit
        // float Bo(const float x, const float rho, const float g, const float sigma) const { return sq(x)*rho*g/sigma; } // Bond number Bo = x^2*rho*g/sigma = [1] no unit
        // float Mo(const float rho, const float delta_rho, const float g, const float sigma, const float mu) const { return g*delta_rho*sq(sq(mu))/(cb(sigma)*sq(rho)); } // Morton number g*delta_rho*mu^4/(sigma^3*rho^2)
        // float Ga(const float x, const float nu, const float g) { return cb(x)*g/sq(nu); } // Galilei number Ga = x^3*g/nu^2 = [1] no unit
        // float Ga(const float x, const float rho, const float nu, const float f) { return cb(x)*f/(sq(nu)*rho); } // Galilei number Ga = x^3*g/nu^2 = [1] no unit
        // float Ma(const float u) const { return u/0.57735027f; } // Mach number Ma = u/c = [1] no unit, c = 1/sqrt(3) is lattice speed of sound
        // float rho_from_p(const float p) const { return 3.0f*p; } // density rho = p/c^2 = 3*p = [kg/(m*s^2)], p is pressure, c = 1/sqrt(3) is lattice speed of sound
        // float rho_laplace(const float sigma, const float R) { return 6.0f*sigma/R; } // sphere laplace pressure p = 2*sigma/R, density rho = 3*p
        // float rho_hydrostatic(const float f, const float h, const float h0) { return 3.0f*f*(h0-h)+1.0f; } // hydrostatic pressure p = rho*g*h+p0 = f*h+1/3, density rho = 3*p, force per volume f = rho*g, rho0 = 1
        // float nu_from_mu(const float mu, const float rho) const { return mu/rho; } // kinematic shear viscosity nu = mu/rho = [m^2/s]
        // float nu_from_tau(const float tau) const { return (tau-0.5f)/3.0f; } // kinematic shear viscosity nu = (tau-1/2)/3 = [m^2/s]
        // float nu_from_Re(const float Re, const float x, const float u) const { return x*u/Re; } // kinematic shear viscosity nu = x*u/Re = [m^2/s]
        // float f_from_F(const float F, const float V) const { return F/V; } // force per volume f = F/V = [kg/(m*s)^2]
        // float f_from_g(const float g, const float rho) const { return rho*g; } // force per volume f = rho*g = [kg/(m*s)^2]
        // float g_from_f(const float f, const float rho) const { return f/rho; } // force per volume f = rho*g = [kg/(m*s)^2]
        // float u_from_Re(const float Re, const float x, const float nu) const { return Re*nu/x; } // velocity u = Re*nu/x, Reynolds number Re = x*u/nu = [1] no unit
        // float u_from_Re(const float Re, const float x, const float mu, const float rho) const { return Re*mu/(x*rho); } // velocity u = Re*nu/x, Reynolds number Re = x*u*rho/mu = [1] no unit
        // float u_from_Ma(const float Ma) const { return 0.57735027f*Ma; } // velocity u = c*Ma, Mach number Ma = u/c = [1] no unit, c = 1/sqrt(3) is lattice speed of sound
        // float u_from_We(const float We, const float x, const float sigma, const float rho) const { return sqrt(We*sigma/(x*rho)); } // velocity u = sqrt(We*sigma/(x*rho)) = [m/s], Weber number We = x*u^2*rho/sigma = [1] no unit
        // float u_from_Fr(const float Fr, const float x, const float g) const { return Fr*sqrt(x*g); } // velocity u = Fr*sqrt(x*g) = [m/s], Froude number Fr = u/sqrt(x*g) = [1] no unit
        // float u_from_Ca(const float Ca, const float sigma, const float nu, const float rho) const { return Ca*sigma/(rho*nu); } // velocity u = Ca*sigma/(rho*nu) = [m/s], Capillary number Ca = u*rho*nu/sigma = [1] no unit
        // float u_from_Ca(const float Ca, const float sigma, const float mu) const { return Ca*sigma/mu; } // velocity u = Ca*sigma/(rho*nu) = [m/s], Capillary number Ca = u*mu/sigma = [1] no unit
        // float u_from_f_Poiseuille_2D(const float f, const float rho, const float nu, const float R) const { return f*sq(R)/(2.0f*rho*nu); } // center velocity in 2D Poiseuille flow in channel u = f/(2*rho*nu)*R^2, force per volume f = 2*u*rho*nu/R^2
        // float u_from_f_Poiseuille_3D(const float f, const float rho, const float nu, const float R) const { return f*sq(R)/(4.0f*rho*nu); } // center velocity in 3D Poiseuille flow in channel u = f/(4*rho*nu)*R^2, force per volume f = 4*u*rho*nu/R^2
        // float u_from_f_Poiseuille_2D(const float Q, const float R) const { return 0.75f*Q/R; } // center velocity in 2D Poiseuille flow in channel u = 3/4*Q/R, force per volume f = 2*u*rho*nu/R^2, flow rate Q = 2/3*f*R^3/(rho*nu) = 4/3*u*R
        // //float u_from_f_Poiseuille_3D(const float Q, const float R) const { return 2.0f/pif*Q/sq(R); } // center velocity in 3D Poiseuille flow in channel u = 2/pi*Q/R^2, force per volume f = 4*u*rho*nu/R^2, flow rate Q = pi/8*f*R^4/(rho*nu) = pi/2*R^2*u
        // float f_from_u_Poiseuille_2D(const float u, const float rho, const float nu, const float R) const { return 2.0f*u*rho*nu/sq(R); } // force per volume for 2D Poiseuille flow f = 2*u*rho*nu/R^2, u is center velocity, R is channel radius
        // float f_from_u_Poiseuille_3D(const float u, const float rho, const float nu, const float R) const { return 4.0f*u*rho*nu/sq(R); } // force per volume for 3D Poiseuille flow in cylinder f = 4*u*rho*nu/R^2, u is center velocity, R is cylinder radius
       
        // // other conversions in SI units (can be called before set_m_kg_s(...);)
        // float si_Re(const float Re) const { return Re; } // Reynolds number Re = x*u/nu = [1] no unit
        // float si_Re(const float si_x, const float si_u, const float si_nu) const { return si_x*si_u/si_nu; } // Reynolds number Re = x*u/nu = [1] no unit
        // float si_Re(const float si_x, const float si_u, const float si_mu, const float si_rho) const { return si_x*si_u*si_rho/si_mu; } // Reynolds number Re = x*u*rho/mu = [1] no unit
        // float si_We(const float si_x, const float si_u, const float si_rho, const float si_sigma) const { return si_x*sq(si_u)*si_rho/si_sigma; } // Weber number We = x*u^2*rho/sigma = [1] no unit
        // float si_Fr(const float si_x, const float si_u, const float si_g) const { return si_u/sqrt(si_x*si_g); } // Froude number Fr = u/sqrt(x*g) = [1] no unit
        // float si_Ca(const float si_u, const float si_mu, const float si_sigma) const { return si_u*si_mu/si_sigma; } // Capillary number Ca = u*mu/sigma = [1] no unit
        // float si_Ca(const float si_u, const float si_rho, const float si_nu, const float si_sigma) const { return si_u*si_rho*si_nu/si_sigma; } // Capillary number Ca = u*rho*nu/sigma = [1] no unit
        // float si_Bo(const float si_x, const float si_rho, const float si_g, const float si_sigma) const { return sq(si_x)*si_rho*si_g/si_sigma; } // Bond number Bo = x^2*rho*g/sigma = [1] no unit
        // float si_Mo(const float si_rho, const float si_delta_rho, const float si_g, const float si_sigma, const float si_mu) const { return Mo(si_rho, si_delta_rho, si_g, si_sigma, si_mu); } // Morton number g*delta_rho*mu^4/(sigma^3*rho^2)
        // float si_Ga(const float si_x, const float si_nu, const float si_g) { return cb(si_x)*si_g/sq(si_nu); } // Galilei number Ga = x^3*g/nu^2 = [1] no unit
        // float si_Ga(const float si_x, const float si_rho, const float si_nu, const float si_f) { return cb(si_x)*si_f/(sq(si_nu)*si_rho); } // Galilei number Ga = x^3*g/nu^2 = [1] no unit
        // float si_nu_from_si_mu(const float si_mu, const float si_rho) const { return si_mu/si_rho; } // kinematic shear viscosity nu = mu/rho = [m^2/s]
        // float si_nu_from_si_Re(const float si_Re, const float si_x, const float si_u) const { return si_x*si_u/si_Re; } // kinematic shear viscosity nu = x*u/Re = [m^2/s]
        // float si_mu_from_si_nu(const float si_nu, const float si_rho) const { return si_nu*si_rho; } // dynamic viscosity mu = nu*rho = [kg/(m*s)]
        // float si_f_from_si_g(const float si_g, const float si_rho) const { return si_rho*si_g; } // force per volume f = rho*g = [kg/(m*s)^2]
        // float si_g_from_si_f(const float si_f, const float si_rho) const { return si_f/si_rho; } // force per volume f = rho*g = [kg/(m*s)^2]
        // float si_u_from_si_Re(const float si_Re, const float si_x, const float si_nu) const { return si_Re*si_nu/si_x; } // velocity u = Re*nu/x, Reynolds number Re = x*u/nu = [1] no unit
        // float si_u_from_si_Re(const float si_Re, const float si_x, const float si_mu, const float si_rho) const { return si_Re*si_mu/(si_x*si_rho); } // velocity u = Re*nu/x, Reynolds number Re = x*u*rho/mu = [1] no unit
        // float si_u_from_si_We(const float si_We, const float si_x, const float si_sigma, const float si_rho) const { return sqrt(si_We*si_sigma/(si_x*si_rho)); } // velocity u = sqrt(We*sigma/(x*rho)) = [m/s], Weber number We = x*u^2*rho/sigma = [1] no unit
        // float si_u_from_si_Fr(const float si_Fr, const float si_x, const float si_g) const { return si_Fr*sqrt(si_x*si_g); } // velocity u = Fr*sqrt(x*g) = [m/s], Froude number Fr = u/sqrt(x*g) = [1] no unit
        // float si_u_from_si_h(const float si_h, const float si_g) const { return sqrt(2.0f*si_h*si_g); } // free fall drop height h = 1/2*g*t^2 -> t = sqrt(2*h/g), u = a*t = g*sqrt(2*h/g) = sqrt(2*h*g) = [m/s]
    };

    extern Units_Conv units;
#endif