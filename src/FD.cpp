#include "math_util.hpp"
#include "defines.hpp"

// Finite-Difference First Order Upwind Scheme
double fd_fuw(double stc_l, double stc_c, double stc_r, double dx, double vel, short type_l, short type_r)
{
    double res = 0.0;
    
    if (vel > 0){
        if (type_l == TYPE_S || type_l == TYPE_A || type_l == TYPE_FS || type_l == TYPE_P || type_l == TYPE_O_E|| type_l == TYPE_I_E)
            res = 0.0;
        else
            res = (stc_c - stc_l) / dx;
    }
    else{
        if (type_r == TYPE_S || type_r == TYPE_A || type_r == TYPE_FS || type_r == TYPE_P || type_r == TYPE_O_E|| type_r == TYPE_I_E)
            res = 0.0;
        else
            res = (stc_r - stc_c) / dx;
    }

    return res;
}

double fd_uw_2(double stc_0, double stc_1, double stc_2, double dx, double vel)
{
    double res = 0.0;

    if(vel >= 0)
        res = (3.0*stc_0 - 4.0*stc_1 + 1.0*stc_2) / (2*dx);
    else
        res = (-3.0*stc_0 + 4.0*stc_1 - 1.0*stc_2) / (2*dx);

    return res;
}

double fd_uw(double stc_0, double stc_1, double stc_2, double dx, double vel)
{
    double res = 0.0;

    if(vel >= 0)
        res = (3.0*stc_0 - 4*stc_1 + stc_2) / (2*dx);
    else
        res = (-3.0*stc_0 + 4*stc_1 - stc_2) / (2*dx);

    return res;
}

double fd_central(double stc_l, double stc_c, double stc_r, double dx, double vel, short type_l, short type_r)
{
    double res = 0.0;
    
    if (type_l == TYPE_S || type_l == TYPE_A || type_l == TYPE_FS || type_l == TYPE_P || type_l == TYPE_O_E|| type_l == TYPE_I_E)
        res = (stc_r - stc_c) / (dx);
    else if (type_r == TYPE_S || type_r == TYPE_A || type_r == TYPE_FS || type_r == TYPE_P || type_r == TYPE_O_E|| type_r == TYPE_I_E)
        res = (stc_c - stc_l) / (dx);
    else
        res = (stc_r - stc_l) / (2*dx);
    

    return res;
}