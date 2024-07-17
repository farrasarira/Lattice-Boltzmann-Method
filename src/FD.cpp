#include "math_util.hpp"
#include "defines.hpp"

// Finite-Difference First Order Upwind Scheme
double fd_fuw(double stc_l, double stc_c, double stc_r, double dx, double vel, short type_l, short type_r)
{
    double res = 0.0;
    
    if (vel > 0){
        if (type_l == TYPE_S || type_l == TYPE_A || type_l == TYPE_FS || type_l == TYPE_P)
            res = 0.0;
        else
            res = (stc_c - stc_l) / dx;
    }
    else{
        if (type_r == TYPE_S || type_r == TYPE_A || type_r == TYPE_FS || type_r == TYPE_P)
            res = 0.0;
        else
            res = (stc_r - stc_c) / dx;
    }

    return res;
}

double fd_central(double stc_l, double stc_c, double stc_r, double dx, double vel, short type_l, short type_r)
{
    double res = 0.0;
    
    if (type_l == TYPE_S || type_l == TYPE_A || type_l == TYPE_FS || type_l == TYPE_P)
        res = (stc_r - stc_c) / (dx);
    else if (type_r == TYPE_S || type_r == TYPE_A || type_r == TYPE_FS || type_r == TYPE_P)
        res = (stc_c - stc_l) / (dx);
    else
        res = (stc_r - stc_l) / (2*dx);
    

    return res;
}