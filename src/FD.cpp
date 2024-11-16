#include "math_util.hpp"
#include "defines.hpp"

// Finite-Difference First Order Upwind Scheme
double fd_fuw(double stc_l, double stc_c, double stc_r, double dx, double vel, short type_l, short type_r)
{
    double res = 0.0;
    
    if (vel > 0){
        if (type_l == TYPE_S || type_l == TYPE_A || type_l == TYPE_FS || type_l == TYPE_P || type_l == TYPE_O || type_l == TYPE_O_C || type_l == TYPE_I || type_l == TYPE_I_C)
            res = 0.0;
        else
            res = (stc_c - stc_l) / dx;
    }
    else{
        if (type_r == TYPE_S || type_r == TYPE_A || type_r == TYPE_FS || type_r == TYPE_P || type_r == TYPE_O ||  type_r == TYPE_O_C || type_r == TYPE_I || type_r == TYPE_I_C)
            res = 0.0;
        else
            res = (stc_r - stc_c) / dx;
    }

    return res;
}

double fd_central_2der(double stc_l, double stc_c, double stc_r, double dx, short type_l, short type_r)
{
    double res = 0.0;
    
    if (type_l == TYPE_A || type_l == TYPE_FS || type_l == TYPE_S || type_l == TYPE_P)
        res = 0.0;
    else if (type_r == TYPE_A || type_r == TYPE_FS || type_r == TYPE_S || type_r == TYPE_P)
        res = 0.0;
    else
        res = (stc_r - 2.0*stc_c + stc_l) / (dx*dx);
    
    return res;
}


double fd_laplace(double stc_c, double stc_x_l, double stc_x_r, double stc_y_l, double stc_y_r, double stc_z_l, double stc_z_r, double dx, double dy, double dz, short type_x_l, short type_x_r, short type_y_l, short type_y_r, short type_z_l, short type_z_r)
{
    double res = 0.0;

    double d2phi_dx2 = fd_central_2der(stc_x_l, stc_c, stc_x_r, dx, type_x_l, type_x_r);
    double d2phi_dy2 = fd_central_2der(stc_y_l, stc_c, stc_y_r, dx, type_y_l, type_y_r);
    double d2phi_dz2 = fd_central_2der(stc_z_l, stc_c, stc_z_r, dx, type_z_l, type_z_r);
    res = d2phi_dx2 + d2phi_dy2 + d2phi_dz2; 

    return res;
}

double fd_uw2(double stc_0, double stc_1, double stc_2, double dx, double vel)
{
    double res = 0.0;

    if(vel < 0)
        res = (3.0*stc_0 - 4*stc_1 + stc_2) / (2*dx);
    else
        res = (-3.0*stc_0 + 4*stc_1 - stc_2) / (2*dx);

    return res;
}

double fd_uw(double stc_0, double stc_1, double stc_2, double dx, double vel)
{
    double res = 0.0;

    if(vel < 0)
        res = (stc_0 - stc_1) / (dx);
    else
        res = (-stc_0 + stc_1) / (dx);

    return res;
}

// double fd_central(double stc_l, double stc_c, double stc_r, double dx, double vel, short type_l, short type_r)
// {
//     double res = 0.0;
    
//     if (type_l == TYPE_S || type_l == TYPE_A || type_l == TYPE_FS || type_l == TYPE_P || type_l == TYPE_O || type_l == TYPE_O_C)
//         if (type_r == TYPE_P)
//             res = 0.0;
//         else
//             res = (stc_r - stc_c) / (dx);
//     else if (type_r == TYPE_S || type_r == TYPE_A || type_r == TYPE_FS || type_r == TYPE_P || type_r == TYPE_O || type_r == TYPE_O_C)
//         if (type_l == TYPE_P)
//             res = 0.0;
//         else
//             res = (stc_c - stc_l) / (dx);
//     else
//         res = (stc_r - stc_l) / (2*dx);
    

//     return res;
// }