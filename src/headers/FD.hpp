
#ifndef FD_H

    #define FD_H
    
    // Finite-Difference First Order Upwind Scheme
    double fd_fuw(double stc_l, double stc_c, double stc_r, double dx, double vel, short type_l, short type_r);
    double fd_central_2der(double stc_l, double stc_c, double stc_r, double dx, short type_l, short type_r);
    double fd_laplace(double stc_c, double stc_x_l, double stc_x_r, double stc_y_l, double stc_y_r, double stc_z_l, double stc_z_r, double dx, double dy, double dz, short type_x_l, short type_x_r, short type_y_l, short type_y_r, short type_z_l, short type_z_r);
    double fd_uw2(double stc_0, double stc_1, double stc_2, double dx, double vel);
    double fd_uw(double stc_0, double stc_1, double stc_2, double dx, double vel);

    // double fd_central(double stc_l, double stc_c, double stc_r, double dx, double vel, short type_l, short type_r);
    inline double fd_central(double stc_l, double stc_c, double stc_r, double dx, double vel, short type_l, short type_r)
    {
        double res = 0.0;
        
        if (type_l == TYPE_S || type_l == TYPE_A || type_l == TYPE_FS || type_l == TYPE_P || type_l == TYPE_O || type_l == TYPE_O_C)
            if (type_r == TYPE_P)
                res = 0.0;
            else
                res = (stc_r - stc_c) / (dx);
        else if (type_r == TYPE_S || type_r == TYPE_A || type_r == TYPE_FS || type_r == TYPE_P || type_r == TYPE_O || type_r == TYPE_O_C)
            if (type_l == TYPE_P)
                res = 0.0;
            else
                res = (stc_c - stc_l) / (dx);
        else
            res = (stc_r - stc_l) / (2*dx);
        

        return res;
    }

#endif