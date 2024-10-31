
#ifndef FD_H

    #define FD_H
    
    // Finite-Difference First Order Upwind Scheme
    double fd_fuw(double stc_l, double stc_c, double stc_r, double dx, double vel, short type_l, short type_r);
    double fd_central(double stc_l, double stc_c, double stc_r, double dx, double vel, short type_l, short type_r);
    double fd_central_2der(double stc_l, double stc_c, double stc_r, double dx, short type_l, short type_r);
    double fd_laplace(double stc_c, double stc_x_l, double stc_x_r, double stc_y_l, double stc_y_r, double stc_z_l, double stc_z_r, double dx, double dy, double dz, short type_x_l, short type_x_r, short type_y_l, short type_y_r, short type_z_l, short type_z_r);
    double fd_uw2(double stc_0, double stc_1, double stc_2, double dx, double vel);
    double fd_uw(double stc_0, double stc_1, double stc_2, double dx, double vel);

#endif