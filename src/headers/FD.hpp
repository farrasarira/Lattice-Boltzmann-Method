
#ifndef FD_H

    #define FD_H
    
    // Finite-Difference First Order Upwind Scheme
    double fd_fuw(double stc_l, double stc_c, double stc_r, double dx, double vel, short type_l, short type_r);
    double fd_uw_2(double stc_0, double stc_1, double stc_2, double dx, double vel);
    double fd_central(double stc_l, double stc_c, double stc_r, double dx, double vel, short type_l, short type_r);
    double fd_uw(double stc_0, double stc_1, double stc_2, double dx, double vel);

#endif