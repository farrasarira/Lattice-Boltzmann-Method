#ifndef MATH_UTIL
    #define MATH_UTIL

    double dotproduct_Vec3(double a[], double b[]);
    double sq(double a);
    double cb(double a);
    double v_sqr(double u, double v, double w);
    double v_mag(double u, double v, double w);

    double limiterVanleer(double sL, double sR);
    double limiterMinmod(double sL, double sR);
    double limiterMaxmod(double sL, double sR);
    double limiterSuperbee(double sL, double sR);
    double limiterMC(double sL, double sR);

    float smootherstep(float, float, float);
    double smooth(double left, double right, double x, double center, double alpha);
    double smooth2D(double inside, double outside, double x, double y, double core_rad, double alpha);

    double calc_ratio_slopes(double stc_n, double stc_c, double stc_p);
    double limiterMinmod(double theta);
    double limiterVanleer(double theta);
    double FD_limiterMinmod(double stc_n, double stc_c, double stc_p, double dx);
    double FD_limiterVanleer(double stc_n, double stc_c, double stc_p, double dx);

    bool onlyOne_1(double a, double b, double c);

#endif