#ifndef MATH_UTIL
    #define MATH_UTIL

    double dotproduct_Vec3(double a[], double b[]);
    double sq(double a);
    double cb(double a);
    double v_sqr(double u, double v, double w);

    double limiterVanleer(double sL, double sR);
    double limiterMinmod(double sL, double sR);
    double limiterMaxmod(double sL, double sR);
    double limiterSuperbee(double sL, double sR);
    double limiterMC(double sL, double sR);

    float smootherstep(float, float, float);
    double smooth(double left, double right, double x, double center, double alpha);

#endif