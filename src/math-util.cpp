
#include"math_util.hpp"
#include<math.h>
#include<algorithm>

double dotproduct_Vec3(double a[], double b[])
{
    double product = 0.0;
    for(int i = 0; i < 3; ++i)
    {
        product += a[i] * b[i];
    }
    return product;
}

double sq(double a)
{
    return a*a;
}

double cb(double a)
{
    return a*a*a;
}

double v_sqr(double u, double v, double w)
{
    return u*u+v*v+w*w;
}

double limiterVanleer(double sL, double sR)
{
    double slope=(abs(sR)*sL+sR*abs(sL))/(abs(sL)+abs(sR));
    if (isfinite(slope) == false) slope = 0;
    return slope;
}

double limiterMinmod(double sL, double sR)
{
    double slope = 0.0*sL; 
    slope =(signbit(sL)==signbit(sR))*signbit(sL)*std::min(abs(sL),abs(sR));
    return slope;
}

double limiterMaxmod(double sL, double sR)
{
    double slope = 0.0*sL; 
    slope =(signbit(sL)==signbit(sR))*signbit(sL)*std::max(abs(sL),abs(sR));
    return slope;
}

double limiterSuperbee(double sL, double sR)
{
    return limiterMaxmod(limiterMinmod(2.0*sL, sR) , limiterMinmod(sL, 2.0*sR));
}

double limiterMinmodM(double sL, double sM, double sR)
{
    double slope = 0.0*sL;
    slope = ((signbit(sL) == signbit(sM)) && (signbit(sM) == signbit(sR))) * signbit(sL) * std::min({abs(sL), abs(sM), abs(sR)});
    return slope;
}

double limiterMC(double sL, double sR)
{
    return limiterMinmodM(2.0*sL, 0.5*(sL+sR), 2.0*sR);
}

float clamp(float x, float lowerlimit = 0.0f, float upperlimit = 1.0f) 
{
  if (x < lowerlimit) return lowerlimit;
  if (x > upperlimit) return upperlimit;
  return x;
}

float smootherstep(float x, float edge0 = 0.0f, float edge1 = 1.0f) 
{
  // Scale, and clamp x to 0..1 range
  x = clamp((x - edge0) / (edge1 - edge0));

  return x * x * x * (3.0f * x * (2.0f * x - 5.0f) + 10.0f);
}

double smooth(double left, double right, double x, double center, double alpha) {
    return left +  0.5 * (1 + std::tanh(alpha * (x - center))) * (right-left);
}

// calculate ratio of consecutive slopes (theta)
double calc_ratio_slopes(double stc_n, double stc_c, double stc_p)
{
    if (stc_p == stc_c)
        return 0;
    else
        return (stc_c - stc_n)/(stc_p - stc_c);
}

double limiterMinmod(double theta)
{
    return std::max(0.0, std::min(1.0, theta));
}

double limiterVanleer(double theta)
{
    return (theta + abs(theta)) / (1 + abs(theta)); 
}

double FD_limiterMinmod(double stc_n, double stc_c, double stc_p, double dx)
{
    double theta = calc_ratio_slopes(stc_n, stc_c, stc_p);
    double phi = limiterMinmod(theta);
    return phi * (stc_p - stc_c) / dx;
}

double FD_limiterVanleer(double stc_n, double stc_c, double stc_p, double dx)
{
    double theta = calc_ratio_slopes(stc_n, stc_c, stc_p);
    double phi = limiterVanleer(theta);
    return -1.0 * phi * (stc_p - stc_c) / dx;
}