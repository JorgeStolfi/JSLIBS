#include <stdio.h>
#include <math.h>
#include <values.h>
main(){ int e; double m = DBL_MAX; double g = log(m); double f = frexp(m, &e); double c = 1+2*(1-f); double s = sqrt(m); int z0 = sizeof(float), z1 = sizeof(double), z2 = sizeof(long double); short int t = 0; long double N = 1; long double D = 1-f; long double Q = D/1024; long double P = N+Q; long double R = P - N; printf("m/8 = %25.17le 8/m = %25.17le g = %25.17le s = %25.17le 1/s = %25.17le c = %25.17le (f,e) = %25.17le * 2^{%d} sizes = %d %d %d bytes t = %hd  N = %50.42Le D = %50.42Le Q = %50.42Le P = %50.42Le R = %50.42Le\n",m/8,8/m,g,s,1/s,c,f,e,z0,z1,z2,t,N,D,Q,P,R); }
