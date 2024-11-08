// See test_minu_problems.h
// Last edited on 2024-11-08 11:35:43 by stolfi

#define _GNU_SOURCE
#include <values.h>
#include <math.h>

#include <test_minu.h>

#include <test_minu_problems.h>

#define ABS(x) ((x) < 0 ? -(x) : (x))
#define MAX(x,y) ((x) > (y) ? (x) : (y))
#define MIN(x,y) ((x) < (y) ? (x) : (y))

double Log2(double x);

#define QuadR (+1.2)

void QuadEval(void *prb, double x, double *fx, double *dfx)
{
  double u = x - QuadR;
  *fx = u*u;
  *dfx = 2.0 * u;
}

double QuadError (void *prb, double x, double fx)
{
  double dist = ABS(x - QuadR);
  return Log2(dist/(2.0*((Problem *)prb)->tol));
}

#define BentR (+1.2)

void BentEval (void *prb, double x, double *fx, double *dfx)
{
  double a = (1.0 + x * x);
  double aDx = 2.0 * x;
  double u = x - BentR;
  double b = u * u;
  double bDx = 2.0 * u;
  *fx = a * b;
  *dfx = a * bDx + aDx * b;
}

double BentError (void *prb, double x, double fx)
{
  double dist = ABS(x - BentR);
  return Log2(dist/(2.0*((Problem *)prb)->tol));
}

#define CornR (+1.2)

void CornEval (void *prb, double x, double *fx, double *dfx)
{
  double u = x - CornR;
  *fx = 0.5 * u + ABS(u);
  *dfx = (u < 0.0 ? -0.5 : 1.5);
}

double CornError (void *prb, double x, double fx)
{
  double dist = ABS(x - CornR);
  return Log2(dist/(2.0*((Problem *)prb)->tol));
}

#define WavyLambda (2.0)
#define WavyR (+1.2)

void WavyEval (void *prb, double x, double *fx, double *dfx)
{
  double u = x - WavyR;
  double y = u * u;
  double yDx = 2.0*u;

  double W = 2.0 * M_PI / WavyLambda;
  double A = 0.1;
  double phase = W * u;
  double phaseDx = W;
  double m = 1.0 + A * sin(phase);
  double mDx = A * cos(phase) * phaseDx;
  *fx = m * y;
  *dfx = m * yDx + mDx * y;
}

double WavyError (void *prb, double x, double fx)
{
  double dist = ABS(x - WavyR);
  return Log2(dist/(2.0*((Problem *)prb)->tol));
}

#define BiquR (+3.5)
#define BiquS (-2.5)

void BiquEval (void *prb, double x, double *fx, double *dfx)
{
  double u = x - BiquR;
  double v = x - BiquS;
  double z = u * v;
  double zDx = u + v;
  double y = z * z;
  double yDx = 2.0 * z * zDx;
  *fx = y;
  *dfx = yDx;
}

double BiquError (void *prb, double x, double fx)
{
  double d1 = ABS(x - BiquR);
  double d2 = ABS(x - BiquS);
  double dist = MIN(d1, d2);
  return Log2(dist/(2.0*((Problem *)prb)->tol));
}

#define HoleAX (+3.5) 
#define HoleAW (0.05)

#define HoleBX (+5.5)  
#define HoleBW (0.10)

#define HoleCX (-4.5)  
#define HoleCW (0.20)

void HoleEval (void *prb, double x, double *fx, double *dfx)
{
  double ra = (x - HoleAX)/HoleAW;
  double raDx = 1.0/HoleAW;
  double sa = 1.0/(1.0 + ra*ra);
  double ta = 1.0 - sa;
  double taDx = 2.0 * ra * raDx / (sa * sa);

  double rb = (x - HoleBX)/HoleBW;
  double rbDx = 1.0/HoleBW;
  double sb = 1.0/(1.0 + rb*rb);
  double tb = 1.0 - sb;
  double tbDx = 2.0 * rb * rbDx / (sb * sb);

  double rc = (x - HoleCX)/HoleCW;
  double rcDx = 1.0/HoleCW;
  double sc = 1.0/(1.0 + rc*rc);
  double tc = 1.0 - sc;
  double tcDx = 2.0 * rc * rcDx / (sc * sc);
  *fx = ta * tb * tc;
  *dfx = taDx * tb * tc + ta * tbDx * tc + ta * tb * tcDx;
}

double HoleError (void *prb, double x, double fx)
{
  double da = ABS(x - HoleAX);
  double db = ABS(x - HoleBX);
  double dc = ABS(x - HoleCX);
  double dist = MIN(da, MIN(db, dc));
  return Log2(dist/(2.0*((Problem *)prb)->tol));
}

double Log2(double x)
/* 
  Log base 2 of x, clamped at smallpos */
{
  double xx = MAX(MINDOUBLE, ABS(x));
  return log(xx)/log(2.0);
}

