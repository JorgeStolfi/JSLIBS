/* See sym_eigen.h */
/* Last edited on 2024-11-20 22:49:44 by stolfi */

#include <stdlib.h>
#include <stdio.h>
#include <stdint.h>
#include <assert.h>
#include <math.h>

#include <affirm.h>
#include <bool.h>

#include <sym_eigen.h>

#define Pr fprintf
#define Er stderr
  /* Shoter names. */

void sym_eigen_tridiagonalize
  ( uint32_t n,
    double A[],
    double d[],
    double e[],
    double R[]
  )
  {
    /* The algorithm is based on the EISPACK routines "tred1.f" and
      "tred2.f" by Burton S. Garbow, Argonne National Laboratory
      (aug/1983); originally from the ALGOL procedures {tred1} and
      {tred2}, Num. Math. 11, 181-195(1968) by Martin, Reinsch, and
      Wilkinson. See Handbook for Auto. Comp., vol.II-Linear Algebra,
      212-226(1971). Re-implemented in C by Jorge Stolfi, Unicamp
      (dec/2002). */
      
    demand(R != NULL, "{R} must be non-null");

    double f, g, h, hh, scale;
    for (int32_t i = 0; i < n; i++)
     { uint32_t rowi = n*i;
       uint32_t rown1 = n*(n-1);
       d[i] = A[rown1+i]; 
       for (int32_t j = i; j < n; j++) 
         { uint32_t rowj = n*j; R[rowi+j] = A[rowj+i]; }
     }
    for (int32_t ii = (int32_t)n-1; ii >= 2; ii--)
      { assert(ii >= 1);
        uint32_t i = (uint32_t)ii;
        h = 0.0;
        /* Scale row: */
        scale = 0.0;
        for (int32_t k = 0; k < i; k++) { scale = scale + fabs(d[k]); }
        if (scale == 0.0)
          { e[i] = d[i-1];
            for (int32_t j = 0; j < i; j++)
              { double *Rji = &(R[n*j+i]);
                double *Rji1 = &(R[n*j+(i-1)]);
                d[j] = *Rji1; 
                R[n*i+j] = 0.0;
                *Rji = 0.0;
              }
          }
        else
          { for (int32_t k = 0; k < i; k++) { d[k] /= scale; h = h + d[k]*d[k]; }
            f = d[i-1];
            g = - copysign(sqrt(h), f);
            e[i] = scale * g;
            h = h - f * g;
            d[i-1] = f - g;
            /* form A*u: */
            for (int32_t j = 0; j < i; j++) { e[j] = 0.0; }
            for (int32_t j = 0; j < i; j++)
              { uint32_t rowj = n*j;
                double *Rjj = &(R[rowj+j]);
                f = d[j];
                g = e[j] + (*Rjj) * f;
                for (int32_t k = j+1; k < i; k++)
                  { double Rjk = R[n*j+k];
                    g += Rjk*d[k]; e[k] += Rjk*f;
                  }
                e[j] = g;
              }
            /* Compute {p}: */
            f = 0.0;
            for (int32_t j = 0; j < i; j++) { e[j] /=  h; f += e[j]*d[j]; }
            hh = f / (h + h);
            /* Compute {q}: */
            for (int32_t j = 0; j < i; j++) { e[j] -= hh*d[j]; }
            /* Compute the reduced {A}: */
            for (int32_t j = 0; j < i; j++)
              { uint32_t rowj = n*j;
                f = d[j]; g = e[j];
                for (int32_t k = j; k < i; k++)
                  { double *Rjk = &(R[rowj+k]);
                    (*Rjk) -= (f*e[k] + g*d[k]);
                  }
              }
          }
      }
      { uint32_t rown1 = n*(n-1);
        uint32_t row0 = n*0;
        uint32_t row1 = n*1;
        e[1] = d[0];
        /* Accumulate of transformation matrices: */
        d[0] = R[row0+0]; R[row0+1] = 0.0; R[row1+0] = 0.0;
        d[1] = 0.0;
        for (int32_t i = 1; i < n; i++)
          { uint32_t rowi = n*i;
            uint32_t rowi1 = n*(i-1);
            R[rowi1+n-1] = R[rowi1+i-1]; R[rowi1+i-1] = 1.0;
            h = d[i];
            if (h != 0.0)
              { for (int32_t k = 0; k < i; k++) { d[k] = R[rowi+k] / h; }
                for (int32_t j = 0; j < i; j++)
                  { uint32_t rowj = n*j;
                    g = 0.0;
                    for (int32_t k = 0; k < i; k++) { g += R[rowi+k]*R[rowj+k]; }
                    for (int32_t k = 0; k < i; k++) { R[rowj+k] -= g*d[k]; }
                  }
              }
            for (int32_t k = 0; k < i; k++) { R[rowi+k] = 0.0; }
          }
        for (int32_t i = 0; i < n; i++) 
          { uint32_t rowi = n*i;
            d[i] = R[rowi+n-1]; R[rowi+n-1] = 0.0;
          }
        R[rown1+(n-1)] = 1.0;
      }
    e[0] = 0.0;
  }

void sym_eigen_trid_eigen(uint32_t n, double *d, double *e, double *R, uint32_t *p, uint32_t absrt)
  {
    /* The algorithm is based on the EISPACK routines "tql1.f" and
      "tql2.f" by Burton S. Garbow, Argonne National Laboratory
      (aug/1983); orignally from the ALGOL procedures {tql1} and
      {tql2}, Num. Math. 11, 293-306(1968) by Bowdler, Martin,
      Reinsch, and Wilkinson. See Handbook for Auto. Comp., vol.II -
      Linear Algebra, 227-240(1971). Re-implemented in C by Jorge
      Stolfi, Unicamp (dec/2002). */
    
    bool_t debug = TRUE;
    
    if (debug) { Pr(Er, "  > --- %s  n = %d ---\n", __FUNCTION__, n); }
    
    double f, magn;
    uint32_t maxiter = 30;
    
    if (n == 1) { *p = 1; return; }
    /* Downshift {e} one slot (so {e[i] = T[i+1,i]}), and clear {e[n-1]}: */
    for (int32_t i = 1; i < n; i++) { e[i-1] = e[i]; }
    e[n-1] = 0.0;

    f = 0.0;
    magn = 0.0;
    int32_t L = 0;
    while(L < n)
      { if (debug) { Pr(Er, "    ...... L = %d  f = %24.16e  magn = %24.16e ......\n", L, f, magn); }
        uint32_t niter = 0;
        { double w = fabs(d[L]) + fabs(e[L]); if (magn < w) { magn = w; } }
        double test = magn + fabs(e[L]);
        if (debug) { Pr(Er, "      magn = %24.16e  e[L] = %24.16f  test1 = %24.16f\n", magn, e[L], test); }
        if (test != magn)
          { /* Can't be at the sentinal {e[n-1]}. */
            assert(L < n-1); 
            /* Must clear {e[L]} */
            /* look for small sub-diagonal element {e[t]} in cols {L+1..n-1}; */
            /* Note that it will always find {e[n-1] == 0} as last resort. */
            int32_t t = L;
            do { t++; } while ((magn + fabs(e[t])) != magn);
            if (debug) { Pr(Er, "      (0) t = %d\n", t); }
            if (debug) { Pr(Er, "      magn = %24.16e  e[t] = %24.16f  test = %24.16f\n", magn, e[t], test); }
            assert((L >= 0) && (t > L) && (t < n));
            /* Try to anihilate {e[L]}: */
            do
              { /* Compute shift: */
                int32_t L1 = L+1;
                double g = d[L];
                double pp = (d[L1] - g)/(2.0*e[L]);
                double r = hypot(pp, 1.0);
                double s2 = 0, c3 = 0;
                d[L] = e[L]/(pp + copysign(r, pp));
                d[L1] = e[L]*(pp + copysign(r, pp));
                double h = g - d[L];
                if (debug) { Pr(Er, "        (1) h =  %24.16e  d[%d] = %24.16e  d[%d] = %24.16e\n", h, L, d[L], L1, d[L1]); }
                for (int32_t i = L + 2; i < n; i++) { d[i] -= h; }
                f = f + h;
                /* QL transformation: */
                pp = d[t];
                double c2 = 1.0;
                double c = 1.0; 
                double s = 0.0;
                double dL1 = d[L1]; 
                double eL1 = e[L1];
                if (debug) { Pr(Er, "        (2) pp = %24.16e  e[%d] = %24.16e\n", pp, L1, e[L1]); }
                /* Note {int} not {uint} in case loop ends at -1. */
                for (int32_t i = ((int32_t)t)-1; i >= L; i--)
                  { int32_t rowi = (int32_t)n*i;
                    int32_t rowi1 = (int32_t)n*(i+1);
                    if (debug) { Pr(Er, "          (3) i = %d  rowi = %d rowi1 = %d\n", i, rowi, rowi1); }
                    c3 = c2; 
                    c2 = c; s2 = s;
                    g = c*e[i];
                    h = c*pp;
                    r = hypot(pp, e[i]);
                    e[i+1] = s*r;
                    s = e[i]/r; c = pp/r;
                    pp = c*d[i] - s*g;
                    d[i+1] = h + s*(c*g + s*d[i]);
                    if (debug) { Pr(Er, "          (4) i = %d  d[%d] = %24.16e\n", i, i+1, d[i+1]); }
                    /* Update eigenvectors: */
                    for (int32_t j = 0; j < n; j++)
                      { double R0 = R[rowi + j], R1 = R[rowi1 + j];
                        R[rowi  + j] = c*R0 - s*R1;
                        R[rowi1 + j] = s*R0 + c*R1;
                      }
                  }
                if (debug) { Pr(Er, "        (5) eL1 = %24.16e  e[%d] = %24.16f  dL1 = %24.16f\n", eL1, L, e[L], dL1); }
                pp = - s*s2*c3*eL1*e[L]/dL1;
                e[L] = s*pp; 
                d[L] = c*pp;
                if (debug) { Pr(Er, "        (6) e[%d] = %24.16f  d[%d] = %24.16f\n", L, e[L], L, d[L]); }
                niter++;
                test = magn + fabs(e[L]);
                if (debug) { Pr(Er, "        (7) magn = %24.16f  test = %24.16f\n", magn, test); }
              }
            while ((niter < maxiter) && (test > magn));
          }
        if (test != magn) { /* Failed to converge: */ break; }
        d[L] = d[L] + f;
        if (debug) { Pr(Er, "      (8) d[%d] = %24.16f\n", L, d[L]); }
        L++;
      }
    /* We did what we could: */
    assert(L >= 0);
    *p = (uint32_t)L;
    if (debug) { Pr(Er, "      (9) *p = %d\n", *p); }
    /* Sort the eigenvalues {d[0..L-1]}, carrying along their eigenvectors: */
    for (int32_t i = 0; i < L; i++)
      { /* Find {(i+1)}th smallest eigenvalue {dk = d[k]} */
        uint32_t k = i; 
        double dk = d[i];
        if (absrt) { dk = fabs(dk); }
        for (int32_t j = i+1; j < L; j++) 
          { double dj = d[j];
            if (absrt) { dj = fabs(d[j]); }
            if (dj < dk) { k = j; dk = dj; }
          }
        /* Swap with {d[i], R[i,*]} if necessary: */
        if (k != i)
          { d[k] = d[i]; d[i] = dk;
            uint32_t rowi = n*i;
            uint32_t rowk = n*k;
            for (int32_t j = 0; j < n; j++)
              { double Rij = R[rowi+j]; R[rowi+j] = R[rowk+j]; R[rowk+j] = Rij; }
          }
      }
  }
