/* See sym_eigen.h */
/* Last edited on 2021-06-09 20:36:04 by jstolfi */

#define _GNU_SOURCE
#include <stdlib.h>
#include <stdio.h>
#include <stdint.h>
#include <math.h>

#include <affirm.h>

#include <sym_eigen.h>

void syei_tridiagonalize(int32_t n, double *A, double *d, double *e, double *R)
  {
    /* The algorithm is based on the EISPACK routines "tred1.f" and
      "tred2.f" by Burton S. Garbow, Argonne National Laboratory
      (aug/1983); originally from the ALGOL procedures {tred1} and
      {tred2}, Num. Math. 11, 181-195(1968) by Martin, Reinsch, and
      Wilkinson. See Handbook for Auto. Comp., vol.II-Linear Algebra,
      212-226(1971). Re-implemented in C by Jorge Stolfi, Unicamp
      (dec/2002). */

    int32_t i, j, k;
    double f, g, h, hh, scale;
    for (i = 0; i < n; i++)
     { int32_t rowi = n*i;
       int32_t rown1 = n*(n-1);
       d[i] = A[rown1+i]; 
       if (R == NULL)
         { A[rown1+i] = A[rowi+i]; }
       else
         { for (j = i; j < n; j++) 
             { int32_t rowj = n*j; R[rowi+j] = A[rowj+i]; }
         }
     }
    for (i = n-1; i >= 2; i--)
      { h = 0.0;
        /* Scale row: */
        scale = 0.0;
        for (k = 0; k < i; k++) { scale = scale + fabs(d[k]); }
        if (scale == 0.0)
          { e[i] = (R == NULL ? 0.0 : d[i-1]);
            for (j = 0; j < i; j++)
              { double *Rji = (R == NULL ? &(A[n*i+j]) : &(R[n*j+i]));
                double *Rji1 = (R == NULL ? &(A[n*(i-1)+j]) : &(R[n*j+(i-1)]));
                d[j] = *Rji1; 
                if (R == NULL) { *Rji1 = *Rji; } else { R[n*i+j] = 0.0; }
                *Rji = 0.0;
              }
          }
        else
          { int32_t rowi = n*i;
            int32_t rowi1 = n*(i-1);
            for (k = 0; k < i; k++) { d[k] /= scale; h = h + d[k]*d[k]; }
            f = d[i-1];
            g = - copysign(sqrt(h), f);
            e[i] = scale * g;
            h = h - f * g;
            d[i-1] = f - g;
            /* form A*u: */
            for (j = 0; j < i; j++) { e[j] = 0.0; }
            for (j = 0; j < i; j++)
              { int32_t rowj = n*j;
                double *Rjj = (R == NULL ? &(A[rowj+j]) : &(R[rowj+j]));
                f = d[j];
                if (R != NULL) { R[rowi+j] = f; }
                g = e[j] + (*Rjj) * f;
                for (k = j+1; k < i; k++)
                  { double Rjk = (R == NULL ? A[n*k+j] : R[n*j+k]);
                    g += Rjk*d[k]; e[k] += Rjk*f;
                  }
                e[j] = g;
              }
            /* Compute {p}: */
            f = 0.0;
            for (j = 0; j < i; j++) { e[j] /=  h; f += e[j]*d[j]; }
            hh = f / (h + h);
            /* Compute {q}: */
            for (j = 0; j < i; j++) { e[j] -= hh*d[j]; }
            /* Compute the reduced {A}: */
            for (j = 0; j < i; j++)
              { int32_t rowj = n*j;
                f = d[j]; g = e[j];
                for (k = j; k < i; k++)
                  { double *Rjk = (R == NULL ? &(A[n*k+j]) : &(R[rowj+k]));
                    (*Rjk) -= (f*e[k] + g*d[k]);
                  }
                if (R != NULL) { d[j] = R[rowj+i-1]; R[rowj+i] = 0.0; }
              }
            if (R == NULL)
              { for (j = 0; j < i; j++)
                  { f = d[j];
                    d[j] = A[rowi1+j]; A[rowi1+j] = A[rowi+j]; A[rowi+j] = f*scale;
                  }
              }
          }
        if (R != NULL) { d[i] = h; }
      }
    if (R == NULL)
      { int32_t row0 = n*0;
        int32_t row1 = n*1;
        e[1] = -d[0];
        d[0] = A[row0+0]; A[row0+0] = A[row1+0]; A[row1+0] = -2.0*e[1];
      }
    else
      { int32_t rown1 = n*(n-1);
        int32_t row0 = n*0;
        int32_t row1 = n*1;
        e[1] = d[0];
        /* Accumulate of transformation matrices: */
        d[0] = R[row0+0]; R[row0+1] = 0.0; R[row1+0] = 0.0;
        d[1] = 0.0;
        for (i = 1; i < n; i++)
          { int32_t rowi = n*i;
            int32_t rowi1 = n*(i-1);
            R[rowi1+n-1] = R[rowi1+i-1]; R[rowi1+i-1] = 1.0;
            h = d[i];
            if (h != 0.0)
              { for (k = 0; k < i; k++) { d[k] = R[rowi+k] / h; }
                for (j = 0; j < i; j++)
                  { int32_t rowj = n*j;
                    g = 0.0;
                    for (k = 0; k < i; k++) { g += R[rowi+k]*R[rowj+k]; }
                    for (k = 0; k < i; k++) { R[rowj+k] -= g*d[k]; }
                  }
              }
            for (k = 0; k < i; k++) { R[rowi+k] = 0.0; }
          }
        for (i = 0; i < n; i++) 
          { int32_t rowi = n*i;
            d[i] = R[rowi+n-1]; R[rowi+n-1] = 0.0;
          }
        R[rown1+(n-1)] = 1.0;
      }
    e[0] = 0.0;
  }

void syei_trid_eigen(int32_t n, double *d, double *e, double *R, int32_t *p, int32_t absrt)
  {
    /* The algorithm is based on the EISPACK routines "tql1.f" and
      "tql2.f" by Burton S. Garbow, Argonne National Laboratory
      (aug/1983); orignally from the ALGOL procedures {tql1} and
      {tql2}, Num. Math. 11, 293-306(1968) by Bowdler, Martin,
      Reinsch, and Wilkinson. See Handbook for Auto. Comp., vol.II -
      Linear Algebra, 227-240(1971). Re-implemented in C by Jorge
      Stolfi, Unicamp (dec/2002). */
    
    int32_t i, j, t, L;
    double f, magn;
    int32_t maxiter = 30;
    
    if (n == 1) { *p = 1; return; }
    /* Downshift {e} one slot (so {e[i] = T[i+1,i]}), and clear {e[n-1]}: */
    for (i = 1; i < n; i++) { e[i-1] = e[i]; }
    e[n-1] = 0.0;

    f = 0.0;
    magn = 0.0;
    for (L = 0; L < n; L++)
      { int32_t niter = 0;
        double test;
        { double w = fabs(d[L]) + fabs(e[L]); if (magn < w) { magn = w; } }
        test = magn + fabs(e[L]);
        if (test != magn)
          { /* Must clear {e[L]} */
            affirm(L < n-1, "missed sentinel e[n-1]");
            /* look for small sub-diagonal element {e[t]} in cols {L+1..n-1}; */
            /* Note that it will always find {e[n-1] == 0} as last resort. */
            t = L;
            do { t++; test = magn + fabs(e[t]); } while (test != magn);
            /* Try to anihilate {e[L]}: */
            do
              { /* Compute shift: */
                int32_t L1 = L+1;
                double g = d[L];
                double pp = (d[L1] - g)/(2.0*e[L]);
                double r = hypot(pp, 1.0);
                double h, dL1, eL1;
                double c, s, c2, s2 = 0, c3 = 0;
                d[L] = e[L]/(pp + copysign(r, pp));
                d[L1] = e[L]*(pp + copysign(r, pp));
                h = g - d[L];
                for (i = L + 2; i < n; i++) { d[i] -= h; }
                f = f + h;
                /* QL transformation: */
                pp = d[t];
                c2 = 1.0;
                c = 1.0; s = 0.0;
                dL1 = d[L1]; eL1 = e[L1];
                for (i = t-1; i >= L; i--)
                  { int32_t rowi = n*i;
                    int32_t rowi1 = n*(i+1);
                    c3 = c2; 
                    c2 = c; s2 = s;
                    g = c*e[i];
                    h = c*pp;
                    r = hypot(pp, e[i]);
                    e[i+1] = s*r;
                    s = e[i]/r; c = pp/r;
                    pp = c*d[i] - s*g;
                    d[i+1] = h + s*(c*g + s*d[i]);
                    if (R != NULL)
                      { /* Update eigenvectors: */
                        for (j = 0; j < n; j++)
                          { double R0 = R[rowi + j], R1 = R[rowi1 + j];
                            R[rowi  + j] = c*R0 - s*R1;
                            R[rowi1 + j] = s*R0 + c*R1;
                          }
                      }
                  }
                pp = - s*s2*c3*eL1*e[L]/dL1;
                e[L] = s*pp; d[L] = c*pp;
                niter++;
                test = magn + fabs(e[L]);
              }
            while ((niter < maxiter) && (test > magn));
          }
        if (test != magn) { /* Failed to converge: */ break; }
        d[L] = d[L] + f;
      }
    /* We did what we could: */
    *p = L;
    /* Sort the eigenvalues {d[0..L-1]}, carrying along their eigenvectors: */
    for (i = 0; i < L; i++)
      { /* Find {(i+1)}th smallest eigenvalue {dk = d[k]} */
        int32_t k = i; 
        double dk = d[i];
        if (absrt) { dk = fabs(dk); }
        for (j = i+1; j < L; j++) 
          { double dj = d[j];
            if (absrt) { dj = fabs(d[j]); }
            if (dj < dk) { k = j; dk = dj; }
          }
        /* Swap with {d[i], R[i,*]} if necessary: */
        if (k != i)
          { d[k] = d[i]; d[i] = dk;
            if (R != NULL)
              { int32_t rowi = n*i;
                int32_t rowk = n*k;
                for (j = 0; j < n; j++)
                  { double Rij = R[rowi+j]; R[rowi+j] = R[rowk+j]; R[rowk+j] = Rij; }
              }
          }
      }
  }
