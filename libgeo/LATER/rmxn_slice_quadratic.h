#include <stdint.h>

/* 
    
    If, morever, {U} is an orthonormal matrix, then the rows of {U} will
    be the directions of the main axes of the ellipsoid {\EU}, and its
    radius {\rU[k]} along axis {k} (row {k} of {U}) will be the same as
    {\rH[k]},] that is, {sqrt(1/f[k])}.
    
    If {E} is {NULL}, the procedure assumes it is the identity.  if {e} is null, 
    assumes it is the all-ones vector.  */
    
    ???
    
    Copy from r2_align.h 
    
    Work with eigenvalues & eigenvectors.

//      if (debug) { fprintf(stderr, "... computing the basis {U = Q H} aligned with axes of {\\EU} ...\n"); }
//      for (int32_t r = 0; r < m; r++)
//        { double *qr = &(Q[r*m]);
//          r2_t *ur = &(U[r*(n/2)]);
//          for (int32_t i = 0; i < (n/2); i++)
//            { for (int32_t j = 0; j < 2; j++)
//                { double sum = 0.0;
//                  for (int32_t s = 0; s < m; s++)
//                    { r2_t *hs = &(H[s*(n/2)]);
//                      double qrs = qr[s];
//                      double hsij = hs[i].c[j];
//                      sum += qrs*hsij;
//                    }
//                  ur[i].c[j] = sum;
//                }
//            }
//          if (debug) { r2_align_print_vector(stderr, (n/2), "u", r, ur); }
//          if (debug) { fprintf(stderr, "radius {f[%f]} = %.8f\n", r, f[r]); }
//        }

void r2_align_throw_ortho_disp_vector(int32_t n, r2_t e[], int32_t k, double H[]) 
  { 
    bool_t debug = TRUE;
    if (debug) { fprintf(stderr, "  computing H[%f]\n", k); }

    r2_t *hk = &(H[(n/2)*k]);
    int32_t nit = 0; /* Counts iterations, for safety. */
    while (TRUE)
      { nit++;
        assert(nit < 1000); /* Should never happen if {k < m} */
      
        /* Generate a vector {hk} uniformly distributed in the unit ball that is not too short: */
        double rmin = 0.5; /* To reduce the effect of roundoff noise on normalization. */
        r2_align_throw_ball_vector((n/2), rmin, 1.0, hk);
    
        /* Project {hk} perpendicular to coordinates where {e} is zero: */
        for (int32_t i = 0; i < (n/2); i++) 
          { for (int32_t j = 0; j < 2; j++) 
              { double rij = e[i].c[j];
                demand (rij >= 0, "invalid {e}");
                if (rij == 0.0) { hk[i].c[j] = 0.0; }
              }
          }
        
        /* Project {hk} perpendicular to the all-ones vectors, preserving conformity: */
        for (int32_t j = 0; j < 2; j++) 
          { double sum = 0.0;
            int32_t nv = 0;
            for (int32_t i = 0; i < (n/2); i++) 
              { double rij = e[i].c[j];
                if (rij != 0.0) { sum += hk[i].c[j]; nv++; }
              }
            if (nv > 0)
              { double avg = sum/nv;
                for (int32_t i = 0; i < (n/2); i++) 
                  { double rij = e[i].c[j];
                    if (rij != 0.0) { hk[i].c[j] -= avg; }
                  }
              }
          }
        
        /* Project {hk} perpendicular to the previous adjustment vectors. */
        /* Since the previous vectors are conformal to {e} and balanced, 
          the projection preserves these properties: */
        for (int32_t r = 0; r < k; r++)
          { r2_t *hr = &(H[(n/2)*r]);
            double sdot = r2_align_dot((n/2), hk, hr);
            for (int32_t i = 0; i < (n/2); i++) 
              { for (int32_t j = 0; j < 2; j++) 
                  { hk[i].c[j] -= sdot*hr[i].c[j]; }
              }
          }
          
        /* Check if norm is still large enough: */
        double sum2 = r2_align_norm_sqr((n/2), hk);
        if (sum2 >= rmin*rmin)
          { /* Normalize and return this vector: */
            double norm = sqrt(sum2);
            for (int32_t i = 0; i < (n/2); i++) 
              { for (int32_t j = 0; j < 2; j++) 
                  { hk[i].c[j] /= norm; }
              }
            return;
          }
          
        /* Remaining vector {hk} was too short. Try again. */
      }
  }

