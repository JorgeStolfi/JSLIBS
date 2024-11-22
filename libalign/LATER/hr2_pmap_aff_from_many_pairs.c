/* See hr2_pmap_aff_from_point_pairs.h */
/* Last edited on 2024-11-15 13:00:54 by stolfi */ 

/* Based on HR2.m3 created 1994-05-04 by J. Stolfi. */

#define _GNU_SOURCE
#include <stdint.h>
#include <math.h>
#include <assert.h>

#include <r2.h>
#include <r3.h>
#include <rn.h>
#include <r3x3.h>
#include <r2x2.h>
#include <affirm.h>
#include <sign.h>
#include <sign_get.h>
#include <hr2.h>

#include <hr2_pmap_aff_from_point_pairs.h>

#define NH 3
  /* Number of homogeneous coordinates in a point. */

#define NC 2
  /* Number of Cartesian coordinates in a point. */

hr2_pmap_t hr2_pmap_aff_from_point_pairs(int32_t np, r2_t p1[], r2_t p2[], double w[])
  {
    double debug = FALSE;
    
    if (debug) { fprintf(stderr, "--- computing the affine matrix ---\n"); }
    
    r3x3_t A;
    r3x3_ident(&A);
    if (np > 0)
      { /* Compute the barycenters of {p1} and {p2}: */
        r2_t bar1; r2_barycenter(np, p1, w, &bar1);
        if (debug) { r2_gen_print(stderr, &bar1, "%10.4f", "  bar1 = ( ", " ", " )\n"); }
        demand(r2_is_finite(&bar1), "barycenter of {p1} is undefined");

        r2_t bar2; r2_barycenter(np, p2, w, &bar2);
        if (debug) { r2_gen_print(stderr, &bar2, "%10.4f", "  bar2 = ( ", " ", " )\n"); }
        demand(r2_is_finite(&bar2), "barycenter of {p2} is undefined");

        r2_t d; /* Displacement vector */
        
        if (np == 1)
          { /* Just translate {bar1} to {bar2}: */
            r2_sub(&bar2, &bar1, &d);
          }
        else
          { /* Determine the linear map matrix {L = A.c[1..2][1..2]}: */
            r2x2_t L;
            if (np == 2)
              { /* Compute {L} by composing rotation & scale matrices: */
                r2_t q1; r2_sub(&(p1[0]), &bar1, &q1);
                r2x2_t R1; r2x2_rot_and_scale(&q1, &R1);

                r2_t q2; r2_sub(&(p2[0]), &bar2, &q2);
                r2x2_t R2; r2x2_rot_and_scale(&q2, &R2);

                r2x2_inv(&R1, &R1); r2x2_mul(&R1, &R2, &L);
              }
            else if (np == 3)
              { /* Compute {L} by least squares: */
                ???
                r2x2_t E; r2x2_zero(&E); /* Moment matrix. */
                r2x2_t P; r2x2_zero(&P); /* Projection matrix. */
                for (int32_t k = 0; k < np; k++)
                  { double wk = (w == NULL ? 1.0 : w[k]);
                    /* Reduce points relative to barycenter: */
                    r2_t q1k, q2k;
                    r2_sub(&(p1[k]), &bar1, &q1k);
                    r2_sub(&(p2[k]), &bar2, &q2k);
                    /* Accumulate moments and projections: */
                    for (int32_t i = 0; i < 2; i ++)
                      { for (int32_t j = 0; j < 2; j++)
                          { E.c[i][j] += wk*q1k.c[i]*q1k.c[j];
                            P.c[i][j] += wk*q1k.c[i]*q2k.c[j];
                          }
                      }
                  }
                r2x2_t Z; r2x2_inv(&E, &Z);
                r2x2_mul(&Z, &P, &L);
            else
              { /* Compute {L} by least squares: */

                r2x2_t E; r2x2_zero(&E); /* Moment matrix. */
                r2x2_t P; r2x2_zero(&P); /* Projection matrix. */
                for (int32_t k = 0; k < np; k++)
                  { double wk = (w == NULL ? 1.0 : w[k]);
                    /* Reduce points relative to barycenter: */
                    r2_t q1k, q2k;
                    r2_sub(&(p1[k]), &bar1, &q1k);
                    r2_sub(&(p2[k]), &bar2, &q2k);
                    /* Accumulate moments and projections: */
                    for (int32_t i = 0; i < 2; i ++)
                      { for (int32_t j = 0; j < 2; j++)
                          { E.c[i][j] += wk*q1k.c[i]*q1k.c[j];
                            P.c[i][j] += wk*q1k.c[i]*q2k.c[j];
                          }
                      }
                  }
                r2x2_t Z; r2x2_inv(&E, &Z);
                r2x2_mul(&Z, &P, &L);
              }
              
            /* Store linear part into {A}: */
            A.c[1][1] = L.c[0][0];
            A.c[1][2] = L.c[0][1];
            A.c[2][1] = L.c[1][0];
            A.c[2][2] = L.c[1][1];
            
            /* Compute the displacement taking {L} into account: */
            r2_t v1; r2x2_map_row(&bar1, &L, &v1);
            r2_sub(&bar2, &v1, &d);
          }

        /* Store the displacement vector {d}: */
        A.c[0][1] = d.c[0];
        A.c[0][2] = d.c[1];
      }

    if (debug) 
      { fprintf(stderr, "  matrix:\n");
        r3x3_gen_print(stderr, &(A), "%13.6e", "", "\n", "\n", "    [ ", " ", " ]");
      }
      
    hr2_pmap_t M;
    M.dir = A;
    r3x3_inv(&A, &(M.inv));
    return M;
  }

