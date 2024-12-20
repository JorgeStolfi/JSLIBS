/* See {hr2test_tools.h}. */
/* Last edited on 2024-12-05 10:19:47 by stolfi */

#include <stdlib.h>
#include <stdint.h>
#include <stdio.h>
#include <assert.h>
#include <math.h>

#include <jsrandom.h>
#include <affirm.h>

#include <rn.h>
#include <r3x3.h>
#include <r3.h>
#include <r2x2.h>
#include <r2.h>
#include <hr2.h>
#include <hr2_pmap.h>
#include <hr2_pmap_throw_by_type.h>
#include <hr2_test_tools.h>
#include <hr2_pmap_test_tools.h>
#include <hr2_pmap_from_many_pairs.h>

#include <hr2_pmap_opt_test_tools.h>

#define NH 3
  /* Number of homogeneous coordinates in a point. */

#define NC 2
  /* Number of Cartesian coordinates in a point. */

#define ht_IDENTITY    hr2_pmap_type_IDENTITY
#define ht_TRANSLATION hr2_pmap_type_TRANSLATION
#define ht_CONGRUENCE  hr2_pmap_type_CONGRUENCE
#define ht_SIMILARITY  hr2_pmap_type_SIMILARITY
#define ht_AFFINE      hr2_pmap_type_AFFINE
#define ht_GENERIC     hr2_pmap_type_GENERIC
  /* Shorter names. */

void hr2_pmap_opt_test_choose_r2_point_pairs
  ( hr2_pmap_type_t type,
    sign_t sgn,
    bool_t tight,
    bool_t ident, 
    double pert,
    bool_t verbose,
    int32_t *np_P,
    r2_t **p1_P,
    r2_t **p2_P,
    double **w_P,
    hr2_pmap_t *M_P
  )
  { 
    char *xtype = hr2_pmap_type_to_string(type);
    if (verbose) { fprintf(stderr, "    > --- %s type = %s sgn = %+d tight = %c ident = %c --- \n", __FUNCTION__, xtype, sgn, "FT"[tight], "FT"[ident]); }
    demand((sgn == -1) || (sgn == +1), "invalid {sgn}");
    int32_t nr = hr2_pmap_from_many_pairs_required_rank(type);
    int32_t np = (tight? nr : nr + 10);
    r2_t *p1 = talloc(np, r2_t);
    r2_t *p2 = talloc(np, r2_t);
    double *w = (drandom() < 0.2 ? NULL : talloc(np, double));
    
    /* Get maps {M1,M2} with the right type and sign: */
    hr2_pmap_t M1 = hr2_pmap_throw_by_type((ident ? ht_IDENTITY : type), +1);
    hr2_pmap_t M2 = hr2_pmap_throw_by_type((ident ? ht_IDENTITY : type), sgn);
    
    /* Generate points {p1[0..np-1]} and weights {w[0..np-1]}: */
    for (uint32_t kp = 0;  kp < np; kp++)
      { /* Fill {p1} with {nr} well-spaced points and some random points: */
        r2_t q;
        if (kp < nr)
          { double x = ((kp&1) != 0 ? -1.0 : +1.0);
            double y = ((kp&2) != 0 ? -1.0 : +1.0);
            q = (r2_t){{ x, y }};
          }
        else
          { q = (r2_t){{ 10.0*(2*drandom()-1),  10.0*(2*drandom()-1) }}; }
        /* Set {p1,p2} to the images of {q} by {M1,M2} plus noise: */
        p1[kp] = hr2_pmap_r2_point(&(q), &M1);
        p2[kp] = hr2_pmap_r2_point(&(q), &M2);
        if ((! ident) && ((! tight) || (kp > nr)))
          { /* Add random noise, small enough to preserve the handedness: */
            r2_t e1; r2_throw_ball(&e1); r2_mix_in(pert, &e1, &(p1[kp]));
            r2_t e2; r2_throw_ball(&e2); r2_mix_in(pert, &e2, &(p2[kp]));
            /* Set the weight to random, less than 1: */
            if (w != NULL) { w[kp] = fmax(0.0, 1.01*(drandom()*drandom()) - 0.01); }
          }
        else
          { /* Set the weight to random, large, less than 1: */
            if (w != NULL) { w[kp] = 0.7 + 0.3*drandom(); }
          }
        if (verbose)
          { fprintf(stderr, "      %3d", kp);
            r2_gen_print(stderr, &(p1[kp]), "%+10.5f", " p1 = ( ", " ", " )");
            r2_gen_print(stderr, &(p2[kp]), "%+10.5f", " p2 = ( ", " ", " )");
            if (w != NULL) { fprintf(stderr, " w = %10.8f", w[kp]); }
            fprintf(stderr, "\n");
          }
      }
    if ((type != ht_AFFINE) || (type == ht_GENERIC))
      { /* Print table of distance ratios: */
        fprintf(stderr, "      rel point distances\n");
        for (uint32_t kp = 0;  kp < np; kp++)
          { fprintf(stderr, "      %3d", kp);
            for (uint32_t jp = 0;  jp < kp; jp++)
              { double D1kj = r2_dist(&(p1[kp]), &(p1[jp]));
                double D2kj = r2_dist(&(p2[kp]), &(p2[jp]));
                double rat = D1kj/(D2kj + 1.0e-200);
                double tar = D2kj/(D1kj + 1.0e-200);
                fprintf(stderr, " %8.6f", 0.5*(rat - tar));
              }
            fprintf(stderr, "\n");
          }
      }

    (*np_P) = np;
    (*p1_P) = p1;
    (*p2_P) = p2;
    (*w_P) = w;
    (*M_P) = hr2_pmap_inv_compose(&M1, &M2);
    if (verbose) { fprintf(stderr, "    < --- %s ---\n", __FUNCTION__); }
  }

void hr2_pmap_opt_test_print_map(char *name, hr2_pmap_t *M, double f2M)
  {
    fprintf(stderr, "  %s", name);
    if (! isnan(f2M)) { fprintf(stderr, " ( f2 = %20.12f = %12.4e )", f2M, f2M); }
    fprintf(stderr, " =\n");
    hr2_pmap_gen_print(stderr, M, "%12.7f", "    ", "[ ","  "," ]\n    ", "[ "," "," ]", "\n");
  }
