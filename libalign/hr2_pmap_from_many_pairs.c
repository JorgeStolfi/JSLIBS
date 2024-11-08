/* See {hr2_pmap_from_many_pairs.h}. */
/* Last edited on 2024-11-01 06:10:08 by stolfi */

#define _GNU_SOURCE
#include <math.h>
#include <stdint.h>
#include <limits.h>
#include <assert.h>
#include <string.h>
#include <limits.h>
#include <math.h>
 
#include <bool.h>
#include <sign.h>
#include <r2.h>
#include <r3x3.h>
#include <affirm.h>
#include <hr2_pmap_opt.h>

#include <hr2_pmap_from_many_pairs.h>
#include <hr2_pmap_from_many_pairs_aux.h>

hr2_pmap_t hr2_pmap_from_many_pairs
  ( hr2_pmap_type_t type_req,
    int32_t np,
    r2_t p1[], 
    r2_t p2[], 
    double w[],
    int32_t maxIter,
    double maxErr,
    bool_t verbose
  )
  {
    bool_t debug = TRUE;
    if (debug) fprintf(stderr, "  > begin {hr2_pmap_from_many_pairs}\n");

    int32_t nr = -1;  /* Frame pair size {0..3}. */
    int32_t ixr[4];   /* The frame pair has indices {ixr[0..nr-1]}. */
    hr2_pmap_from_many_pairs_try_find_frame_pair(type_req, np, p1, p2, w, &nr, ixr);

    /* Downgrade the type if there are not enough points: */
    hr2_pmap_type_t type_eff = type_req;
    if (nr == 0)
      { type_eff = hr2_pmap_type_IDENTITY; }
    else if (nr == 1)
      { if (type_eff > hr2_pmap_type_TRANSLATION) { type_eff = hr2_pmap_type_TRANSLATION; } }
    else if (nr == 2)
      { if (type_eff > hr2_pmap_type_SIMILARITY) { type_eff = hr2_pmap_type_SIMILARITY; } }
    else if (nr == 3)
      { if (type_eff > hr2_pmap_type_AFFINE) { type_eff = hr2_pmap_type_AFFINE; } }
   
    /* For similarity and congruence, we must try twice, varying the {flip}: */
    /* For general projective, we must try four times, varying the sign of two weights: */
    hr2_pmap_t M_best;     /* Best map found. */
    double Fx_best = +INF; /* Its discrepancy: */
    int32_t nsgn = ((type_eff == hr2_pmap_type_IDENTITY) || (type_eff == hr2_pmap_type_TRANSLATION) ? 1 : 2);
    int32_t nclass = ((type_eff == hr2_pmap_type_GENERIC) ? 4 : 1);
    for (int32_t isgn = 0; isgn < nsgn; isgn++)
      { sign_t sgn = (isgn == 1 ? -1 : +1);
        for (int32_t class = 0; class < nclass; class++)
          { hr2_pmap_t M = hr2_pmap_from_many_pairs_initial(type_eff, np, p1, p2, w, nr, ixr, sgn, class);
            double Fx = hr2_pmap_mismatch_sqr(&M, np, p1, p2, w);
            if ((type_eff != hr2_pmap_type_IDENTITY) && (type_eff != hr2_pmap_type_TRANSLATION))
              { /* Needs iterative optimization. */
                if (verbose) { fprintf(stderr, "optimizing...\n"); }
                hr2_pmap_from_many_pairs_optimize(type_eff, sgn, np, p1, p2, w, maxIter, maxErr, &M, &Fx, verbose);
              }
            if (Fx < Fx_best) { Fx_best = Fx; M_best = M; }
          }
      }

    if (verbose) { fprintf(stderr, "final rms error = %13.6f\n", Fx_best); }
    return M_best;
  }


int32_t hr2_pmap_from_many_pairs_required_rank(hr2_pmap_type_t type)
  {
    switch(type)
      { 
        case hr2_pmap_type_IDENTITY:
          return 0;
        case hr2_pmap_type_TRANSLATION:
          return 1;
        case hr2_pmap_type_CONGRUENCE:
        case hr2_pmap_type_SIMILARITY:
          return 2;
        case hr2_pmap_type_AFFINE:
          return 3;
        case hr2_pmap_type_GENERIC:
          return 4;
        case hr2_pmap_type_NONE:
          demand(FALSE, "invalid projective map type");
        default:
          demand(FALSE, "unimplemented map type");
      }
  }
