/* See {hr2_pmap_from_many_pairs_aux.h}. */
/* Last edited on 2024-09-16 18:33:55 by stolfi */

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

void hr2_pmap_from_many_pairs_try_find_frame_pair
  ( hr2_pmap_type_t type,
    int32_t np,
    r2_t p1[],
    r2_t p2[],
    double w[],
    int32_t *nt_P,                /* OUT */
    int32_t ixr[]                 /* OUT */
  )
  {
    int32_t nr_req = hr2_pmap_from_many_pairs_required_rank(type);
    int32_t nr = 0; /* For now. */
    
    typedef double pair_score_proc_t(r2_t *p1k, r2_t *p2k);
      /* Type of a procedure that evaluates a point pair for inclusion
        in the frame pair. If the score is zero or negative, the point
        pair is not good. Otherwise, the higher the score, the
        better. */

    auto int32_t find_best_pair(pair_score_proc_t *score);
      /* Returns the index {kp} in {0..np-1} that is not any of
        the indices {ixr[0..nr-1]}, and
        {w[kp]*score(p1[kp],p2[kp])} is maximum and positive.  If there is no such pair, 
        returns {-1}. */
    
    auto double bitriarea(r2_t *r1, r2_t *r2, r2_t *s1, r2_t *s2, r2_t *t1, r2_t *t2);
      /* Returns the minimum of the absolute values of the areas of the triangles
        {r1,s1,t1} and {r2,s2,t2}. */
    
    bool_t maybe_more = FALSE;
    while (maybe_more && (nr < np) && (nr < nr_req))
      { int32_t kp_best = -1; /* Best index for {ixr[nr]}. */
        switch(nr)
          { 
            case 0:
              { /* Try to find a pair with nonzero weight, as far as possible from the center: */
                /* Compute the barycenters of the point sets: */
                r2_t bar1, bar2;
                r2_barycenter(np, p1, w, &bar1);
                r2_barycenter(np, p2, w, &bar2);
                
                if (r2_is_finite(&bar1) && r2_is_finite(&bar1))
                  { /* There must be at least one point pair with positive weight. */
                    /* Try to find the pair with max weighted dist square from barycenters: */
                    auto double dist_from_center(r2_t *p1k, r2_t *p2k);
                    kp_best = find_best_pair(&dist_from_center);
                    
                    double dist_from_center(r2_t *p1k, r2_t *p2k)
                      { double Q1 = r2_dist_sqr(p1k, &bar1);
                        double Q2 = r2_dist_sqr(p2k, &bar2);
                        /* Add a fudge term so that the pair is accepted even if it is the barycenter: */
                        double Q = fmin(Q1, Q2) + 1.0e-200;
                        return Q;
                      }
                  }
                else
                  { kp_best = -1; }
                break;
              }
            case 1:
              { /* Try to find a pair with nonzero weight, as far as possible from the previous pair: */
                r2_t *a1 = &(p1[ixr[0]]);
                r2_t *a2 = &(p2[ixr[0]]);

                auto double dist_from_prev_pair(r2_t *p1k, r2_t *p2k);

                kp_best = find_best_pair(&dist_from_prev_pair);
                break;

                double dist_from_prev_pair(r2_t *p1k, r2_t *p2k)
                  { double Q1 = r2_dist_sqr(p1k, a1);
                    double Q2 = r2_dist_sqr(p2k, a2);
                    double Q = fmin(Q1, Q2);
                    return Q;
                  }
              }
            case 2:
              { /* Try to find a third pair not collinear with the prev two: */
                r2_t *a1 = &(p1[ixr[0]]);
                r2_t *a2 = &(p2[ixr[0]]);
                r2_t *b1 = &(p1[ixr[1]]);
                r2_t *b2 = &(p2[ixr[1]]);

                auto double dist_from_prev_line(r2_t *p1k, r2_t *p2k);

                kp_best = find_best_pair(&dist_from_prev_line);
                break;

                double dist_from_prev_line(r2_t *p1k, r2_t *p2k)
                  { double Q = bitriarea(a1,a2, b1,b2, p1k,p2k);
                    return Q;
                  }
              }
            case 3:
              { /* Try to find a ourth pair not collinear with any two of the prev three: */
                r2_t *a1 = &(p1[ixr[0]]);
                r2_t *a2 = &(p2[ixr[0]]);
                r2_t *b1 = &(p1[ixr[1]]);
                r2_t *b2 = &(p2[ixr[1]]);
                r2_t *c1 = &(p1[ixr[2]]);
                r2_t *c2 = &(p2[ixr[2]]);

                auto double dist_from_prev_three_lines(r2_t *p1k, r2_t *p2k);

                kp_best = find_best_pair(&dist_from_prev_three_lines);
                break;

                double dist_from_prev_three_lines(r2_t *p1k, r2_t *p2k)
                  { double Qab = bitriarea(a1, a2, b1, b2, p1k, p2k);
                    double Qbc = bitriarea(b1, b2, c1, c2, p1k, p2k);
                    double Qca = bitriarea(c1, c2, a1, a2, p1k, p2k);
                    double Q = fmin(Qab, fmin(Qbc, Qca));
                    return Q;
                  }
              }
            default:
              /* Should not happen: */
              assert(FALSE);
          }
        if (kp_best >= 0)
          { assert((kp_best >= 0) && (kp_best < np));
            ixr[nr] = kp_best;
            nr++;
          }
        else
          { maybe_more = FALSE; }
      }

    int32_t find_best_pair(pair_score_proc_t *score)
      { 
        double wQ_max = -INF;
        int32_t kp_best = -1;
        for (int32_t kp = 0; kp < np; kp++)
          { /* Skip pairs already selected: */
            bool_t chosen = FALSE;
            for (int32_t ir = 0; ir < nr; ir++) { if (kp == ixr[ir]) { chosen = TRUE; } }
            if (! chosen)
              { double wk = (w == NULL ? 1.0 : w[kp]);
                if (wk != 0) 
                  { r2_t *p1k = &(p1[kp]);
                    r2_t *p2k = &(p2[kp]);
                    double wQ = wk*score(p1k, p2k);
                    if ((wQ > 0) && (wQ > wQ_max)) { kp_best = kp; wQ_max = wQ; }
                  }
              }
          }
        return kp_best;
      }

    double bitriarea(r2_t *r1, r2_t *r2, r2_t *s1, r2_t *s2, r2_t *t1, r2_t *t2)
      { r2_t tr1; r2_sub(r1, t1, &tr1);
        r2_t ts1; r2_sub(s1, t1, &ts1);
        double Q1 = fabs(r2_det(&tr1, &ts1));

        r2_t tr2; r2_sub(r2, t2, &tr2);
        r2_t ts2; r2_sub(s2, t2, &ts2);
        double Q2 = fabs(r2_det(&tr2, &ts2));

        double Q = fmin(Q1, Q2);
        return Q;
      }
  }
  
hr2_pmap_t hr2_pmap_from_many_pairs_initial
  ( hr2_pmap_type_t type,
    int32_t np,
    r2_t p1[],
    r2_t p2[],
    double w[],
    int32_t nr,
    int32_t ixr[],
    sign_t sgn,
    int32_t class
  )
  {
    hr2_pmap_t M;
    if (type == hr2_pmap_type_IDENTITY)
      { assert(nr == 0);
        assert(class == 0);
        M = hr2_pmap_identity();
      }
    else if (type == hr2_pmap_type_TRANSLATION)
      { assert(nr == 1);
        assert (class == 0);
        /* Compute the barycenters of the point sets: */
        r2_t bar1, bar2;
        r2_barycenter(np, p1, w, &bar1);
        r2_barycenter(np, p2, w, &bar2);
        /* If the barycenter is NAN or infinite, there must be bad points: */
        demand((r2_is_finite(&bar1) && r2_is_finite(&bar1)), "invalid points");
        r2_t disp; r2_sub(&bar2, &bar1, &disp);
        M = hr2_pmap_translation(&disp);
      }
    else 
      { assert(nr >= 2);
        hr2_pmap_t M1, M2;
        r2_t *a1 = &(p1[ixr[0]]);
        r2_t *a2 = &(p2[ixr[0]]);
        r2_t *b1 = &(p1[ixr[1]]);
        r2_t *b2 = &(p2[ixr[1]]);
        if (type == hr2_pmap_type_CONGRUENCE)
          { assert(nr == 2);
            assert((class == 0) || (class == 1));

            r2_t mid1; r2_mix(0.5, a1, 0.5, b1, &mid1);
            r2_t dir1; r2_sub(b1, a1, &dir1); (void)r2_dir(&dir1, &dir1);
            M1 = hr2_pmap_congruence_from_point_and_dir(&mid1, &dir1, +1);
            
            r2_t mid2; r2_mix(0.5, a2, 0.5, b2, &mid2);
            r2_t dir2; r2_sub(b2, a2, &dir2); (void)r2_dir(&dir2, &dir2);
            M2 = hr2_pmap_congruence_from_point_and_dir(&mid2, &dir2, +1);
          }
        else if (type == hr2_pmap_type_SIMILARITY)
          { assert(nr == 2);
            assert((class == 0) || (class == 1));
            M1 = hr2_pmap_similarity_from_two_points(a1, b1, +1);
            M2 = hr2_pmap_similarity_from_two_points(a2, b2, +1);
          }
        else
          { r2_t *c1 = &(p1[ixr[2]]);
            r2_t *c2 = &(p2[ixr[2]]);
            if (sgn < 0) 
              { /* Swap {b2,c2} to reverse the orientation of the map: */
                r2_t *tmp = b2; b2 = c2; c2 = tmp;
              }
            if (type == hr2_pmap_type_AFFINE)
              { assert(nr == 3);
                M1 = hr2_pmap_aff_from_three_points(a1, b1, c1);
                M2 = hr2_pmap_aff_from_three_points(a2, b2, c2);
              }
            else if (type == hr2_pmap_type_GENERIC)
              { assert(nr == 4);
                r2_t *d1 = &(p1[ixr[3]]);
                r2_t *d2 = &(p2[ixr[3]]);
                M1 = hr2_pmap_from_four_r2_points(a1, b1, c1, d1);
                M2 = hr2_pmap_from_four_r2_points(a2, b2, c2, d2);
                assert((class >= 0) && (class < 4));
                if (class != 0)
                  { /* Apply the hither/yonder options from {class}: */
                    hr2_pmap_t K = hr2_pmap_r2_from_class(class);
                    M1 = hr2_pmap_compose(&K, &M1);
                  }
              }
            else
              { demand(FALSE, "invalid map type"); }
          }
          
        hr2_pmap_t N1 = hr2_pmap_inv(&M1);
        M = hr2_pmap_compose(&N1, &M2);
      }
    if (sgn*r3x3_det(&(M.dir)) < 0)
      { /* Pre-compose with a {XY} swap: */
        hr2_pmap_t S = hr2_pmap_xy_swap();
        M = hr2_pmap_compose(&S, &M);
      }
    return M;
  }

void hr2_pmap_from_many_pairs_optimize
  ( hr2_pmap_type_t type,
    sign_t sgn,
    int32_t np,
    r2_t p1[],
    r2_t p2[],
    double w[],
    int32_t maxIter,
    double maxErr,
    hr2_pmap_t *M,     /* (IN/OUT) The projective map to adjust. */
    double *f2M_P,     /* (IN/OUT) Goal function value for {*A}. */
    bool_t verbose
  )
  {
    auto double goalf(hr2_pmap_t *M);
      /* Computes the mean squared distance between the positions of the
         mapped points {p1*M.dir} and {p2}, and {p2*M.inv} and {p1}. */

    double f2M = goalf(M);

    if (verbose) { fprintf(stderr, "initial rms error = %13.6f\n", sqrt(f2M)); }

    double maxMod = 1.0e6; /* !!! FIX THIS !!! */
    
    hr2_pmap_opt_quadratic(type, sgn, &goalf, maxIter, maxErr, maxMod, M, &f2M, verbose);

    if (verbose) { fprintf(stderr, "final rms error = %13.6f\n", sqrt(f2M)); }

    return;
    
    double goalf(hr2_pmap_t *M)
     { 
       double f2 = hr2_pmap_mismatch_sqr(M, np, p1, p2, w);
       return f2;
     }
  }
