/* See {hr2_pmap_from_many_pairs_aux.h}. */
/* Last edited on 2024-11-21 21:17:35 by stolfi */

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
#include <hr2_pmap_special_opt.h>

#include <hr2_pmap_from_many_pairs.h>
#include <hr2_pmap_from_many_pairs_aux.h>

void hr2_pmap_from_many_pairs_choose_yrad
  ( hr2_pmap_type_t type,
    sign_t sgn,
    int32_t np,
    r2_t p1[],
    r2_t p2[],
    int32_t n,
    double yctr[],
    double yrad[]
  );
  /* Chooses a suitable parameter {yrad[0..n-1]} so that 
    the optimum is in the box {yctr ± yrad},
    and changing each encoding element {yctr[k]} by {eps*yrad[k]} 
    will have about the same effect on the mismatch squared. */

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
        for (uint32_t kp = 0;  kp < np; kp++)
          { /* Skip pairs already selected: */
            bool_t chosen = FALSE;
            for (uint32_t ir = 0;  ir < nr; ir++) { if (kp == ixr[ir]) { chosen = TRUE; } }
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
        M = hr2_pmap_translation_from_disp(&disp);
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
                M1 = hr2_pmap_affine_from_three_points(a1, b1, c1);
                M2 = hr2_pmap_affine_from_three_points(a2, b2, c2);
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
                    hr2_pmap_t K = hr2_pmap_r2_from_sign_class(class);
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

void hr2_pmap_from_many_pairs_choose_yrad
  ( hr2_pmap_type_t type,
    sign_t sgn,
    int32_t np,
    r2_t p1[],
    r2_t p2[],
    int32_t n,
    double yrad[]
  );

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
        
    auto bool_t ok_pred(hr2_pmap_t *A, double fA);
      /* Checks whether {A} is acceptable. Currently, that
        means {fA} is less than {maxErr}. */

    double f2M = goalf(M);

    if (verbose) { fprintf(stderr, "initial rms error = %13.6f\n", sqrt(f2M)); }

    double yrad[n];
    hr2_pmap_from_many_pairs_choose_yrad(type, sgn, np, p1, p2, n, yrad);
    
    hr2_pmap_special_opt_quadratic
      ( type, sgn, &goalf, &ok_pred, 
        yrad, maxIter, M, &f2M, verbose
      );

    if (verbose) { fprintf(stderr, "final rms error = %13.6f\n", sqrt(f2M)); }

    return;
    
    double goalf(hr2_pmap_t *M)
     { double f2 = hr2_pmap_mismatch_sqr(M, np, p1, p2, w);
       return f2;
     }
     
    bool_t ok_pred(hr2_pmap_t *A, double fA)
      { return fA < maxErr; }
      
  }

void hr2_pmap_from_many_pairs_choose_yrad
  ( hr2_pmap_type_t type,
    sign_t sgn,
    int32_t np,
    r2_t p1[],
    r2_t p2[],
    double w[],
    int32_t n,
    double yctr[],
    double yrad[]
  )
  {
    /* !!! Improve !!! */
    
    /* Compute the barycenters of the point clouds: */
    r2_t pBar1; r2_barycenter(np, p1, w, &bar1);
    r2_t pBar2; r2_barycenter(np, p2, w, &bar2);
    double dBar = r2_dist(&pBar1, &pBar2);
    
    /* Compute the RMS radii of the clouds: */
    double pRad1 = sqrt(r2_mean_dist_sqr(np, p1, w, &pBar1));
    double pRad2 = sqrt(r2_mean_dist_sqr(np, p2, w, &pBar2));
    double pRad = fmax(pRad1, pRad2);
    
    /* Compute the overall magnification either way: */
    double pMag = fmax(pRad1/pRad2, pRad2/pRad2);
    
    /* Let {D} be the typical distance between {M(p1[kp])} and {p2[kp]}
      or {p1[kp]} and {M^{-1}(p2[kp])}.  The following numbers
      are the approximate derivative of {D} with respect to 
      each map encoding coordinate {y[0..n-1]} for each type of
      encoding coordinate: */
      
    double fTrn = 1.0;       /* Translation. */
    double fRot = pRad;      /* Rotation. */
    double fMag = pMag;      /* Primary {u} vector length. */
    double fShr = pMag*PMag; /* Relative {v} vector coords. */
    
    /* Let {z[0..n-1]} be the optimization coordinates. The encoding
      coordinates will be {yctr[k] + z[k]*yrad[k]} where {yctr[k]} is
      the encoding of the initial guess, hopefully close to the optimum
      map. To make sure that the optimum map {M} is reached when the
      optimization variables {z[0..n-1]} vary over {[-1 _ +1]}, the
      radii {yrad[0..n-1]} must be at least the following, depending on
      the type of encoding coordinate: */
      
    double mTrn = 2*dBar;      /* Translation. */
    double mRot = 3*M_PI;      /* Rotation. */
    double mMag = 2*pMag;      /* Primary {u} vector length. */
    double mShr = 50.0;        /* Relative {v} vector coords. */
    
    /* !!! The {mShr} value is a hack. !!! */
    /* !!! Should use the eigenvectors of point clouds !!! */
    
    /* The values of {yrad[0..n-1]} must be proportional to 
      {fTrn}, {fRot}, etc, by some factor {srad},
      but at least {mTrn}, {mRot}, etc. */
    
    double srad;
    switch (type)
      { case ht_IDENTITY:
          /* Nothing to do: */
          break;
        case ht_TRANSLATION:
          srad = mTrn/fTrn;
          yrad[0] = srad*fTrn; yrad[1] = srad*fTrn;
          break;
        case ht_CONGRUENCE:
          srad = fmax(mTrn/fTrn, mRot/fRot);
          yrad[0] = srad*fRot;
          yrad[1] = srad*fTrn; yrad[2] = srad*fTrn;
          break;
        case ht_SIMILARITY:
          srad = fmax(fmax(mTrn/fTrn, mRot/fRot), mMag/fMag);
          yrad[0] = sRad*fMag;
          yrad[1] = sRad*fRot;
          yrad[2] = srad*fTrn; yrad[3] = srad*fTrn;
          break;
        case ht_AFFINE:
          srad = fmax(fmax(mTrn/fTrn, mRot/fRot), fmax(mMag/fMag, mShr/fShr));
          yrad[0] = sRad*fMag;
          yrad[1] = sRad*fRot;
          yrad[2] = sRad*fShr; yrad[3] = sRad*fShr;
          yrad[4] = sRad*fTrn; yrad[5] = sRad*fTrn;
          break;
        case ht_GENERIC:
          break;
        case ht_NONE:
        default:
          demand(FALSE, "invalid type");
      }
          
  }
