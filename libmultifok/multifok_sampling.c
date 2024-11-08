/* See {multifok_sampling.h}. */
/* Last edited on 2024-10-29 23:41:06 by stolfi */

#define _GNU_SOURCE
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <math.h>
#include <string.h>
#include <assert.h>

#include <affirm.h>
#include <r2.h>
#include <wt_table.h>
#include <wt_table_hann.h>

#include <multifok_sampling.h>

void multifok_sampling_generate_2D_hann_samples
  ( int32_t H,
    double R,
    int32_t *N_P,
    r2_t **v_P,
    double **w_P
  );
  /* Generates an array {v[0..N-1]} of 2D vectors and an array {w[0..N-1]} of
    weights, where {N = (2*H+1)^2}.  
    
    The vectors {v[0..N-1]} are a regular orthogonal grid of size {2*H+1}
    by {2*H+1}, symmetric about {(0,0)}, with spacing {1/(H+1)} along each 
    coordinate. Thus the max abs value of each coordinate is {H/(H+1)}.
    
    The weights are a 2D windowing function that is the product of a 1D
    Hann (raised cosine) window function on each coordinate. The vectors
    and weights form a partition or unit if replicated over the plane
    with stride 1.0 along each coordinate.
    
    The vectors and weights are sorted do that {|v[k]|} is increasing
    with {k}. Thus {v[0]} is always {(0,0)} and {w[0]} is 1. */

void multifok_sampling_generate_2D_hann_samples
  ( int32_t H,
    double R,
    int32_t *N_P,
    r2_t **v_P,
    double **w_P
  )
  {   
    demand(H >= 0, "invalid {H}");
    bool_t debug = FALSE;
    
    /* Allocate and fill the final arrays: */
    int32_t W = 2*H + 1;
    int32_t N = W*W;
    if (debug) { fprintf(stderr, "generating %d vectors total (H = %d)\n", N, H); }
    r2_t *v = talloc(N, r2_t);
    double *w = talloc(N, double);
    
    /* Define element 0: */
    int32_t ko = 0; /* Elements already filled: */
    v[ko] = (r2_t){{ 0,0 }};
    w[ko] = 1.0;
    ko++;

    if (H > 0)
      { /* Generate a 1D Hann window weight table {wHan[0..W-1]}: */
        double wHan[W];
        wt_table_hann_fill(W, 0.0, wHan, NULL);
        for (int32_t j = 0; j < W; j++)
          { fprintf(stderr, "  wHan[%d] = %10.8f\n", j, wHan[j]); }
        /* Paranoia check of partition-of-unity property: */
        assert(wHan[H] == 1.0);
        for (int32_t j = 0; j < H; j++)
          { assert(fabs(wHan[j] + wHan[j+H+1] - 1.0) < 1.0e-13); }

        /* Generate the 2D samples and weights for the first quadrant: */
        int32_t Q = (N - 1)/4;
        if (debug) { fprintf(stderr, "generating %d vectors in first quadrant\n", Q); }
        r2_t *vq = talloc(Q, r2_t);
        double *wq = talloc(Q, double);
        int32_t kq = 0; /* Sample index in {0..NS-1}. */
        for (int32_t ix = 1; ix <= +H; ix++)
          { double vx = R*((double)ix)/(H+1);
            for (int32_t iy = 0; iy <= +H; iy++)
              { assert(kq < Q);
                double vy = R*((double)iy)/(H+1);
                vq[kq] = (r2_t){{ vx, vy }};
                double wk = wHan[ix+H]*wHan[iy+H];
                wq[kq] = wk;
                kq++;
              }
          }
        assert(kq == Q);

        /* Sort 1st quadrant items by increasing distance: */
        int32_t it[Q]; /* A permutation of the indices {0..Q-1}. */
        for (int32_t j = 0; j < Q; j++) { it[j] = j; }
        auto int comp_weight(const void *a, const void *b);
        qsort(it, Q, sizeof(int32_t), &comp_weight);

        int comp_weight(const void *a, const void *b)
          { int32_t ia = *((int32_t *)a); double da = r2_norm(&(vq[ia]));
            int32_t ib = *((int32_t *)b); double db = r2_norm(&(vq[ib]));
            if (da < db)
              { return -1; }
            else if (da > db)
              { return +1; }
            else
              { return 0; }
          }
          
        /* Paranoia: */
        for (int32_t j = 1; j < Q; j++) 
          { if (debug) 
              { double rj = r2_norm(&(vq[it[j]]));
                fprintf(stderr, "    it[%4d] = %4d vq[%4d] = ", j, it[j], it[j]);
                r2_gen_print(stderr, &(vq[it[j]]), "%+14.8f", " ( ", " ", " )");
                fprintf(stderr, "  norm = %15.9f wq[%4d] = %10.8f\n", rj, it[j], wq[it[j]]);
              }
            assert(r2_norm(&(vq[it[j]])) >= r2_norm(&(vq[it[j-1]])));
          }
          
        /* Replicate the quadrant items into {v[1..N-1],w[1..N-1]: */

        for (int32_t j = 0; j < Q; j++)
          { double wj = wq[it[j]]; 
            double xj = vq[it[j]].c[0];
            double yj = vq[it[j]].c[1]; 
            assert (ko + 4 <= N);
            w[ko] = wj; v[ko] = (r2_t){{ +xj, +yj }}; ko++;
            w[ko] = wj; v[ko] = (r2_t){{ -yj, +xj }}; ko++;
            w[ko] = wj; v[ko] = (r2_t){{ -xj, -yj }}; ko++;
            w[ko] = wj; v[ko] = (r2_t){{ +yj, -xj }}; ko++;
          }

        free(wq);
        free(vq);
      }
    assert(ko == N);

    /* Return results: */
    (*N_P) = N;
    (*v_P) = v;
    (*w_P) = w;
  }

void multifok_sampling_choose_pixel_sampling_points_and_weights
  ( int32_t HS,
    int32_t *NS_P,
    r2_t **uSmp_P,
    double **wSmp_P,
    bool_t verbose
  )
  { /* This implementation chooses the subsampling points 
      in a regular orthogonal grid of {2*HS+1} by {2*HS+1},]
      and 2D Hann (raised cosine) weights. */

    int32_t NS;
    r2_t *uSmp;
    double *wSmp;
    multifok_sampling_generate_2D_hann_samples(HS, 1.0, &NS, &uSmp, &wSmp);

    if (verbose)
      { fprintf(stderr, "      generated %d sampling points and weights:\n", NS);
        for (int32_t ks = 0; ks < NS; ks++)
          { fprintf(stderr, "       %4d uSmp = ", ks);
            r2_gen_print(stderr, &(uSmp[ks]), "%+9.6f", "( ", " ", " )");
            fprintf(stderr, " wSmp = %12.10f\n", wSmp[ks]);
          }
      }

    (*NS_P) = NS;
    (*uSmp_P) = uSmp;
    (*wSmp_P) = wSmp;
  }
   
void multifok_sampling_choose_ray_tilts_and_weights
  ( int32_t KR, 
    int32_t NS,
    int32_t *NR_P,
    r2_t **tRay_P, 
    double **wRay_P,
    bool_t verbose
  )
  {
    /* This implementation chooses the ray tilts in a regular
      orthogonal grid of {2*HR+1} by {2*HR+1}, and 2D Hann 
      (raised cosine) weights; where {HR} is determined
      so that {NR} is 1 if {KR==0}, or at least {KR*NS}
      if {KR>0}. */
    
    int32_t HR;
    if (KR == 0)
      { if (verbose) { fprintf(stderr, "using one vertical ray for all sampling points\n"); }
        HR = 0; 
      }
    else
      { demand(KR >= 0, "invalid {KR}");
        demand(NS >= 1, "invalid {NS}");
        demand((NS % 2) == 1, "invalid {NS}");

        /* Factor {NS} into {SS^2*MS} where {MS} is squarefree: */
        int32_t SS = 1, MS = NS;
        int32_t d = 3;
        while (TRUE)
          { int32_t d2 = d*d;
            if (MS < d2) 
              { break; }
            else if ((MS % d2) == 0)
              { MS /= d2; SS *= d; }
            else
              { d += 2; }
          }
        if (verbose) { fprintf(stderr, "factored NS = %d into *%d^2*%d\n", NS, SS, MS); }
        assert(SS*SS*MS == NS);
        assert((SS*MS % 2) == 1);
        /* Now choose {KR_new} so that {KR_new*MS} is an odd square and {KR_new >= KR}: */
        int32_t TR = (int32_t)ceil((sqrt(((double)KR)/((double)MS)) - 1.0)/2);
        int32_t LR_loc = (2*TR+1)*MS*SS;
        int32_t NR_loc = LR_loc*LR_loc;
        int32_t KR_new = NR_loc/NS;
        HR = (LR_loc - 1)/2;
        if (verbose) 
          { if (KR_new != KR)
              { fprintf(stderr, "rays per sampling point {KR} rounded up from %d to %d\n", KR, KR_new); }
          }
        KR = KR_new;
        assert(KR >= 1);
      }
    int32_t LR = 2*HR + 1;
    int32_t NR = LR*LR;
    fprintf(stderr, "generating %d × %d = %d total rays (HR = %d)\n", LR, LR, NR, HR);
    if (NR != 1) { assert((NR % NS) == 0); }
        
    int32_t NR_gen;
    r2_t *tRay;
    double *wRay;
    multifok_sampling_generate_2D_hann_samples(HR, 1.0, &NR_gen, &tRay, &wRay);
    assert(NR_gen == NR);
    
    if (NR > 1)
      { /* Normalize the tilts to unit RMS radius: */
        double sum_w = 0;
        double sum_w_r2 = 0;
        for (int32_t ir = 0; ir < NR; ir++)
          { sum_w += wRay[ir];
            double r2 = r2_norm_sqr(&(tRay[ir])); 
            sum_w_r2 += wRay[ir]*r2;
          }
        assert(sum_w >= 1.0);
        double r_avg = sqrt(sum_w_r2/sum_w);
        assert(r_avg > 0.0);
        for (int32_t ir = 0; ir < NR; ir++) 
          { tRay[ir].c[0] /= r_avg;
            tRay[ir].c[1] /= r_avg;
          }
      }

    if (verbose)
      { fprintf(stderr, "      generated %d ray tilts and weights:\n", NR);
        for (int32_t ir = 0; ir < NR; ir++)
          { fprintf(stderr, "      ray %4d tRay = ", ir);
            r2_gen_print(stderr, &(tRay[ir]), "%+9.6f", "( ", " ", " )");
            fprintf(stderr, " wRay = %12.10f\n", wRay[ir]);
          }
      }

    /* Return results: */
    (*NR_P) = NR;
    (*tRay_P) = tRay;
    (*wRay_P) = wRay;
    
    return;
  }
