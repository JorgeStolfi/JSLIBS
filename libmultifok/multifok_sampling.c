/* See {multifok_sampling.h}. */
/* Last edited on 2024-12-15 22:34:44 by stolfi */

#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <math.h>
#include <string.h>
#include <assert.h>

#include <affirm.h>
#include <i2.h>
#include <r2.h>
#include <wt_table.h>
#include <wt_table_hann.h>

#include <multifok_sampling.h>

void multifok_sampling_choose_pixel_sampoints_and_weights
  ( uint32_t HS, 
    uint32_t *NS_P,
    i2_t **iSmp_P, 
    double **wSmp_P,
    bool_t verbose
  );
  /* Generates an array {iSmp[0..NS-1]} of 2D integer vectors and an array {wSmp[0..NS-1]} of
    weights, where {NS = (2*HS+1)^2}.  
    
    The vectors {iSmp[0..NS-1]} are a subset of the regular orthogonal
    grid of size {2*HS+1} by {2*HS+1}, symmetric about {(0,0)}. Thus
    each coordinate will be in {-HS..+HS}.
    
    The weights are a 2D windowing function that is the product of a 1D
    Hann (raised cosine) window function on each coordinate. The vectors
    and weights form a partition or unit if replicated over the plane
    with stride 1.0 along each coordinate.
    
    The vectors and weights are sorted so that {|iSmp[k]|} is increasing
    with {k}. Thus {iSmp[0]} is always {(0,0)} and {wSmp[0]} is 1.
    
    Returns {NS,uSMP,wSmp} in {*NS_P,*iSmp_P,*wSmp_P}. */

void multifok_sampling_get_grid_points_in_quadrant
  ( uint32_t HP,
    i2_t **pq_P,
    uint32_t *NPQ_P,
    bool_t verbose
  );
  /* Allocates and fills a list {pq[0..NPQ-1]} of the {NPQ=HP*(HP+1)}
    integer grid points in the first quadrant with coordinates up to
    {HP}; specifically, {(ix,iy)} with {ix} in {1..HP} and {iy} in
    {0..HP}, sorted by distance from the origin. */

void multifok_sampling_sort_points_by_norm(uint32_t NP, i2_t p[]);
  /* Sorts the list of points {p[0..NP-1]} by increasing distance from origin.
    Specifically, according to {multifok_sampling_compare_points}. */ 
      
void multifok_sampling_check_point_order(uint32_t NP, i2_t p[]);    
  /* Checks whether the points {p[0..NP-1]} are sorted in incrasing order
    of distance from the origin. */
  
int32_t multifok_sampling_compare_points(const void *a, const void *b);
  /* Suitable for sorting routines like {qsort}. Assumes
    {a} and {b} point to integer pair ({i2_t}) records {pa} and {pb},
    and returns {-1} or {+1} if {pa} should come before
    or after {pb}, respectively; or 0 if they are the same pair.
    Compares {pa} with {pb} by distance from origin.
    Breaks ties by {min(|x|,|y|)}, then by {|x|+|y|}, 
    then by lex order. */

double *multifok_sampling_hann_weight_table(uint32_t HW, bool_t verbose);
  /* Allocates and fills a 1-dim table {w[0..NW-1]} of Hann (raised cosine)
    weights, where {NW = 2*HW + 1}. 
    
    Element {w[HW]} will be 1.0, and {w[HW-i]} will be equal to {w[HW+i]}.
    The table will have the partition-of-unit property with stride {HW+1};
    that is, {w[HW+1+i]+w[i]=1.0} for {i} in {0..HW-1}. */

void multifok_sampling_print_points_and_weights
  ( uint32_t NP,
    i2_t p[],
    char *pName,
    double step,
    double w[],
    char *wName
  );
  /* Prints {p0..NP-1]} and {w[0..NP-1]} to {stderr}. */

/* IMPLEMENTATIONS */

multifok_sampling_t *multifok_sampling_choose(uint32_t HS, uint32_t KR, bool_t verbose)
  { 
    multifok_sampling_t *samp = talloc(1, multifok_sampling_t);
    multifok_sampling_choose_pixel_sampoints_and_weights
      (HS, &(samp->NS), &(samp->iSmp), &(samp->wSmp), verbose);
    samp->step = 1.0/(double)(HS+1); /* So sampoints almost reach the center of next pixel. */
    if (verbose) 
      { fprintf(stderr, "      generated %d sampling points and weights:\n", samp->NS);
        multifok_sampling_print_points_and_weights(samp->NS, samp->iSmp, "iSmp", samp->step, samp->wSmp, "wSmp");
      }
    samp->KR = KR;
    return samp;
  }
  
void multifok_sampling_free(multifok_sampling_t *samp)
  { free(samp->iSmp);
    free(samp->wSmp);
    free(samp);
  }

void multifok_sampling_choose_pixel_sampoints_and_weights
  ( uint32_t HS,
    uint32_t *NS_P,
    i2_t **iSmp_P,
    double **wSmp_P,
    bool_t verbose
  )
  { /* This implementation chooses the subsampling points as a 
     regular orthogonal grid of {2*HS+1} by {2*HS+1}, 
     and 2D Hann (raised cosine) weights.  This ensures the
     partition of unity property for the sampling points. */

    demand(HS >= 0, "invalid {HS}");
    bool_t debug = TRUE; if (debug) { verbose = TRUE; }
    
    /* Allocate and fill the temporary arrays: */
    uint32_t NW = 2*HS + 1;
    uint32_t NS = NW*NW;
    if (debug) { fprintf(stderr, "generating %d × %d = %d sampling vectors total (HS = %d)\n", NW, NW, NS, HS); }

    i2_t *iSmp = talloc(NS, i2_t);
    double *wSmp = talloc(NS, double);

    /* Define element 0: */
    uint32_t ko = 0; /* Elements already filled: */
    iSmp[ko] = (i2_t){{ 0,0 }};
    wSmp[ko] = 1.0;
    ko++;
    
    if (HS > 0)
      { double *wHan = multifok_sampling_hann_weight_table(HS, verbose);
        uint32_t NSQ;
        i2_t *pq;
        multifok_sampling_get_grid_points_in_quadrant(HS, &pq, &NSQ, verbose);
        assert(NSQ > 0);
        assert(pq != NULL);
        assert(NS == 1 + 4*NSQ);

        /* Replicate the quadrant points into {iSmp[1..NS_max-1],wSmp[1..NS_max-1]: */
        for (uint32_t j = 0;  j < NSQ; j++)
          { int32_t xj = pq[j].c[0]; assert(abs(xj) <= HS);
            int32_t yj = pq[j].c[1]; assert(abs(yj) <= HS);
            double wj = wHan[xj+(int32_t)HS]*wHan[yj+(int32_t)HS]; 
            assert (ko + 4 <= NS);
            wSmp[ko] = wj; iSmp[ko] = (i2_t){{ +xj, +yj }}; ko++;
            wSmp[ko] = wj; iSmp[ko] = (i2_t){{ -yj, +xj }}; ko++;
            wSmp[ko] = wj; iSmp[ko] = (i2_t){{ -xj, -yj }}; ko++;
            wSmp[ko] = wj; iSmp[ko] = (i2_t){{ +yj, -xj }}; ko++;
          }
        free(pq);
        free(wHan);
      }
    assert(ko == NS);

    (*NS_P) = NS;
    (*iSmp_P) = iSmp;
    (*wSmp_P) = wSmp;
  }
     
void multifok_sampling_get_grid_points_in_quadrant
  ( uint32_t HP,
    i2_t **pq_P,
    uint32_t *NPQ_P,
    bool_t verbose
  )
  { 
    demand(HP >= 0, "invalid {HP}");
    bool_t debug = TRUE;
    
    /* Allocate and fill the temporary arrays: */
    uint32_t NW = 2*HP + 1;

    /* Generate the 2D samples and weights for the first quadrant: */
    uint32_t NPQ = (NW*NW - 1)/4;
    if (debug) { fprintf(stderr, "  generating %d grid points in first quadrant\n", NPQ); }
    i2_t *pq = talloc(NPQ, i2_t);

    uint32_t kq = 0; /* First-quadrant points generated so far. */
    for (int32_t ix = 1;  ix <= +HP; ix++)
      { for (int32_t iy = 0;  iy <= +HP; iy++)
          { assert(kq < NPQ);
            pq[kq] = (i2_t){{ ix, iy }};
            kq++;
          }
      }
    assert(kq == NPQ);
    
    multifok_sampling_sort_points_by_norm(NPQ, pq);
    
    (*pq_P) = pq;
    (*NPQ_P) = NPQ;
  } 

void multifok_sampling_sort_points_by_norm(uint32_t NP, i2_t p[])
  { qsort(p, NP, sizeof(i2_t), &multifok_sampling_compare_points); }
    
void multifok_sampling_check_point_order(uint32_t NP, i2_t p[])
  { bool_t debug = FALSE;
    for (uint32_t j = 0;  j < NP; j++) 
      { if (debug)
          { double rj = sqrt((double)i2_norm_sqr(&(p[j])));
            fprintf(stderr, "        p[%4d] = ", j);
            i2_gen_print(stderr, &(p[j]), "%+8d", " ( ", " ", " )");
            fprintf(stderr, "  norm = %15.9f\n", rj);
          }
        if (j > 0) { assert(multifok_sampling_compare_points(&(p[j-1]), &(p[j])) == -1); }
      }
  } 

int32_t multifok_sampling_compare_points(const void *a, const void *b)
  { i2_t *pa = ((i2_t *)a); uint64_t r2a = i2_norm_sqr(pa);
    i2_t *pb = ((i2_t *)b); uint64_t r2b = i2_norm_sqr(pb);
    if (r2a < r2b)
      { return -1; }
    else if (r2a > r2b)
      { return +1; }
    else
      { int32_t xa = pa->c[0], ya =  pa->c[1];
        int32_t xb = pb->c[0], yb =  pb->c[1];
        /* Try to break the tie by {min(|x|,|y|)} */
        int32_t ma = (int32_t)imin(abs(xa),abs(ya));
        int32_t mb = (int32_t)imin(abs(xb),abs(yb));
        if (ma < mb) 
          { return -1; }
        else if (ma > mb) 
          { return +1; }
        else
          { /* Try to break the tie by {|x|+|y|} */
            int32_t sa = abs(xa) + abs(ya);
            int32_t sb = abs(xb) + abs(yb);
            if (sa < sb) 
              { return -1; }
            else if (sa > sb) 
              { return +1; }
            else
              { /* Try to break the tie by lex order: */
                if (xa < xb)
                  { return -1; }
                else if (xa > xb)
                  { return +1; }
                else if (ya < yb)
                  { return -1; }
                else if (ya > yb)
                  { return +1; }
                else
                  { /* The pairs are equal: */
                    return 0;
                  }
              }
          }
      }
  }

void multifok_sampling_print_points_and_weights
  ( uint32_t NP,
    i2_t p[],
    char *pName,
    double step,
    double w[],
    char *wName
  )
  { for (uint32_t k = 0;  k < NP; k++)
      { fprintf(stderr, "       %4d %s = ", k, pName);
        i2_gen_print(stderr, &(p[k]), "%+3d", "( ", " ", " )");
        r2_t rp = (r2_t){{ p[k].c[0]*step,  p[k].c[1]*step }};
        r2_gen_print(stderr, &rp, "%+8.5f", " = ( ", " ", " )");
        if (w != NULL) { fprintf(stderr, " %s = %12.10f", wName, w[k]); }
        fprintf(stderr, "\n");
      }
  }
  
double* multifok_sampling_hann_weight_table(uint32_t HW, bool_t verbose)
  {
    uint32_t NW = 2*HW + 1;
    double *wHan = talloc(NW, double);
    wt_table_hann_fill(NW, 0.0, wHan, NULL);
    if (verbose)
      { for (uint32_t j = 0;  j < NW; j++)
          { fprintf(stderr, "  wHan[%d] = %10.8f\n", j, wHan[j]); }
      }
    assert(wHan[HW] == 1.0);
    /* Paranoia check of partition-of-unity property: */
    for (uint32_t j = 0;  j < HW; j++)
      { assert(fabs(wHan[j] + wHan[j+HW+1] - 1.0) < 1.0e-13); }
    return wHan;
  }
