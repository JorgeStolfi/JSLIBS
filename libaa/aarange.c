/* See aarange.h */
/* Last edited on 2024-12-31 01:14:55 by stolfi */

#include <stdlib.h>
#include <stdint.h>

#include <flt.h>
#include <aa.h>

#include <aarange.h>

/* INTERNAL PROTOS */

void aa_2d_range_get_pairs
  ( AAP x, 
    AAP y, 
    AATermCount *nv,
    double *xv,
    double *yv
  );
  /* Gathers the non-zero pairs of matching coefficients from {x} and
    {y}, and normalizes them so that either {y} is positive, or {y} is
    zero and {x} is positive. The variable {nv} is set to the number
    of gathered pairs, which are stored in {xv[0..nv-1],yv[0..nv-1]}.
    Assumes that there is enough space in {xv,yv}. */
    
void aa_2d_range_sort_pairs
  ( AATermCount nv,
    double *xv,
    double *yv
  );
  /* Assumes that {xv[0..nv-1],yv[0..nv-1]} are non-zero vectors
    with arguments between 0 (inclusive) and {\pi} (exclusive).
    Sorts them in order of increasing argument. */

void aa_2d_range_compute_vertices
  ( Float x0,
    Float y0,
    AATermCount nv,
    double *xv,
    double *yv
  );
  /* Assumes that {xv[0..nv-1],yv[0..nv-1]} are non-zero coeff pairs
    with arguments between 0 (inclusive) and {\pi} (exclusive),
    ordered by increasing argument. Computes the vertices
    {xv[0..2*nv-1],yv[0..2*nv-1]} of the joint range. */


/* IMPLEMENTATIONS */

void aa_2d_range_get_pairs
  ( AAP x, 
    AAP y, 
    AATermCount *nv,
    double *xv,
    double *yv
  )
  {
    AATermP xp = (AATermP) (x + 1); AATermCount xn = (x->nterms);
    AATermP yp = (AATermP) (y + 1); AATermCount yn = (y->nterms);
    AATermCount n = 0;
    Float xt, yt;
    
    /* Merge term sequences, gathering nonzero pairs: */
    /* fprintf(stderr, "gathering pairs...\n"); */
    while ((xn > 0) || (yn > 0))
      { 
        if ((yn == 0) || ((xn != 0) && (xp->id < yp->id)))
          { xt = xp->coef; yt = 0;
            xp++; xn--; 
          }
        else if ((xn == 0) || ((yn != 0) && (yp->id < xp->id)))
          { xt = 0; yt = yp->coef; 
            yp++; yn--; 
          }
        else /* (xn != 0) && (yn != 0) && (yp->id == xp->id) */
          { xt = xp->coef; 
            yt = yp->coef;
            xp++; xn--; yp++; yn--; 
          }
        if ((xt != 0) || (yt != 0))
          {
            if ((yt < 0) || ((yt == 0) && (xt < 0))) { xt = -xt; yt = -yt; }
            /* fprintf(stderr, "  + %ld = (%f,%f)\n", n, xt, yt); */
            xv[n] = (double)xt; yv[n] = (double)yt;
            n++;
          }
      }
    (*nv) = n;
  }

#define ANG_GREATER(xa,ya,xb,yb) ((xb)*(ya) - (xa)*(yb) > 0)

void aa_2d_range_sort_pairs
  ( AATermCount nv,
    double *xv,
    double *yv
  )
  {
    /* Should use a faster sort... */
    /* fprintf(stderr, "sorting pairs\n"); */
    for (int32_t i = 1; i < nv; i++)
      { int32_t j = i, k = j-1;
        double xi = xv[i], yi = yv[i];
        while ((k >= 0) && ANG_GREATER(xi, yi, xv[k], yv[k]))
          { xv[j] = xv[k]; yv[j] = yv[k]; j = k; k--; }
        if (j < i) { xv[j] = xi; yv[j] = yi; }
      }
  }

void aa_2d_range_compute_vertices
  ( Float x0,
    Float y0,
    AATermCount nv,
    double *xv,
    double *yv
  )
  {
    double xs = 0.0, ys = 0.0;
    
    /* fprintf(stderr, "computing vertices (nv = %ld)\n", nv); */
    /* Compute the sum of all vertices: */
    for (int32_t i = 0; i < nv; i++)
      { xs += xv[i]; ys += yv[i]; }
    /* Fill the second half with the partial sums: */
    for (int32_t i = 0; i < nv; i++) 
      { xv[(int32_t)nv+i] = xs; xs -= 2*xv[i];
        yv[(int32_t)nv+i] = ys; ys -= 2*yv[i]; 
      }
    /* Fill the first half with the negated vertices: */
    for (int32_t i = 0; i < nv; i++) 
      { xv[i] = -xv[(int32_t)nv+i]; yv[i] = -yv[(int32_t)nv+i]; }
    /* Shift everything by the center: */
    for (int32_t i = 0; i < 2*nv; i++) 
      { xv[i] += (double)x0; yv[i] += (double)y0; }
  }

void aa_2d_range
  ( AAP x, 
    AAP y, 
    AATermCount *nv,
    double *xv,
    double *yv
  )
  {
    aa_2d_range_get_pairs(x, y, nv, xv, yv);
    aa_2d_range_sort_pairs(*nv, xv, yv);
    aa_2d_range_compute_vertices(x->center, y->center, *nv, xv, yv);
    (*nv) = 2*(*nv);
  }    
