/* See {float_image_geostereo.h}. */
/* Last edited on 2024-11-23 05:54:48 by stolfi */

#include <stdio.h>
#include <stdint.h>
#include <assert.h>
#include <limits.h>
#include <string.h>
#include <math.h>
 
#include <bool.h>
#include <affirm.h>
#include <ix_reduce.h>
#include <float_image.h>
#include <float_image_interpolate.h>

#include <float_image_geostereo.h>

void float_image_geostereo_single_pixel_best
  ( float_image_t *f1,  /* Image 1. */
    float_image_t *f2,  /* Image 2. */
    int32_t x,          /* Central horiz position. */
    int32_t y,          /* Current row index in image. */
    int32_t nwx,        /* Window width. */
    int32_t nwy,        /* Window height. */
    double wt[],        /* Window weights. */
    double dmin,        /* Minimum signed displacement (pixels). */
    double dmax,        /* Maximum signed displacement (pixels). */
    int32_t ncands,     /* Number of best displacements to keep. */
    double dbest[],     /* Best {ncands} displacements. */
    double sbest[],     /* Score (squared mismatch) for {dbest[0..ncands-1]}}. */
    bool_t debug,       /* If true, debugs computations. */
    double smp1[],      /* (WORK) Buffer for image 1 window samples. */
    double smp2[]       /* (WORK) Buffer for image 2 window samples. */
  )
  {
    if (debug) { fprintf(stderr, "  === pixel_best x = %d  y = %d ===\n", x, y); }
            
    /* Initialize the queue: */
    for (uint32_t k = 0;  k < ncands; k++) { dbest[k] = -INF; sbest[k] = +INF; }
    /* Enumerate displacements multiple of 1/3: */
    int32_t idmin = (int32_t)ceil(3.0*dmin*0.99999999);
    int32_t idmax = (int32_t)floor(3.0*dmax*1.00000001);
    double dprev = NAN;    /* Previous candidate displacement. */
    double sprev = +INF;   /* Its score. */
    bool_t mprev = FALSE;  /* True iff {sprev} was strictly smaller than the previous score. */
    for (int32_t id = idmin; id <= idmax; id++)
      { double d = ((double)id)/3.0; /* Candidate displacement. */
        if ((d >= dmin) && (d <= dmax))
          { if (debug) { fprintf(stderr, "  d = %10.3f", d); }
            double s = float_image_geostereo_single_disp_score(f1,f2,x,y,d,nwx,nwy,wt,debug,smp1,smp2);
            if (debug) { fprintf(stderr, "  score = %10.6f", s); }
            /* If there are two or more successive cands with equal score, take the first one */ 
            if (mprev && (sprev <= s))
              { /* Displacement {dprev} is a local minimum of score: */
                if (debug) { fprintf(stderr, " (min"); }
                bool_t rank = float_image_geostereo_queue_insert(dprev, sprev, ncands, dbest, sbest);
                if (debug && (rank < ncands)) { fprintf(stderr, " rank %d", rank); }
                if (debug) { fprintf(stderr, ")"); }
              }
            if (debug) { fprintf(stderr, "\n"); }
            /* Prepare for next displacement: */
            mprev = s < sprev;
            dprev = d;
            sprev = s;
          }
      }
    
    /* Check the last displacement for candidacy: */
    if (mprev) { float_image_geostereo_queue_insert(dprev, sprev, ncands, dbest, sbest); }
    if (debug) { float_image_geostereo_queue_dump(ncands, dbest, sbest); }
  }
  
double float_image_geostereo_single_disp_score
  ( float_image_t *f1,  /* Image 1. */
    float_image_t *f2,  /* Image 2. */
    uint32_t x,         /* Column index in the map domain. */
    uint32_t y,         /* Row index in the map domain. */
    double d,           /* Parallax displacement (pixels). */
    uint32_t nwx,       /* Window width. */
    uint32_t nwy,       /* Window height. */
    double wt[],        /* Window weights. */
    bool_t debug,       /* TRUE to debug the computations. */
    double smp1[],      /* (WORK) Buffer for window samples of image 1. */
    double smp2[]       /* (WORK) Buffer for window samples of image 2. */
  )
  { 
    /* Check image sizes: */
    int32_t NC, NX, NY;
    float_image_get_size(f1 ,&NC, &NX, &NY);
    float_image_check_size(f2, NC, NX, NY);

    /* Window parameters: */
    int32_t npix = nwx*nwy;
    demand((npix % 2) == 1, "window sizes must be odd");

    float_image_geostereo_get_samples(f1, x, y, +d, nwx, nwy, smp1);
    float_image_geostereo_get_samples(f2, x, y, -d, nwx, nwy, smp2);

    if (debug) { fprintf(stderr, "    x = %d  y = %d  d = %.2f\n", x, y, d); }
    if (debug) 
      { fprintf(stderr, "    extracted samples:\n");
        float_image_geostereo_debug_window(smp1, nwx, nwy, NC); 
        float_image_geostereo_debug_window(smp2, nwx, nwy, NC);
      }

    /* Compute discrepancy score, ignoring missing samples: */
    double score = float_image_geostereo_compute_score(smp1, smp2, nwx, nwy, wt, NC);
    if (debug) { fprintf(stderr, "    score = %7.4f\n", score); }
    
    return score;
  }

double float_image_geostereo_compute_score
  ( double smp1[], 
    double smp2[], 
    uint32_t nwx, 
    uint32_t nwy, 
    double wt[],
    uint32_t NC
  )
  { 
    /* Compute discrepancy, ignoring missing samples: */
    double sum_we2 = 0.0;
    double sum_w = 1.0e-200; /* To prevent divide by zero if all are {NAN}. */
    int32_t icxy = 0;  /* Sequential sample index. */
    int32_t kw = 0; /* Index of next elem of weight table. */
    for (uint32_t iy = 0;  iy < nwy; iy++) 
      { for (uint32_t ix = 0;  ix < nwx; ix++) 
          { double w = wt[kw];
            for (uint32_t ic = 0;  ic < NC; ic++) 
              { /* Get the two samples: */
                double s1 = smp1[icxy];
                bool_t ok1 = (!isnan(s1)) && (fabs(s1) != INF);
                double s2 = smp2[icxy];
                bool_t ok2 = (!isnan(s2)) && (fabs(s2) != INF);
                if (ok1 && ok2)
                  { double e = s2 - s1;
                    /* Accumulate: */
                    sum_we2 += w*e*e;
                    sum_w += w;
                  }
                icxy++;
              }
            kw++;
          }
      }
    /* Compute the mean diff squared: */
    double score = sum_we2/sum_w;
    return score;
  }

void float_image_geostereo_get_samples
  ( float_image_t *f,   /* Pixel row buffer for image 1. */
    uint32_t x,         /* Reference column index. */
    uint32_t y,         /* Reference row index. */
    double d,           /* Horizontal displacement (pixels, fractional). */
    uint32_t nwx,       /* Window width. */
    uint32_t nwy,       /* Window height. */
    double smp[]         /* (OUT) Window sample buffer. */
  )
  { 
    double ctrx = (double)x + d;
    double ctry = (double)y;
    int32_t hx = nwx/2;
    int32_t hy = nwy/2;
    double dx = 1.0;
    double dy = 1.0;
    ix_reduce_mode_t red = ix_reduce_mode_SINGLE; /* Surround image with {NAN} pixels. */
    int32_t order = 1; /* Interpolation continuity order (C1). */
    /* double undef = NAN; */ /* Value for pixels outside the domain. */
    float_image_interpolate_grid_pixels(f, ctrx, hx, dx, ctry, hy, dy, order, red, smp);
  }

void float_image_geostereo_debug_window
  ( double smp[],
    uint32_t nwx,       /* Window width. */
    uint32_t nwy,       /* Window height. */
    uint32_t NC         /* Channel count. */
  )
  {
    int32_t k = 0; /* Next element in {smp} array. */
    fprintf(stderr, "\n");
    for (uint32_t y = 0;  y < nwy; y++)
      { fprintf(stderr, "    ");
        for (uint32_t x = 0;  x < nwx; x++)
          { fprintf(stderr, " ");
            for (uint32_t c = 0;  c < NC; c++) 
              { fprintf(stderr, " %7.3f", smp[k]); k++; }
          }
        fprintf(stderr, "\n");
      }
    fprintf(stderr, "\n");
  }

int32_t float_image_geostereo_queue_insert(double d, double s, int32_t nq, double dq[], double sq[])
  { 
    if ((nq == 0) || (s >= sq[nq-1]))
      { return nq; }
    else
      { int32_t rk = nq-1; /* Tentative rank of {s} in {sq}. */
        /* Look for proper rank, pushing up larger entries: */
        while((rk > 0) && (s < sq[rk-1]))
          { dq[rk] = dq[rk-1]; sq[rk] = sq[rk-1]; rk--; }
        dq[rk] = d; sq[rk] = s;
        return rk;
      }
  }

void float_image_geostereo_queue_dump(int32_t nq, double dq[], double sq[])
  { fprintf(stderr, "  best candidates:\n");
    for (uint32_t rk = 0;  rk < nq; rk++)
      { double d = dq[rk];
        double s = sq[rk];
        if (! isnan(d))
          { fprintf(stderr, "  d[%02d] = %10.3f s[%02d] = %10.6f\n", rk, d, rk, s); } 
      }
  }
