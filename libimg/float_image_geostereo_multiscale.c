/* See {float_image_geostereo_multiscale.h}. */
/* Last edited on 2017-06-25 04:55:23 by stolfilocal */

#include <assert.h>
#include <limits.h>
#include <string.h>
#include <math.h>
 
#include <bool.h>
#include <float_image.h>
#include <float_image_mscale.h>
#include <float_image_geostereo_uniscale.h>

#include <float_image_geostereo_multiscale.h>

#define ALPHA (0.5)
  /* Relative weight of higher-scale scores in total score. */

float float_image_geostereo_interpolate(float sa, float sb, float sc, float sd, double fd);
  /* Given four consecutive samples {sa,sb,sc,sd}, returns the interpolated 
    image value at a point {fd} of the way between samples {sb} and {sc}.
    Must be called with {fd} in {[0_1]} only. */

void float_image_geostereo_normalize_samples(float w[], int npix, int NC);
  /* Independently normalizes each channel of the given samples {w[0..NC*npix-1]} to have
    mean 0 and unit variance. Ignores {NAN} samples. If all samples
    are equal, sets them all to 0. */

void float_image_geostereo_multiscale_displacement_map
  ( float_image_t *f1,  /* Image 1. */
    float_image_t *f2,  /* Image 2. */
    int nscales,        /* Number of scales to consider (0 = uniscale). */
    int ncands,         /* Number of candidates to keep. */
    int rx,             /* Window half-width. */
    int ry,             /* Window half-height. */
    int dmin,           /* Minimum signed displacement (pixels). */
    int dmax,           /* Maximum signed displacement (pixels). */
    float_image_t **fd, /* (OUT) Dispmap image. */
    float_image_t **fs  /* (OUT) Scoremap image. */
  )
  {
    int fNX = (int)f1->sz[1];
    int fNY = (int)f1->sz[2];
    if (nscales <= 0)
      { float_image_geostereo_uniscale_displacement_map
          ( f1, f2, 
            /* ncands: */ ncands, 
            /* rx,ry: */ rx,ry,  
            /* dmin,dmax: */ dmin, dmax,
            fd, fs
          );
      }
    else
      { /* Scale images to half-size: */
        int gNX = (fNX+1)/2;
        int gNY = (fNY+1)/2;
        float_image_t *g1 = float_image_mscale_shrink(f1, NULL, gNX, gNY, 1, 1, 3);
        float_image_t *g2 = float_image_mscale_shrink(f2, NULL, gNX, gNY, 1, 1, 3);
        float_image_t *gd;  /* Displacement map. */
        float_image_t *gs;  /* Score map. */
        
        /* Compute dispmaps for {g1,g2}, with twice as many cands: */
        float_image_geostereo_multiscale_displacement_map
          ( g1, g2,
            /* nscales: */ nscales-1, 
            /* ncands: */ 2*ncands,
            /* rx,ry: */ (rx-1)/2, (ry-1)/2,
            /* dmin,dmax: */ dmin/2, (dmax+1)/2,
            &gd, &gs
          );
        float_image_free(g1); float_image_free(g2);
        
        /* Translate displacements, refine, keep {ncands} best ones: */
        /* !!! Check for 0.5 offset !!! */
        (*fd) = float_image_new(ncands, fNX, fNY);
        (*fs) = float_image_new(ncands, fNX, fNY);
        float_image_geostereo_refine_and_prune_displacement_map
          ( gd, gs, 
            f1, f2,
            /* rx,ry: */ (rx-1)/2, (ry-1)/2,
            /* dmin,dmax: */ dmin, dmax,
            (*fd), (*fs)
          );
        float_image_free(gd); float_image_free(gs);
      }
  }

#define HDEBUG 146
#define VDEBUG 055
  /* Print debugging information for this pixel. */
  
void float_image_geostereo_uniscale_displacement_map
  ( float_image_t *f1,  /* Image 1. */
    float_image_t *f2,  /* Image 2. */
    int ncands,         /* Number of candidates to keep. */
    int rx,             /* Window half-width. */
    int ry,             /* Window half-height. */
    int dmin,           /* Minimum signed displacement (pixels). */
    int dmax,           /* Maximum signed displacement (pixels). */
    float_image_t **fd, /* (OUT) Dispmap image. */
    float_image_t **fs  /* (OUT) Scoremap image. */
  )
  {
    int NC = (int)f1->sz[0];
    int NX = (int)f1->sz[1];
    int NY = (int)f1->sz[2];
    
    assert(((int)f2->sz[0]) == NC);
    assert(((int)f2->sz[1]) == NX);
    assert(((int)f2->sz[2]) == NY);

    int wx = 2*rx + 1, wy = 2*ry + 1;
    int x, y;
    
    /* Window buffers: */
    int maxsamp = wx*wy*NC;
    float w1[maxsamp], w2[maxsamp];

    /* Allocate dispmap and scoremap: */
    (*fd) = float_image_new(ncands, NX, NY);
    (*fs) = float_image_new(ncands, NX, NY);

    /* Best displacements for {f1,f2} and their scores: */
    float disp[ncands];
    float scor[ncands];

    /* Compute dispmap/scoremap: */
    for (y = 0; y < NY; ++y)
      { for (x = 0; x < NX; ++x)
         { int j, id;
           int debug = ((x == HDEBUG) && (y == VDEBUG));
           /* Initialize {disp,scor} with nulls: */
           for (j = 0; j < ncands; j++) { scor[j] = INFINITY; disp[j] = 0; }
           /* Compute the {ncands} best int disps {disp,scor} for this pix: */
           for (id = 3*dmin; id <= 3*dmax; id += 3)
             { int dbest; double sbest;
               if (debug) { fprintf(stderr, "  id = %d\n", id); }
               float_image_geostereo_local_match(f1,f2, x,y, id-1,id+1, rx,ry, &dbest,&sbest, w1,w2);
               if (debug) { fprintf(stderr, "  d = %d s = %7.4f", dbest, sbest); }
               float_image_geostereo_insert_disp(dbest, (float)sbest, disp, scor, ncands);
               if (debug) 
                 { fprintf(stderr, " dbest = %.0f sbest = %7.4f\n", disp[0], scor[0]); }
             }
           float_image_set_pixel((*fd), x, y, disp);
           float_image_set_pixel((*fs), x, y, scor);
         }
      }
  }

void float_image_geostereo_refine_and_prune_displacement_map
  ( float_image_t *gd,  /* Displacement map for halfsize images. */
    float_image_t *gs,  /* Score map for halfsize images. */
    float_image_t *f1,  /* Full-size image 1. */
    float_image_t *f2,  /* Full-size image 2. */
    int rx,             /* Window half-width. */
    int ry,             /* Window half-height. */
    int dmin,           /* Minimum signed displacement (pixels). */
    int dmax,           /* Maximum signed displacement (pixels). */
    float_image_t *fd,  /* (OUT) Dispmap for full-size images. */
    float_image_t *fs   /* (OUT) Scoremap for full-size images. */
  )
  {
    int NC = (int)f1->sz[0];

    int fNX = (int)f1->sz[1];
    int fNY = (int)f1->sz[2];

    int gNX = (int)gd->sz[1];
    int gNY = (int)gd->sz[2];
    
    int gncands = (int)gd->sz[0];
    int fncands = (int)fd->sz[0];
    
    int wx = 2*rx + 1, wy = 2*ry + 1;
    int fx, fy;

    /* Window buffers: */
    int maxsamp = wx*wy*NC;
    float w1[maxsamp], w2[maxsamp];
    
    /* Best displacement for {f1,f2} and their scores: */
    float fdisp[fncands];
    float fscor[fncands];

    /* Best displacements for {g1,g2} and their scores: */
    float gdisp[gncands];
    float gscor[gncands];

    fprintf(stderr, "expanding from %d×%d to %d×%d...\n", gNX,gNY,fNX,fNY); 
    for (fy = 0; fy < fNY; ++fy)
      { for (fx = 0; fx < fNX; ++fx)
          { int gx = fx/2, gy = fy/2;
            assert((gy >= 0) && (gy < gNY));
            assert((gx >= 0) && (gx < gNX));
            float_image_get_pixel(gd, gx, gy, gdisp);
            float_image_get_pixel(gs, gx, gy, gscor);
            int fj, gj;
            /* Initialize {fdisp,fscor} with nulls: */
            for (fj = 0; fj < fncands; fj++) { fscor[fj] = INFINITY; fdisp[fj] = 0; }
            /* Copy and adjust best {fncands} disps found at higher level: */
            for (gj = 0; gj < fncands; gj++)
              { int id = 2*(int)(gdisp[gj]);
                int dbest;
                double sbest;
                float_image_geostereo_local_match(f1,f2, fx,fy, id, id+1, rx,ry, &dbest,&sbest, w1,w2);
                double stot = (1.0-ALPHA)*sbest + ALPHA*gscor[gj];
                float_image_geostereo_insert_disp(dbest, (float)stot, fdisp, fscor, fncands);
              }
            float_image_set_pixel(fd, fx, fy, fdisp);
            float_image_set_pixel(fs, fx, fy, fscor);
          }
      }
  }

/* Smoothing weights ({WT0} = central pixel, {WT1} = adjacent): */
#define WT0 (0.425)
#define WT1 (0.075)

/*
#define WTERR (2.0*sqrt(WT1*WT1*WT1*WT1 + 2.0*WT0*WT0*WT1*WT1 + WT0*WT0*WT0*WT0))
*/

void float_image_geostereo_local_match
  ( float_image_t *f1, /* Image 1. */
    float_image_t *f2, /* Image 2. */
    int x,             /* Central column (origin for displacement). */
    int y,             /* Current row index in image. */
    int dmin,          /* Min displacement, in 1/3 pixels. */
    int dmax,          /* Max displacement, in 1/3 pixels. */
    int rx,            /* Window half-width. */
    int ry,            /* Window half-height. */
    int *dbest,        /* Adjusted displacement, in 1/3 pixels. */
    double *sbest,     /* Score (squared mismatch) for {dbest}. */
    float *w1,         /* Buffer for image 1 window samples. */
    float *w2          /* Buffer for image 2 window samples. */
  )
  { 
    int NC = (int)f1->sz[0];
    int NX = (int)f1->sz[1];
    int npix = (2*rx+1)*(2*ry+1);
    int nsamp = npix*NC;
    int d;
    int debug = ((x == HDEBUG) && (y == VDEBUG));
    (*sbest) = INFINITY;
    for (d = dmin; d <= dmax; d++)
      { double s = 0.0; 
        int nok = 0, i;
        if (debug) { fprintf(stderr, "    d = %d\n", d); }
        float_image_geostereo_get_samples(f1, x, y, +d, rx, ry, w1);
        float_image_geostereo_get_samples(f2, x, y, -d, rx, ry, w2);
        if (debug) 
          { float_image_geostereo_debug_window(w1, rx, ry, NX); 
            float_image_geostereo_debug_window(w2, rx, ry, NX);
          }
        /* Normalize samples for zero mean and unit variance: */
        float_image_geostereo_normalize_samples(w1, npix, NX);
        float_image_geostereo_normalize_samples(w2, npix, NX);
        if (debug) 
          { float_image_geostereo_debug_window(w1, rx, ry, NX); 
            float_image_geostereo_debug_window(w2, rx, ry, NX);
          }
        /* Compute discrepancy, ignoring missing samples: */
        for (i = 0; i < nsamp; i++) 
          { double a = w1[i], b = w2[i];
            if ((!isnan(a)) && (fabs(a) != INF) && (! isnan(b)) && (fabs(b) != INF))
              { double e = a - b;
                s += e*e; 
                nok ++;
              }
          }
        /* Adjust discrepancy to account for missing samples: */
        if (nok > 0) { s *= (double)nsamp/(double)nok; }
        /* Save best discrepancy: */
        if (s < (*sbest)) { (*dbest) = d; (*sbest) = s; }
        if (debug) { fprintf(stderr, "    s = %7.4f\n", s); }
      }
  }

void float_image_geostereo_get_samples
  ( float_image_t *f,  /* Pixel row buffer for image 1. */
    int x,             /* Central column (origin for displacement). */
    int y,             /* Current row index in image. */
    int d,             /* The displacement, in 1/3 pixels. */
    int rx,            /* Window half-width. */
    int ry,            /* Window half-height. */
    float *w           /* (OUT) Window sample buffer. */
  )
  { 
    int NC = (int)f->sz[0];
    int NX = (int)f->sz[1];
    int NY = (int)f->sz[2];

    int rd = (d + 3000000) % 3; /* Fractional displacement in 1/3 pixs. */
    int id = (d - rd)/3;        /* Displacement in whole pixels. */
    int nw = 0, ic, iy, ix;
    /* Index range for interpolation (rel to {id}): */
    int ixlo = -1;
    int ixhi = (rd == 0 ? +1 : +2);
    for (iy = y-ry; iy <= y+ry; iy++)
      { bool_t iyok = ((iy >= 0) || (iy < NY));
        for (ix = x+id-rx; ix <= x+id+rx; ix++)
          { bool_t ixok = ((ix+ixlo >= 0) || (ix+ixhi < NX));
            for (ic = 0; ic < NC; ic++)
              { if (! (iyok && ixok)) 
                  { /* Out of bounds (partially or totally): */
                    w[nw] = NAN;
                  }
                else if (rd == 0)
                  { /* Exact index: */
                    w[nw] = float_image_get_sample(f, ic, ix, iy);
                  }
                else
                  { /* Fractional index, float_image_geostereo_interpolate: */
                    float sa = float_image_get_sample(f, ic, ix - 1, iy);
                    float sb = float_image_get_sample(f, ic, ix + 0, iy);
                    float sc = float_image_get_sample(f, ic, ix + 1, iy);
                    float sd = float_image_get_sample(f, ic, ix + 2, iy);
                    w[nw] = float_image_geostereo_interpolate(sa, sb, sc, sd, rd);
                  }
                nw++;
              }
          }
      }
  }

void float_image_geostereo_normalize_samples(float *w, int npix, int NC)
  { int i, c;
    for (c = 0; c < NC; c++)
      { double s; int nok;
        /* Shift so that mean of valid pixels is 0: */
        s = 0.0; nok = 0;
        for (i = c; i < npix; i += NC) 
          { double wi = w[i]; if (! isnan(wi)) { s += wi; nok++; } }
        if (nok == 0) { return; }
        s /= (double)nok;
        for (i = c; i < npix; i += NC)
          { if (! isnan(w[i])) { w[i] = (float)(w[i] - s); } }
        /* Scale so that variance of valid pixels is 1: */
        s = 0.0;
        for (i = c; i < npix; i += NC)
          { double wi = w[i]; if (! isnan(wi)) { s += wi*wi; } }
        s = sqrt(s/(double)nok);
        if (s == 0.0) { return; }
        for (i = c; i < npix; i += NC) 
          { if (! isnan(w[i])) { w[i] = (float)(w[i]/s); } }
      }
  }

float float_image_geostereo_interpolate(float sa, float sb, float sc, float sd, int rd)
  {
    if (rd == 1)
      { return (float)(- 0.1334*sa + 0.8169*sb + 0.3902*sc - 0.0507*sd); }
    else if (rd == 2)
      { return (float)(- 0.0507*sa + 0.3902*sb + 0.8169*sc - 0.1334*sd); }
    assert(0);
  }

void float_image_geostereo_insert_disp(int d, float s, float *disp, float *scor, int nd)
  { int j;
    /* Should not be called twice for the same {d} in the same pixel. */
    if (s < scor[nd-1])
      { j = nd;
        while((j > 0) && (scor[j-1] > s))
          { disp[j] = disp[j-1]; scor[j] = scor[j-1]; j--; }
        disp[j] = (float)d; scor[j] = s;
      }
  }

void float_image_geostereo_debug_window(float *w, int rx, int ry, int NC)
  {
    int c, x, y, k;
    k = 0;
    fprintf(stderr, "\n");
    for (y = -ry; y <= ry; y++)
      { fprintf(stderr, "    ");
        for (x = -rx; x <= rx; x++)
          { fprintf(stderr, " ");
            for (c = 0; c < NC; c++) 
              { fprintf(stderr, " %7.3f", w[k]); k++; }
          }
        fprintf(stderr, "\n");
      }
    fprintf(stderr, "\n");
  }
