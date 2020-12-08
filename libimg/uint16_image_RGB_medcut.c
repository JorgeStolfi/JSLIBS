/* See ppm_medcutx.h
** Last edited on 2017-06-21 00:32:28 by stolfilocal
**
** Copyright (C) 1989, 1991 by Jef Poskanzer. See note at end of file.
*/

#include <affirm.h>
#include <jspnm.h>

#include <uint16_image_RGB_medcut.h>
#include <uint16_image_RGB_hist.h>
#include <uint16_image_RGB_table.h>

/* Options for the median cut algorithm: */

/* How to choose the axis of greatest variation: */

#define LARGE_LUM
/* #define LARGE_NORM */

/* How to choose the representative color in each box: */

#define REP_AVERAGE_PIXELS
/* #define REP_CENTER_BOX */
/* #define REP_AVERAGE_COLORS */

#define vmin(X,Y) ((X)<(Y) ? (X) : (Y))
#define vmax(X,Y) ((X)>(Y) ? (X) : (Y))

/* Internal types */

typedef struct 
  { int lo, hi;       /* Index range in color histogram */
    float spread;     /* Box size */
    int longaxis;     /* Axis of maximum spread (0,1,2) for (R,G,B) */
    float cR, cG, cB; /* Box centroid */
    ppm_pixel_t rep;  /* Pixel representative */
  } ppm_box;
  
typedef ppm_box* ppm_box_vector;

typedef int (*qcomparefn)(const void *, const void *) ;

typedef int (*matchfn)(long R, long G, long B, uint16_image_RGB_hist_vector cm, int ncolors);
  /* A function that finds the best match to {(R,G,B)} in the colormap {cm}.
    Note that {(R,G,B)} may be somewhat outside of the color cube. */

/* Internal prototypes */

void ppm_set_box(
    ppm_box_vector b,
    uint16_image_RGB_hist_vector ch,
    int lo, int hi,
    uint16_t maxval
  );
  /* Sets {lo} and {hi}, and computes {center}, {rep}, {longaxis}, 
    and {spread}. */

void ppm_box_center(
    uint16_image_RGB_hist_vector ch,
    int lo, int hi,
    uint16_t maxval,
    float *cRp, float *cGp, float *cBp,
    ppm_pixel_t *repp
  );
  /*
    Computes a "central" rgb value for this box,
    and rounds it to a "representative" pixel.
  */

void ppm_box_spread(
    uint16_image_RGB_hist_vector ch,
    int lo,
    int hi,
    ppm_pixel_t rep,
    uint16_t maxval,
    int *laxp,
    float *spreadp
  );
  /* Compute the box spread {*spreadp} and largest 
    dimension {*laxp}, relative to the appointed representative {rep}. */

int uint16_image_RGB_red_compare(const uint16_image_RGB_hist_vector ch1, const uint16_image_RGB_hist_vector ch2);

int uint16_image_RGB_green_compare(const uint16_image_RGB_hist_vector ch1, const uint16_image_RGB_hist_vector ch2);

int uint16_image_RGB_blue_compare(const uint16_image_RGB_hist_vector ch1, const uint16_image_RGB_hist_vector ch2);

int ppm_spread_compare(const ppm_box *b1, const ppm_box *b2);

ppm_box_vector new_ppm_box_vector(int n);

uint16_image_RGB_hist_vector uint16_image_RGB_hist_vector_new(int n);

/* Implementation */

uint16_image_RGB_hist_vector uint16_image_RGB_median_cut
  ( uint16_image_RGB_hist_vector ch,
    int colors, 
    uint16_t maxval,
    int *newcolorsp
  )
  {
    uint16_image_RGB_hist_vector cm; /* The colormap */
    ppm_box_vector bv;   /* The box tree */
    register int i;
    int boxes;

    bv = new_ppm_box_vector(*newcolorsp);

    /* Set up the initial box. */
    
    ppm_set_box(&(bv[0]), ch, 0, colors, maxval);
    boxes = 1;

    /* Main loop: split boxes until we have enough. */
    while (boxes < (*newcolorsp))
      { 
        register int lo, hi;
        int lax;
        float ctrv;
        qcomparefn cmp;
        
        /* Box 0 should have the largest spread. */
        if (bv[0].spread <= 0.0) break; /* exact */
        lo = bv[0].lo;
        hi = bv[0].hi;
        if (lo >= hi) pnm_error("bad ppm_box spread");
        lax = bv[0].longaxis;
        
        /* Sort box colors in the direction of largest spread: */
        switch(lax)
          { 
            case 0: cmp = (qcomparefn) uint16_image_RGB_red_compare;   ctrv = bv[0].cR; break;
            case 1: cmp = (qcomparefn) uint16_image_RGB_green_compare; ctrv = bv[0].cG; break;
            case 2: cmp = (qcomparefn) uint16_image_RGB_blue_compare;  ctrv = bv[0].cB; break;
            default: pnm_error("huh?"); cmp = NULL; ctrv = 0;
          }
        qsort((void*)&(ch[lo]), hi-lo+1, sizeof(struct uint16_image_RGB_hist_item), cmp); 
        
        /* Split colors in two groups at centroid: */
        for (i = lo; i <= hi; ++i)
          { uint16_t v = ch[i].color.c[lax];
            if ((float)v >= ctrv) break;
          }
        if ((i > hi) || (i <= lo)) pnm_error("bad centroid in pgm box");
        ppm_set_box(&(bv[0]), ch, lo, i-1, maxval);
        ppm_set_box(&(bv[boxes]), ch, i, hi, maxval);
        ++boxes;

        /* Sort to bring the worst boxes to the top. */
        qsort((void*) bv, boxes, sizeof(ppm_box), (qcomparefn) ppm_spread_compare);
      }

    /* Ok, we've got enough boxes.  Collect their centroids: */
    cm = uint16_image_RGB_hist_vector_new(boxes);
    for (i = 0; i < boxes; ++i) 
      { cm[i].color = bv[i].rep; }
    (*newcolorsp) = boxes;
    free(bv);
    return (cm);
  }

void ppm_set_box
  ( ppm_box_vector b,
    uint16_image_RGB_hist_vector ch,
    int lo, int hi,
    uint16_t maxval
  )
  {
    b->lo = lo;
    b->hi = hi;
    ppm_box_center(ch, lo, hi, maxval, &(b->cR), &(b->cG), &(b->cB), &(b->rep));
    ppm_box_spread(ch, lo, hi, b->rep, maxval, &(b->longaxis), &(b->spread));
  }

void ppm_box_center
  ( uint16_image_RGB_hist_vector ch,
    int lo, int hi,
    uint16_t maxval,
    float *cRp, float *cGp, float *cBp,
    ppm_pixel_t *repp
  )
  {
    double sumv[3];
    uint16_t maxv[3];
    uint16_t minv[3];
    double sumw = 0;
    ppm_pixel_t rep;
    int i, c;
    
    for (c = 0; c < 3; c++)
      { maxv[c] = 0; minv[c] = maxval; sumv[c] = 0; }
    
    for (i = lo; i <= hi; ++i)
      { double fw = (double)ch[i].value;
        for (c = 0; c < 3; c++)
          { uint16_t v = ch[i].color.c[c];
            minv[c] = vmin(minv[c], v);
            maxv[c] = vmax(maxv[c], v);
            sumv[c] += fw*(double)v; 
          }
        sumw += fw; 
      }
    if (sumw <= 0.0) pnm_error("empty uint16_image_GRAY_box");

    float rgb[3];
    for (c = 0; c < 3; c++)
      { if (minv[c] >= maxv[c])
          { rgb[c] = minv[c]; 
            rep.c[c] = minv[c]; }
        else
          { double fv = sumv[c] / sumw; 
            int64_t zv = (int64_t)(fv + 0.5);
            if (zv > maxv[c]) { zv = maxv[c]; }
            if (zv < minv[c]) { zv = minv[c]; }
            rgb[c] = (float)fv; 
            rep.c[c] = (uint16_t)zv;
          }
      }
    (*cRp) = rgb[0];
    (*cGp) = rgb[1];
    (*cBp) = rgb[2];
    (*repp) = rep;
    /*  
      pm_message("[%d..%d] = [%d-%d]x[%d-%d]x[%d-%d] ctr=(%.1f,%.1f,%.1f) rep=(%d,%d,%d)",
        lo, hi,
        minr, maxr, ming, maxg, minb, maxb,
        (*cRp), (*cGp), (*cBp),
        PPM_GETR(rep), PPM_GETG(rep), PPM_GETB(rep)
      );
    */
  }

void ppm_box_spread
  ( uint16_image_RGB_hist_vector ch,
    int lo,
    int hi,
    ppm_pixel_t rep,
    uint16_t maxval,
    int *laxp,
    float *spreadp
  )
  {
    float fv, fw;
    float sumRR = 0.0, sumGG = 0.0, sumBB = 0.0;
    uint16_t v;
    long w;
    register int i;
    float frepR = (float)rep.c[0];
    float frepG = (float)rep.c[1];
    float frepB = (float)rep.c[2];
    for (i = lo; i <= hi; ++i)
      {
        w = ch[i].value;
        fw = (float)w;
        
        v = ch[i].color.c[0];
        fv = (float)v - frepR;
        sumRR += fw*fv*fv;
        
        v = ch[i].color.c[1];
        fv = (float)v - frepG;
        sumGG += fw*fv*fv;
        
        v = ch[i].color.c[2];
        fv = (float)v - frepB;
        sumBB += fw*fv*fv;
      }
    if ((sumRR >= sumGG) && (sumRR >= sumBB))
      { (*laxp) = 0; 
        (*spreadp) = sumRR;
      }
    else if ((sumGG >= sumRR) && (sumRR >= sumBB))
      { (*laxp) = 1; 
        (*spreadp) = sumGG;
      }
    else /* if ((sumbb >= sumgg) && (sumbb >= sumrr)) */
      { (*laxp) = 2; 
        (*spreadp) = sumBB;
      }
  }

int uint16_image_RGB_red_compare(uint16_image_RGB_hist_vector ch1, uint16_image_RGB_hist_vector ch2)
  {
    return (int) ch1->color.c[0] - (int) ch2->color.c[0];
  }

int uint16_image_RGB_green_compare(uint16_image_RGB_hist_vector ch1, uint16_image_RGB_hist_vector ch2)
  {
    return (int) ch1->color.c[1] - (int) ch2->color.c[1];
  }

int uint16_image_RGB_blue_compare(uint16_image_RGB_hist_vector ch1, uint16_image_RGB_hist_vector ch2)
  {
    return (int) ch1->color.c[2] - (int) ch2->color.c[2];
  }

int ppm_spread_compare(const ppm_box *b1, const ppm_box *b2)
  {
    return ((int)
      (b1->spread > b2->spread ? -1 : 
      (b1->spread < b2->spread ?  1 :
      0))
    );
  }

ppm_box_vector new_ppm_box_vector(int n)
  { ppm_box_vector bv = (ppm_box_vector)pnm_malloc(n*sizeof(ppm_box));
    return(bv);
  }

uint16_image_RGB_hist_vector uint16_image_RGB_hist_vector_new(int n)
  { uint16_image_RGB_hist_vector ch = (uint16_image_RGB_hist_vector)pnm_malloc(n*sizeof(struct uint16_image_RGB_hist_item));
    return(ch);
  }

/* Copyright (C) 1989, 1991 by Jef Poskanzer.
**
** Permission to use, copy, modify, and distribute this software and its
** documentation for any purpose and without fee is hereby granted, provided
** that the above copyright notice appear in all copies and that both that
** copyright notice and this permission notice appear in supporting
** documentation.  This software is provided "as is" without express or
** implied warranty.
*/

