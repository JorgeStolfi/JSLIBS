/* See uint16_image_GRAY_medcut.h
** Last edited on 2013-10-21 00:32:14 by stolfilocal
**
** Copyright (C) 1989, 1991 by Jef Poskanzer. See note at end of file.
*/

#include <bool.h>
#include <affirm.h>
#include <jspnm.h>
#include <uint16_image_GRAY_medcut.h>

/* Internal types */

typedef struct 
  { uint16_t lo;       /* Lowest gray in box */
    uint16_t hi;       /* Highest gray in box */
    float center;         /* The ideal box centroid */
    uint16_t rep;      /* The box representative */
    float spread;         /* Total quantization error for entries in box */
  } uint16_image_GRAY_box;

typedef uint16_image_GRAY_box* uint16_image_GRAY_box_vector;

/* Internal prototypes */

typedef int32_t (*qcomparefn)(const void *, const void *);

extern void 
qsort (void *items, size_t nitems, size_t itemsize, qcomparefn cmp);
        
static void
uint16_image_GRAY_set_box(
    uint16_image_GRAY_box_vector b,
    long *gh,
    uint16_t lo, uint16_t hi,
    uint16_t maxval
  );

static void
uint16_image_GRAY_box_center(
    long *gh,
    uint16_t lo, uint16_t hi,
    uint16_t maxval,
    float *centerp,
    uint16_t *repp
  );

static float
uint16_image_GRAY_box_spread(
    long *gh,
    uint16_t lo,
    uint16_t hi,
    uint16_t rep,
    uint16_t maxval
  );

static int32_t
uint16_image_GRAY_spread_compare(const uint16_image_GRAY_box *b1, const uint16_image_GRAY_box *b2);

static int32_t
uint16_image_GRAY_sample_compare(const uint16_t *g1, const uint16_t *g2);

static uint16_image_GRAY_box_vector
uint16_image_GRAY_box_vector_new(int32_t n);

static uint16_t*
uint16_image_GRAY_pixel_vector_new(int32_t n);

/* Implementation */

extern pgm_pixel_vector
uint16_image_GRAY_median_cut(
    long *gh,      /* Pixel count indexed by gray level */
    uint16_t maxval,   /* maxval of pixel values, also size of "gh" */
    uint32_t *newgraysp /* In: desired number of grays, out: number chosen */
  )
  /*
    Here is the fun part, the median-cut graymap generator.  This is
    based on Paul Heckbert's paper "Color Image Quantization for Frame
    Buffer Display", SIGGRAPH '82 Proceedings, page 297.
  */
  {
    pgm_pixel_vector gm;     /* The graymap */
    uint16_image_GRAY_box_vector bv;  /* The box list */
    register int32_t bi;
    register uint16_t g, lo, hi;
    float ctr;
    uint32_t boxes;
    
    bv = uint16_image_GRAY_box_vector_new((int32_t)*newgraysp);

    /* Set up the initial box. */
    
    uint16_image_GRAY_set_box(&(bv[0]), gh, 0, maxval, maxval);

    boxes = 1;

    /* Main loop: split boxes until we have enough. */
    while (boxes < (*newgraysp))
      {
        /* Box bv[0] should have the largest spread: */
        if (bv[0].spread <= 0.0) break; /* exact graymap! */
        lo = bv[0].lo;
        hi = bv[0].hi;
        if (lo >= hi) pnm_error("bad uint16_image_GRAY_box spread");
        ctr = bv[0].center;
        /* Split box at center: */
        if (ctr < (float)lo) pnm_error("bad uint16_image_GRAY_box center");
        if (ctr >= (float)hi) pnm_error("bad uint16_image_GRAY_box center");
        g = (uint16_t)ctr;
        uint16_t g1 = (uint16_t)(g+1);
        uint16_image_GRAY_set_box(&(bv[0]),     gh, lo,  g,  maxval);
        uint16_image_GRAY_set_box(&(bv[boxes]), gh, g1,  hi, maxval);
        ++boxes;

        /* Sort to bring the worst boxes to the top. */
        qsort((void*) bv, boxes, sizeof(uint16_image_GRAY_box), (qcomparefn) uint16_image_GRAY_spread_compare);
      }

    /* Ok, we've got enough boxes.  Collect their centroids: */
    gm = uint16_image_GRAY_pixel_vector_new((int32_t)boxes);
    for (bi = 0; bi < boxes; ++bi)
      { gm[bi] = bv[bi].rep; }
    qsort((void*)gm, boxes, sizeof(uint16_t), (qcomparefn) uint16_image_GRAY_sample_compare);
    (*newgraysp) = boxes;
    free(bv);
    return (gm);
  }

static void
uint16_image_GRAY_set_box(
    uint16_image_GRAY_box_vector b,
    long *gh,
    uint16_t lo, uint16_t hi,
    uint16_t maxval
  )
  /*
    Sets lo and hi, and computes center, rep, spread
  */
  {
    b->lo = lo;
    b->hi = hi;
    uint16_image_GRAY_box_center(gh, lo, hi, maxval, &(b->center), &(b->rep));
    b->spread = uint16_image_GRAY_box_spread(gh, lo, hi, b->rep, maxval);
  }

static void
uint16_image_GRAY_box_center(
    long *gh,
    uint16_t lo, uint16_t hi,
    uint16_t maxval,
    float *centerp,
    uint16_t *repp
  )
  /*
    Compute a "central" gray level for this box,
    and rounds it to a "representative" level.
  */
  {
    long w;
    float fg, fw;
    float sumg = 0.0, sumw = 0.0;
    uint16_t maxg = 0;
    uint16_t ming = maxval;
    uint16_t g;
#ifdef USE_LUM
    float offset = 0.05 * ((float)maxval);
#endif /*USE_LUM*/
    for (g = lo; g <= hi; ++g)
      { w = gh[g];
        if (w > 0)
          { 
            if (g < ming) ming = g;
            if (g > maxg) maxg = g;
            fg = (float)(g);
            fw = (float)(w);
#ifdef USE_LUM
            fg = fg/(fg + offset);        
#endif /*USE_LUM*/
            sumg += fw*fg; sumw += fw;
          }
      }
    if (sumw <= 0.0) pnm_error("empty uint16_image_GRAY_box");
    fg = sumg / sumw;
#ifdef USE_LUM
    fg = (offset*fg)/(1.0 - fg);
#endif /*USE_LUM*/
    g = (uint16_t)(fg);
    if (g > maxg) g = maxg;
    if (g < ming) g = ming;
    (*centerp) = fg;
    (*repp) = g;
  }

static float
uint16_image_GRAY_box_spread(
    long *gh,
    uint16_t lo,
    uint16_t hi,
    uint16_t rep,
    uint16_t maxval
  )
  /* 
    Compute the box spread, relative to the 
    appointed representative.
  */
  { 
    float fg, fw;
    float sumgg = 0.0;
    uint16_t g;
    long w;
    float frep = (float)rep;
#ifdef USE_LUM
    float offset = 0.05 * ((float)maxval);
    frep = frep/(frep + offset);
#endif /*USE_LUM*/
    for (g = lo; g <= hi; ++g)
      { w = gh[g];
        if (w > 0)
          { 
            fg = (float)(g);
            fw = (float)(w);
#ifdef USE_LUM
            fg = fg/(fg + offset);        
#endif /*USE_LUM*/
            fg = fg - frep;
            sumgg += fw*(fg*fg);
          }
      }
    return(sumgg);
  }

static int32_t
uint16_image_GRAY_spread_compare(const uint16_image_GRAY_box *b1, const uint16_image_GRAY_box *b2)
  {
    return ((int32_t)
      (b1->spread > b2->spread ? -1 : 
      (b1->spread < b2->spread ?  1 :
      0))
    );
  }

static int32_t
uint16_image_GRAY_sample_compare(const uint16_t *g1, const uint16_t *g2)
  {
    return ((int32_t) (*g1) - (*g2));
  }

static uint16_image_GRAY_box_vector
uint16_image_GRAY_box_vector_new(int32_t n)
  { uint16_image_GRAY_box_vector bv = talloc(n, uint16_image_GRAY_box);
    return(bv);
  }

static uint16_t*
uint16_image_GRAY_pixel_vector_new(int32_t n)
  { uint16_t *gv = talloc(n, uint16_t);
    return(gv);
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

