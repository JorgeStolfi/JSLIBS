/* See {uint16_color_tree.h} */
/* Last edited on 2024-12-27 07:26:49 by stolfi */

/* Copyright (C) 1989, 1991 by Jef Poskanzer. See note at end of file. */

#include <affirm.h>
#include <jspnm.h>
#include <vec.h>

#include <uint16_color_table.h>

#include <uint16_color_tree.h>

/* Options for the median cut algorithm: */

/* How to choose the axis of greatest variation: */

#define LARGE_LUM
/* #define LARGE_NORM */

/* How to choose the representative color in each box: */

#define REP_AVERAGE_PIXELS
/* #define REP_CENTER_BOX */
/* #define REP_AVERAGE_COLORS */

typdef struct uct_queue_entry_t 
  { uint16_t *min;
    uint16_t *max;
    uint32_t cmax;
    uint16_color_tree_node_t *nds;
  } uct_queue_entry_t;
  /* Type of an entry in the queue used by the algorithm. Each entry
    {qe} entry contains a set of colors, represented by list {qe.nds} of
    (incomplete) tree nodes; a bounding box {qe.min[0..chns-1]} and
    {qe.max[0..chsn-1]} for those colors; and the index {qe.cmax} of the
    channel where the difference {max[c]-min[c]} is maximum.
    
    The tree nodes in the list {nds} are linked serially by the 
    {nd.sub[0]} fields; the fields {nd.sub[1]} are {NULL}, and the 
    other fields {nd.lim[0..1]}, {nc.csplit} are undefined. */


typedef struct ppm_box_t
  { int32_t lo, hi;       /* Index range in color histogram */
    float spread;         /* Box size */
    int32_t longaxis;     /* Axis of maximum spread (0,1,2) for (R,G,B) */
    float cR, cG, cB;     /* Box centroid */
    ppm_pixel_t rep;      /* Pixel representative */
  } ppm_box_t;
  
vec_typedef(ppm_box_vec_t, ppm_box_vec, ppm_box_t);

typedef int32_t qcomparefn_t(const void *, const void *);

/* Internal prototypes */

void ppm_set_box(
    ppm_box_t *b,
    uint16_color_table_t *ch,
    int32_t lo,
    int32_t hi,
    uint16_t maxval
  );
  /* Sets {lo} and {hi}, and computes {center}, {rep}, {longaxis}, 
    and {spread}. */

void ppm_box_center
  ( uint16_color_table_t *ch,
    int32_t lo, 
    int32_t hi,
    uint16_t maxval,
    float *cRp, float *cGp, float *cBp,
    ppm_pixel_t *repp
  );
  /* Computes a "central" rgb value for this box,
    and rounds it to a "representative" pixel. */

void ppm_box_spread
  ( uint16_color_table_t *ch,
    int32_t lo,
    int32_t hi,
    ppm_pixel_t rep,
    uint16_t maxval,
    int32_t *laxp,
    float *spreadp
  );
  /* Compute the box spread {*spreadp} and largest 
    dimension {*laxp}, relative to the appointed representative {rep}. */

int32_t uint16_image_RGB_red_compare(const uint16_color_table_node_t *ch1, const uint16_color_table_node_t *ch2);

int32_t uint16_image_RGB_green_compare(const uint16_color_table_node_t *ch1, const uint16_color_table_node_t *ch2);

int32_t uint16_image_RGB_blue_compare(const uint16_color_table_node_t *ch1, const uint16_color_table_node_t *ch2);

int32_t ppm_spread_compare(const ppm_box_t *b1, const ppm_box_t *b2);

ppm_box_vec_t new_ppm_box_vec_t(int32_t n);

uint16_color_table_t uint16_color_table_t_new(int32_t n);

/* Implementation */

uint16_color_tree_t *uint16_color_tree_build
  ( uint32_t NN,                       /* Number of input colors. */
    uint32_t chns,                     /* Number of channels. */
    uint16_t *smp,                     /* Lienarized color samples */
    uint16_t maxval,                   /* Max sample value. */
    uint32_t maxColors                 /* Max desired colors. */
  )
  //uint16_color_table_t *ch,
  //uint32_t colors, 
  //uint16_t maxval,
  //uint32_t *NH_P
  {
    /* Build a list {all} of tree nodes, linked by {sub[0]}, with one node for each 
      input color.  Also get the overall bounding box {min,max}: */
    uint16_color_tree_node_t *all = NULL;
    uint16_t *min = talloc(chns, uint16_t);
    uint16_t *max = talloc(chns, uint16_t);
    for (int32_t c = 0; c < chns; c++) { min[c] = maxval; max[c] = 0; }
    uint16_t *p = smp;
    for (int32_t i = 0; i < NN; i++)
      { /* Prepend a new node to {all}: */
        uint16_color_tree_node_t *nd = talloc(1, uint16_color_tree_node_t);
        nd->clr = p;
        nd->sub[0] = all;
        nd->sub[1] = NULL;
        nd->csplit = 0; /* Meaningless for now. */
        nd->lim[0] = 0; /* Meaningless for now. */
        nd->lim[1] = 0; /* Meaningless for now. */
        all = nd;
        /* Update {min,max}: */
        for (int32_t c = 0; c < chns; c++)
          { uint16_t smp = (*p);
            demand(smp < maxval, "invalid sample in input array");
            if (smp < min[c]) { min[c] = smp; }
            if (smp > max[c]) { max[c] = smp; }
            p++;
          }
      }
      
    uct_queue_entry_t *queue = talloc(maxColors, uct_queue_entry_t);

    /* The queue is an array of {uct_queue_entry_t} records
      organized as a heap structure, sorted by the spread {sprd(qe)} of the set
      {qe.nds}, where {sprd(qe)} is the difference {qe.max[qe.cmax] -
      qe.min[qe.cmax]}. The children of queue entry
      {queue[i]} are {queue[2*i+1]} and {queue[2*i+2]}, if they exist,
      and their spread does not exceed that of their parent. So the
      root at {queue[0]} has the largest spread among all To.
    
      In order to save allocations, all the {min} and {max} vectors of those queue entries
      are packed in two flat arrays {all_min} and {all_max}. */
      
    uint16_t *all_min = talloc(maxColors*chns, uint16_t);
    uint16_t *all_max = talloc(maxColors*chns, uint16_t);
    uint32_t NQ = 0; /* Number of entries in the queue. */
    
    void add_entry(uint16_color_tree_node_t *nds);
      /* Appends to the queue, at position {NQ}, a new entry {qe} consisting of the set
        of tree nodes {nds}.  Sets the pointers {qe.min} and {qe.max}
        to the next {chns} free entries of {all_min} and {all_max},
        assumed to begin at position {NQ*chns}, and computes the min and max.
        Sets {qe.split.
        Increments {NQ}. */
    
    /* Create the first entry: */
    add_entry(all);
    
    /* Main loop: split NQ until we have enough. */
    uint32_t NT = 0; /* Number of nodes selected for the tree. */
    while (NT < maxColors)
      { 
        /* Take the first entry as the next tree node: */
        uct_queue_entry_t *q0 = pop_queue();
        
        /* Sort the nodes in {q0.nds} by the coordinate of the largest spread: */
        
        
        /* Sort box colors in the direction of largest spread: */
        float ctrv;
        switch(lax)
          { 
            case 0: cmp = (qcomparefn_t*) uint16_image_RGB_red_compare;   ctrv = queue.e[0].cR; break;
            case 1: cmp = (qcomparefn_t*) uint16_image_RGB_green_compare; ctrv = queue.e[0].cG; break;
            case 2: cmp = (qcomparefn_t*) uint16_image_RGB_blue_compare;  ctrv = queue.e[0].cB; break;
            default: pnm_error("huh?"); cmp = NULL; ctrv = 0;
          }
        qsort((void*)&(ch->e[lo]), (size_t)(hi-lo+1), sizeof(uint16_color_table_node_t), cmp); 
        
        /* Compute the 
        
        
        /* Box 0 should have the largest spread. */
        if (queue.e[0].spread <= 0.0) break; /* exact */
        int32_t lo = queue.e[0].lo;
        int32_t hi = queue.e[0].hi;
        if (lo >= hi) pnm_error("bad ppm_box spread");
        int32_t lax = queue.e[0].longaxis;
        /* Split colors in two groups at centroid: */
        int32_t imd = -1;
        for (int32_t i = lo; i <= hi; ++i)
          { uint16_t v = ch->e[i].color.c[lax];
            if ((float)v >= ctrv) {imd = i; break; }
          }
        if ((imd > hi) || (imd <= lo)) pnm_error("bad centroid in pgm box");
        ppm_set_box(&(queue.e[0]), ch, lo, imd-1, maxval);
        ppm_set_box(&(queue.e[NQ]), ch, imd, hi, maxval);
        ++NQ;

        /* Sort to bring the worst NQ to the top. */
        qsort((void*)queue.e, NQ, sizeof(ppm_box_t), (qcomparefn_t*)ppm_spread_compare);
      }

    /* Ok, we've got enough NQ.  Collect their centroids for the colormap: */
    uint16_color_table_t cm = uint16_color_table_new((vec_size_t)NQ);
    for (int32_t i = 0; i < NQ; ++i) 
      { cm.e[i].color = queue.e[i].rep; }
    (*NH_P) = NQ;
    free(queue.e);
    return cm;
  }

void ppm_set_box
  ( ppm_box_t *b,
    uint16_color_table_t *ch,
    int32_t lo, 
    int32_t hi,
    uint16_t maxval
  )
  {
    b->lo = lo;
    b->hi = hi;
    ppm_box_center(ch, lo, hi, maxval, &(b->cR), &(b->cG), &(b->cB), &(b->rep));
    ppm_box_sprd(ch, lo, hi, b->rep, maxval, &(b->longaxis), &(b->spread));
  }

void ppm_box_center
  ( uint16_color_table_t *ch,
    int32_t lo, int32_t hi,
    uint16_t maxval,
    float *cRp, float *cGp, float *cBp,
    ppm_pixel_t *repp
  )
  {
    double sumv[3];
    uint16_t maxv[3];
    uint16_t minv[3];
    double sumw = 0;
    
    for (int32_t c = 0; c < 3; c++)
      { maxv[c] = 0; minv[c] = maxval; sumv[c] = 0; }
    
    for (int32_t i = lo; i <= hi; ++i)
      { double fw = (double)ch->e[i].value;
        for (int32_t c = 0; c < 3; c++)
          { uint16_t v = ch->e[i].color.c[c];
            minv[c] = vmin(minv[c], v);
            maxv[c] = vmax(maxv[c], v);
            sumv[c] += fw*(double)v; 
          }
        sumw += fw; 
      }
    if (sumw <= 0.0) pnm_error("empty uint16_image_GRAY_box");

    ppm_pixel_t rep;
    float rgb[3];
    for (int32_t c = 0; c < 3; c++)
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
  }

void ppm_box_spread
  ( uint16_color_table_t *ch,
    int32_t lo,
    int32_t hi,
    ppm_pixel_t rep,
    uint16_t maxval,
    int32_t *laxp,
    float *spreadp
  )
  {
    float fv, fw;
    float sumRR = 0.0, sumGG = 0.0, sumBB = 0.0;
    uint16_t v;
    float frepR = (float)rep.c[0];
    float frepG = (float)rep.c[1];
    float frepB = (float)rep.c[2];
    long w;
    for (int32_t i = lo; i <= hi; ++i)
      {
        w = ch->e[i].value;
        fw = (float)w;
        
        v = ch->e[i].color.c[0];
        fv = (float)v - frepR;
        sumRR += fw*fv*fv;
        
        v = ch->e[i].color.c[1];
        fv = (float)v - frepG;
        sumGG += fw*fv*fv;
        
        v = ch->e[i].color.c[2];
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

int32_t uint16_image_RGB_red_compare(const uint16_color_table_node_t *ch1, const uint16_color_table_node_t *ch2)
  {
    return (int32_t) ch1->color.c[0] - (int32_t) ch2->color.c[0];
  }

int32_t uint16_image_RGB_green_compare(const uint16_color_table_node_t *ch1, const uint16_color_table_node_t *ch2)
  {
    return (int32_t) ch1->color.c[1] - (int32_t) ch2->color.c[1];
  }

int32_t uint16_image_RGB_blue_compare(const uint16_color_table_node_t *ch1, const uint16_color_table_node_t *ch2)
  {
    return (int32_t) ch1->color.c[2] - (int32_t) ch2->color.c[2];
  }

int32_t ppm_spread_compare(const ppm_box_t *b1, const ppm_box_t *b2)
  {
    return ((int32_t)
      (b1->spread > b2->spread ? -1 : 
      (b1->spread < b2->spread ?  1 :
      0))
    );
  }
  
vec_typeimpl(ppm_box_vec_t, ppm_box_vec, ppm_box_t);

/* Copyright (C) 1989, 1991 by Jef Poskanzer.
**
** Permission to use, copy, modify, and distribute this software and its
** documentation for any purpose and without fee is hereby granted, provided
** that the above copyright notice appear in all copies and that both that
** copyright notice and this permission notice appear in supporting
** documentation.  This software is provided "as is" without express or
** implied warranty.
*/

