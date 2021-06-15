/* See ppv_brush.h */
/* Last edited on 2021-06-13 13:36:22 by jstolfi */

#define _GNU_SOURCE
#include <stdio.h>
#include <stdint.h>
#include <stdlib.h>
#include <math.h>
#include <assert.h>

#include <bool.h>
#include <jsmath.h>
#include <indexing.h>
#include <indexing_io.h>
#include <affirm.h>
#include <ppv_array.h>

#include <ppv_brush.h>

typedef struct ppv_brush_t
  { ppv_dim_t d;            /* Dimension of brush. */
    ppv_sample_count_t nv;  /* Total count of voxels in brush. */
    ppv_index_t *ixlo;     /* Min voxel index along each axis. */
    ppv_index_t *ixhi;     /* Max voxel index along each axis. */
    ppv_index_t *index;     /* Indices of those voxels, linearized. */
  } ppv_brush_t;
  /* A record {b} of this type describes a {d}-dimensional digital brush
    as a set of {V(b)} of voxel indices, where {d = b.d}.  
  
    The brush voxel indices are relative to central voxel of the brush.
    The minimum value of {ix[ax]}, for each axis {ax} in {0..d-1}, is
    {b.ixlo[ax]}, typically negative. The maximum value is
    {b.ixhi[ax]}, typically positive.
    
    The value {b.nv} is the total count of voxels in the brush. Element
    {b.index} is a vector that contains the indices of those voxels,
    linearized. More precisely, {V(b)} is all index tuples {ix} such
    that {ix[ax]} is {b.index[k*d+ax]}, for each {ax} in {1..d-1} and
    for some {k} in {0..b.nv-1}.
    
    The voxel with indices {ix = {0,0,0...}} is the brush
    center, which corresponds to the center of the neighborhood in local
    operations. */ 

void ppv_brush_free(ppv_brush_t *b)
  { free(b->ixlo);
    free(b->ixhi);
    free(b->index);
    free(b);
  }

bool_t ppv_brush_enum(ppv_index_op_t op, ppv_brush_t *b)
  {
    ppv_index_t ix[b->d];
    for(ppv_sample_count_t k = 0;  k < b->nv; k++)
      { ppv_sample_count_t ka0 = k*b->d;
        for (ppv_axis_t ax = 0; ax < b->d; ax++) { ix[ax] = b->index[ka0 + ax]; }
        bool_t stop = op(ix);
        if (stop) { return TRUE; }
      }
    return FALSE;
  }

ppv_brush_t *ppv_brush_make_ball(ppv_dim_t d, double rad)
  { 
    bool_t debug = TRUE;
    
    ppv_brush_t *b = notnull(malloc(sizeof(ppv_brush_t)), "no mem");
    b->d = d;
    b->ixlo = notnull(malloc(d*sizeof(ppv_index_t)), "no mem");
    b->ixhi = notnull(malloc(d*sizeof(ppv_index_t)), "no mem");
    
    demand(rad >= 0, "brush radius is negative");
    ppv_size_t irad  = (ppv_size_t)floor(fabs(rad)); /* Radius of digital ball. */
    if (debug) { fprintf(stderr, "ball brush voxel indices range in {%+ld..%+ld}\n", -irad, +irad); }
    ppv_size_t sz = 2*irad + 1; /* Number of slices in brush along all axes. */
    
    ppv_sample_count_t nvmax = 1;
    ppv_size_t size[d]; /* Number of potential axis indices in each axis. */
    for (ppv_axis_t ax =0; ax < d; ax++)
       { b->ixlo[ax] = -irad;
         b->ixhi[ax] = +irad;
         size[ax] = sz;
         demand(nvmax <= ppv_MAX_SAMPLES/sz, "brush radius too big");
         nvmax = nvmax*sz;
       }
    b->index = notnull(malloc(nvmax*d*sizeof(ppv_index_t)), "no mem");
    
    ppv_sample_count_t nvtot = 0; /* Total number of voxels in ball. */
    
    auto bool_t add_voxel(const ppv_index_t ix[], ppv_pos_t pA, ppv_pos_t pB, ppv_pos_t pC);
      /* If the center of the voxel with index tuple {ix+b.ixlo} is within 
        distance {rad} of the voxel {{0,0,.. 0}}, adds it to {b}. 
        The positions {pA,pB,pC} are ignored.  Always returns {FALSE}. */
    
    (void)ix_enum(add_voxel, d, size, ix_order_L, FALSE, 0,NULL, 0,NULL, 0,NULL);
    if (debug) { fprintf(stderr, "ball brush has %lu voxels\n", nvtot); }
    b->index = realloc(b->index, nvtot*d*sizeof(ppv_index_t));
    b->nv = nvtot;
    return b;

    /* Internal procs: */

    bool_t add_voxel(const ppv_index_t ix[], ppv_pos_t pA, ppv_pos_t pB, ppv_pos_t pC)
      { /* The following computation should be exact for reasonable {rad}: */
        double dist2 = 0; /* Dist squared from center of voxel {ix+b.ixlo} to center of origin voxxel. */
        for (ppv_axis_t ax = 0; ax < d; ax++)
          { ppv_index_t ixsk = ix[ax] + b->ixlo[ax];
            dist2 += ((double)(ixsk*ixsk));
          }
        if (dist2 <= rad*rad)
          { /* Voxel is in ball: */
            ppv_sample_count_t k0 = nvtot*d;
            for (ppv_axis_t ax = 0; ax < d; ax++)
              { b->index[k0 + ax] = ix[ax] + b->ixlo[ax]; }
            nvtot++;
          }
        return FALSE;
      }
  }

ppv_dim_t ppv_brush_dimension(ppv_brush_t *b)
  { return b->d; }

ppv_sample_count_t ppv_brush_voxel_count(ppv_brush_t *b)
  { return b->nv; }
  
void ppv_brush_index_ranges(ppv_brush_t *b, ppv_index_t ixlo[], ppv_index_t ixhi[])
  { for (ppv_axis_t ax = 0; ax < b->d; ax++) 
      { ixlo[ax] = b->ixlo[ax]; ixhi[ax] = b->ixhi[ax]; }
  }
