/* See {group_sort_uint32.h}. */
/* Last edited on 2024-11-22 03:28:32 by stolfi */

#include <stdint.h>
#include <stdlib.h>
#include <assert.h>

#include <affirm.h>

#include <group_sort_uint32.h>

void group_sort_uint32
  ( uint32_t ni,        /* Number of integers to sort. */
    uint32_t ng,        /* Number of groups. */
    int32_t group[],    /* Group assigned to each integer. */
    uint32_t isort[],   /* (OUT) Integers {0..ni-1} sorted by group. */
    int32_t gstart[]    /* (OUT) Beg and end of each group in {isort}. */
  )
  {
    /* Determine the number of integers {gsize[ig]} in each group {ig}: */
    uint32_t *gsize = notnull(malloc(ng*sizeof(uint32_t)), "no mem"); /* Num of integers in each group. */
    for (uint32_t ig = 0;  ig < ng; ig++) { gsize[ig] = 0; }
    for (uint32_t it = 0;  it < ni; it++)
      { int32_t ig = group[it];
        assert(ig < ng);
        gsize[ig]++; 
      }
      
    /* Compute the start {gstart[ig]} of each group in the {ix} list: */
    gstart[0] = 0;
    for (uint32_t ig = 0;  ig < ng; ig++)
      { gstart[ig+1] = gstart[ig] + (int32_t)gsize[ig]; }
    
    /* Distribute each integer {0..ni-1} into the proper group {ig}, decrementing {gsize[ig]}. */
    for (uint32_t it = 0;  it < ni; it++)
      { /* Get the group of integer {it}: */
        int32_t ig = group[it];
        /* At this point, {gsize[ig]} is how many faces are still missing in group {ig}. */
        /* Store integer {it} in group {ig}, decrement {gsize}: */
        assert((0 < gsize[ig]) && (gsize[ig] <= gstart[ig+1]));
        int32_t k = gstart[ig+1] - (int32_t)gsize[ig]; /* First free {ix} slot in group {ig}. */
        assert(k >= gstart[ig]);
        assert(k < ni);
        isort[k] = (uint32_t)it;
        gsize[ig]--;
      }
      
    /* Paranoia check (should be in tests directory): */
    for (uint32_t ig = 0;  ig < ng; ig++) 
      { assert(gsize[ig] == 0);
        for (int32_t k = gstart[ig]; k < gstart[ig+1]; k++) 
          { assert(group[isort[k]] == ig); }
      }
    assert(gstart[ng] == ni);
    
    /* Cleanup: */
    free(gsize);
  }
