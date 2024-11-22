/* See {group_sort_uint32.h}. */
/* Last edited on 2024-11-22 02:19:56 by stolfi */

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
    uint32_t ig;
    for (ig = 0; ig < ng; ig++) { gsize[ig] = 0; }
    uint32_t it;
    for (it = 0; it < ni; it++)
      { ig = group[it];
        assert(ig < ng);
        gsize[ig]++; 
      }
      
    /* Compute the start {gstart[ig]} of each group in the {ix} list: */
    gstart[0] = 0;
    for (ig = 0; ig < ng; ig++)
      { gstart[ig+1] = gstart[ig] + gsize[ig]; }
    
    /* Distribute each integer {0..ni-1} into the proper group {ig}, decrementing {gsize[ig]}. */
    for (it = 0; it < ni; it++)
      { /* Get the group of integer {it}: */
        ig = group[it];
        /* At this point, {gsize[ig]} is how many faces are still missing in group {ig}. */
        /* Store integer {it} in group {ig}, decrement {gsize}: */
        assert((0 < gsize[ig]) && (gsize[ig] <= gstart[ig+1]));
        uint32_t k = gstart[ig+1] - gsize[ig]; /* First free {ix} slot in group {ig}. */
        assert(k >= gstart[ig]);
        assert(k < ni);
        isort[k] = it;
        gsize[ig]--;
      }
      
    /* Paranoia check (should be in tests directory): */
    for (ig = 0; ig < ng; ig++) 
      { assert(gsize[ig] == 0);
        for (int32_t k = gstart[ig]; k < gstart[ig+1]; k++) 
          { assert(group[isort[k]] == ig); }
      }
    assert(gstart[ng] == ni);
    
    /* Cleanup: */
    free(gsize);
  }
