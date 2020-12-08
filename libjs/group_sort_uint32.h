/* Group sorting of {uint32_t} integers */
/* Last edited on 2016-04-17 10:00:59 by stolfilocal */

#ifndef group_sort_uint32_H
#define group_sort_uint32_H

#define _GNU_SOURCE
#include <stdint.h>

void group_sort_uint32
  ( uint32_t ni,        /* Number of integers to sort. */
    uint32_t ng,        /* Number of groups. */
    uint32_t group[],   /* Group assigned to each integer. */
    uint32_t isort[],   /* (OUT) Integers {0..ni-1} sorted by group. */
    uint32_t gstart[]   /* (OUT) Beg and end of each group in {isort}. */
  );
  /* Sorts the integers {0..ni-1} into {ng} groups assigned to 
    them by the client. 
    
    Assumes that each group is identified by a /group index/ in
    {0..ng-1}; specifically, that integer {i} belongs to the group with
    index {group[i]}.
    
    The procedure stores the integers {0..ni-1} into {isort[0..ni-1]},
    sorted by increasing group index. It also sets {gstart[0..ng]} to
    indicate the beginning and end of each group in {f}. Namely, for any
    group index {ig} in {0..ng-1}, the integers in group {ig} are
    {isort[k0..k1-1]}, where {k0=gstart[ig]} and {k1=gstart[ig+1]}. Note
    that {gstart} must have {ng+1} elements, not {ng}. */

#endif
