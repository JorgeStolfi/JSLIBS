/* See {kdtom.h}. */
/* Last edited on 2021-06-24 01:36:00 by jstolfi */

#define _GNU_SOURCE
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <assert.h>

#include <bool.h>
#include <ppv_array.h>
#include <affirm.h>

#include <kdtom.h>
#include <kdtom_split.h>
#include <kdtom_array.h>
#include <kdtom_const.h>
#include <kdtom_ixmap.h>

ppv_sample_t kdtom_get_sample(kdtom_t *T, ppv_index_t ix[])
  {
    /* Check index bounds.  Early warning is good. */
    for (ppv_axis_t ax = 0; ax < T->d; ax++)
      { demand((0 <= ix[ax]) && (ix[ax]< T->size[ax]), "invalid voxel index"); }
      
    /* Switch to proper handler: */
    switch (T->kind)
      { 
        case kdtom_kind_ARRAY:
          return kdtom_array_get_sample((kdtom_array_t *)T, ix);
        case kdtom_kind_CONST:
          return kdtom_const_get_sample((kdtom_const_t *)T, ix);
        case kdtom_kind_SPLIT:
          return kdtom_split_get_sample((kdtom_split_t *)T, ix);
        case kdtom_kind_IXMAP:
          return kdtom_ixmap_get_sample((kdtom_ixmap_t *)T, ix);
        default:
          assert(FALSE);
      }
  }

kdtom_t *kdtom_alloc(ppv_dim_t d, size_t rec_bytes, char **pendP)
  {
    size_t hdr_bytes = sizeof(kdtom_t);        /* Bytesize of the {kdom_t} fixed fields. */
    size_t szv_bytes = d * sizeof(ppv_size_t); /* Bytesize of the {size} vector. */
    size_t tot_bytes = rec_bytes + szv_bytes;  /* Total bytes to allocate. */
    
    kdtom_t *T = (kdtom_t *)notnull(malloc(tot_bytes), "no mem");
    char *pend = ((char *)T) + hdr_bytes; /* Freespace pointer. */
    
    /* Set the fields {d,size}: */
    T->d = d;
    T->size = (ppv_size_t *)pend;  pend += szv_bytes;

    /* Return results: */
    if ((pendP != NULL) && (rec_bytes > hdr_bytes)) { (*pendP) = pend; }
    return T;
  }
