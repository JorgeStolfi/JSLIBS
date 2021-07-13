/* See {kdtom.h}. */
/* Last edited on 2021-07-12 22:08:54 by jstolfi */

#define _GNU_SOURCE
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <assert.h>

#include <bool.h>
#include <ppv_array.h>
#include <jsmath.h>
#include <affirm.h>

#include <kdtom.h>
#include <kdtom_split.h>
#include <kdtom_array.h>
#include <kdtom_const.h>

bool_t kdtom_box_is_empty(ppv_dim_t d, ppv_size_t  size[])
  { for (ppv_axis_t k = 0; k < d; k++)
      { if (size[k] == 0) { return TRUE; } }
    return FALSE;
  }

bool_t kdtom_has_empty_core(kdtom_t  *T)
  { return (T->size[0] == 0); }

ppv_sample_t kdtom_get_sample(kdtom_t *T, ppv_index_t ix[])
  {
    ppv_dim_t d = T->d;
    
    /* Check index bounds, compute relative index {dx}: */
    bool_t inside_core = TRUE;
    ppv_index_t dx[d];
    for (ppv_axis_t k = 0; k < d; k++)
      { ppv_index_t dxk = ix[k] - T->ixlo[k];
        inside_core = inside_core && (0 <= dxk) && (dxk < T->size[k]);
        assert((T->size[k] == 0) == (T->size[0] == 0)); /* Empty domain invariant. */
        dx[k] = dxk;
      }
    
    ppv_sample_t smp;
    if (! inside_core)
      { smp = T->fill; }
    else
      { switch (T->kind)
          { case kdtom_kind_ARRAY:
              smp = kdtom_array_get_core_sample((kdtom_array_t *)T, dx);
              break;
            case kdtom_kind_CONST:
              smp = kdtom_const_get_core_sample((kdtom_const_t *)T, dx);
              break;
            case kdtom_kind_SPLIT:
              smp = kdtom_split_get_core_sample((kdtom_split_t *)T, dx);
              break;
            default:
              assert(FALSE);
          }
      }
    assert(smp <= T->maxsmp);
    return smp;
  }

size_t kdtom_bytesize(kdtom_t *T, bool_t total)
  {
    size_t bytes = 0;
    switch (T->kind)
      { case kdtom_kind_ARRAY:
          bytes = kdtom_array_bytesize((kdtom_array_t *)T, total);
          break;
        case kdtom_kind_CONST:
          bytes = kdtom_const_bytesize((kdtom_const_t *)T);
          break;
        case kdtom_kind_SPLIT:
          bytes = kdtom_split_bytesize((kdtom_split_t *)T, total);
          break;
        default:
          assert(FALSE);
      }
    assert(bytes > 0);
    return bytes;
  }

kdtom_t *kdtom_clip(kdtom_t *T, ppv_index_t ixlo[], ppv_size_t size[])
  {
    ppv_dim_t d = T->d;
    bool_t empty = FALSE; /* Is the clip box empty? */
    for (ppv_axis_t k = 0; k < d; k++) { if (size[k] == 0) { empty = TRUE; break; } }
    
    bool_t empty_core = empty; /* Is the core {S.K} empty?  */
    ppv_index_t ixlo_B[d]; /* Low corner of {T.DK} clipped. */
    ppv_size_t size_B[d];  /* Size of {T.DK} clipped. */
    if (! empty)
      { for (ppv_axis_t k = 0; k < d; k++)
          { ppv_index_t ixlok = (ixlo == NULL ? 0 : ixlo[k]);
            ixlo_B[k] = (ppv_index_t)imax(ixlok, T->ixlo[k]);
            ppv_index_t ixhik_B = (ppv_index_t)imin(ixlok+size[k], T->ixlo[k]+T->size[k]);
            if (ixlo_B[k] >= ixhik_B) { empty_core = TRUE; break; }
            size_B[k] = (ppv_size_t)(ixhik_B - ixlo_B[k]);
          }
      }
    
    kdtom_t *S = NULL;
    if (empty_core)
      { S = (kdtom_t *)kdtom_const_make(d, T->maxsmp, T->fill, NULL, NULL, T->fill); }
    else
      { switch (T->kind)
          { case kdtom_kind_ARRAY:
              S = (kdtom_t *)kdtom_array_clip((kdtom_array_t *)T, ixlo_B, size_B);
              break;
            case kdtom_kind_CONST:
              S = (kdtom_t *)kdtom_const_clip((kdtom_const_t *)T, ixlo_B, size_B);
              break;
            case kdtom_kind_SPLIT:
              S = (kdtom_t *)kdtom_split_clip((kdtom_split_t *)T, ixlo_B, size_B);
              break;
            default:
              assert(FALSE);
          }
      }
    assert(S != NULL);
    return S;
  }

char *kdtom_alloc_internal_vector(kdtom_t *T, size_t sztot, ppv_dim_t d, size_t elsz, char **pendP)
  {
    demand((d > 0) && (d <= ppv_MAX_DIM), "invalid dimension {d}");
    char *pend = (*pendP);
    
    /* Check where the address {pend} is properly aligned for the element size: */
    demand(((uint64_t)pend) % elsz == 0, "{*pendP} sync error");
    
    /* Check for enough space, alocate the vector, and update {*pendP}: */
    size_t vec_bytes = d*elsz;
    assert (sztot - (pend - (char*)T) >= vec_bytes);
    char *vaddr = pend; pend += vec_bytes;
    (*pendP) = pend;
    
    return vaddr;
  }

void kdtom_node_init
  ( kdtom_t *T, 
    kdtom_kind_t kind,
    ppv_dim_t d, 
    ppv_sample_t maxsmp,
    ppv_sample_t fill,
    ppv_index_t ixlo[], 
    ppv_size_t size[], 
    size_t sztot,
    char **pendP
  )
  {
    demand((d > 0) && (d <= ppv_MAX_DIM), "invalid dimension {d}");
    demand(maxsmp <= ppv_MAX_SAMPLE_VAL, "invalid {maxsmp}");
    demand((kind >= kdtom_kind_FIRST) && (kind <= kdtom_kind_LAST), "invalid {kind}");
    demand(fill <= maxsmp, "invalid {fill}");
    bool_t empty = (size == NULL);
    if (! empty)
      { for(ppv_axis_t k = 0; k < d; k++) 
          { if (size[k] == 0) { empty = TRUE; break; } }
      }
    if (!empty)
      { for (ppv_axis_t k = 0; k < d; k++) 
          { demand(size[k] <= ppv_MAX_SIZE, "invalid {size}");
            if (ixlo != NULL)
              { ppv_index_t ixmin = - ppv_MAX_INDEX;
                ppv_index_t ixmax = + ppv_MAX_INDEX - size[k];
                demand((ixlo[k] >= ixmin) && (ixlo[k] <= ixmax), "invalid {ixlo}");
              }
          }
      }
    
    /* Check for clobbering of the fixed head fields: */
    size_t fix_bytes = sizeof(kdtom_t);        /* Fixed fields, incl. {ixlo,size} pointers. */
    demand((*pendP) - (char *)T >= fix_bytes, "bad {*pendP}");
    
    /* Allocate the {T.ixlo} vector: */
    size_t elsz = sizeof(ppv_index_t);
    T->ixlo = (ppv_index_t *)kdtom_alloc_internal_vector(T, sztot, d, elsz, pendP);
    
    /* Allocate the {T.size} vector: */
    elsz = sizeof(ppv_size_t);
    T->size = (ppv_size_t *)kdtom_alloc_internal_vector(T, sztot, d, elsz, pendP);
    
    /* Stores the given fields: */
    T->kind = kind;
    T->d = d;
    T->maxsmp = maxsmp;
    T->fill = fill;
    for (ppv_axis_t k = 0; k < d; k++) 
      { if (empty)
          { T->ixlo[k] = 0; T->size[k] = 0; }
        else
          { T->ixlo[k] = (ixlo == NULL ? 0 : ixlo[k]); T->size[k] = size[k]; }
      }
  }

void kdtom_translate(kdtom_t  *T, ppv_index_t dx[])
  {
    bool_t empty = kdtom_has_empty_core(T);
    for (ppv_axis_t k = 0; k < T->d; k++) 
      { if (empty)
          { assert(T->ixlo[k] == 0);}
        else
          { T->ixlo[k] += dx[k]; }
      }
  }

bool_t kdtom_is_all_fill(kdtom_t *T, ppv_sample_t fill)
  {
    if (T->kind != kdtom_kind_CONST) { return FALSE; }
    kdtom_const_t *Tc = (kdtom_const_t *)T;
    return kdtom_const_is_all_fill(Tc, fill);
  }

kdtom_t *kdtom_join_nodes
  ( ppv_size_t size[], 
    ppv_axis_t ax,
    kdtom_t *T0, 
    ppv_size_t sz0, 
    kdtom_t *T1, 
    ppv_size_t sz1
  )
  {
    if (T0->kind != kdtom_kind_CONST) { return NULL; }
    if (T1->kind != kdtom_kind_CONST) { return NULL; }
    kdtom_const_t *Tc0 = (kdtom_const_t *)T0;
    kdtom_const_t *Tc1 = (kdtom_const_t *)T1;
    kdtom_const_t *Tc = kdtom_const_join_nodes(size, ax, Tc0, sz0,Tc1, sz1);
    return (kdtom_t *)Tc;
  }

void kdtom_intersect_boxes
  ( ppv_dim_t d,
    ppv_index_t ixlo_A[], 
    ppv_size_t  size_A[],
    ppv_index_t ixlo_B[], 
    ppv_size_t  size_B[],
    ppv_index_t ixlo_R[], 
    ppv_size_t  size_R[]
  )
  {
    bool_t empty_A = kdtom_box_is_empty(d, size_A);
    bool_t empty_B = kdtom_box_is_empty(d, size_B);
    
    /* Compute the new core domain {ixlo_R,size_R} and {empty_R}: */
    bool_t empty_R = FALSE;
    for (ppv_axis_t k = 0; k < d;  k++)
      { /* Get the index range  {ixlok_A..ixhik_A} of box{A} on axis {k}: */
        ppv_index_t ixlok_A = (ppv_size_t)(empty_A ? 0 : ixlo_A[k]);
        ppv_index_t ixhik_A = (ppv_size_t)(empty_A ? -1 : ixlok_A + size_A[k] - 1);
                
        /* Get the index range {ixlok_B..ixhik_B} of  box {B} on axis {k}: */
        ppv_index_t ixlok_B = (ppv_size_t)(empty_B ? 0 : ixlo_B[k]);
        ppv_index_t ixhik_B = (ppv_size_t)(empty_B ? -1 : ixlo_B[k] + size_B[k] - 1);
        
        /* Get the expected clipped core box {ixlo_R[k],size_R[k]}: */
        ppv_index_t ixlok_R = (ppv_index_t)imax(ixlok_A, ixlok_B);
        ppv_index_t ixhik_R = (ppv_index_t)imin(ixhik_A, ixhik_B);
        if (ixlok_R > ixhik_R) { empty_R = TRUE; break; }
        ixlo_R[k] = ixlok_R; 
        size_R[k] = (ppv_size_t)(ixhik_R - ixlok_R + 1);
      }
      
    if (empty_R) { for (ppv_axis_t k = 0; k < d;  k++) { ixlo_R[k] = 0; size_R[k] = 0; } }
    return;
  }
   
