/* See {kdtom.h}. */
/* Last edited on 2021-07-11 18:09:13 by jstolfi */

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
            case kdtom_kind_CONST:
              smp = kdtom_const_get_core_sample((kdtom_const_t *)T, dx);
            case kdtom_kind_SPLIT:
              smp = kdtom_split_get_core_sample((kdtom_split_t *)T, dx);
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
        case kdtom_kind_CONST:
          bytes = kdtom_const_bytesize((kdtom_const_t *)T);
        case kdtom_kind_SPLIT:
          bytes = kdtom_split_bytesize((kdtom_split_t *)T, total);
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
            case kdtom_kind_CONST:
              S = (kdtom_t *)kdtom_const_clip((kdtom_const_t *)T, ixlo_B, size_B);
            case kdtom_kind_SPLIT:
              S = (kdtom_t *)kdtom_split_clip((kdtom_split_t *)T, ixlo_B, size_B);
            default:
              assert(FALSE);
          }
      }
    assert(S != NULL);
    return S;
  }

char *kdtom_alloc_internal_vector(kdtom_t *T, size_t sztot, ppv_dim_t d, size_t elsz, char **pendP)
  {
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
    demand(maxsmp <= ppv_MAX_SAMPLE_VAL, "invalid {maxsmp}");
    demand(d <= ppv_MAX_DIM, "invalid {d}");
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
            demand((ixlo[k] >= -ppv_MAX_INDEX) && (ixlo[k] <= ppv_MAX_INDEX - size[k]), "invalid {ixlo}");
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
          { T->ixlo[k] = ixlo[k]; T->size[k] = size[k]; }
      }
  }

void kdtom_translate(kdtom_t  *T, ppv_index_t dx[])
  {
    for (ppv_axis_t k = 0; k < T->d; k++) { T->ixlo[k] += dx[k]; }
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
