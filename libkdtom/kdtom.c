/* See {kdtom.h}. */
/* Last edited on 2021-07-19 05:06:26 by jstolfi */

#define _GNU_SOURCE
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <assert.h>

#include <bool.h>
#include <ppv_array.h>
#include <jsmath.h>
#include <affirm.h>
#include <ixbox.h>

#include <kdtom.h>
#include <kdtom_split.h>
#include <kdtom_array.h>
#include <kdtom_const.h>
#include <kdtom_test.h>

#define ixFMT ppv_index_t_FMT
#define smpFMT ppv_sample_t_FMT
#define szFMT ppv_size_t_FMT

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

kdtom_t *kdtom_clone(kdtom_t *T)
  { 
    kdtom_t *S = NULL;
    switch (T->kind)
      { case kdtom_kind_ARRAY:
          S = (kdtom_t *)kdtom_array_clone((kdtom_array_t *)T);
          break;
        case kdtom_kind_CONST:
          S = (kdtom_t *)kdtom_const_clone((kdtom_const_t *)T);
          break;
        case kdtom_kind_SPLIT:
          S = (kdtom_t *)kdtom_split_clone((kdtom_split_t *)T);
          break;
        default:
          assert(FALSE);
      }
    assert(S != NULL);
    return S;
  }

kdtom_t *kdtom_clip_core(kdtom_t *T, ppv_index_t ixlo[], ppv_size_t size[])
  {
    bool_t debug = TRUE;
    #define dbtag "{kdtom_clip_core}: "
    
    ppv_dim_t d = T->d;
    
    ppv_index_t ixlo_new[d]; /* Low corner of {T.DK} clipped. */
    ppv_size_t size_new[d];  /* Size of {T.DK} clipped. */
    if (debug) { ixbox_print(stderr, d, dbtag "T.DK    =  [ ", T->ixlo, T->size, " ]\n"); }
    if (debug) { ixbox_print(stderr, d, dbtag "clip box = [ ", ixlo, size, " ]\n"); }
    ixbox_intersect(d, T->ixlo, T->size, ixlo, size, ixlo_new, size_new);
    if (debug) { ixbox_print(stderr, d, dbtag "clipped =  [ ", ixlo_new, size_new, " ]\n"); }
    
    kdtom_t *S = NULL;
    if (ixbox_equal(d, T->ixlo, T->size, ixlo_new, size_new))
      { /* Core did not change: */
        S = T;
      }
    else if (size_new[0] == 0)
      { /* Result has an empty core: */
        S = (kdtom_t *)kdtom_const_make(d, T->maxsmp, T->fill, NULL, NULL, T->fill);
        if (debug) { fprintf(stderr, dbtag "empty core\n"); }
      }
    else 
      { /* Result has a non-empty core: */
        switch (T->kind)
          { case kdtom_kind_ARRAY:
              S = (kdtom_t *)kdtom_array_clip_core((kdtom_array_t *)T, ixlo_new, size_new);
              break;
            case kdtom_kind_CONST:
              S = (kdtom_t *)kdtom_const_clip_core((kdtom_const_t *)T, ixlo_new, size_new);
              break;
            case kdtom_kind_SPLIT:
              S = (kdtom_t *)kdtom_split_clip_core((kdtom_split_t *)T, ixlo_new, size_new);
              break;
            default:
              assert(FALSE);
          }
        if (debug) { ixbox_print(stderr, d, dbtag "S.DK =  [ ", S->ixlo, S->size, " ]\n"); }
      }
    assert(S != NULL);
    return S;
    #undef dbtag
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
    bool_t empty = ixbox_is_empty(d, size);
    if (! empty)
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

void kdtom_translate_one(kdtom_t  *T, ppv_axis_t ax, ppv_index_t dx)
  {
    bool_t empty = kdtom_has_empty_core(T);
    if (empty)
      { assert(T->ixlo[ax] == 0);}
    else
      { T->ixlo[ax] += dx; }
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
    ppv_size_t size0, 
    kdtom_t *T1, 
    ppv_size_t size1
  )
  {
    if (T0->kind != kdtom_kind_CONST) { return NULL; }
    if (T1->kind != kdtom_kind_CONST) { return NULL; }
    kdtom_const_t *Tc0 = (kdtom_const_t *)T0;
    kdtom_const_t *Tc1 = (kdtom_const_t *)T1;
    kdtom_const_t *Tc = kdtom_const_join_nodes(size, ax, Tc0, size0, Tc1, size1);
    return (kdtom_t *)Tc;
  }

void kdtom_print_node(FILE *wr, int32_t ind, char *name, kdtom_t *T, bool_t rec)
  { 
    fprintf(wr, "%*s", ind, "");
    if (name != NULL) { fprintf(wr, "%s = ", name); }
    if (T == NULL)
      { fputs("NULL", wr); }
    else
      { fputs("{", wr);
        fprintf(wr, " .d = %d", T->d);
        fprintf(wr, " .maxsmp = " smpFMT, T->maxsmp);
        fprintf(wr, " .fill = " smpFMT, T->fill);
        ixbox_print(wr, T->d, " .KD = [", T->ixlo, T->size, " ]");
        switch (T->kind)
          { case kdtom_kind_ARRAY:
              kdtom_array_print_fields(wr, (kdtom_array_t *)T);
              fputs(" }\n", wr);
              break;
            case kdtom_kind_CONST:
              kdtom_const_print_fields(wr, (kdtom_const_t *)T);
              fputs(" }\n", wr);
              break;
            case kdtom_kind_SPLIT:
              kdtom_split_print_fields(wr, (kdtom_split_t *)T);
              fputs(" }\n", wr);
              if (rec)
                { kdtom_split_print_subtrees(wr, ind+2, (kdtom_split_t *)T); }
              break;
            default:
              assert(FALSE);
          }
      }

    return;
  }

