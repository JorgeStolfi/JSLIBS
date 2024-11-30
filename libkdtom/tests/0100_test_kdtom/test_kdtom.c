#define PROG_NAME "test_kdtom"
#define PROG_DESC "Test the basics of {kdtom.h}"
#define PROG_VERS "1.0"

#define tkdt_C_COPYRIGHT \
  "Copyright © 2021 by the State University of Campinas (UNICAMP)"

/* Last edited on 2021-07-16 06:29:08 by jstolfi */

#define _GNU_SOURCE
#include <stdio.h>
#include <assert.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include <math.h>

#include <bool.h>
#include <affirm.h>
#include <jsfile.h>
#include <jsmath.h>
#include <jsrandom.h>
#include <ppv_array.h>

#include <kdtom.h>
#include <kdtom_test.h>

/* INTERNAL PROTOTYPES */

#define kdtom_kind_DUMMY kdtom_kind_CONST
  /* Just for testing. */
  
#define smpFMT ppv_sample_t_FMT
  
kdtom_t *tkdt_make_dummy_node(ppv_dim_t d, ppv_sample_t maxsmp, bool_t empty);
size_t tkdt_dummy_node_bytesize(ppv_dim_t d);

void tkdt_do_tests(ppv_dim_t d, ppv_sample_t maxsmp, bool_t empty);

int32_t main(int32_t argc,char** argv);

typedef struct tkdt_dummy_t 
  { kdtom_t h;
    ppv_sample_t *bogus;
    ppv_nbits_t silly;
    ppv_word_t *elmo;
  } tkdt_dummy_t;
  /* A dummy {kdtom_t} node variant.  The {bogus} is a vector allocated
    in the record itself. The {elmo} pointer is supposed to point somewhere
    outside or be {NULL}. */

/* IMPLEMENTATIONS */

int32_t main(int32_t argc, char** argv)
  {
    int32_t nms = 33;
    ppv_sample_t *ms = kdtom_test_pick_max_samples(nms);
    
    for (ppv_dim_t d = 1; d <= 6; d++)
      { for (uint32_t k = 0;  k < nms; k++)
          { for (uint32_t ety = 0;  ety <= 1; ety++)
              { tkdt_do_tests(d, ms[k], (bool_t)ety); }
          }
      }
    
    fprintf(stderr, "done.\n");
      
    return 0;
  }
  
void tkdt_do_tests(ppv_dim_t d, ppv_sample_t maxsmp, bool_t empty)
  { 
    fprintf(stderr, "=== testing d = %u maxsmp = " smpFMT " empty = %c ===\n", d, maxsmp, "FT"[empty]);
    
    /* Create a dummy variant of {kdtom_t}: */
    kdtom_t *T = tkdt_make_dummy_node(d, maxsmp, empty);
    
    ppv_index_t ixlo[d];
    for (ppv_axis_t k = 0; k < d;  k++) { ixlo[k] = 18 + 10*k; }
    kdtom_translate(T, ixlo);
    for (ppv_axis_t k = 0; k < d;  k++) 
      { demand((T->ixlo[k] == (empty ? 0 : ixlo[k])), "{kdtom_translate} error"); }

    ppv_index_t dx[d];
    for (ppv_axis_t k = 0; k < d;  k++) { dx[k] = 400 -10*k; }
    kdtom_translate(T, dx);
    for (ppv_axis_t k = 0; k < d;  k++) 
      { demand((T->ixlo[k] == (empty ? 0 : ixlo[k] + dx[k])), "{kdtom_translate} error"); }

    return;
  }
    
kdtom_t *tkdt_make_dummy_node(ppv_dim_t d, ppv_sample_t maxsmp, bool_t empty)
  { 
    ppv_sample_t fill = (ppv_sample_t)(maxsmp <= 1 ? maxsmp : uint64_abrandom(1, maxsmp-1));
    assert(fill <= maxsmp);
    
    /* Allocate the full record including all internal vectors: */
    size_t tot_bytes = tkdt_dummy_node_bytesize(d); /* Total node bytesize including all vectors. */
    tkdt_dummy_t *T = (tkdt_dummy_t *)notnull(malloc(tot_bytes), "no mem");
    
    /* Set {pend} to the free space inside the record, past all fixed fields: */
    size_t fix_bytes = sizeof(tkdt_dummy_t);  /* Size of fixed felds incl. those of {h}. */
    char *pend = addrsync(((char*)T) + fix_bytes, 8);
    
    /* Initialize the {h} fields, including the {T.h.size} vector: */
    ppv_size_t size[d];
    if (empty)
      { for (ppv_axis_t k = 0; k < d; k++) { size[k] = 0; } }
    else
      { ppv_sample_count_t npos = 10000;
        ppv_choose_test_size(d, npos, size);
      }
    kdtom_node_init((kdtom_t *)T, kdtom_kind_DUMMY, d, maxsmp, fill, NULL, size, tot_bytes, &pend);
    
    /* Allocate and fill the {T.bogus} vector: */
    size_t elsz = sizeof(ppv_sample_t);
    pend = addrsync(pend, 8);
    T->bogus = (ppv_sample_t *)kdtom_alloc_internal_vector((kdtom_t *)T, tot_bytes, d, elsz, &pend);
    for (ppv_axis_t k = 0; k < d; k++) { T->bogus[k] = k; }
    
    /* Initialize all remaining fields: */
    T->silly = 23;
    T->elmo = NULL;

    demand(T->h.kind == kdtom_kind_DUMMY, "{T.kind} not set correctly");
    demand(T->h.d == d, "{T.d} not set correctly");
    demand(T->h.maxsmp == maxsmp, "{T.maxsmp} not set correctly");
    demand(T->h.fill == fill, "{T.fill} not set correctly");
    for (ppv_axis_t k = 0; k < d;  k++) 
      { demand(T->h.ixlo[k] == 0, "{T.ixlo} not set correctly");
        demand(T->h.size[k] == size[k], "{T.size} not set correctly");
      }
    return (kdtom_t *)T;
  }
   
size_t tkdt_dummy_node_bytesize(ppv_dim_t d)
  {
    size_t fixf_bytes = sizeof(tkdt_dummy_t);     /* Fixed fields incl those of head part {h}. */
    size_t tot_bytes = iroundup(fixf_bytes, 8);    /* Account for address sync. */

    size_t sizv_bytes = d * sizeof(ppv_size_t);    /* Bytesize for {h.size} vector. */
    tot_bytes += iroundup(sizv_bytes, 8);          /* Paranoia, account for address sync. */

    size_t step_bytes = d * sizeof(ppv_step_t);    /* Bytesize of the {step} vector. */
    tot_bytes += iroundup(step_bytes, 8);          /* Paranoia, account for address sync. */

    size_t bogus_bytes = d * sizeof(ppv_sample_t); /* Bytesize of the {bogus} vector. */
    tot_bytes += iroundup(bogus_bytes, 8);         /* Paranoia, account for address sync. */

    return tot_bytes;
  }

