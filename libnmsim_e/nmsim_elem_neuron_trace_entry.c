/* See {nmsim_elem_neuron_trace_entry.h} */
/* Last edited on 2019-05-28 14:39:36 by jstolfi */

#define _GNU_SOURCE
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <assert.h>
#include <math.h>

#include <bool.h>
#include <affirm.h>
#include <nget.h>
#include <fget.h>
#include <filefmt.h>

#include <nmsim_basic.h>
#include <nmsim_write.h>
#include <nmsim_read.h>
#include <nmsim_compare.h>

#include <nmsim_elem_neuron_trace_entry.h>

void nmsim_elem_neuron_trace_entry_clear(nmsim_elem_neuron_trace_entry_t *ts)
  {
    (*ts) = (nmsim_elem_neuron_trace_entry_t)
      { .V = NAN, .age = -1,
        .M = NAN, .H = NAN, 
        .X = FALSE, .I = NAN, .J = NAN
      };
  }

bool_t nmsim_elem_neuron_trace_entry_is_undefined(nmsim_elem_neuron_trace_entry_t *ts)
  {
    if (ts == NULL)
      { return TRUE; }
    else if (isnan(ts->V))
      { /* Every field must be undefined: */
        assert(ts->age == -1);
        assert(isnan(ts->M) && isnan(ts->H));
        assert(! ts->X);
        assert(isnan(ts->I) && isnan(ts->J));
        return TRUE;
      }
    else
      { /* Every state field must be valid, but the evolution fields need not be: */
        assert(ts->age >= 0);
        assert((!isnan(ts->M)) && (!isnan(ts->H)));
        return FALSE;
      }
  }

void nmsim_elem_neuron_trace_entry_write(FILE *wr, nmsim_elem_neuron_trace_entry_t *ts)
  {
    fputs(" ", wr);  nmsim_write_double_value(wr, ts->V, nmsim_write_VIJ_PREC, TRUE,  FALSE, FALSE);
    fprintf(wr, " %10ld", ts->age);
    fputs("  ", wr); nmsim_write_double_value(wr, ts->M, nmsim_write_MH_PREC,  FALSE, TRUE,  TRUE );
    fputs(" ", wr);  nmsim_write_double_value(wr, ts->H, nmsim_write_MH_PREC,  FALSE, TRUE,  TRUE );
    fprintf(wr, " %1d", (int32_t)(ts->X));
    fputs("  ", wr); nmsim_write_double_value(wr, ts->I, nmsim_write_VIJ_PREC, TRUE,  TRUE,  FALSE);
    fputs(" ", wr);  nmsim_write_double_value(wr, ts->J, nmsim_write_VIJ_PREC, TRUE,  TRUE,  FALSE);
  }
    
void nmsim_elem_neuron_trace_entry_read(FILE *rd, nmsim_elem_neuron_trace_entry_t *ts)
  {
    double VIJmax = 10000.0; /* A max double value for {V,I,J}. */
    double MHmax = 10000.0; /* A max double value for {M,H}. */
    fget_skip_spaces(rd);
    ts->V = nmsim_read_double_value(rd, "V", -VIJmax, +VIJmax);
    ts->age = (nmsim_step_count_t)nmsim_read_int64_value(rd, "age", -1, nmsim_step_count_MAX);
    ts->M = nmsim_read_double_value(rd, "M", 0.0, MHmax);
    ts->H = nmsim_read_double_value(rd, "H", 0.0, MHmax);
    ts->X = (nmsim_read_int64_value(rd, "X", 0, 1) == 1);
    ts->I = nmsim_read_double_value(rd, "I", -VIJmax, VIJmax);
    ts->J = nmsim_read_double_value(rd, "J", -VIJmax, VIJmax);
  }

void nmsim_elem_neuron_trace_entry_compare
  ( nmsim_elem_neuron_trace_entry_t *ts_read, 
    nmsim_elem_neuron_trace_entry_t *ts_orig
  )
  {
    /* Check for undefined (or {NULL}) entries: */
    bool_t undef_read = nmsim_elem_neuron_trace_entry_is_undefined(ts_read);
    bool_t undef_orig = nmsim_elem_neuron_trace_entry_is_undefined(ts_orig);
    if (undef_read != undef_orig)
      { fprintf(stderr, "  ** only one of the traces is defined");
        if (undef_read)
          { fprintf(stderr, " -- entry {ts_orig}\n"); }
        else
          { fprintf(stderr, " -- entry {ts_read}\n"); }
        demand(FALSE, "aborted");
      }
    if (! (undef_read | undef_orig))
      { /* Compare entries: */ 
        nmsim_compare_double_param("V",    ts_read->V,          ts_orig->V,    nmsim_write_VIJ_PREC, FALSE, FALSE);
        nmsim_compare_int64_param ("age",  ts_read->age,        ts_orig->age);
        nmsim_compare_double_param("M",    ts_read->M,          ts_orig->M,    nmsim_write_MH_PREC, TRUE,  TRUE );
        nmsim_compare_double_param("H",    ts_read->H,          ts_orig->H,    nmsim_write_MH_PREC, TRUE,  TRUE );
        nmsim_compare_int64_param ("X",    (int64_t)ts_read->X, (int64_t)ts_orig->X);
        nmsim_compare_double_param("I",    ts_read->I,          ts_orig->I,    nmsim_write_VIJ_PREC, TRUE,  FALSE);
        nmsim_compare_double_param("J",    ts_read->J,          ts_orig->J,    nmsim_write_VIJ_PREC, TRUE,  FALSE);
      }
  }
