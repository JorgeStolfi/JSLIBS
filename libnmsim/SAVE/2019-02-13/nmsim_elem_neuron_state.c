/* See {nmsim_elem_neuron_state.h} */
/* Last edited on 2019-02-13 18:28:32 by jstolfi */

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

#include <nmsim_elem_neuron_state.h>

void nmsim_elem_neuron_state_write(FILE *wr, nmsim_time_t t, nmsim_elem_neuron_state_t *st)
  {
    fprintf(wr, "%10ld  %+10.6f %10ld", t, st->V, st->age);
    fprintf(wr, "  %10.7f %10.7f", st->M, st->H);
    fprintf(wr, " %1d", (int32_t)(st->X));
    fprintf(wr, "  %+10.6f %+10.6f", st->I, st->J);
    fputs("\n", wr);
  }
    
void nmsim_elem_neuron_state_read
  ( FILE *rd, 
    nmsim_time_t t,
    nmsim_elem_neuron_state_t *st
  )
  {
    double dmax = 10000.0; /* A max double value for {V,M,H,I,J}. */
    fget_skip_formatting_chars(rd);
    (void)read_int64_value(rd, "time", t, t);
    st->V = read_double_value(rd, "V", -dmax, +dmax);
    st->age = (nmsim_step_count_t)read_int64_value(rd, "age", 0, nmsim_step_count_MAX);
    st->M = read_double_value(rd, "M", 0.0, dmax);
    st->H = read_double_value(rd, "H", 0.0, dmax);
    st->X = (read_int64_value(rd, "X", 0, 1) == 1);
    st->I = read_double_value(rd, "I", -dmax, dmax);
    st->J = read_double_value(rd, "J", -dmax, dmax);
  }
