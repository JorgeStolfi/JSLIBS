/* See {nmsim_elem_neuron_state.h} */
/* Last edited on 2019-05-28 14:42:15 by jstolfi */

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

#include <nmsim_elem_neuron_state.h>

void nmsim_elem_neuron_state_write(FILE *wr, nmsim_elem_neuron_state_t *st)
  {
    fputs(" ", wr); nmsim_write_double_value(wr, st->V, nmsim_write_VIJ_PREC, TRUE, FALSE, FALSE);
    fprintf(wr, " %ld", st->age);
    fputs("\n", wr);
  }
    
void nmsim_elem_neuron_state_read(FILE *rd, nmsim_elem_neuron_state_t *st)
  {
    double Vmax = 10000.0; /* A max double value for {V}. */
    fget_skip_spaces(rd);
    st->V = nmsim_read_double_value(rd, "V", -Vmax, +Vmax);
    st->age = (nmsim_step_count_t)nmsim_read_int64_value(rd, "age", -1, nmsim_step_count_MAX);
  }
