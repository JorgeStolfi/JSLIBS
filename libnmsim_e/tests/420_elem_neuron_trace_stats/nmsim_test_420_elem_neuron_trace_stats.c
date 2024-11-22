#define PROG_NAME "nmsim_test_420_elem_neuron_trace_stats"
#define PROG_DESC "tests of {limnmism} neuron-level network simulation"
#define PROG_VERS "1.0"

/* Last edited on 2022-10-20 06:32:55 by stolfi */ 

#define PROG_COPYRIGHT \
  "Copyright © 2020  State University of Campinas (UNICAMP)"
  
#define PROG_AUTH \
  "Created by J. Stolfi on 2020-12-15"
  
#define PROG_HIST
  
#define PROG_HELP \
  "  " PROG_NAME ""

#define _GNU_SOURCE
#include <stdio.h>
#include <stdint.h>
#include <string.h>
#include <stdlib.h>
#include <assert.h>
#include <math.h>

#include <bool.h>
#include <argparser.h>
#include <jsfile.h>
#include <jsrandom.h>
#include <jsmath.h>
#include <affirm.h>

#include <nmsim_basic.h>
#include <nmsim_firing_func.h>
#include <nmsim_class_neuron.h>
#include <nmsim_class_synapse.h>
#include <nmsim_group_neuron.h>
#include <nmsim_group_synapse.h>
#include <nmsim_elem_neuron.h>
#include <nmsim_elem_synapse.h>

#include <nmsim_class_net.h>
#include <nmsim_class_net_throw.h>

#include <nmsim_group_net.h>
#include <nmsim_group_net_throw.h>

#include <nmsim_elem_net.h>
#include <nmsim_elem_net_throw.h>
#include <nmsim_elem_net_sim_stats.h>

#include <nmsim_elem_neuron_trace.h>
#include <nmsim_elem_neuron_trace_entry.h>
#include <nmsim_elem_neuron_trace_stats.h>

void nmsim_test_elem_neuron_trace_stats(int32_t nne);
  /* Tests the simulation statistics tools with a hypothetical neuron of 
    {nne} neurons.  The network is actually not even created.
    
    Output files will be called "out/sim_{TAG}_stats.txt" where {TAG} is a string formed 
    from the parameters. */

void nmsim_test_elem_neuron_trace_stats_write(char *prefix, nmsim_elem_net_sim_stats_t *S);
  /* Writes the statistics in {S} with name "{prefix}.txt". */

int32_t main(int32_t argc, char **argv);

/* IMPLEMENTATIONS: */

int32_t main(int32_t argc, char **argv)
  { 
    nmsim_test_elem_neuron_trace_stats(21);
    return 0;
  }
    
void nmsim_test_elem_neuron_trace_stats(int32_t nne)
  {
    /* We need a sufficient number of neurons: */
    demand(nne >= 10, "bad neuron count");

    /* Choose nominal simulation time parameters: */
    nmsim_time_t nSteps = 10000;  /* Pretend to simulate from {t=0} to {t=nSteps}. */
    
    /* Allocate the trace structure, choose neuron to trace: */
    nmsim_elem_neuron_ix_t ineLo = int32_abrandom(nne/4, nne/2); /* Frst neuron to trace. */
    nmsim_elem_neuron_ix_t ineHi = int32_abrandom(ineLo, 4*nne/5); /* Last neuron to trace. */
    assert((ineLo > 0) && (ineLo <= ineHi) && (ineHi < nne-1));  
    nmsim_time_t trtLo = nSteps/4;    /* First time to trace. */
    nmsim_time_t trtHi = 3*nSteps/4;  /* Last time to trace. */
    fprintf(stderr, "neurons %d .. %d  times %ld .. %ld\n", ineLo, ineHi, trtLo, trtHi);
    
    nmsim_elem_neuron_trace_t *trne = nmsim_elem_neuron_trace_new(ineLo, ineHi, trtLo, trtHi);

    /* Pretend to simulate: */
    srandom(19501129);
    for (nmsim_time_t t = 0; t <= nSteps; t++)
      { /* Generate random values for state variables: */
        double V =    10.0 + 5.0*dgaussrand();
        nmsim_step_count_t age =  (nmsim_step_count_t)(t % 100);
        double M =    1.0 +0.5*dgaussrand();
        double H =    2.0 +3.0*dgaussrand();
        /* Gather statistics of state variables: */
        nmsim_elem_neuron_trace_set_V_age_M_H(trne, t, V, (double)age, M, H);
        if (t < nSteps) { 
          /* Generate random values for evolution variables: */
          bool_t X =    (drandom() < 0.10);
          double I =    10.0 + 5.0*dgaussrand();
          double J =    20.0 + 10.0*dgaussrand();
          /* Gather statistics of evolution variables: */
          nmsim_elem_neuron_trace_set_X_I_J(trne, t, (double)X, I, J); 
        }
      }

    /* Compute the statistics of that trace: */
    nmsim_elem_net_sim_stats_t *S = nmsim_elem_net_sim_stats_new(ineLo, ineHi, trtLo, trtHi);
    nmsim_elem_neuron_trace_stats_compute(trne, S);

    /* Write the statistics: */
    char *prefix = jsprintf("out/sim_%06d", nne);
    nmsim_test_elem_neuron_trace_stats_write(prefix, S);
    
    free(trne);
    free(S);
    free(prefix);
  }
  
void nmsim_test_elem_neuron_trace_stats_write(char *prefix, nmsim_elem_net_sim_stats_t *S)
  { char *fname = NULL;
    char *fname = jsprintf("%s_ne%010d--%010d_stats.txt", prefix, S->ineLo, S->ineHi);
    FILE *wr = open_write(fname, TRUE);
    nmsim_elem_net_sim_stats_write(wr, S);
    fclose(wr);
    free(fname);
  }
         
