#define PROG_NAME "nmsim_test_330_elem_neuron_trace"
#define PROG_DESC "tests of {limnmism} single neuron simulation trace"
#define PROG_VERS "1.0"

/* Last edited on 2022-10-20 06:33:30 by stolfi */ 

#define PROG_COPYRIGHT \
  "Copyright © 2019  State University of Campinas (UNICAMP)"
  
#define PROG_AUTH \
  "Created by J. Stolfi on 2019-01-11"
  
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
#include <nmsim_test.h>

#include <nmsim_class_neuron.h>
#include <nmsim_group_neuron.h>
#include <nmsim_elem_neuron.h>

#include <nmsim_elem_neuron_trace.h>

void nmsim_test_elem_neuron_trace
  ( char *tag,
    nmsim_time_t tLo_s, 
    nmsim_time_t tHi_s,
    nmsim_time_t tLo_t, 
    nmsim_time_t tHi_t
  );
  /* Tests a neuron element trace structure that can save neuron states for time
    {t} and step data for the step from {t} to {t+1}, where {t} is in
    {tLo_t .. tHi_t}; that is, with {tHi_t - tLo_t + 1} entries.  
    
    The simulation is run from {tLo_s}
    to {tHi_s} inclusive, that is {tHi_s - tLo_s} steps.
    
    Output file will be called "out/trace_{tag}.txt". */

void nmsim_test_elem_neuron_trace_write(char *fname, nmsim_elem_neuron_trace_t *etrace);
  /* Writes the neuron element trace {etrace} to a file with name {fname}. */

void nmsim_test_elem_neuron_trace_read(char *fname, nmsim_elem_neuron_trace_t *etrace);
  /* Reads a neuron element trace from file {fname} and compares
    its contents to that of the given {etrace}. */

int32_t main(int32_t argc, char **argv);

/* IMPLEMENTATIONS: */

int32_t main(int32_t argc, char **argv)
  { 
    nmsim_test_elem_neuron_trace("A", 418, 418, 400, 450);
    nmsim_test_elem_neuron_trace("B",  10,  20,   7,  17);
    nmsim_test_elem_neuron_trace("C", -10,  10,  -7,  20);
    nmsim_test_elem_neuron_trace("C", -30,  30, -17,  27);
    return 0;
  }
    
void nmsim_test_elem_neuron_trace
  ( char *tag,
    nmsim_time_t tLo_s, 
    nmsim_time_t tHi_s,
    nmsim_time_t tLo_t, 
    nmsim_time_t tHi_t
  )
  {
    fprintf(stderr, "simulation times %5ld .. %5ld ( %5ld steps   )\n", tLo_s, tHi_s, tHi_s - tLo_s);
    fprintf(stderr, "tracing times    %5ld .. %5ld ( %5ld entries )\n", tLo_t, tHi_t, tHi_t - tLo_t + 1);
    
    /* Create a neuron class: */
    nmsim_class_neuron_count_t nnc = 10; /* Arbitrary number of neuron classes in network. */
    nmsim_class_neuron_t *nclass = nmsim_class_neuron_throw();
    
    /* Pick an index for the neuron: */
    nmsim_group_neuron_count_t nng = 2*nnc; /* Arbitrary number of neuron groups in network. */ 
    nmsim_elem_neuron_count_t nne = 10*nng; /* Arbitrary num of neuron elems in network. */
      
    /* Create the filename: */
    char *fname = NULL;
    asprintf(&fname, "out/trace_%s.txt", tag);
      
    /* Allocate neuron element trace structure: */
    nmsim_elem_neuron_count_t nne_tr_MAX = 4; /* Max number of neurons in monitored set. */
    nmsim_elem_neuron_ix_t ineLo = (nmsim_elem_neuron_ix_t)int64_abrandom(0, nne-1);
    nmsim_elem_neuron_ix_t ineHi = (nmsim_elem_neuron_ix_t)(ineLo + int64_abrandom(0, nne_tr_MAX) - 1);
    if (ineHi >= nne) { ineHi = nne-1; }
    assert((ineLo >= 0) && (ineLo <= ineHi) && (ineHi < nne));
    bool_t single = (ineLo == ineHi);
    nmsim_elem_neuron_trace_t *trne = nmsim_elem_neuron_trace_new(ineLo, ineHi, tLo_t, tHi_t);
    
    /* Initialize the neuron ages, potentials, and modulators: */
    double V, age, M, H;
    nmsim_step_count_t age1;
    nmsim_class_neuron_throw_state(nclass, &(V), &(age1)); 
    age = (double)age1;
    M = nmsim_class_neuron_compute_M(nclass, age1);
    H = nmsim_class_neuron_compute_H(nclass, age1);

    /* Fake a simulation for a range of times that overlaps {tLo..tHi}: */
    nmsim_time_t t = tLo_s;
    while (t < tHi_s)
      { /* Store into {trne} the state {V,age,M,H} at time {t}, if applicable: */
        nmsim_elem_neuron_trace_set_V_age_M_H(trne, t, V, age, M, H);
        /* Define firing indicator for step {t} to {t+1}: */
        double X = (single ? (double)(drandom() < 0.05) : drandom());
        /* Define sinusoidal {I}, random {J} for step between {t} and {t+1}: */
        double freq = 0.07;
        double I = 10.0*sin(0.3 + freq*((double)t));
        double J = I + 0.5*dgaussrand();
        nmsim_elem_neuron_trace_set_X_I_J(trne, t, X, I, J);
        /* Update {V,age,M,H}: */
        V = X*nclass->V_R + (1-X)*((V + I)*M*nclass->c_B + 0.5*dgaussrand());
        age = X*0 + (1-X)*(age + 1);
        M = X*nclass->M_R + (1-X)*(1 - (1 - M)*nclass->M_mu);
        H = X*nclass->H_R + (1-X)*(1 - (1 - H)*nclass->H_mu);
        /* Advance the time: */ 
        t++;
      }
    /* Store into {trne} the state {V,age,M,H} at final simulated time {t}, if applicable: */
    nmsim_elem_neuron_trace_set_V_age_M_H(trne, t, V, age, M, H);     
      
    /* Write the traces: */
    nmsim_test_elem_neuron_trace_write(fname, trne);
    
    /* Read tem back and compare: */
    nmsim_test_elem_neuron_trace_read(fname, trne);
    
    free(fname);
    nmsim_elem_neuron_trace_free(trne);
  }

void nmsim_test_elem_neuron_trace_write(char *fname, nmsim_elem_neuron_trace_t *trne)
  { 
    FILE *wr = open_write(fname, TRUE);
    nmsim_elem_neuron_trace_write(wr, trne);
    fclose(wr);
  }

void nmsim_test_elem_neuron_trace_read(char *fname, nmsim_elem_neuron_trace_t *trne)
  {
    FILE *rd = open_read(fname, TRUE);
    nmsim_elem_neuron_trace_t *trne_read = nmsim_elem_neuron_trace_read(rd);
    fclose(rd);
    nmsim_elem_neuron_trace_compare(trne_read, trne);
    nmsim_elem_neuron_trace_free(trne_read);
  }
