#define PROG_NAME "nmsim_test_330_elem_neuron_trace"
#define PROG_DESC "tests of {limnmism} single neuron simulation trace"
#define PROG_VERS "1.0"

/* Last edited on 2019-03-28 22:08:27 by jstolfi */ 

#define PROG_COPYRIGHT \
  "Copyright © 2019  State University of Campinas (UNICAMP)"
  
#define PROG_AUTH \
  "Created by J. Stolfi on 2019-01-11"
  
#define PROG_HIST
  
#define PROG_HELP \
  "  " PROG_NAME ""

#define _GNU_SOURCE
#include <stdio.h>
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
    nmsim_time_t tlo_s, 
    nmsim_time_t thi_s,
    nmsim_time_t tlo_t, 
    nmsim_time_t thi_t
  );
  /* Tests a neuron element trace structure that can save neuron states for time
    {t} and step data for the step from {t} to {t+1}, where {t} is in
    {tlo_t .. thi_t}; that is, with {thi_t - tlo_t + 1} entries.  
    
    The simulation is run from {tlo_s}
    to {thi_s} inclusive, that is {thi_s - tlo_s} steps.
    
    Output file will be called "out/trace_{tag}.txt". */

void nmsim_test_elem_neuron_trace_write(char *fname, nmsim_elem_neuron_trace_t *etrace);
  /* Writes the neuron element trace {etrace} to a file with name {fname}. */

void nmsim_test_elem_neuron_trace_read(char *fname, nmsim_elem_neuron_trace_t *etrace);
  /* Reads a neuron element trace from file {fname} and compares
    its contents to that of the given {etrace}. */

int main(int argc, char **argv);

/* IMPLEMENTATIONS: */

int main(int argc, char **argv)
  { 
    nmsim_test_elem_neuron_trace("A", 418, 418, 400, 450);
    nmsim_test_elem_neuron_trace("B",  10,  20,   7,  17);
    nmsim_test_elem_neuron_trace("C", -10,  10,  -7,  20);
    nmsim_test_elem_neuron_trace("C", -30,  30, -17,  27);
    return 0;
  }
    
void nmsim_test_elem_neuron_trace
  ( char *tag,
    nmsim_time_t tlo_s, 
    nmsim_time_t thi_s,
    nmsim_time_t tlo_t, 
    nmsim_time_t thi_t
  )
  {
    fprintf(stderr, "simulation times %5ld .. %5ld ( %5ld steps   )\n", tlo_s, thi_s, thi_s - tlo_s);
    fprintf(stderr, "tracing times    %5ld .. %5ld ( %5ld entries )\n", tlo_t, thi_t, thi_t - tlo_t + 1);
    
    /* Create a neuron class: */
    nmsim_class_neuron_count_t nnc = 10; /* Arbitrary number of neuron classes in network. */
    nmsim_class_neuron_t *nclass = nmsim_class_neuron_throw();
    
    /* Pick an index for the neuron: */
    nmsim_group_neuron_count_t nng = 2*nnc; /* Arbitrary number of neuron groups in network. */ 
    nmsim_elem_neuron_count_t nne = 10*nng; /* Arbitrary num of neuron elems in network. */
    nmsim_elem_neuron_ix_t ine = (nmsim_elem_neuron_ix_t)int64_abrandom(0, nne-1);
      
    /* Create the filename: */
    char *fname = NULL;
    asprintf(&fname, "out/trace_%s.txt", tag);
      
    /* Allocate neuron element trace structure: */
    nmsim_elem_neuron_trace_t *trne = nmsim_elem_neuron_trace_new(ine, tlo_t, thi_t);
    
    /* Full neuron state for simulation: */
    nmsim_elem_neuron_trace_entry_t ts;
    nmsim_elem_neuron_trace_entry_clear(&ts);
    
    /* Initialize the neuron ages, potentials, and modulators: */
    nmsim_class_neuron_throw_state(nclass, &(ts.V), &(ts.age));
    ts.M = nmsim_class_neuron_compute_M(nclass, ts.age);
    ts.H = nmsim_class_neuron_compute_H(nclass, ts.age);

    /* Fake a simulation for a range of times that overlaps {tlo..thi}: */
    nmsim_time_t t = tlo_s;
    while (t < thi_s)
      { 
        /* Store into {trne} the state {V,age,M,H} at time {t}, if applicable: */
        nmsim_elem_neuron_trace_set_V_age_M_H(trne, t, ts.V, ts.age, ts.M, ts.H);
        /* Define firing indicator for step {t} to {t+1}: */
        ts.X = (drandom() < 0.05);
        /* Define sinusoidal {I}, random {J} for step between {t} and {t+1}: */
        double freq = 0.07;
        ts.I = 10.0*sin(0.3 + freq*((double)t));
        ts.J = ts.I + 0.5*dgaussrand();
        nmsim_elem_neuron_trace_set_X_I_J(trne, t, ts.X, ts.I, ts.J);
        /* Update {V,age,M,H}: */
        if (ts.X)
          { /* Reset: */
            ts.V = nclass->V_R;
            ts.age = 0;
            ts.M = nclass->M_R;
            ts.H = nclass->H_R;
          }
        else
          { /* Integrate: */
            ts.V = (ts.V + ts.I)*ts.M*nclass->c_B + 0.5*dgaussrand();
            ts.age = ts.age + 1;
            ts.M = 1 - (1 - ts.M)*nclass->M_mu;
            ts.H = 1 - (1 - ts.H)*nclass->H_mu;
          }
        /* Advance the time: */ 
        t++;
      }
    /* Store into {trne} the state {V,age,M,H} at final simulated time {t}, if applicable: */
    nmsim_elem_neuron_trace_set_V_age_M_H(trne, t, ts.V, ts.age, ts.M, ts.H);     
      
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
