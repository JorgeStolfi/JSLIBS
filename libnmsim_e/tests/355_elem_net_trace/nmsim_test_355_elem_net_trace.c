#define PROG_NAME "nmsim_test_355_elem_net_trace"
#define PROG_DESC "tests of {limnmism} neuron-level simulation trace"
#define PROG_VERS "1.0"

/* Last edited on 2022-10-20 06:48:41 by stolfi */ 

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

#include <nmsim_elem_neuron.h>

#include <nmsim_elem_neuron_trace.h>
#include <nmsim_elem_net_trace.h>

void nmsim_test_elem_net_trace
  ( char *tag,
    nmsim_elem_neuron_count_t nne,
    nmsim_elem_neuron_count_t tne
  );
  /* Tests an element-level neuron net trace structure, that
    traces {tne}  out of {nne} neurons.
    Output file will be called "out/trace_{tag}_n{NN}.txt"
    where {NN} is the neuron index. */

void nmsim_test_elem_net_trace_write(char *prefix, nmsim_elem_net_trace_t *etrace);
  /* Writes the elem-level traces in {etrace} with names "{prefix}_n{NN}.txt"
    where {NN} is the neuron index. */

void nmsim_test_elem_net_trace_read(char *prefix, nmsim_elem_net_trace_t *etrace);
  /* Reads the elem-level traces from files "{prefix}_n{NN}.txt",
    where {NN} is the neuron numbers listed in {etrace}, and compares
    their contents to those in the given {etrace}. */

int32_t main(int32_t argc, char **argv);

/* IMPLEMENTATIONS: */

int32_t main(int32_t argc, char **argv)
  { 
    nmsim_test_elem_net_trace("A",1,1);
    nmsim_test_elem_net_trace("B",20,3);
    nmsim_test_elem_net_trace("C",20,3);
    return 0;
  }
    
void nmsim_test_elem_net_trace
  ( char *tag,
    nmsim_elem_neuron_count_t nne,
    nmsim_elem_neuron_count_t tne
  )
  {
    /* Choose simulation time parameters: */
    nmsim_time_t nSteps = 100;          /* Simulate from {t=0} to {t=nSteps}. */
    nmsim_time_t tLo = 418;            /* First time to save in traces. */
    nmsim_time_t tHi = tLo + nSteps;  /* Last time to save in traces. */
    
    /* Create the filename prefix: */
    char *prefix = jsprintf("out/trace_%s", tag);
      
    /* Allocate element-level trace structure, choose neurons to trace: */
    fprintf(stderr, "creating elem net trace structure...\n");
    nmsim_elem_net_trace_t *etrace = nmsim_elem_net_trace_throw(nne, tLo, tHi, tne);
    for(nmsim_elem_neuron_ix_t k = 0; k < tne; k++)
      { fprintf(stderr, "  entry[%2d]", k);
        nmsim_elem_neuron_trace_t *trnek = etrace->trne[k];
        if (trnek != NULL)
          { fprintf(stderr, " neurons %d..%d", trnek->ineLo, trnek->ineHi);
            fprintf(stderr, " time range %ld .. %ld\n", trnek->tLo, trnek->tHi); 
          }
        else
          { fprintf(stderr, " NULL\n"); }
      }
    
    /* Allocate work arrays: */
    fprintf(stderr, "allocating pseudo-simulation state variables...\n");
    double *V = nmsim_test_NAN_vector(nne);
    nmsim_step_count_t *age = notnull(malloc(nne*sizeof(nmsim_step_count_t)), "no mem");
    double *M = nmsim_test_NAN_vector(nne);
    double *H = nmsim_test_NAN_vector(nne);
    bool_t *X = notnull(malloc(nne*sizeof(bool_t)), "no mem");
    double *I = nmsim_test_NAN_vector(nne);
    double *J = nmsim_test_NAN_vector(nne);
    
    /* Initialize the neuron ages, potentials, and modulators: */
    fprintf(stderr, "initializing state variables...\n");
    for (nmsim_elem_neuron_ix_t ine = 0; ine < nne; ine++) 
      { V[ine] = dabrandom(-60.0, -20.0);
        age[ine] = (nmsim_step_count_t)int64_abrandom(0, 100);
        M[ine] = dabrandom(0.5, 2.0);
        H[ine] = dabrandom(0.5, 2.0);
      }

    /* Fill trace: */
    fprintf(stderr, "faking the simulation...\n");
    for (nmsim_time_t t = tLo-10; t < tHi-10; t++)
      { nmsim_elem_net_trace_set_V_age_M_H(etrace, t, nne, V, age, M, H);
        for (nmsim_elem_neuron_ix_t ine = 0; ine < nne; ine++) 
          { /* Define firing indicator: */
            X[ine] = (drandom() < 0.05);
            /* Define sinusoidal {I}, random {J} for step between  {t}: */
            double freq = 0.2 + 0.02*((double)ine);
            I[ine] = 10.0*sin(0.3 + freq*((double)t));
            J[ine] = I[ine] + 0.5*dgaussrand();
          }
        nmsim_elem_net_trace_set_X_I_J(etrace, t, nne, X, I, J);
        for (nmsim_elem_neuron_ix_t ine = 0; ine < nne; ine++) 
          { /* Update {V,age,M,H}: */
            if (X[ine])
               { /* Reset: */
                 V[ine] = -40.0;
                 age[ine] = 0;
                 M[ine] = 0.5;
                 H[ine] = 2.0;
               }
             else
               { /* Integrate: */
                 V[ine] = (V[ine] + 30)*0.8 - 30;
                 age[ine] = age[ine] + 1;
                 M[ine] = 1 - (1 - M[ine])*0.75;
                 H[ine] = 1 - (1 - H[ine])*0.75;
               }
          }
        nmsim_elem_net_trace_set_V_age_M_H(etrace, t, nne, V, age, M, H);
      }
      
    /* Write the traces: */
    nmsim_test_elem_net_trace_write(prefix, etrace);
    
    /* Read tem back and compare: */
    nmsim_test_elem_net_trace_read(prefix, etrace);
    
    free(V);
    free(age);
    free(X);
    free(M); free(H);
    free(I); free(J);
    free(prefix);
    nmsim_elem_net_trace_free(etrace);
  }

void nmsim_test_elem_net_trace_write(char *prefix, nmsim_elem_net_trace_t *etrace)
  { 
    nmsim_elem_net_trace_write(prefix, etrace);
  }

void nmsim_test_elem_net_trace_read(char *prefix, nmsim_elem_net_trace_t *etrace)
  {
    nmsim_elem_neuron_count_t tne = etrace->tne;
    for (nmsim_elem_neuron_ix_t k = 0; k < tne; k++)
      { nmsim_elem_neuron_trace_t *trnek = etrace->trne[k];
        nmsim_elem_neuron_ix_t ineLo = trnek->ineLo;
        nmsim_elem_neuron_ix_t ineHi = trnek->ineHi;
        fprintf(stderr, "reading and checking trace trne[%d] of neurons %d..%d\n", k, ineLo, ineHi);
        char *fname = jsprintf("%s_ne%010d--%010d_trace.txt", prefix, ineLo, ineHi);
        FILE *rd = open_read(fname, TRUE);
        nmsim_elem_neuron_trace_t *trnek_read = nmsim_elem_neuron_trace_read(rd);
        fclose(rd);
        nmsim_elem_neuron_trace_compare(trnek_read, trnek);
        nmsim_elem_neuron_trace_free(trnek_read);
        free(fname);
      }
  }
