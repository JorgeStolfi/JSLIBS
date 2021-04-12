#define PROG_NAME "nmsim_test_410_elem_net_sim_stats"
#define PROG_DESC "tests of {limnmism} neuron-level network simulation"
#define PROG_VERS "1.0"

/* Last edited on 2020-12-17 11:08:38 by jstolfi */ 

#define PROG_COPYRIGHT \
  "Copyright © 2020  State University of Campinas (UNICAMP)"
  
#define PROG_AUTH \
  "Created by J. Stolfi on 2020-12-15"
  
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

void nmsim_test_elem_net_sim_stats(int32_t nne);
  /* Tests the simulation statistics tools with a hypothetical neuron of 
    {nne} neurons.  The network is actually not even created.
    
    Output files will be called "out/sim_{TAG}_stats.txt" where {TAG} is a string formed 
    from the parameters. */

void nmsim_test_elem_net_sim_stats_write(char *prefix, nmsim_elem_net_sim_stats_t *S);
  /* Writes the statistics in {S} with name "{prefix}.txt". */

int main(int argc, char **argv);

/* IMPLEMENTATIONS: */

int main(int argc, char **argv)
  { 
    nmsim_test_elem_net_sim_stats(21);
    return 0;
  }
    
void nmsim_test_elem_net_sim_stats(int32_t nne)
  {
    /* We need a sufficient number of neurons: */
    demand(nne >= 10, "bad neuron count");

    /* Choose nominal simulation time parameters: */
    nmsim_time_t nSteps = 10000;  /* Pretend to simulate from {t=0} to {t=nSteps}. */
    
    /* Allocate the statistics structure, choose (proper) subset of neurons to statize: */
    nmsim_elem_neuron_ix_t ineLo = (nmsim_elem_neuron_ix_t)imax(nne/4, nne/2-3); /* First neuron to consider. */
    nmsim_elem_neuron_ix_t ineHi = (nmsim_elem_neuron_ix_t)imin(3*nne/4, nne/2+3); /* Last neuron to consider. */
    assert((ineLo > 0) && (ineLo < ineHi) && (ineHi < nne-1));  
    nmsim_time_t sttLo = nSteps/4;    /* First time to statize. */
    nmsim_time_t sttHi = 3*nSteps/4;  /* Last time to statize. */

    nmsim_elem_net_sim_stats_t *S = nmsim_elem_net_sim_stats_new(ineLo, ineHi, sttLo, sttHi);
    nmsim_elem_net_sim_stats_initialize(S);
    
    /* Allocate work arrays: */
    double *V = notnull(malloc(nne*sizeof(double)), "no mem");
    nmsim_step_count_t *age = notnull(malloc(nne*sizeof(nmsim_step_count_t)), "no mem");
    bool_t *X = notnull(malloc(nne*sizeof(bool_t)), "no mem");
    double *M = notnull(malloc(nne*sizeof(double)), "no mem");
    double *H = notnull(malloc(nne*sizeof(double)), "no mem");
    double *I = notnull(malloc(nne*sizeof(double)), "no mem");
    double *J = notnull(malloc(nne*sizeof(double)), "no mem");
    
    /* Pretend to simulate: */
    srandom(19501129);
    for (nmsim_time_t t = 0; t < nSteps; t++)
      { /* Fill the state vectors with random numbers: */
        for (nmsim_elem_neuron_ix_t ine = 0; ine < nne; ine++) 
          { V[ine] =    10.0 + 5.0*dgaussrand();
            age[ine] =  (t + ine) % 100;
            M[ine] =    1.0 +0.5*dgaussrand();
            H[ine] =    2.0 +3.0*dgaussrand();
            X[ine] =    (drandom() < 0.10);
            I[ine] =    10.0 + 5.0*dgaussrand();
            J[ine] =    20.0 + 10.0*dgaussrand();
          }
        /* Gather statistics: */
        nmsim_elem_net_sim_stats_accumulate_V_age_M_H(S, t, nne, V, age, M, H);
        nmsim_elem_net_sim_stats_accumulate_VF_AF_X_I_J(S, t, nne, V, age, X, I, J);
      }
      
    nmsim_elem_net_sim_stats_finalize(S);

    /* Write the statistics: */
    char *prefix = NULL;
    asprintf(&prefix, "out/sim_%06dne", nne);
    nmsim_test_elem_net_sim_stats_write(prefix, S);
    
    free(V);
    free(age);
    free(X);
    free(M); free(H);
    free(I); free(J);
    free(prefix);
    free(S);
  }
  
void nmsim_test_elem_net_sim_stats_write(char *prefix, nmsim_elem_net_sim_stats_t *S)
  { char *fname = NULL;
    asprintf(&fname, "%s_ne%010d--%010d_stats.txt", prefix, S->ineLo, S->ineHi);
    FILE *wr = open_write(fname, TRUE);
    nmsim_elem_net_sim_stats_write(wr, S);
    fclose(wr);
    free(fname);
  }
         
