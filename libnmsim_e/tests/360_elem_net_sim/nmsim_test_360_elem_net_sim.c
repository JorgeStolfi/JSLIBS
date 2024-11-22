#define PROG_NAME "nmsim_test_360_elem_net_sim"
#define PROG_DESC "tests of {limnmism} neuron-level network simulation"
#define PROG_VERS "1.0"

/* Last edited on 2024-11-20 05:28:30 by stolfi */ 

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
#include <nmsim_elem_net_trace.h>

#include <nmsim_elem_net_sim.h>

void nmsim_test_elem_net_sim
  ( int32_t nnc,
    int32_t nsc,
    int32_t nng, 
    int32_t nsg, 
    int32_t nne,
    int32_t nse
  );
  /* Tests a neuron net with {nnc} neuron classes, {nsc} synapse classes,
    {nng} neuron groups, {nsg} neuron groups, {nne} neurons, and {nse} synapses. 
    Each neuron class will have about {nng/nnc} groups, and each neuron group
    will have about {nne/nng} neurons; and similarly for synapses.
    
    Output files will be called "out/sim_{TAG}_*.txt" where {TAG} is a string formed 
    from the parameters.   
    
    A trace of some neurons will be written to files
    "out/sim_{TAG}_elem_n{INE}.txt", where {INE} is the neuron index
    {i}, zero-padded to 10 digits. If the number of neurons is not
    small, only a small randomly selected subset will be traced.
    
    !!! IMPLEMENT THIS: A trace of every neuron group, with lines
    of the form "{t} {bV_k[t]} {bage_k[t]} {bX_k[t]} {bI_k[t]} {bJ_k[t]} {bM_k[t]} {bH_k[t]}", 
    will be written to files "out/sim_{TAG}_grp_n{ING}.txt", where {ING} is the neuron group index
    {k}, zero-padded to 10 digits. All the "b" quantities above are averages
    over all neurons in the group. !!!
    
   A description of the simulated network is written to file "out/sim_{TAG}_elem_net.txt". */

void nmsim_test_elem_net_write(char *prefix, nmsim_elem_net_t *enet, double timeStep);
  /* Writes a description of the network to file "{prefix}_elem_net.txt". */

void nmsim_test_elem_net_trace_write(char *prefix, nmsim_elem_net_trace_t *etrace);
  /* Writes the elem-level traces in {etrace} with names "{prefix}_elem_n{NN}.txt"
    where {NN} is the neuron index. */

int32_t main(int32_t argc, char **argv);

/* IMPLEMENTATIONS: */

int32_t main(int32_t argc, char **argv)
  { 
    /* nmsim_test_elem_net_sim(1,1,1,1,1,1); */
    nmsim_test_elem_net_sim(1,1,3,12,21,120);
    return 0;
  }
    
void nmsim_test_elem_net_sim
  ( int32_t nnc,
    int32_t nsc,
    int32_t nng, 
    int32_t nsg, 
    int32_t nne,
    int32_t nse
  )
  {
    demand((nnc >= 1) && (nsc >= 1), "bad neuron/synapse class count");
    nmsim_class_net_t *cnet = nmsim_class_net_throw(nnc, nsc);
    
    /* Compute number of neurons per neuron group: */
    demand(nng >= 1, "bad neuron group count");
    demand((nne % nng) == 0, "bad neuron elem count");
    
    /* Compute the average indegree of each neuron: */
    demand(nsg >= 1, "bad synapse group count");
    demand((nse % nsg) == 0, "bad synapse elem count");

    /* Choose simulation time parameters: */
    nmsim_time_t nSteps = 10000; /* Simulate from {t=0} to {t=nSteps}. */
    nmsim_time_t tLo = 0;         /* First time to save in traces. */
    nmsim_time_t tHi = nSteps-1;  /* Last time to save in traces. */
    double timeStep = 1.0; /* Nominal time step (ms). */
    
    /* Create the filename prefix: */
    char *prefix = jsprintf(
        "out/sim_%04dnc_%04dsc_%06dng_%06dsg_%06dne_%06dse", 
        nnc, nsc, nng, nsg, nne, nse
      );
   
    /* Create the group-levek network: */
    nmsim_group_net_t *gnet = nmsim_group_net_throw(cnet, nng, nsg, nne, nse);
    
    /* Create and write out the element_level network: */
    nmsim_elem_net_t *enet = nmsim_elem_net_throw(gnet);
    nmsim_test_elem_net_write(prefix, enet, timeStep); 
      
    /* Allocate element-level trace structure, choose neurons to trace: */
    nmsim_elem_neuron_count_t tne_max = 5; /* Max neurons to trace. */
    nmsim_elem_neuron_count_t tne = (nmsim_elem_neuron_count_t)imin(tne_max, nne);
    nmsim_elem_net_trace_t *etrace = nmsim_elem_net_trace_throw(nne, tLo, tHi, tne);
      
    /* Allocate and initialize per-group neurons state statistics structure: */
    nmsim_elem_net_sim_group_stats_t *gstats = NULL;
    
    /* Allocate work arrays: */
    double *V = notnull(malloc(nne*sizeof(double)), "no mem");
    nmsim_step_count_t *age = notnull(malloc(nne*sizeof(nmsim_step_count_t)), "no mem");
    bool_t *X = notnull(malloc(nne*sizeof(bool_t)), "no mem");
    double *M = notnull(malloc(nne*sizeof(double)), "no mem");
    double *H = notnull(malloc(nne*sizeof(double)), "no mem");
    double *I = notnull(malloc(nne*sizeof(double)), "no mem");
    double *J = notnull(malloc(nne*sizeof(double)), "no mem");
    
    /* Initialize the neuron ages, potentials, and modulators: */
    for (nmsim_elem_neuron_ix_t ine = 0; ine < nne; ine++) 
      { nmsim_group_neuron_ix_t ing = enet->neu[ine].ing; /* Neuron group index. */
        nmsim_class_neuron_ix_t inc = gnet->ngrp[ing].inc; /* Neuron class index. */
        nmsim_class_neuron_throw_state(cnet->nclass[inc], &(V[ine]), &(age[ine]));
      }
    nmsim_elem_net_sim_compute_modulators(enet, 0, age, M, H);

    /* Simulate: */
    for (nmsim_time_t t = 0; t < nSteps; t++)
      { /* Define the external inputs (zero for now): */
        for (nmsim_elem_neuron_ix_t ine = 0; ine < nne; ine++) { I[ine] = 0.0; } 
        /* Apply the evolution equations: */
        nmsim_elem_net_sim_step(enet, t, V, age, M, H, X, I, J, etrace, gstats);
      }
      
    /* Write the traces: */
    nmsim_test_elem_net_trace_write(prefix, etrace);
    
    free(V);
    free(age);
    free(X);
    free(M); free(H);
    free(I); free(J);
    free(prefix);
    nmsim_elem_net_trace_free(etrace);
  }
  
void nmsim_test_elem_net_write(char *prefix, nmsim_elem_net_t *enet, double timeStep)
  { char *fname = NULL;
    char *fname = jsprintf("%s_elem_net.txt", prefix);
    FILE *wr = open_write(fname, TRUE);
    nmsim_elem_net_write(wr, enet, timeStep);
    free(fname);
  }

void nmsim_test_elem_net_trace_write(char *prefix, nmsim_elem_net_trace_t *etrace)
  { char *epref = NULL;
    char *epref = jsprintf("%s_elem", prefix);
    nmsim_elem_net_trace_write(epref, etrace);
    free(epref);
  }
         
