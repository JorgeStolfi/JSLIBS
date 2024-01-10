#define PROG_NAME "nmsim_test_052_elem_net_sim"
#define PROG_DESC "tests of {limnmism} neuron-level network simulation"
#define PROG_VERS "1.0"

/* Last edited on 2019-02-13 17:47:39 by jstolfi */ 

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
#include <nmsim_firing_func.h>
#include <nmsim_class_neuron.h>
#include <nmsim_class_synapse.h>
#include <nmsim_group_neuron.h>
#include <nmsim_group_synapse.h>
#include <nmsim_elem_neuron.h>
#include <nmsim_elem_synapse.h>

#include <nmsim_class_net.h>
#include <nmsim_group_net.h>
#include <nmsim_elem_net.h>
#include <nmsim_elem_net_trace.h>

#include <nmsim_elem_net_sim.h>

void nmsim_elem_net_test
  ( int32_t nnc,
    int32_t nsc,
    int32_t nng, 
    int32_t nsg, 
    int32_t nne,
    int32_t nse,
    int32_t nnp
  );
  /* Tests a neuron net with {nnc} neuron classes, {nsc} synapse classes,
    {nng} neuron groups, {nsg} neuron groups, {nne} neurons, and {nse} synapses. 
    Each neuron class will have about {nng/nnc} groups, and each neuron group
    will have about {nne/nng} neurons; and similarly for synapses.
    
    Output files will be called "out/sim_{TAG}_*.txt" where {TAG} is a string formed 
    from the parameters.  A trace of {nnp} randomly selected neurons 
    will be written to files "out/sim_{TAG}_elem_n{INE}.txt", where {INE} is the neuron index
    {i}, zero-padded to 10 digits.  
    
    A trace of every neuron group, with lines
    of the form "{t} {bV_k[t]} {bage_k[t]} {bX_k[t]} {bI_k[t]} {bJ_k[t]} {bM_k[t]} {bH_k[t]}", 
    will be written to files "out/sim_{TAG}_grp_n{ING}.txt", where {ING} is the neuron group index
    {k}, zero-padded to 10 digits. All the "b" quantities above are averages
    over all neurons in the group.
    
    A description of the simulated network is written to file "out/sim_{TAG}_elem_net.txt". */

nmsim_elem_net_trace_t *nmsim_elem_net_test_trace_new
  ( nmsim_elem_net_t *enet, 
    nmsim_time_t tini, 
    nmsim_time_t tfin, 
    nmsim_elem_neuron_count_t nnp
  );
  /* Creates a trace data structure for of {nnp} randomly selected neurons
    in the network {enet}, for discrete times {tini..tfin}. */

void nmsim_elem_net_test_trace_write(char *prefix, nmsim_elem_net_trace_t *etrace);
  /* Writes the elem-level traces in {etrace} with names "{prefix}_elem_n{NN}.txt"
    where {NN} is the neuron index. */

int main(int argc, char **argv);

/* IMPLEMENTATIONS: */

int main(int argc, char **argv)
  { 
    nmsim_elem_net_test(1,1,1,1,1,1,1);
    nmsim_elem_net_test(1,1,2,1,20,100,3);
    return 0;
  }
    
void nmsim_elem_net_test
  ( int32_t nnc,
    int32_t nsc,
    int32_t nng, 
    int32_t nsg, 
    int32_t nne,
    int32_t nse,
    int32_t nnp
  )
  {
    demand((nnc >= 1) && (nsc >= 1), "bad neuron/synapse class count");
    nmsim_class_net_t *cnet = nmsim_class_net_throw(nnc, nsc);
    
    /* Compute number of neurons per neuron group: */
    demand(nng >= 1, "bad neuron group count");
    demand((nne % nng) == 0, "bad neuron elem count");
    int32_t nne_g = nne/nng; /* Num neurons per neuron group. */
    assert(nne_g > 0);
    
    /* Compute the average indegree of each neuron: */
    demand(nsg >= 1, "bad synapse group count");
    demand((nse % nsg) == 0, "bad synapse elem count");
    int32_t nse_g = nse/nsg; /* Desired num synapses per synapse group. */
    assert(nse_g > 0);
    double K = ((double)nse_g)/((double)nne_g); /* Desired avg indegre. */
   
    /* Create the group-levek network: */
    nmsim_group_net_t *gnet = nmsim_group_net_throw(cnet, nng, nsg, nne_g, nne_g, K, K);
    
    /* Create the element_level network: */
    nmsim_elem_net_t *enet = nmsim_elem_net_throw(gnet);

    /* Choose simulation time parameters: */
    nmsim_time_t nSteps = 10000; /* Simulate from {t=0} to {t=nSteps}. */
    nmsim_time_t tini = 0;         /* First time to save in traces. */
    nmsim_time_t tfin = nSteps-1;  /* Last time to save in traces. */
    
    /* Create the filename prefix: */
    char *prefix = NULL;
    asprintf
      ( &prefix, 
        "out/sim_%04dnc_%04dsc_%06dng_%06dsg_%06dne_%06dse", 
        nnc, nsc, nng, nsg, nne, nse
      );
      
    /* Allocate element-level trace structure, choose neurons to trace: */
    nmsim_elem_net_trace_t *etrace = nmsim_elem_net_test_trace_new(enet, tini, tfin, nnp);
    
    /* Allocate work arrays: */
    double *V = notnull(malloc(nne*sizeof(double)), "no mem");
    nmsim_step_count_t *age = notnull(malloc(nne*sizeof(nmsim_step_count_t)), "no mem");
    bool_t *X = notnull(malloc(nne*sizeof(bool_t)), "no mem");
    double *M = notnull(malloc(nne*sizeof(double)), "no mem");
    double *H = notnull(malloc(nne*sizeof(double)), "no mem");
    double *I = notnull(malloc(nne*sizeof(double)), "no mem");
    double *J = notnull(malloc(nne*sizeof(double)), "no mem");
    
    /* Initialize the neuron ages and potentials: */
    for (nmsim_elem_neuron_ix_t ine = 0; ine < nne; ine++) 
      { nmsim_group_neuron_ix_t ing = enet->neu.e[ine].ing; /* Neuron group index. */
        nmsim_class_neuron_ix_t inc = gnet->ngrp.e[ing].inc; /* Neuron class index. */
        nmsim_class_neuron_throw_state(cnet->nclass.e[inc], &(V[ine]), &(age[ine]));
      }
    /* Simulate: */
    for (nmsim_time_t t = 0; t < nSteps; t++)
      { /* Define the external inputs (zero for now): */
        for (nmsim_elem_neuron_ix_t ine = 0; ine < nne; ine++) { I[ine] = 0.0; } 
        /* Apply the evolution equations: */
        nmsim_elem_net_sim_step(enet, t, I, V, age, X, H, M, J, etrace);
      }
      
    /* Write the traces: */
    nmsim_elem_net_test_trace_write(prefix, etrace);
    
    free(V);
    free(age);
    free(X);
    free(M); free(H);
    free(I); free(J);
    free(prefix);
    nmsim_elem_net_trace_free(etrace);
  }

nmsim_elem_net_trace_t *nmsim_elem_net_test_trace_new
  ( nmsim_elem_net_t *enet, 
    nmsim_time_t tini, 
    nmsim_time_t tfin, 
    nmsim_elem_neuron_count_t nnp
  )
  {
    nmsim_group_net_t *gnet = enet->gnet; /* Group-level network description. */
    /* nmsim_class_net_t *cnet = gnet->cnet; */ /* Class-level network description. */

    nmsim_elem_neuron_count_t nne = enet->nne; /* Num of neurons in network. */
    demand((nnp >= 0) & (nnp <= nne), "invalid monitored neuron count");
    
    /* Create trace structure: */
    nmsim_elem_net_trace_t *etrace = nmsim_elem_net_trace_new(tini, tfin, nnp);
    
    /* Select {nnp} raandom neurons out of {nne}, put in {perm[0..nnp-1]}: */
    nmsim_elem_neuron_ix_t *perm = notnull(malloc(nne*sizeof(nmsim_elem_neuron_ix_t)), "no mem");
    for (nmsim_elem_neuron_ix_t ine = 0; ine < nne; ine++) { perm[ine] = ine; }
    if (nnp < nne)
      { for (nmsim_elem_neuron_ix_t i = 0; i < nnp; i++) 
          { nmsim_elem_neuron_ix_t j = (nmsim_elem_neuron_ix_t)int64_abrandom(i, nne-1); 
            nmsim_elem_neuron_ix_t k = perm[i]; perm[i] = perm[j]; perm[j] = k;
          }
      }
      
    /* Create trace structures for those neurons: */
    for (nmsim_elem_neuron_ix_t i = 0; i < nnp; i++) 
      { nmsim_elem_neuron_ix_t ine = perm[i]; /* Neuron elem index. */
        nmsim_group_neuron_ix_t ing = enet->neu.e[ine].ing; /* Neuron group index. */
        nmsim_class_neuron_ix_t inc = gnet->ngrp.e[ing].inc; /* Neuron class index. */
        nmsim_elem_neuron_trace_t *trne = nmsim_elem_neuron_trace_new(ine, ing, inc, tini, tfin);
        etrace->trne[i] = trne;
      }
      
    free(perm);
    return etrace;
  }

void nmsim_elem_net_test_trace_write(char *prefix, nmsim_elem_net_trace_t *etrace)
  { char *epref = NULL;
    asprintf(&epref, "%s_elem", prefix);
    nmsim_elem_net_trace_write(epref, etrace);
    free(epref);
  }
         
