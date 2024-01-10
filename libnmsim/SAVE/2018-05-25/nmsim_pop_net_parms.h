#ifndef nmsim_pop_net_parms_H
#define nmsim_pop_net_parms_H
 
/* Population-level modeling of networks of Galves-Löcherbach neurons. */
/* Last edited on 2018-04-11 10:47:20 by jstolfi */

#define _GNU_SOURCE

#include <nmsim_basic.h>
#include <nmsim_neuron_parms.h>
#include <nmsim_bundle_parms.h>

typedef struct nmsim_pop_net_parms_t
  { nmsim_pop_count_t np;         /* Number of populations. */
    nmsim_bundle_count_t nb;      /* Number of synaptic bundles between populations. */
    /* These arrays are indexed with population index {0..np-1}: */
    nmsim_neuron_parms_t **ppop;  /* Neuron parameters in each population. */
    nmsim_neuron_count_t *nn;     /* Number of neurons in each population. */
    nmsim_pop_count_t *deg_in;    /* Number of input bundles in each population. */
    nmsim_pop_count_t *deg_ot;    /* Number of output bundles in each population. */
    /* These parameters are indexed with the connection bundle index {0..nb-1}: */
    nmsim_bundle_parms_t **pbun;  /* Connections between populations. */
  } nmsim_pop_net_parms_t;
  /* Statistical description of a network with {np} homogeneous populations of neurons,
    with specified ensemble parameters for neurons in each population
    and ensemble parameters for the connections between each pair of populations.
    
    The parameters of population {p} is {.ppop[p]}, for {p} in {0..np-1}.
    The population is assumed to have {.nn[r]} neurons.
    
    The parameters of the bundle of synapses {b}, including its origin and 
    destination populations, are in {.pbun[b]}, for {b} in {0..nb-1}. */
    
nmsim_pop_net_parms_t *nmsim_pop_net_parms_new
  ( nmsim_pop_count_t np,
    nmsim_bundle_count_t nb
  );
  /* Allocates an {nmsim_pop_net_t} structure with {np} populations
    neurons and {nb} synaptic bundles.  Allocates the tables {.ppop[0..np-1]},
    {.nn[0..np-1]}, {.deg_in[0..np-1]}, {.deg_ot[0..np-1]},
    and {.pbun[0..nb-1]}, but sets all the entres to {NULL}. */

#endif
