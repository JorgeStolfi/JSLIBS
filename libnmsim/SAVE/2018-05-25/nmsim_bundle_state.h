#ifndef nmsim_bundle_state_H
#define nmsim_bundle_state_H

/* State of a bundle of connections two populations. */ 
/* Last edited on 2018-04-11 10:35:12 by jstolfi */

#define _GNU_SOURCE
#include <stdint.h>

#include <nmsim_neuron_parms.h>
#include <nmsim_bundle_parms.h>

typedef struct nmsim_bundle_state_t
  { 
  } nmsim_bundle_state_t;
  /* State of a set of synapses from a homogeneous population {A}
    of neurons to another homogeneous population {B} (possibly the
    same as {A}). */

#endif
