#ifndef nmsim_elem_neuron_state_H
#define nmsim_elem_neuron_state_H
 
/* Minimal state of a Galves-Löcherbach neuron. */
/* Last edited on 2019-03-21 10:37:56 by jstolfi */

#define _GNU_SOURCE
#include <stdint.h>
#include <stdio.h>

#include <bool.h>

#include <nmsim_basic.h>
#include <nmsim_write.h>
#include <nmsim_read.h>
#include <nmsim_compare.h>

typedef struct nmsim_elem_neuron_state_t
  { 
    double V;               /* Membrane potential at {t}. */
    nmsim_step_count_t age; /* Firing age at {t}(whole time steps since last firing). */
  } nmsim_elem_neuron_state_t;
  /* In some contexts, the values {V=NAN} and {age=-1} may be used to indicate 
    an "undefined" state.  Otherwise, {V} must be a valid potential in millivolts,
    and {age} must be non-negative. */
  
void nmsim_elem_neuron_state_write(FILE *wr, nmsim_elem_neuron_state_t *st);
  /* Writes {st} to file {wr}, in the format described by 
    {nmsim_elem_neuron_state_read_INFO} below. */
    
void nmsim_elem_neuron_state_read(FILE *rd, nmsim_elem_neuron_state_t *st);
  /* Reads a neuron state from {rd}, in the format described by 
    {nmsim_elem_neuron_state_read_INFO} below, and stores it into {*st}. 
    
    The data is assumed to be all in one line. Also reads the final
    end-of-line. Skips lines that are totally blank. Fails if the line
    is malformed or runs into end-of-file.  May return an undefined 
    state, with {V=NAN} and {age=-1}.  */
    
#define nmsim_elem_neuron_state_read_INFO \
  "Each state is written in a single line of the form\n" \
  "      \"{V} {AGE}\"\n" \
  "where {V} is the membrane potential and {age} is the firing" \
  " age.  A negative age indicates that the state is unknown or undefined."

#endif
