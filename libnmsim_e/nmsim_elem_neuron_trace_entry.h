#ifndef nmsim_elem_neuron_trace_entry_H
#define nmsim_elem_neuron_trace_entry_H
 
/* Full state and step evolution data of a Galves-Löcherbach neuron. */
/* Last edited on 2020-12-06 19:24:21 by jstolfi */

#define _GNU_SOURCE
#include <stdint.h>
#include <stdio.h>

#include <bool.h>

#include <nmsim_basic.h>
#include <nmsim_write.h>
#include <nmsim_read.h>
#include <nmsim_compare.h>

typedef struct nmsim_elem_neuron_trace_entry_t
  { 
    double V;               /* Membrane potential at {t}. */
    nmsim_step_count_t age; /* Firing age at {t} (whole time steps since last firing). */
    double M;               /* Potential decay modulator at time {t}. */
    double H;               /* Output synapse strength modulator at time {t}. */
    bool_t X;               /* Firing indicator between {t} and {t+1}. */
    double I;               /* External input between {t} and {t+1} (mV). */
    double J;               /* Total input between {t} and {t+1} (mV). */
  } nmsim_elem_neuron_trace_entry_t;
  /* Full state of a neuron element. */
  
void nmsim_elem_neuron_trace_entry_clear(nmsim_elem_neuron_trace_entry_t *ts);
  /* Sets all {double} fields of {ts} to {NAN}, the firing age {ts.age} to {-1}, and the firing
    indicator {ts.X} with {FALSE}. */
    
bool_t nmsim_elem_neuron_trace_entry_is_undefined(nmsim_elem_neuron_trace_entry_t *ts);
  /* Returns true if {ts} is {NULL}, or if all the {double} fields
    are {NAN}, the firing age {ts.age} is {-1}, and the firing indicator {ts.X} is {FALSE}.
    Returns false {ts} is not NULL, none of {ts.V,ts.M,ts.H} is {NAN}, and the firing age 
    is non-negative. Fails in all other cases.  
    
    Note that the entries {ts.I} and {ts.J} may be {NAN} even when {ts.V,ts.M,ts.H} are not {NAN}. */
  
void nmsim_elem_neuron_trace_entry_write(FILE *wr, nmsim_elem_neuron_trace_entry_t *ts);
  /* Writes {ts} to file {wr}, in the format described by 
    {nmsim_elem_neuron_trace_entry_read_INFO} below.  Does NOT write
    a final newline. */
    
void nmsim_elem_neuron_trace_entry_read(FILE *rd, nmsim_elem_neuron_trace_entry_t *ts);
  /* Reads a neuron state from {rd}, in the format described by 
    {nmsim_elem_neuron_trace_entry_read_INFO} below, and stores it into {*ts}. 
    
    The data is assumed to be all in one line. Also reads the final
    end-of-line. Skips lines that are totally blank. Fails if the line
    is malformed or runs into end-of-file. */
    
#define nmsim_elem_neuron_trace_entry_read_INFO \
  "Each entry contains the following fields, all in the same line:\n" \
  "      \"{V} {AGE} {M} {H} {X} {I} {J}\"\n" \
  " where {V,AGE,M,H} are the membrane potential, firing age, potential decay" \
  " modulator, and output modulator at some time {t}, {X,I,J} are the" \
  " firing indicator, external input, and total input between" \
  " time {t} and {t+1}.\n" \
  "\n" \
  "   The firing indicator is a Boolean (0 or 1), and all other fields are" \
  " floats.  A negative {age} indicates an unknown or undefined neuron state."

/* TESTING AND DEBUGGING */

void nmsim_elem_neuron_trace_entry_compare
  ( nmsim_elem_neuron_trace_entry_t *ts_read, 
    nmsim_elem_neuron_trace_entry_t *ts_orig
  );
  /* Compares the elem-level neuron trace entry {ts_read} read from a file
    with the expected trace entry {ts_orig}. Aborts if there are discrepancies (beyond
    roundoff tolerance). In particular, if one of them is {NULL}, the other must be {NULL}
    or totally undefined. */

#endif
