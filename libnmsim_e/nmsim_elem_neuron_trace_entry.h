#ifndef nmsim_elem_neuron_trace_entry_H
#define nmsim_elem_neuron_trace_entry_H
 
/* Full state and step evolution data of a Galves-Löcherbach neuron. */
/* Last edited on 2020-12-16 23:05:50 by jstolfi */

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
    double V;     /* Membrane potential at {t}. */
    double age;   /* Firing age at {t} (whole time steps since last firing). */
    double M;     /* Potential decay modulator at time {t}. */
    double H;     /* Output synapse strength modulator at time {t}. */
    double X;     /* Firing indicator between {t} and {t+1}. */
    double I;     /* External input between {t} and {t+1} (mV). */
    double J;     /* Total input between {t} and {t+1} (mV). */
  } nmsim_elem_neuron_trace_entry_t;
  /* State and activity variables of a neuron element, or average 
    of variables of two or more neurons. */
  
void nmsim_elem_neuron_trace_entry_clear(nmsim_elem_neuron_trace_entry_t *ts);
  /* Sets all fields of {ts} to {NAN}. */
    
bool_t nmsim_elem_neuron_trace_entry_is_undefined(nmsim_elem_neuron_trace_entry_t *ts);
  /* Returns true if {ts} is {NULL}, or if all the fields
    are {NAN}.  Returns false {ts} is not NULL and none of its fields is {NAN}. Fails in all other cases.  
    
    Note that the entries {ts.X,ts.I,ts.J} may be {NAN} even when {ts.V,ts.age,ts.M,ts.H} are not {NAN}. */
  
void nmsim_elem_neuron_trace_entry_write(FILE *wr, nmsim_elem_neuron_trace_entry_t *ts, bool_t single);
  /* Writes {ts} to file {wr}, in the format described by 
    {nmsim_elem_neuron_trace_entry_read_INFO} below.  Does NOT write
    a final newline. 
    
    If {single} is true, assumes that the data refers to a single neuron.  It implies 
    that {ts.age} must be a non-negative integer, and {ts.X} must be 0 or 1. */
    
void nmsim_elem_neuron_trace_entry_read(FILE *rd, nmsim_elem_neuron_trace_entry_t *ts);
  /* Reads a neuron state from {rd}, in the format described by 
    {nmsim_elem_neuron_trace_entry_read_INFO} below, and stores it into {*ts}. 
    
    The data is assumed to be all in one line. Also reads the final
    end-of-line. Skips lines that are totally blank. Fails if the line
    is malformed or runs into end-of-file. */
    
#define nmsim_elem_neuron_trace_entry_read_INFO \
  "Each entry contains the following fields, all in the same line:\n" \
  "      \"{V} {age} {M} {H} {X} {I} {J}\"\n" \
  " where {V,AGE,M,H} are the membrane potential, firing age, potential decay" \
  " modulator, and output modulator at some time {t}, {X,I,J} are the" \
  " firing indicator, external input, and total input between" \
  " time {t} and {t+1}.\n" \
  "\n" \
  "   For a single neuron, the firing indicator {X} is a Boolean" \
  " value (0 or 1), the {age} is an integer value. Otherwise, {X} is a" \
  " fraction between 0 and 1, and {age} is a non-negative fraction. All other fields " \
  " can be fractions. Some of the values may be \"nan\"."

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
