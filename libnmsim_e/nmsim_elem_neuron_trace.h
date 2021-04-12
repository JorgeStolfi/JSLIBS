#ifndef nmsim_elem_neuron_trace_H
#define nmsim_elem_neuron_trace_H
 
/* Trace of a single neuron in an elem-level simulation of a Galves-Löcherbach net. */
/* Last edited on 2020-12-17 02:23:58 by jstolfi */

#define _GNU_SOURCE
#include <stdint.h>
#include <stdio.h>

#include <bool.h>

#include <nmsim_basic.h>
#include <nmsim_write.h>
#include <nmsim_read.h>
#include <nmsim_compare.h>
#include <nmsim_class_neuron.h>
#include <nmsim_group_neuron.h>
#include <nmsim_elem_neuron.h>
#include <nmsim_elem_neuron_trace_entry.h>

typedef struct nmsim_elem_neuron_trace_t
  { nmsim_elem_neuron_ix_t ineLo;        /* Index of the first neuron in this trace. */
    nmsim_elem_neuron_ix_t ineHi;        /* Index of the last neuron in this trace. */
    nmsim_time_t tLo;                     /* Discrete time of first trace entry. */
    nmsim_time_t tHi;                     /* Discrete time of last trace entry. */
    nmsim_elem_neuron_trace_entry_t *ts;  /* Full states and evolution data of the neuron. */
  } nmsim_elem_neuron_trace_t;
  /* A data structure that records the state and activity data of a selected neuron, or a set of 
    consecutive neurons, with indices {ineLo..ineHi}, over the discrete times {tLo..tHi} inclusive, in
    an elem-level simulation of a GL network.
    
    The state of the neuron at discrete time {t} and the evolution
    data from {t} to {t+1} are to be stored in {ts[t-tLo]}. Note that
    the array {ts} must be allocated with at least {tHi-tLo+1}
    entries. */
   
#define nmsim_elem_neuron_trace_entries_MAX (((int64_t)1) << 24)
  /* Max number of trace entries (times), for memory allocation purposes. */
    
nmsim_elem_neuron_trace_t *nmsim_elem_neuron_trace_new
  ( nmsim_elem_neuron_ix_t ineLo,   /* Index of the first neuron of the set to trace. */
    nmsim_elem_neuron_ix_t ineHi,   /* Index of the last neuron of the set to trace. */
    nmsim_time_t tLo,                /* Discrete time of first trace entry. */
    nmsim_time_t tHi                 /* Discrete time of last trace entry. */
  );
  /* Allocates a new {nmsim_elem_neuron_trace_t} record on the heap
    for neurons {ineLo..ineHi} with space for state variables at 
    times {tLo..tHi}.  All entries are cleared with 
    {nmsim_elem_neuron_trace_entry_clear}. */
 
void nmsim_elem_neuron_trace_free(nmsim_elem_neuron_trace_t *trne);
  /* Releases the storage of the record {trne}, including internal tables. */

/* TRACE SAVING PROCEDURES 
  
  These procedures receive variables {V,X,age,M,H,I,J} for a given time {t}
  or the time step from {t} to {t+1}.  They store those values
  in the corresponding fields of the corresponding entry of the
  neuron trace {trne}, if it exists.
  
  If {trne} covers more than one neuron, the variables stored should be
  averages over those neurons. */

void nmsim_elem_neuron_trace_set_V_age_M_H
  ( nmsim_elem_neuron_trace_t *trne,
    nmsim_time_t t,  /* Discrete time. */
    double V,        /* Membrane potential of each neuron at {t}. */
    double age,      /* Firing age of each neuron at time {t}. */
    double M,        /* Recharge modulator of each neuron at time {t}. */
    double H         /* Output modulator of each neuron at time {t}. */
  );
  /* If {trne} is not {NULL}, and {t} is in {trne.tLo .. trne.tHi},
    stores the given values in the entry of {trne} that corresponds to time {t}, 
    namely {trne.ts[t-trne.tLo]}. If {trne} refers to more than one 
    neuron, the given values must be averages over those neurons. */

void nmsim_elem_neuron_trace_set_X_I_J
  ( nmsim_elem_neuron_trace_t *trne,
    nmsim_time_t t,              /* Time at start of step. */
    double X,                    /* Firing indicators between {t} and {t+1}. */
    double I,                    /* External input of each neuron from {t} to {t+1} (mV). */
    double J                     /* Total input of each neuron from {t} to {t+1} (mV). */
  );
  /* If {trne} is not {NULL}, and {t} is in {trne.tLo .. trne.tHi-1},
    stores the given values in the entry of {trne} that corresponds to the step from
    {t}, to {t+1}, namely {trne.ts[t-trne.tLo]}. If {trne} refers to more than one 
    neuron, the given values must be averages over those neurons. */

/* INPUT-OUTPUT */

void nmsim_elem_neuron_trace_write(FILE *wr, nmsim_elem_neuron_trace_t *trne);
  /* Writes the neuron trace {trne} to file {wr}, in the format described by 
    {nmsim_elem_neuron_trace_read_INFO} below.  Does not write any initial and 
    final entries with are completely undefined. */
   
nmsim_elem_neuron_trace_t *nmsim_elem_neuron_trace_read(FILE *rd);
  /* Reads a neuron trace {trne} from file {rd}, in the format described by 
    {nmsim_elem_neuron_trace_read_INFO} below. */
    
#define nmsim_elem_neuron_trace_read_INFO \
  "  The file begins with a line \n" \
  "      \"begin " nmsim_elem_neuron_trace_FILE_TYPE " format of " nmsim_elem_neuron_trace_VERSION "\"\n" \
  "  Then follow a set of lines in the format {NAME} = {VALUE}, where {NAME}" \
  " is \"first_neuron\" \"last_neuron\", " \
  " \"initial_time\", and \"final time\".  Then" \
  " there follows one line for each discrete time {t} in the monitored interval.  Each" \
  " line has the time {t}, followed by the state and evolution fields as described below.\n" \
  "\n" \
  "TRACE ENTRY FORMAT\n" \
  "  " nmsim_elem_neuron_trace_entry_read_INFO ""

#define nmsim_elem_neuron_trace_FILE_TYPE "nmsim_elem_neuron_trace"
    
#define nmsim_elem_neuron_trace_VERSION "2020-12-16"

/* TESTING AND DEBUGGING */

void nmsim_elem_neuron_trace_compare
  ( nmsim_elem_neuron_trace_t *trne_read, 
    nmsim_elem_neuron_trace_t *trne_orig
  );
  /* Compares the elem-level neuron trace {trne_read} read from a file
    with the expected trace {trne_orig}.  Aborts if there are discrepancies, beyond
    roundoff tolerance.
    
    The two traces may have different time ranges, but any entries that are
    missing (outside of the time range) in one trace must be either
    missing or undefined () in the other. */

#endif
