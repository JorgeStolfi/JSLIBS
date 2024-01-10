#ifndef nmsim_elem_neuron_state_H
#define nmsim_elem_neuron_state_H
 
/* Full state of a Galves-Löcherbach neuron. */
/* Last edited on 2019-02-13 18:18:31 by jstolfi */

#define _GNU_SOURCE
#include <stdint.h>
#include <stdio.h>

#include <bool.h>

#include <nmsim_basic.h>

typedef struct nmsim_elem_neuron_state_t
  { 
    double V;               /* Membrane potential at {t}. */
    nmsim_step_count_t age; /* Firing age at {t}(whole time steps since last firing). */
    double M;               /* Potential decay modulator at time {t}. */
    double H;               /* Output synapse strength modulator at time {t}. */
    bool_t X;               /* Firing indicator between {t} and {t+1}. */
    double I;               /* External input between {t} and {t+1}. */
    double J;               /* Total input between {t} and {t+1}. */
  } nmsim_elem_neuron_state_t;
  
void nmsim_elem_neuron_state_write(FILE *wr, nmsim_time_t t, nmsim_elem_neuron_state_t *st);
  /* Writes {t} and {st} to file {wr}, in the format described by 
    {nmsim_elem_neuron_state_read_INFO} below. */
    
void nmsim_elem_neuron_state_read
  ( FILE *rd, 
    nmsim_time_t t,
    nmsim_elem_neuron_state_t *st
  );
  /* Reads a neuron state from {rd}, in the format described by 
    {nmsim_elem_neuron_state_read_INFO} below, and stores it into {*st}.
    Checks whether the time {st.t} read from the file is indeed {t}a. 
    
    The data is assumed to be all in one line. Also reads the final
    end-of-line. Skips lines that are totally blank. Fails if the line
    is malformed or runs into end-of-file. */
    
#define nmsim_elem_neuron_state_read_INFO \
  "Each state is written in a single line of the form\n" \
  "      \"{t} {V} {AGE} {M} {H} {X} {I} {J}\"\n" \
  " where {t} is the discrete sampling time," \
  " {V,AGE,M,H} are the membrane potential, firing age, potential decay" \
  " modulator, and output modulator at time {t}, {X,I,J} are the" \
  " firing indicator, external input, and total input between" \
  " time {t} and {t+1}.\n" \
  "   The discrete time {t} is an integer, the firing indicator" \
  " is a Boolean (0 or 1), and all other fields are floats."

#endif
