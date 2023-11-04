#ifndef neuromat_H
#define neuromat_H

/* NeuroMat tools. */
/* Last edited on 2023-10-21 21:48:21 by stolfi */

#define _GNU_SOURCE
#include <stdio.h>
#include <stdint.h>

#include <bool.h>

void neuromat_eeg_get_channel_names(char *capType, int32_t nv, char *evnames[], int32_t *ne_P, char ***chname_P);
  /* Returns the number and names of EEG channels in the cap with the given {capType}.
  
    Assumes that the channels comprise {ne} electrode voltages, with standard names,
    that depend on {capType}; followed by {nv} event channels named {evname[0..nv-1]}.
  
    Returns in {*ne_P} the number {ne} of electrodes and in {*chname_P} 
    a vectro {chname[0..ne+nv-1]} with the standard names of the {ne}
    electrodes followed by the given event channel names {evnames[0..ne-1]}. 
    
    Currently {capType} must be "R20" (the 20-electrode cap
    used in the early INDC experiments), "R128" (the 128-electrode cap
    used in later INDC experiments), "R129" (the 128-electrodes plus the
    voltage reference electrode "CZ"), or "FN3" (the three virtual electrodes
    used by Fernando Najman in the 2023 "renewal points" paper). In all cases,  (If {nv} is
    zero then {evname} may be null.)
    
    All the strings as well as the vector {chname} are newly allocated by the procedure.
    In particular, the strings {evname[0..nv-1]} are copied to new storage. */

int32_t neuromat_eeg_find_channel_by_name(char *name, int32_t ic_start, int32_t ic_end, char *chname[], bool_t die);
  /* Returns the index {i} in {ic_start..ic_end} such that the string {chname[i]} 
    is equal to the string {name}.  If there is no such {i}, returns {-1} if {die}
    is false, or fails with message if {die} id true. */

/* PULSE IDENTIFICATION */
    
int32_t neuromat_eeg_locate_pulses
  ( int32_t nt,              /* Number of data frames. */
    int32_t nc,              /* Number of data channels. */
    double **val,        /* The EEG dataset ({nt} frames with {nc} channels each). */
    int32_t ict,             /* Index of trigger channel. */
    double vmin,         /* "Low" trigger channel value. */
    double vmax,         /* "High" trigger channel value. */
    int32_t np_max,          /* Max expected number of pulses in file. */
    int32_t it_pulse_ini[],  /* (OUT) Index of first frame in each pulse. */
    int32_t it_pulse_fin[]   /* (OUT) Index of last frame in each pulse. */
  );
  /* Locates the trigger pulses in an EEG
    dataset{val[0..nt-1][0..nc-1]}, which is assumed to have {nt} data
    frames each with {nc} data channels.  Returns the number {np}
    of pulses found.  Fails if there are more than {np_max} pulses.
  
    Assumes that channel {val[0..nt-1][ict]} is a trigger signal
    that can only take values {vmin} or {vmax}. Each trigger pulse is
    defined as a consecutive run of frames where that channel has value
    {vmax}, bounded on both sides by frames where it has value {vmin}.
    The dataset must not begin or end in mid-pulse; that is, the samples
    {val[0][ict]} and {val[nt-1][ict]} must be {vmin}.
  
    Stores into {it_pulse_ini[jp]} and {it_pulse_fin[jp]} the indices of
    the first and last frames of pulse number {jp} (in {0..np-1}).
    These arrays should have at least {np_max} elements. */ 

void neuromat_eeg_report_pulse(FILE *wr, char *pre, int32_t ic, char *name, int32_t it_ini, int32_t it_fin, double fsmp, char *suf);
  /* Reports in {wr} a pulse in channel {ic} spanning frames {it_ini..it_fin}.
    Assumes that the channel is called {name} and that the samplig rate is {fsmp} (Hz).
    The procedure rints only one line and does not print a final newline.
    If {pre} is not NULL, it is printed before everything; if {suf} is not NULL,
    it is printed after everything. */

#endif
