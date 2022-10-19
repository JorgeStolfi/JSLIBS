#ifndef neuromat_H
#define neuromat_H

/* NeuroMat tools. */
/* Last edited on 2021-08-28 02:04:04 by stolfi */

#define _GNU_SOURCE
#include <stdio.h>
#include <stdint.h>

#include <bool.h>

char **neuromat_eeg_get_channel_names(int32_t ne, int32_t nv, char *evname[]);
  /* Returns a vector {chname[0..nc-1]} with the names of the 
    {nc=ne+nv} channels for {ne}-electrode EEG datasets. 
    
    Currently accepts only {ne=20} (the 20-electrode setup), {ne=128}
    (the 128-electrode setup), or {ne=129} (the 128-electrodes plus the
    voltage reference electrode "CZ"). In all cases, assumes that the
    channels comprise {ne} electrode voltages, with standard names,
    followed by {nv} event channels named {evname[0..nv-1]}. (If {nv} is
    zero then {evname} may be null.)
    
    All the strings as well as the vector are newly allocated by the procedure.
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
