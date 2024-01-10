#ifndef neuromat_H
#define neuromat_H

/* NeuroMat tools. */
/* Last edited on 2013-12-02 03:12:57 by stolfilocal */

#define _GNU_SOURCE
#include <stdio.h>

#include <bool.h>

void neuromat_eeg_get_channel_names(int ne, int nv, char **evname, int *ncP, char ***chnamesP);
  /* Returns in {*chnamesP} a vector {chnames[0..nc-1]} of {nc} strings with 
    the names of the channels for {ne}-electrode EEG datasets. 
    The number {nc} of channels is returned in {*ncP}.
    
    Currently accepts only {ne=20} (the 20-electrode setup) of {ne=128}
    (the 128-electrode setup). In both cases, assumes that the channels comprise
    {ne} electrode voltages, an extra channel (trigger for {ne=20},
    voltage reference for {ne=128}), and finally {nv} event channels
    named {evname[0..nv-1]}.  (If {nv} is zero then {evname} may be null.)
    
    All the strings as well as the vector are newly allocated by the procedure.
    In particular, the strings {evname[0..nv-1]} are copied to new storage. */

void neuromat_eeg_get_20_channel_names(int *ncP, char ***chnamesP);
  /* Returns in {*chnamesP} a vector {chnames[0..20]} of 21 strings that
    are the names of the raw channels in the 20-electrode experiment.
    The number of channels (129) is returned in {*ncP}.
    
    The channels comprise 20 electrode voltages and one trigger channel.
    Note that the C-language indices ("C" in the table below) range from
    0 to 19 while the Fortran/Matlab indices ("F") range from 1 to 20.
    
      |    C  F  Name    C  F  Name    C  F  Name    C  F  Name  
      |   -- -- -----   -- -- -----   -- -- -----   -- -- -----
      |    0  1  "F7"    5  6 "C3"    10 11 "T6"    15 16 "O2"
      |    1  2  "T3"    6  7 "P3"    11 12 "Fp2"   16 17 "Fz"
      |    2  3  "T5"    7  8 "O1"    12 13 "F4"    17 18 "Cz"
      |    3  4  "Fp1"   8  9 "F8"    13 14 "C4"    18 19 "Pz"
      |    4  5  "F3"    9 10 "T4"    14 15 "P4"    19 20 "Oz"
      
    The trigger channel has C-index 20 (F-index 21) and name "TR". 
    
    All the strings as well as the vector are newly allocated by the procedure.
  */

void neuromat_eeg_get_128_channel_names(int *ncP, char ***chnamesP);
  /* Returns in {*chnamesP} a vector {chnames[0..128]} of 129 strings with 
    the names of the channels for the 128-electrode EEG datasets. 
    The number of channels (129) is returned in {*ncP}.
    
    Specifically, names {chnames[0..127]} are "C1" to "C128" 
    (ie. {chnames[i]} is "C{i+1}"), for {i} in {0..127}.  
    The last channel is assumed to be the reference (ground) electrode,
    so {chnames[128]} is set to "CZ" (name used by the Neuromat team; 
    labeled "VREF" in the electrode diagram?).
    
    All the strings as well as the vector are newly allocated by the procedure. */

int neuromat_eeg_find_channel_by_name(char *name, int ic_start, int ic_end, char *chnames[], bool_t die);
  /* Returns the index {i} in {ic_start..ic_end} such that the string {chnames[i]} 
    is equal to the string {name}.  If there is no such {i}, returns {-1} if {die}
    is false, or fails with message if {die} id true. */

/* PULSE IDENTIFICATION */
    
int neuromat_eeg_locate_pulses
  ( int nt,              /* Number of data frames. */
    int nc,              /* Number of data channels. */
    double **val,        /* The EEG dataset ({nt} frames with {nc} channels each). */
    int ict,             /* Index of trigger channel. */
    double vmin,         /* "Low" trigger channel value. */
    double vmax,         /* "High" trigger channel value. */
    int np_max,          /* Max expected number of pulses in file. */
    int it_pulse_ini[],  /* (OUT) Index of first frame in each pulse. */
    int it_pulse_fin[]   /* (OUT) Index of last frame in each pulse. */
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

void neuromat_eeg_report_pulse(FILE *wr, char *pre, int ic, char *name, int it_ini, int it_fin, double fsmp, char *suf);
  /* Reports in {wr} a pulse in channel {ic} spanning frames {it_ini..it_fin}.
    Assumes that the channel is called {name} and that the samplig rate is {fsmp} (Hz).
    The procedure rints only one line and does not print a final newline.
    If {pre} is not NULL, it is printed before everything; if {suf} is not NULL,
    it is printed after everything. */

#endif
