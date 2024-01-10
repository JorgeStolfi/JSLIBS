#ifndef neuromat_H
#define neuromat_H

/* NeuroMat tools. */
/* Last edited on 2013-05-30 23:58:02 by stolfilocal */

#define _GNU_SOURCE
#include <stdio.h>

char **neuromat_e20_signal_names(void);
  /* Returns a vector of 21 strings that are the names of the data channels in
    the 20-electrode experiment, comprising 20 electrodes and one trigger
    channel.  Note that the C-language indices ("C" in the table below) range
    from 0 to 19 while the Fortran/Matlab indices ("F") range from 1 to 20.
    
      |    C  F  Name    C  F  Name    C  F  Name    C  F  Name  
      |   -- -- -----   -- -- -----   -- -- -----   -- -- -----
      |    0  1  "F7"    5  6 "C3"    10 11 "T6"    15 16 "O2"
      |    1  2  "T3"    6  7 "P3"    11 12 "Fp2"   16 17 "Fz"
      |    2  3  "T5"    7  8 "O1"    12 13 "F4"    17 18 "Cz"
      |    3  4  "Fp1"   8  9 "F8"    13 14 "C4"    18 19 "Pz"
      |    4  5  "F3"    9 10 "T4"    14 15 "P4"    19 20 "Oz"
      
    The trigger channel has index 20 and name "TR". 
    
    The vector of names is newly allocated by the procedure, but the 
    strings are read-only and should not be deallocated.
  */

/* INPUT/OUTPUT */

double **neuromat_read_eeg_signals(FILE *rd, int nskip, int nt, int nc);
  /* Reads an EEG dataset from {rd}. Assumes that each line contains
    {nc} values in "%f" of "%e" format, comprising the simultaneous
    readings of {nc} analog channels (electrodes, triggers, etc.).
    Skips {nskip} lines and reads {nt} lines after that.
    Returns an array {val[0..nt-1][0..nc-1]} where {val[it][i]} is the
    reading of channel {i} at line {skip+it}. 
    
    Skips any lines at the beginning of the file that begin with "#"
    or with an ascii letter (assumed to be header lines). */
    
void neuromat_write_eeg_signals(FILE *wr, int nt, int nc, double **val, int it_ini, int it_fin, int it_step);
  /* Assumes that {val[0..nt-1][0..nc-1]} is an EEG dataset consisting
    of {nc} analog signals sampled at {nt} equally spaced moments. 
    Writes to {wr} the samples {val[it][0..nc-1]} where {it} 
    begins with {it_ini}, increases by {it_step} and does not exceed {it_fin}.
    Requires {0 <= it_ini <= it_fin < nt} and {it-step > 0}.
    
    The values are written in "%e" format with 8 significant digits, separated by 
    blanks. */
      
#endif
