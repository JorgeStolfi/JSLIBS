#ifndef neuromat_eeg_io_H
#define neuromat_eeg_io_H

/* Reading plain-text NeuroMat EEG signals. */
/* Last edited on 2023-10-21 21:46:05 by stolfi */

#define _GNU_SOURCE
#include <stdio.h>
#include <stdint.h>

#include <neuromat_eeg.h>

/* !!! Change the sample type to {float} in order to save memory space? !!! */

int32_t neuromat_eeg_frame_read(FILE *rd, int32_t nc, double frm[], int32_t *nlP, int32_t *nfP);
  /* Read from {rd} the next data frame from an EEG dataset file.
  
    Assumes that each data frame is a line of {rd}, containing {nc}
    values in "%f" of "%e" format, which are the simultaneous readings
    of {nc} analog channels (electrodes, triggers, etc.).
    If {frm} is not null, the values are stored into {frm[0..nc-1]}.
    If {frm} is null, the data frame is read but discarded.
    
    Assumes that the first data frame may be preceded by header or comment
    lines, which are identified by their first character being '#' or an
    ascii letter. 
    
    The result is 1 if the procedure succeeded in reading one more data frame,
    or 0 if it ran into end-of-file (possibly after skipping header and
    comment lines) before the next data frame.
    
    If {nlP} is not null, the procedure assumes that {*nlP} is the
    number of lines already read from {rd} (including header and comment
    lines), and increments {*nlP} with the number of lines actually read
    (including any '#'-comments and header lines).
    
    Similarly if {nfP} is not null, the procedure assumes that {*nfP} is
    the number of data frames already read from {rd} (excluding header
    and comment lines), and increments {*nfP} by the number of data
    frames that it reads (0 or 1). */

double **neuromat_eeg_data_read(FILE *rd, int32_t nskip, int32_t nread, int32_t nc, /*(I/O):*/ int32_t *nlP, /*(OUT):*/ int32_t *ntP);
  /* Reads from {rd} multiple data frames of an EEG dataset. 
    
    The procedure will read but discard the first {nskip} data frames in
    {rd}. Then, if {nread} is positive, it will read at most {nread}
    data frames after that. If {nread} is zero or negative, it will read
    all remaining data frames to the end of the file.
    
    In any case, the procedure returns an array {val[0..nt-1][0..nc-1]}
    where {nt} is the number of data frames actually read and
    {val[it][i]} is the reading of channel {i} at line {skip+it}.
    
    Besides the first {nskip} data frames, the procedure will skip also
    any lines at the beginning of the file that begin with "#" or with
    an ascii letter (assumed to be header lines). lines, which are
    ignored).
    
    If {ntP} is not null, the procedure stores into it the number {nt}
    of data frames actually read (excluding all '#'-comments and header
    lines). */
    
void neuromat_eeg_data_write(FILE *wr, int32_t nt, int32_t nc, double **val, int32_t it_ini, int32_t it_fin, int32_t it_step);
  /* Assumes that {val[0..nt-1][0..nc-1]} is an EEG dataset consisting
    of {nc} analog signals sampled at {nt} equally spaced moments. 
    Writes to {wr} the samples {val[it][0..nc-1]} where {it} 
    begins with {it_ini}, increases by {it_step} and does not exceed {it_fin}.
    Requires {0 <= it_ini <= it_fin < nt} and {it_step > 0}.
    
    The values are written in "%e" format with 8 significant digits, separated by 
    blanks.  Does not write any header or trailer. */
      
void neuromat_eeg_frame_write(FILE *wr, int32_t nc, double val[]);
  /* Assumes that {val[0..nc-1]} is an EEG data frame consisting of
    {nc} analog signal samples. Writes to {wr} the samples
    {val[0..nc-1]}. The values are written in "%e" format with 8
    significant digits, separated by blanks. */
    
void neuromat_eeg_frame_print(FILE *wr, char *pre, int32_t nc, char **chname, double val[], char *sep, char *suf);
  /* Prints to {wr} the samples {val[0..nc-1]}.  The list is preceded by the 
    string {pre}, closed with the string {suf}, and values are separated by {sep};
    if any of these strings is null, it is omitted. 
    Each sample {val[i]} ispreceded by "{chname[i]} = " if {chname} is not null.  */


#endif
