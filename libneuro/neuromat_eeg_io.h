#ifndef neuromat_eeg_io_H
#define neuromat_eeg_io_H

/* Reading plain-text NeuroMat EEG signals. */
/* Last edited on 2023-12-05 23:34:13 by stolfi */

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
    
void neuromat_eeg_data_write
  ( FILE *wr,
    int32_t nt,
    int32_t nc,
    double **val,
    char *fmt,
    int32_t it_ini,
    int32_t it_fin,
    int32_t it_step
  );
  /* Assumes that {val[0..nt-1][0..nc-1]} is an EEG dataset consisting
    of {nt} frames each with {nc} analog signals. 
    
    Writes to {wr} the samples {val[it][0..nc-1]} where {it} 
    begins with {it_ini}, increases by {it_step} and does not exceed {it_fin}.
    Requires {0 <= it_ini <= it_fin < nt} and {it_step > 0}.
    
    Each frame is written with {neuromat_eeg_frame_write(wr,nc,val[it],fmt)}
    Does not write any header or trailer. */
      
void neuromat_eeg_frame_write(FILE *wr, int32_t nc, double val[], char *fmt);
  /* Assumes that {val[0..nc-1]} is an EEG data frame consisting of {nc}
    analog signal samples. Writes to {wr} the samples {val[0..nc-1]}.
    
    Each value is written with the format {fmt}, which must include
    exactly one '%' float {printf} format spec ('e' or 'f'). If {fmt} is
    null or empty, provides the default "%14.8e". An extra blank is always written
    between values. */
    
void neuromat_eeg_frame_print
  ( FILE *wr,
    char *pre,
    int32_t nc,
    char **chname,
    double val[],
    char *fmt,
    char *sep,
    char *suf
  );
  /* Prints to {wr} the samples {val[0..nc-1]}.  The list is preceded by the 
    string {pre}, closed with the string {suf}, and values are separated by {sep};
    if any of these strings is null, it is omitted. 
    
    If {chname} is not null, each sample {val[i]} will be preceded by "{chname[i]} = " if.
    Each value is written with the format {fmt}, which must include
    exactly one '%' float {printf} format spec ('e' or 'f'). If {fmt} is
    null or empty, provides the default "%14.8e". */
    
#define neuromat_eeg_io_FORMAT_OPT_HELP \
  "[ -format {FMT} ]"

#define neuromat_eeg_io_FORMAT_OPT_INFO(default) \
  " -format {FMT}\n" \
  "    This optional argument specifies the format for the data samples in" \
  " output files. The string {FMT} must include exactly" \
  " one {printf}-sryle '%' float format spec ('e' or 'f'). If omitted or empty," \
  " it defaults to default \"" default "\"."
  
#endif
