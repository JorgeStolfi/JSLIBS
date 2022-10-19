#ifndef neuromat_eeg_header_H
#define neuromat_eeg_header_H

/* Tools for reading and writing headers of plain-text NeuroMat EEG datasets. */
/* Last edited on 2021-08-26 12:32:15 by stolfi */

#define _GNU_SOURCE
#include <stdio.h>
#include <stdint.h>
  
/* !!! Add to plain header the date and other fields from raw header? !!! */

/* HEADERS PLAIN TEXT EEG FILES */

/* Plain-text EEG dataset files (and some files derived from them, such as power
  spectra) may have a header containing several lines of the form
  "{name} = {value}".
  
  In any such line, {name} is an identifier that begins with an ascii letter
  and continues with ascii letters, decimal digits, periods '.' or
  underscores '_', without embedded blanks. The {value} is one or more
  tokens separated by blanks. Multiple consecutive blanks are equivalent
  to one blank, the blanks around the "=" are required, and there must
  be no blanks before the name. These lines can be distinguished from
  data lines by starting with a letter and/or by containing a '='
  character.
  
  The {name} should be the name of one field of {neuromat_eeg_header_t}
  structure below, and the value must be compatible with it. Any fields
  of that structure that are needed by an application and are not
  specified in the file header must be provided by some other means.
  
  Depending on the contet, parameter values in {neuromat_eeg_header_t} and its sub-records
  may be unspecified or unknown. This is indicated by the value {INT32_MIN}
  for integer values, {NAN} for doubles, or {NULL} for lists.  */

typedef struct neuromat_eeg_header_t
  { 
    /* Dataset size and content: */
    int32_t nt;            /* Number of sampling times in the EEG dataset. Positive. */
    int32_t nc;            /* Number of data channels (tokens per line) in EEG the dataset. Positive. */

    /* Nature of data: */
    double fsmp;       /* Sampling frequency (Hz) of dataset. */
    int32_t ne;            /* Number of electric measurement channels, in {0..nc}. */
    char *type;        /* Type of experiment, e.g. "B" (bio) or "N" (non-bio). */
    char **chnames;    /* Names of channels, indexed {0..nc-1}. */
    int32_t kfmax;         /* Maximum frequency index (if spectral data), in {0..nt/2}. */
    char *component;   /* Name of a component of original signal, or NULL. */
    
    /* Filtering data: */
    int32_t tdeg;      /* Degree of trend polynomial, or negative if none. */
    int32_t tkeep;     /* Tells whether trend polynomial was kept (1) of suppressed (0). */
    double flo0;       /* Lower cutoff frequency, {-INF}. */
    double flo1;       /* Nominal lowest preserved frequency, {-INF}. */
    double fhi1;       /* Nominal highest preserved frequency, {+INF}. */
    double fhi0;       /* Higher cutoff frequency, {+INF}. */
    int32_t finvert;   /* Filter complementation: 0 for bandpass, 1 for bandkill. */
    double *rebase_wt; /* Electrode weights used for rebasing, or {NULL}. */
    
    /* Data history: */
    struct neuromat_eeg_source_t *orig; /* Source EEG data in original raw recording. */
    
  } neuromat_eeg_header_t;
  /* Notes:
    
    o If the file contains time-domain data, {nt} is the number of time
      frames (lines), and {fsmp} is the number of time frames per
      second. In that case the index {kfmax} is irrelevant.
    
    o If the file contains spectral data, the number of data lines is
      {kfmax+1}, and line with index {kf} (from 0) refers to the the
      Fourier/Hartley components with absolute frequency {fsmp*kf/nt}
      (Hz). Even in this case, {nt} must be specified and must be at least
      {2*kfmax}.
  
    o The parameter {ne} must be at most {nc}. Channels {0..ne-1} in the dataset
      are assumed to be electrode potentials or data derived from them (e.g.
      PCA components, power spectra). The remaining {nc-ne} channels, if any,
      are assumed to be trigger signals or other non-physical data. 
  
    o The {chnames} field is specified in the file by a line
    
        "channel = {chnames[0]} {chnames[1]} ...  {chnames[nc-1]}"
    
      Each channel name in {chnames[0..ne-1]} must be an identifier,
      beginnng with ascii letter and continuing with decimal digits,
      ascii letters, periods or underscores, with no embedded blanks.
      The pointer {chnames} may be {NULL}.
      
    o If the field {rebase_wt} is not {NULL}, it means that the
      electrode potentials were changed to a different reference
      potential, after filtering and trend replacement. Then {rebase_wt}
      should be be a vector of {ne} non-negative weights that were used
      to compute the average potential of each frame, which was
      subtracted from every sample in that frame. If {NULL}, the data is
      assumed to use the original potential reference, such as the "CZ"
      electrode.
  */
  
typedef struct neuromat_eeg_source_t
  { char *file;        /* Name of original raw data file. */
    int32_t nt;        /* Number of frames in original data file. */
    int32_t it_ini;    /* Index of first time frame in original data file (from 0). */
    int32_t it_fin;    /* Index of last time frame in original data file (from 0). */
    double fsmp;       /* Sampling frequency (Hz) of original dataset. */
    int32_t subject;   /* Index of subject (from 1). */
    int32_t run;       /* Index of experimental run (from 1). */
  } neuromat_eeg_source_t;
  /* Describes the raw EEG data this file was ultimately derived from.
    In the header lines, these field names are prefixed with "orig.",
    except that {it_ini,it_fin} are given in a single line: 
    
    "orig.index_range = {it_ini} {it_fin}"
  */
  
neuromat_eeg_header_t *neuromat_eeg_header_new(void);
  /* Allocates a new {neuromat_eeg_header_t} record and all its sub-records
    (such as {orig}).  All fields are set to their
    `null' values. */
   
void neuromat_eeg_header_free(neuromat_eeg_header_t *h);
  /* Reclaims the header record {*h} and all its sub-records
    (such as {h->orig}).  Any strings and vector values
    are assumed to be in `owned' storage and freed too. */

neuromat_eeg_header_t *neuromat_eeg_header_read(FILE *rd, int32_t neDef, double fsmpDef, int32_t *nlP);
  /* Reads zero or more lines from {rd} and parses them as header
    data. Always returns a non-{NULL} pointer {h} with a non-null
    field {h->orig}, both newly allocated, even if no header lines
    are read.  Any strings and vectors are newly allocated, too.
    
    There must be no blanks before the field name in header lines, and
    every header line must end with newline. The procedure stops reading
    from {rd} right before the first line that does not begin with an
    ascii letter [a-zA-Z] (or at EOF).
   
    The procedure ignores any lines that begin with '#'. These lines too
    must end with a newline. Trailing '#'-comments in header lines are
    NOT allowed.
    
    If {nlP} is not null, the procedure assumes that {*nlP} is the
    number of lines already read from {rd}, and increments it with the
    number of header lines and '#'-comment lines read.
    
    If the electrode count {ne} is not specified in the file,
    assumes {h->nc-1} if {h->nc} is specified, else assumes {neDef}.
    If the number of channels is not specified, assumes {h->ne+1}.
    If the sampling frequency is not specified,
    assumes {fsmpDef}. */
   
void neuromat_eeg_header_write(FILE *wr, neuromat_eeg_header_t *h);
  /* Writes to {wr} zero or more header lines with 
    information from {h}.  Fields that have `null' values
    ({INT32_MIN} for integers, {NAN} for [double}s, {NULL} for strings
    and lists) are omitted. */
    
void neuromat_eeg_header_merge(neuromat_eeg_header_t *dst, neuromat_eeg_header_t *src);
void neuromat_eeg_header_merge_orig(neuromat_eeg_source_t *dst, neuromat_eeg_source_t *src);
  /* Compares all corresponding fields of {dst} and {src} and their
    sub-records (such as {src->orig} and {dst->orig}).  

    If a field is defined in both records, the two values must be equal,
    otherwise the procedure fails.

    If a field is undefined in {dst} but defined in {src}, the procedure
    copies its value to {dst}. In that case, string and vector values
    are copied into newly allocated storage. */
    
int32_t neuromat_eeg_header_append_electrode_channel(neuromat_eeg_header_t *h, char *name);
int32_t neuromat_eeg_header_append_marker_channel(neuromat_eeg_header_t *h, char *name);
  /* These procedures insert a new channel with the given {name} in the header,
    at the end of the electrode channels and at the end of all channels, respectively.
    The {h.chnames} table is reallocated and the counts {h.ne,h.nc} are updated,
    as appropriate.  They return the index of the inserted channel, in {0..h.nc-1}
    for the new {h.nc}. */
  
/* LOW-LEVEL PROCEDURES */

void neuromat_eeg_header_read_field_value(FILE *rd, char *name, neuromat_eeg_header_t *h);
  /* Reads from {rd} the value of the header field of {h} with the given {name}, skipping any blanks. 
    Also checks for the for the value read being in valid range. */

void neuromat_eeg_header_write_field_string(FILE *wr, char *pref, char *name, char *value);
void neuromat_eeg_header_write_field_double(FILE *wr, char *pref, char *name, double value, double vmin, double vmax);
void neuromat_eeg_header_write_field_int(FILE *wr, char *pref, char *name, int32_t value, int32_t vmin, int32_t vmax);
void neuromat_eeg_header_write_field_int_range(FILE *wr, char *pref, char *name, int32_t vini, int32_t vfin, int32_t vmin, int32_t vmax);
  /* These procedures write the generic header field as a line
    
      "{pref}{name} = {value(s)}"
      
    They also check that the value is in the valid range {[vmin _
    vmax]}. For integer ranges, the two values {vini} and {vfin} are
    written on the same lines after the "=", and must satisfy
    {vmin<=vini<=vfin<=vmax}. String values must not be empty and must
    not contain any embedded blanks; they are written without any
    quotes, and. String vector elements are written all on the same
    line, without quotes, separated by blanks. */

#endif
