#ifndef neuromat_eeg_channel_stats_H
#define neuromat_eeg_channel_stats_H

/* NeuroMat statistics tools. */
/* Last edited on 2021-09-01 17:51:52 by stolfi */

#define _GNU_SOURCE
#include <stdio.h>
#include <stdint.h>

#include <bool.h>

#include <neuromat_eeg.h>

typedef struct neuromat_eeg_channel_stats_t
  { int32_t num;    /* Number of observations (frames). */
    double twt;     /* Total frame weight. */
    double min;     /* Min value. */
    int32_t icmin;  /* Channel index of min value, or {-1} if none. */
    int32_t itmin;  /* Frame index of min value, or {-1} if none. */
    double max;     /* Max value. */
    int32_t icmax;  /* Channel index of max value, or {-1} if none. */
    int32_t itmax;  /* Frame index of max value, or {-1} if none. */
    double avg;     /* Average (mean). */
    double var;     /* Variance around the mean. */
    double dev;     /* Standard deviation {sqrt(var)}. */
    double msq;     /* Average of squared value. */
    double rms;     /* Root mean square {sqrt(msq)}. */
  } neuromat_eeg_channel_stats_t;
  /* Statistics for a single channel in an EEG dataset.
    
    When the statistics are computed, each frame may be given 
    a non-negative weight.  In that case, {num}, {min} and {max}
    refer only to frames with non-zero weight.  The {avg} and {msq}
    quantities are weighted averages. */
    
/* ALLOCATION AND INITIALIZATION */

neuromat_eeg_channel_stats_t *neuromat_eeg_channel_stats_new(int32_t nc);
  /* Allocates a vector of {nc} channel statistics records. 
    Initializes them all with {neuromat_eeg_channel_stats_clear}. */

void neuromat_eeg_channel_stats_clear(neuromat_eeg_channel_stats_t *st);
  /* Sets the contents of {*st} to an undefined value.
    The fields {.num} and {.twt} are set to zero,
    {.min} is set to {+INF}, {.max} is set to {-INF}, the corresponding 
    indices are set to {-1}  and all other fields are set to {NAN}. */

/* GATHERING STATISTICS

  The following procedures assume that {val[0..nt-1][0..nc-1]} 
  contains the values of {nc} channels sampled at {nt} times.
  
  The array {wt}, if not {NULL}, must have {nt} finite non-negative
  elements, assumed to be frame weights.  That is, frame {it}
  is treated as if it was {wt[it]} identical frames.  If {wt}
  is {NULL}, then the weights are assumed to be all 1. 
  
  The {eps} parameter is assumed to be the standard deviation of
  the measurement error.  It will contribute to {min}, {max},
  {avg} and {msq}.
  
  When applicable, the first {ne} channels are assumed to be
  electrode potentials (in µV), and the remaining {nc-ne} are assumed
  to be event or phase markers (in arbitrary units).
  
  If {nt} is zero, or all weights are zero, {min} and {max} will be
  {+INF} and {-INF}, and the other fields will be {NAN}. */
  
void neuromat_eeg_channel_stats_gather
  ( int32_t nt, 
    int32_t nc, 
    double **val,
    double wt[],
    double eps,
    int32_t ic,
    neuromat_eeg_channel_stats_t *st
  );
  /* Computes the statistics of channel {ic}, which must be in {0..nc-1}, and
    stores them in {*st}. */
 
void neuromat_eeg_channel_stats_gather_all
  ( int32_t nt, 
    int32_t nc, 
    double **val,
    double wt[],
    double eps,
    neuromat_eeg_channel_stats_t st[],
    int32_t ne, 
    neuromat_eeg_channel_stats_t *stg
  );
  /* Computes the statistics of the values {val[0..nt-1][ic]}
    of each channel {ic} in {0..nc-1}, stores them in {st[ic]}. 
    
    Then, if {stg} is not {NULL} merges the statistics of electrode
    channels {st[0..ne-1]} in a single statistics record {*stg}. The
    fields {.num}, {.twt}, and {.min} of {stg} will be the minima of the
    corresponding fields in {st[0..ne-1]}. The fields
    {.max,.var,.dev,.msq,.rms} of {stg} will be the maxima of the
    corresponding fields. The field {stg.avg} will always be {NAN}.
    
    If any field {st[0..ne-1].num} is zero, the fields will be
    {+INF}, {-INF}, or {NAN}, as per above. */

void neuromat_eeg_channel_stats_print
  ( FILE *wr, 
    int32_t indent,
    int32_t ic,
    char *name,
    bool_t pnum,
    neuromat_eeg_channel_stats_t *st
  );
  /* Prints to {wr} one line with the statistics {*st} of channel 
    with index {ic}, assumed to be called {*name}.  The line will
    be indented by {indent} columns.  If {pnum} is true, includes
    the fields {st.num} and {st.twt} in the printout. */

void neuromat_eeg_channel_stats_print_all
  ( FILE *wr, 
    int32_t indent,
    int32_t nc, 
    char *name[],
    bool_t pnum,
    neuromat_eeg_channel_stats_t st[],
    int32_t ne,
    neuromat_eeg_channel_stats_t *stg
  );
  /* Prints to {wr} the statistics {st[0..nc-1]} of the {nc}
    channels.  Assumes that channel {ic} is called {name[ic]}.  
    The other arguments are as in {neuromat_eeg_channel_stats_print}.
    
    At the end, if {stg} is not {NULL}, prints to {wr} the merged
    statistics {*stg} assumed to refer to {ne} electrode channels.*/
    
double neuromat_eeg_channel_stats_avg
  ( int32_t nt, 
    int32_t nc, 
    double **val,
    double wt[],
    int32_t ic
  );
  /* Computes the weighted mean of the values of channel {ic}
    in the sample array {val[0..nt-1][0..nc-1]}, where 
    {nt} is the number of data frames and {ic} is the number of channels. 
    
    If {wt} is not {NULL}, it must be a vector with {nt} finite
    non-negative elements.  The value for frame {it} will be
    weighted with weight {wt[it]}.  If the sum of all weights
    is zero, the procedure returns zero. */

/* COVARIANCE MATRIX */    
  
double *neuromat_eeg_channel_stats_covariance_matrix
  ( int32_t nt,          /* Number of frames. */
    int32_t ne,          /* Number of electrodes. */
    double **val,    /* The samples per frame and electrode. */
    double vshift[], /* Values to subtract from each channel. */
    double wt[]      /* Weight of each frame. */
  );
  /* Returns the covariance matrix of the EEG dataset
    {val[0..nt-1][0..ne-1]}, with the electrode potentials shifted down
    by {vshift[0..ne-1]}.
    
    Assumes that {val[it][ie]} is the value of electrode {ie} at time
    {it}, and subtracts from it {vshift[ic]} before computing the
    covariances.  
    
    If {wt} is not {NULL}, it must be a vector with {nt} finite
    non-negative elements. The values for frame {it} will be weighted
    with weight {wt[it]}. If the sum of all weights is zero, all
    covariances will be zero.
    
    The result is a vector with {ne^2} elements, containing a symetric
    matrix with size {ne} by {ne}, linearized by rows. */

void neuromat_eeg_channel_stats_accum_covariance_matrix
  ( int32_t nt,          /* Number of frames. */
    int32_t ne,          /* Number of electrodes. */
    double **val,    /* The samples per frame and electrode. */
    double vshift[], /* Values to subtract from each channel. */
    double wt[],     /* Weight of each frame. */
    double cnew,     /* Coefficient for newly computed matrix. */
    double cold,     /* Coefficient for old matrix. */
    double *Cv       /* (IN/OUT) Covariance matrix. */
  );
  /* Stores in {Cv} the linear combination
    {cold*Cv + cnew*Cn}, where {Cn} is the covariance matrix of {val}.
  
    The parameter {Cv} should be a vector with {ne^2} elements, containing a
    matrix with size {ne} by {ne}, linearized by rows.  The other parameters
    are as in {neuromat_eeg_channel_stats_covariance_matrix}. */

#endif
