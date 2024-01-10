#ifndef neuromat_eeg_raw_header_H
#define neuromat_eeg_raw_header_H

/* Tools for reading and writing headers of "raw" NeuroMat EEG dataset. */
/* Last edited on 2013-11-12 00:42:07 by stolfilocal */

#define _GNU_SOURCE
#include <stdio.h>
#include <stdint.h>

#include <neuromat_eeg_header.h>

/* RAW EEG FILES */

typedef struct neuromat_eeg_raw_header_t
  { int32_t version; /* Sample data type (2 = {short int}, 4 = {float}). */
    int16_t year;    /* Recording time: year (4-digit). */
    int16_t month;   /* Recording time: month. */
    int16_t day;     /* Recording time: day. */
    int16_t hour;    /* Recording time: hour. */
    int16_t min;     /* Recording time: minute. */
    int16_t sec;     /* Recording time: second. */
    int32_t msec;    /* Recording time: millisecond. */
    int16_t rate;    /* Sampling rate (samples per second). */
    int16_t nc;      /* Number of channels (incl. ref channel). */
    int16_t gain;    /* Board gain (1,2,4,8). */
    int16_t bits;    /* Number of conversion bits. */
    int16_t uVmax;   /* Full-scale range of amplifier in uV. */
    int32_t nt;      /* Number of data frames. */
    int16_t nv;      /* Number of event channels. */
    char **evnames;  /* Names of event channels. */
  } neuromat_eeg_raw_header_t;
  /* Specification for an extra `run phase' marker channel. */
  
neuromat_eeg_raw_header_t *neuromat_eeg_raw_header_read(FILE *rd);
  /* Reads the header from a raw EEG file {rd} and packs it as a {neuromat_eeg_raw_header_t} record. */
  
void neuromat_eeg_raw_header_print(FILE *wr, neuromat_eeg_raw_header_t *hr);
  /* Prints {hr} in readable form to {wr}. */

neuromat_eeg_header_t *neuromat_eeg_raw_header_to_plain_header
  ( neuromat_eeg_raw_header_t *hr,
    char *file,
    int skip,
    int copy,
    int subject,
    int run
  );
  /* Returns a standard NeuroMat dataset header with fields taken from the raw 
    header {hr}, assuming that {skip} frames will be skipped
    and {copy} frames will be copied afterwards.
    Some fields are discarded, others made up. In particular, the 
    channel names are created with {neuromat_eeg_128_channel_names}. */
  
#endif
