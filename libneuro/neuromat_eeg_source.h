#ifndef neuromat_eeg_source_H
#define neuromat_eeg_source_H

/* Tools for reading and writing the data-source lines of NeuroMat EEG dataset headers. */
/* Last edited on 2013-06-05 21:47:02 by stolfilocal */

#define _GNU_SOURCE
#include <stdio.h>

#include <neuromat_eeg_header.h>

/* The type {neuromat_eeg_source_t} is defined in {neuromat_eeg_header.h}. */

neuromat_eeg_source_t *neuromat_eeg_source_new(void);
  /* Allocates a new {neuromat_eeg_source_t} record and all its sub-records,
    if any.  All fields are set to their `null' values. */

void neuromat_eeg_source_free(neuromat_eeg_source_t *ho);
  /* Reclaims the header record {*ho} and all its sub-records, if any.
    Any strings and vector field values are assumed to be in `owned' storage
    and freed too. */

void neuromat_eeg_source_write(FILE *wr, char *pref, neuromat_eeg_source_t *ho);
  /* Writes to {wr} zero or more header lines with 
    information from {ho}.  Field names are prefixed with {pref}, Fields that have `null' values
    ({INT_MIN} for integers, {NAN} for [double}s, {NULL} for strings
    and lists) are omitted. */

void neuromat_eeg_source_read_field(FILE *rd, char *name, neuromat_eeg_source_t *ho);
  /* Reads from {rd} the value of the header field of {ho} with the given {name}, 
    skipping any blanks. Also checks for the value read being in its valid range. */

#endif
