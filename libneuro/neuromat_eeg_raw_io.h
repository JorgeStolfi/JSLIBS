#ifndef neuromat_eeg_raw_io_H
#define neuromat_eeg_raw_io_H

/* Reading plain-text NeuroMat EEG signals. */
/* Last edited on 2023-10-21 21:45:40 by stolfi */

#define _GNU_SOURCE
#include <stdio.h>
#include <stdint.h>

#include <neuromat_eeg.h>

int32_t neuromat_eeg_raw_frame_read(FILE *rd, int32_t version, double unit, int32_t nc, char **chname, double val[]);
  /* Reads a data frame from the raw EEG file {rd}, assumed to consist
    of {nc} channel samples, called {chname[0..nc-1]}.  Converts
    the samples to {double} values and saves  them into {val[0..nc-1]}. 
    Each sample is read with {neuromat_eeg_raw_read_raw_sample}.
    
    Returns 0 if EOF occurs right before the first sample, 1 if
    succeeded. Noisily fails if the file is malformed or EOF occurs in
    the middle of the frame. */
    
/* LOW LEVEL RAW I/O */

double neuromat_eeg_raw_read_raw_sample(FILE *rd, int32_t version, double unit, char *name);
  /* Reads a single sample from the raw EEG file {rd}. The sample is
    assumed to be stored in the file in binary, as a {int16_t} if
    {version=2}, or as a IEEE single-precision {float} if {version=4}.
    In the first case, the integer sample is multiplied by {unit}.
    In either case the result is returned as {double}. */
  
int16_t neuromat_eeg_raw_read_int16(FILE *rd, char *name, int32_t vmin, int32_t vmax);
  /* Reads the next 2 bytes from {rd} as a 16-bit signed integer, in big-endian
    order. */
    
int32_t neuromat_eeg_raw_read_int32(FILE *rd, char *name, int32_t vmin, int32_t vmax);
  /* Reads the next 4 bytes from {rd} as a 32-bit signed integer, in big-endian
    order. */
    
float neuromat_eeg_raw_read_float(FILE *rd, char *name, float vmin, float vmax);
  /* Reads the next 4 bytes from {rd} as a single-precision binary IEEE 
    {float}, in big-endian order. */
    
char *neuromat_eeg_raw_read_event_code(FILE *rd, char *name);
  /* Reads the next 4 bytes from {rd} as an event code; returns as a new string,
    zero-terminated, with all blanks removed. */

void neuromat_eeg_raw_check_int32_val(int32_t res, char *name, int32_t vmin, int32_t vmax); 
void neuromat_eeg_raw_check_float_val(float res, char *name, float vmin, float vmax);
  /* These procedures check whether {vmin <= res <= vmax},
    and noisily fail if not. */

#endif
