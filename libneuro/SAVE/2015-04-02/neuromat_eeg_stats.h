#ifndef neuromat_stats_H
#define neuromat_stats_H

/* NeuroMat statistics tools. */
/* Last edited on 2013-11-30 04:39:19 by stolfilocal */

#define _GNU_SOURCE
#include <stdio.h>

#include <bool.h>

#include <neuromat_eeg.h>
  
void neuromat_eeg_stats_per_channel
  ( int nt, 
    int nc, 
    double **val,
    bool_t zeroMean,
    double vavg[], 
    double vvar[],
    double vmin[], 
    double vmax[], 
    int ne, 
    double *minMinP, 
    double *maxMaxP, 
    double *maxDevP
  );
  /* Assumes {val[0..nt-1][0..nc-1]} has values of {nc} channels sampled
    at {nt} times. 
    
    For {ic} in {0..nc-1}, if {zeroMean} is false, stores into {vavg[ic]} the average of all the
    measurements {val[0..nt-1][ic]} of channel {ic} (which is {NAN} is {nt} is zero).  
    If {zeroMean} is true, sets {vavg[ic]} to zero.
    
    In either case, for each {ic} stores into {vvar[ic]} the variance of
    the samples of hannel {ic} about the mean {vavg[ic]} (which is {NAN} 
    {nt} is zero).
    
    Also saves the min and max values of each channel {ic} in
    {vmin[ic],vmax[ic]}.
    
    Finally returns in {*minMinP} and {*maxMaxP} the minimum and maximum
    values observed in the first {ne} channels, and in {*maxDevP} the
    maximum of the standard deviations {sqrt(vvar[ie])} for those
    channels. */

void neuromat_eeg_stats_per_channel_print
  ( FILE *wr, 
    int nc, 
    char *name[], 
    double vavg[], 
    double vvar[], 
    double vmin[], 
    double vmax[] 
  );
  /* Writes to {wr} the statstics of {nc} channels. Assumes that, for
    {ic} in {0..nc-1}, the string {name[ic]} is the channel name,
    {vavg[ic]} and {vvar[ic]} are its average and variance, {vmin[ix]}
    and {vmax[ic]} are its extremal values. */
  
double *neuromat_eeg_stats_covariance_matrix(int nt, int ne, double **val, double vavg[]);
  /* Returns the covariance matrix of the EEG dataset {val[0..nt-1][0..ne-1]}.
    Assumes that {val[it][ie]} is the value of electrode {ie} at time {it}.
    The result is a symetric matrix with size {ne} by {ne}. */

#endif
