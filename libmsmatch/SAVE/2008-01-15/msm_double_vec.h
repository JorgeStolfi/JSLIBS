#ifndef msm_double_vec_H
#define msm_double_vec_H

/* Tools to manipulate vectors of doubles. */
/* Last edited on 2007-12-21 14:41:28 by stolfi */

#define msm_double_vec_H_COPYRIGHT \
  "Copyright © 2007  by the State University of Campinas (UNICAMP)"

#include <msm_basic.h>
#include <msm_rung.h>

#include <vec.h>

double *msm_double_vec_get_address(double_vec_t *smp, int i, bool_t circ);
  /* Returns the address of element {i} in the vector {smp}.
    If {circ} is TRUE, the index {i} will be reduced modulo
    {smp.nel}.  Otherwise {i} must be in {0..smp.nel-1}. */

double msm_double_vec_interpolate(double_vec_t *smp, double z, bool_t circ);
  /* Returns the value of sequence {smp} at the fractional 
    position {z}. 
    
    If {z} is integer, the result is the element of {smp} with index
    {z}. If {z} is non-integer, the result is obtained by
    interpolating the elements with indices {floor(z)} and {ceil(z)}.
    
    In either case, if {circ} is TRUE, the element indices may have
    any value, and will be reduced modulo {smp.nel}. Otherwise they
    must be in {0..smp.nel-1}. */

double_vec_t msm_double_vec_throw_normal(int ns);
 /* Returns a vector with {ns} random numbers taken from 
   a Gaussian distribution with mean 0 and deviation 1. */

void msm_double_vec_smooth(double_vec_t *smp);
  /* Replaces the sequence {smp} by its cumulative sum,
    then adds to it a linear ramp so as to make it
    periodic. */ 

void msm_double_vec_normalize(double_vec_t *smp);
  /* Applies an affine map to every sample of {smp}, 
    so that they have zero mean and unit standard deviation. */ 

void msm_double_vec_mutate
  ( double_vec_t  *smpo, 
    bool_t circ,
    double mutDev, 
    double delProb,
    double_vec_t *smpd,
    msm_rung_vec_t *gv
  );
  /* Generates a real sequence {smpd} by mutating the sequence
    {smpo}, with noise deviation {mutDev} and deletion/insertion
    probability {delProb}.
    
    The input and output sequences are assumed to be circular
    iff {circ} is TRUE.
  
    The new sequence {smpd} is built one sample at a time. With
    probability {delProb/2}, the next sample of {smpo} is deleted.
    With probabilty {delProb/2}, the average of the previous and next
    sample is inserted before the the next sample. With the remaining
    probability {1-delProb}, the next sample gets copied. In any case,
    each output sample is mixed with a Gaussian-distributed random
    sample, with weights {sqrt(1-mutDev^2)} and {mutDev},
    respectively.
    
    The procedure also returns in {*gv} the list of `true' rungs
    between the two sequences, i.e. the list of pairs {ix,iy} such
    that {smpd[ix]} was copied from {smpo[iy]} (rather than being an
    insertion). Beware that {gv} will be strictly increasing but not
    atomic in general. */

void msm_double_vec_write(FILE *wr, double_vec_t *smp);
  /* Writes to {wr} the real vector {smp}, one element per line. */

void msm_double_vec_write_named(double_vec_t *smp, char *name, char *tag);
  /* Writes the real vector {smp} to a disk file named
    "{name}{tag}.seq", one element per line. */

double_vec_t msm_double_vec_read(FILE *rd);
  /* Reads from {rd} a real vector, one element per line. */ 

double_vec_t msm_double_vec_read_named(char *name, char *tag);
  /* Reads a real vector, one element per line, from a disk 
    file named "{name}{tag}.seq". */ 

#endif
