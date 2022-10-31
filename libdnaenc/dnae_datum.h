#ifndef dnae_datum_H
#define dnae_datum_H

/* Numerical encoding of DNA/RNA bases */
/* Last edited on 2022-10-31 11:22:52 by stolfi */

#define dnae_datum_H_COPYRIGHT \
  "Copyright © 2006  by the State University of Campinas (UNICAMP)"

#define _GNU_SOURCE
#include <stdio.h>
#include <stdint.h>

#include <vec.h>

#include <dnae_sample.h>
  
/* DATUMS */

#define dnae_CHANNELS 3

typedef struct dnae_datum_t { dnae_sample_enc_t c[dnae_CHANNELS]; } dnae_datum_t;
  /* A digitalized signal is a sequence of {dnae_datum_t}
    values, each being a vector with {dnae_CHANNELS} samples. */

#define dnae_datum_ZERO ((dnae_datum_t){0, 0, 0})
  
/* SCALING FACTORS FOR DATUMS */

typedef struct dnae_datum_scale_t { double f[dnae_CHANNELS]; } dnae_datum_scale_t;
  /* Per-channel scale factors for decoded samples of a {dnae_datum_t}. */

#define dnae_NUCLEIC_RAW_SCALE 1.0
  /* Standard deviation of decoded sample values in each channel of a
    {dnae_datum_t} that corresponds to a single nucleotide taken from a
    uniform random distribution. */

dnae_datum_t dnae_datum_mix
  ( double sx, 
    dnae_datum_t *fx, 
    dnae_datum_scale_t *xscale,
    double sy, 
    dnae_datum_t *fy, 
    dnae_datum_scale_t *yscale,
    dnae_datum_scale_t *rscale
  );
  /* Computes the linear combination {sx*fx + sy*fy}, channel by channel.
    The samples are decoded with the scale factors {xscale,yscale}
    and the result is encoded with {rscale}. */
 
/* COMPARING DATUM VALUES LITERALLY */

double dnae_datum_euc_distsq(dnae_datum_t *fx, dnae_datum_scale_t *xscale, dnae_datum_t *fy, dnae_datum_scale_t *yscale);
  /* Returns the Euclidean distance squared between the the datums {*fx} and {*fy},
    decoded using the scales {xscale} and {yscale}, respectively, normalized 
    so that it is 1 between two distinct raw (unfiltered)
    nucleotide vectors, e.g. between {(+1,+1,+1)} ('A') and 
    {(-1,-1,+1)} ('T'), and 0 between any two equal datums,
    raw or not. */
 
/* COMPARING DATUM VALUES AS SURROGATES OF NUCLEOTIDE SEQUENCES */

double dnae_datum_diffsq(dnae_datum_t *fx, dnae_datum_scale_t *xscale, dnae_datum_t *fy, dnae_datum_scale_t *yscale);
  /* Returns an estimate of the average difference squared between the
    raw datums which were averaged into {*fx} and {*fy}. 
    
    The result is a quadratic function of the sample values, normalized
    so that is has extreme values when {fx} and {fy} are raw (unfiltered)
    nucleotide vectors, e.g. {(+1,+1,+1)} ('A') or 
    {(-1,-1,+1)} ('T').  In those cases, the function is 0 if {fx} and {fy}
    are the same, and 1 if they are different.  For other datums, 
    the function ranges between 0 and 1.  Note that if {fx} and {fy} are
    the same but not raw, the distance is not zero. */
    
double dnae_datum_step_diffsq
  ( dnae_datum_t *fx0, 
    dnae_datum_t *fy0, 
    dnae_datum_t *fx1, 
    dnae_datum_t *fy1, 
    dnae_datum_scale_t *xscale, 
    dnae_datum_scale_t *yscale
  );
  /* Average value of {dnae_datum_diffsq(X(t),Y(t))} when {X(t)} and
    {Y(t)} interpolate linearly between {fx0,fy0} and {fx1,fy1}.
    
    Assuming that all arguments lie in the datum simplex, the result
    will range over {[0 _ 1]}. */

double dnae_datum_half_step_diffsq
  ( dnae_datum_t *fx0, 
    dnae_datum_t *fy0, 
    dnae_datum_t *fx1, 
    dnae_datum_t *fy1, 
    dnae_datum_scale_t *xscale, 
    dnae_datum_scale_t *yscale
  );
  /* Integral of {(1-t)*dnae_datum_diffsq(X(t),Y(t))} when {X(t)} and {Y(t)} interpolate
    linearly between {fx0,fy0} and {fx1,fy1}, as {t} ranges over {[0_1]}.
    
    Assuming that all arguments lie in the datum simplex, the result
    will range over {[0.0 _ 0.5]}. 
     
    Informally, the result is the part of {dnae_datum_step_diffsq} that
    comes mostly from the first half of the step. Indeed, the value of
    {dnae_datum_step_diffsq(fx0,fy0,fx1,fy1)} is the sum of the first
    and second half-step values, namely
    {dnae_datum_half_step_diffsq(fx0,fy0,fx1,fy1)} and
    {dnae_datum_half_step_diffsq(fx1,fy1,fx0,fy0)}. */

/* DATUM VECTORS */

vec_typedef(dnae_datum_vec_t,dnae_datum_vec,dnae_datum_t);
/* Vector of {dnae_datum_t}. */

/* DATUM I/O */

void dnae_datum_encoded_write(FILE *wr, dnae_datum_t *d, char *lp, char *sep, char *rp);
  /* Writes the datum {d} (encoded) to file {wr}.
    Samples are preceded by {lp}, separated by {sep}, and followed by {rp}. */

void dnae_datum_decoded_write(FILE *wr, dnae_datum_t *d, dnae_datum_scale_t *dscale, char *lp, char *sep, char *rp);
  /* Writes the datum {d} (decoded with scale {dscale}) to file {wr}.
    Samples are preceded by {lp}, separated by {sep}, and followed by {rp}. */

/* DATUMS FROM DNA/RNA NUCLEOTIDES 

  When used to represent elements of nucleic acid sequences, or
  samples of smoothed versions thereof, the three channels {c[0..2]}
  are {(A-T)+(C-G)}, {(G-C)+(A-T)}, and {(A+T)-(C+G)}; where {A} is
  the local density of 'A' nucleotides at a specified point of the
  sequence, and ditto for {C}, {G} and {T}.
  
  Ordinarily, the datum should be in the convex hull of the four
  /nucleotide datums/ {dnae_datum_from_nucleic_char(b)}
  where {b} is 'A', 'T', 'C', or 'G'. */

void dnae_datum_decoded_from_nucleic_char(char b, int32_t *d);
  /* Stores in {d[0..2]} the numeric representation of the DNA/RNA
    nucleotide character {b}, namely
    
      A    (+1,+1,+1)
      T,U  (-1,-1,+1)
      C    (+1,-1,-1)
      G    (-1,+1,-1)
      
    The result is {(0,0,0)} if {b} is not in [ATCGUatcgu]. */

dnae_datum_t dnae_datum_encoded_from_nucleic_char(char b);
  /* Converts the DNA/RNA nucleotide character {b} to a numeric
    vector, as decribed under {dnae_datum_decoded_from_nucleic_char},
    and encodes it as a {dnae_datum_t}.
    Assumes the scale factor {dnae_NUCLEIC_RAW_SCALE}. */

dnae_datum_vec_t dnae_datum_vec_from_nucleic_string(char *s);
  /* Converts a DNA/RNA nucleotide sequence {c[0..n-1]} to a 
    vector of {dnae_datum_t}s, as explained under
    {dnae_datum_from_nucleic_char}. 
    Assumes the scale factor {dnae_NUCLEIC_RAW_SCALE}. */

void dnae_datum_to_nucleic_densities(dnae_datum_t *d, dnae_datum_scale_t *dscale, double *A, double *T, double *C, double *G);
  /* Returns in {*A,*T,*C,*G} the local densities of the four nuclotides that 
    are implied by the datum {d}. Ordinarily, those four numbers shoudl add to 1.
    The samples of {d} are decoded with the scale factor {dscale}. */

#endif

