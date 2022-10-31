#ifndef dnae_seq_multi_H
#define dnae_seq_multi_H

/* Multiscale filtering of DNA sequences */
/* Last edited on 2022-10-31 15:13:26 by stolfi */

#define dnae_seq_multi_H_COPYRIGHT \
  "Copyright © 2005  by the State University of Campinas (UNICAMP)" \
  " and the Federal Fluminense University (UFF)"

#define _GNU_SOURCE
#include <stdio.h>
#include <stdint.h>

#include <bool.h>
#include <vec.h>

#include <dnae_sample.h>
#include <dnae_datum.h>
#include <dnae_seq.h>

void dnae_seq_multi_free_datums(dnae_seq_t seq[], int32_t maxLevel);
  /* Reclaims the internal storage of sequences {seq[0..maxLevel]} 
    (but not the array {*seq} itself). */
   
void dnae_seq_multi_free(dnae_seq_t *seq[], int32_t maxLevel);
  /* Reclaims the sequences {*(seq[0..maxLevel])} and all
    their internal storage. */
    
void dnae_seq_multi_filter
  ( dnae_seq_t *s, 
    int32_t maxLevel, 
    double_vec_t *wt0, 
    char *wname0,
    double_vec_t *wt1, 
    char *wname1,
    int8_t ek0,
    dnae_seq_t sr[]
  );
  /* Places into {sr[0..maxLevel]} various filtered and resampled versions of the
    sequence {s}. Namely {sr[0]} is a copy of {s} (with 
    separately allocated storage), and {sr[r]} is the filtered and possibly resampled
    version of {sr[r-1]}, for {r} in {1..maxLevel}.  
    
    Specifically, the filter table {wt0} is applied to sequence
    {sr[0]}, and the resulting sequence is up- or downsampled with step {2^ek0} to
    obtain sequence {sr[1]}. From then on, each sequence {sr[r+1]} is
    obtained by filtering the previous sequence with filter table
    {wt1}, and then downsampling the result with step 2. The comment
    strings {wname0} and {wname1} are used in the {cmt} field of the
    resulting sequences, as appropriate.
    
    Therefore, sequence {sr[0]} has the same number of datums as {s},
    and thereafter each sequence {sr[r+1]} will have at most
    {(n+1)/2^{ek0+r}} datums (but usually less, to account for the size
    of the filter tables). Therefore, the sequence {sr[r]} may be empty
    for {r} beyond a certain value. */

void dnae_seq_multi_get_2017_paper_weights
  ( double_vec_t *wt0_P, 
    char **wname0_P, 
    double_vec_t *wt1_P, 
    char **wname1_P,
    bool_t verbose
  );
  /* Stores in {*wt0_P} and {*wt1_P} the initial and incremental 
    weight tables {wt0} and {wt1}, suitable for {dnae_seq_multi_filter}, described in 
    the 2017 paper submitted to /Algorithms for Molecular Biology/
    and to /IEEE Transactions on Computational Biology and Bioinformatics/.
    
    Also stores in {*wname0_P} and {*wname1_P} textual descriptions of the
    tables for documentation.
    
    Level 1 of the multiscale pyramid should be the raw
    nuclotides, mapped to corners of the canonical tetrahedorn,
    filtered with {wt0} and no subsampling.  Each subsequent 
    level should be obtained by applying {wt1} to the previous 
    level and subsampling with stride 2. */

#endif
