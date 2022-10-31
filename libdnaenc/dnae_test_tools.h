#ifndef dnae_test_tools_H
#define dnae_test_tools_H

/* Miscellaneous tools for DNA encoding test programs. */
/* Last edited on 2022-10-31 11:32:53 by stolfi */

#define dnae_test_tools_H_COPYRIGHT \
  "Copyright © 2006  by the State University of Campinas (UNICAMP)"

#define _GNU_SOURCE
#include <stdint.h>

#include <vec.h>

#include <msm_basic.h>
#include <msm_pairing.h>

#include <dnae_seq.h>

#define dnae_test_tools_FONT_SIZE (10.0)
  /* Default plot font size in pt. */

void dnae_test_tools_seq_write_and_plot_named
  ( dnae_seq_t *seq, 
    char *title,
    char *name,
    char *tag,
    bool_t plot,
    double hSize,
    double vSize,
    double fontSize,
    int32_t maxXLabChars
  );
  /* Writes a test sequence to a disk file called "{name}{tag}.txt".
    If {plot} is true, also writes a plot of its sample values to a disk file called
    "{name}{tag}.eps". Assumes that the sequence samples are in the
    vector {*smp}, and that the sequence's descriptor is {*seq}, with
    subsampling factor {den}. 
    
    If [plot} is true, also writes an EPS plot. The parameters {hSize}
    and {vSize} are the total width and height of the figure in mm,
    while the {fontSize} is the nominal font size in points. The
    parameter {maxXLabChars} should be an upper bound on the number of
    characters in the X plot scale labels. If negative, the procedure
    makes a guess.
    
    Assumes that the decoded samples do not exceed 1.0 in absolute value. */

void dnae_test_tools_seq_multi_write_and_plot_named
  ( dnae_seq_t seq[],
    int32_t maxLevel,
    char *title,
    char *name,
    char *tag,
    bool_t plot,
    double hSize,
    double vSize,
    double fontSize,
    int32_t maxXLabChars
  );
  /* Assumes that {seq[i]}, for {i} in {0..maxLevel-1}, is a versions
    of sequence {seq[0]} filtered and downsampled {i} times. Calls
    {dnae_test_tools_seq_write_and_plot_named}
    with arguments {seq[i],titlei,name,tagi,plot,hSize,vSize,fontSize,maxXLabChars}
    on each sequence {seq[i]}; where {titlei} is {title} with " level
    {i}" appended, and {tagi} is {tag} with "-{II}" appended, where
    {II} is the two-digit level {i}.
    
    If [plot} is true, also writes an EPS plot. The parameters {hSize}
    and {vSize} are the total width and height of the figure in mm,
    while the {fontSize} is the nominal font size in points. The
    parameter {maxXLabChars} should be an upper bound on the number of
    characters in the X plot scale labels. If negative, the procedure
    makes a guess, based on the length of {seq[0]}, and uses that value
    for all plots.
    
    Assumes that the decoded samples do not exceed 1.0 in absolute value. */

void dnae_test_tools_make_seq_pair
  ( char *borg, 
    double mutProb, 
    double delProb, 
    dnae_seq_id_t xid, 
    char *xtag,
    dnae_seq_t *xP,
    dnae_seq_id_t yid, 
    char *ytag,
    dnae_seq_t *yP,
    msm_pairing_t **prP,
    char *outName
  );
  /* Generates two test sequences {*xP,*yP} derived from 
    the nuceotide string {borg} by mutations, insertions and deletions.
  
    The sequences will have numeric identifier
    {xid,yid}, name {xtag,ytag}, and appropriate {cmt} strings.
    
    The procedure also returns in {*prP} the list of `true' rungs between
    between the two sequences, i.e. the list of pairs {ix,iy} such that
    {xP[ix]} and {yP[iy]} were both copied from the same letter of 
    {borg}, without mutation.  Beware that {**prP} will be 1-increasing but
    not 1-atomic in general.
    
    The procedure also writes the sequences to files
    "{outName}-{xtag}.bas" and "{outName}-{ytag}.bas".
    
    Each output sequence is generated idependently of the other.
    The precise meaning of the parameters {mutProb} and {delProb} is
    defined by the procedure {dnae_nucleic_string_mutate}. */

void dnae_test_tools_write_generated_sequence(char *outName, char *b, dnae_seq_t *s);
  /* Writes the nucleotide string {b} and its numerically encoded 
    version {s} to disk files named "{outName}-{s.name}.bas" 
    and "{outName}-{s.name}.egs", respectively. */

#endif

