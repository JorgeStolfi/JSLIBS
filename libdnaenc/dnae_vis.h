#ifndef dnae_vis_H
#define dnae_vis_H

/* Filtered DNA sequences */
/* Last edited on 2015-10-29 00:06:35 by stolfilocal */

#define dnae_vis_H_COPYRIGHT \
  "Copyright © 2014  by the State University of Campinas (UNICAMP)" \
  " and the Federal Fluminense University (UFF)"
  
#define _GNU_SOURCE
#include <stdint.h>

#include <bool.h>  
#include <argparser.h>  
#include <r3.h>

#include <dnae_sample.h>
#include <dnae_datum.h>
#include <dnae_seq.h>

/* This interface provides tools to visualize an encoded DNA sequence ({dnae_seq_t}) as a curve in 3D space. */

#define stringify(x) strngf(x)
#define strngf(x) #x

/* SEQUENCE IO AND INTERPOLATION */

dnae_seq_t dnae_vis_seq_read(char *name, int8_t ek);
  /* Reads a filtered DNA sequence from file {name}.
    The {ek} argument must be zero or negative.
    If negative, the input sequence is subsampled 
    with step {2^ek} times the original step.
    See {dnae_seq_interpolate}. */


/* GEOMETRY CONVERSION */

r3_vec_t dnae_vis_datums_to_points(dnae_seq_t *seq, double magnify, r3_t *pvec);
  /* Converts the three channels of each datum of {seq} to a Cartesian coordinate triple, 
  multiplying its coordinates by {magnify} and adding the vector {pvec} to it. */
    
void dnae_vis_choose_perturbations(double perturb, int ns, r3_t pvec[]);
  /* Sets {pvec[0..ns-1]} to {ns} distinct vectors with length {perturb}. */

double dnae_vis_max_datum_euc_distsq(int ns, dnae_seq_t seq[], int k, int dir);
  /* Returns the maxmum Euclidean distance squared between any two datums with same index {k}
    in sequences {seq[0..ns-1]}, as computed by {dnae_datum_euc_distsq}. If {dir} is {+1}, the index {k} is counted
    from the beginnining of each sequence; if {dir} is {-1}, it is counted
    from the end of each sequence. Returns {+INF} if any of the sequences
    has {k} of fewer datums. */

/* HIGHLIGHTING */

void dnae_vis_determine_visible_segments
  ( int ns, 
    dnae_seq_t seq[], 
    bool_t showMatch, 
    bool_t hideMatch, 
    double maxDist, 
    int ini[], 
    int fin[]
  );
  /* Determines the segments of sequences {seq[0..ns-1]} that are to 
    receive special treatment, as described in {dnae_vis_show_hide_match_option_INFO}.
    
    Namely, left {ni} be the number of datums in {seq[i]}. If {showMatch} is
    true, sets {ini[i]} and {fin[i]} to the indices if the first and last datums
    of the segment of {seq[i]} that mathces some other segment. In that case
    datums {ini[i]..fin[i]} are to be shown in full, the rest effaced. If there
    is no matching segment, sets {ini[i] = n}, {fin[i] = -1}.
    
    If {hideMatch} is true, sets {ini[i]} to the index of the last datum of the
    shared prefix (0 if no such datum) and {fin[i]} to the first datum of the
    shared suffix ({ni-1} if no such datum). In that case datums
    {ini[i]..fin[i]} are to be shown in full, the rest effaced. However, if the
    shard prefix and shared suffix of {seq[i]} overlap, sets {ini[i] = n},
    {fin[i] = -1}.
    
    If neither {showMatch} nor {hideMatch}} are true, sets {ini[i] = 0}, {fin[i]
    = ni-1}, just in case. */

int dnae_vis_find_pref_suff_match(int ns, dnae_seq_t seq[], double maxDist, int dir);
  /* Finds the longest prefix (if {dir} is {+1}) or the longest suffix
    (if {dir} is {-1}) of the sequences {seq[0..ns-1]} in which all
    coresponding datums differ by less than {maxDist}.   The distance is
    defined as {sqrt(dnae_datum_euc_distsq)}.  Returns the length (number of
    datums) in the prefix or suffix. If there are no matching datums at
    the selected end, returns {0}. */

void dnae_vis_find_mid_match(int ns, dnae_seq_t seq[], double maxDist, int minSize, int ini[], int fin[]);
  /* Finds the longest segment in each of the sequences {seq[0..ns-1]} with at least {minSize} datums
    that matches a segment of same length some other sequence.  Two segments are considered to match 
    if all coresponding datums differ by less than {maxDist}. The distance is
    defined as {sqrt(dnae_datum_euc_distsq)}. 
    
    Returns in {ini[k]} and {fin[k]} the first and last datum indices of that longest matching segment. 
    If there are no matching datums in some sequence, returns {ini[k] > fin[k]}.  Note that 
    if there are three or more sequences, the longest matching segment in one need not match the 
    longest matching segment in another.  */

/* COMMAND LINE PARSING */

void dnae_vis_parse_seqFile_options(argparser_t *pp, string_vec_t *seqFile,  int_vec_t *texture);
  /* Parses the "-seqFile" arguments as described in {dnae_vis_seqFile_option_HELP} and 
    {dnae_vis_seqFile_option_INFO} below.  The vectors {seqFile} and {texture} are allocated 
    by the procedure itself, so they should be empty on input. */

#define dnae_vis_seqFile_option_HELP \
  "{ -seqFile {SEQ_FILE} [ -texture {TEXTURE_CODE} ] }.."

#define dnae_vis_seqFile_option_INFO \
  "  -seqFile {SEQ_FILE} {TEXTURE_CODE}\n" \
  "    Each occurrence of this flag specifies an additional" \
  " sequence to be visualized.  If the file name is followed" \
  " by \"-texture\", the sequence will be visualized using" \
  " the given {TEXTURE_CODE} (a non-negative integer); otherwise" \
  " the {TEXTURE_CODE} will be set to sequence's index" \
  " in the sequence list, counting from 0."

void dnae_vis_parse_resample_option(argparser_t *pp, int8_t *resample_expP);
  /* Parses the option "-resample" as described by  {dnae_vis_resample_option_HELP}
    below, and the corresponding {*_INFO} string.
    
    The exponent {k} of the resampling step {R_STEP} (a non-positive integer) is
    returned in {*resample_expP}. */

#define dna_vis_RESAMPLE_STEPS_MAX (128)
  /* Max interpolation steps per original step. 
    Must not exceed {2^msm_seq_desc_estep_MAX}. */

#define dnae_vis_resample_option_HELP \
  "[ -resample {R_STEP} ]"
  
#define dnae_vis_resample_option_INFO \
  "  -resample {R_STEP}\n" \
  "    This optional argument specifies that the sequences should" \
  " be interpolated so that successive datums will be connected" \
  " by a poligonal curve instead of a straight rod.   The number {R_STEP} should" \
  " be 1 (meaning no resampling), or a negative power of 2, namely" \
  " 0.5, 0.25, 0.125, 0.0625, etc..  If this parameter omitted, the" \
  " program assumes {R_STEP = 1} (no interpolation)."  

void dnae_vis_parse_render_options
  ( argparser_t *pp, 
    bool_t *showMatchP,
    bool_t *hideMatchP,
    double *maxDistP, 
    double *magnifyP,
    double *perturbP
  );
  /* Parses the options "-showMatch", "-hideMatch",
    "-magnify", and "-perturb", as described by
    {dnae_vis_show_hide_match_option_HELP}, {dnae_vis_magnify_option_HELP}, 
    and {dnae_vis_perturb_option_HELP}  below, and the corresponding {*_INFO} strings.
    
    The parameter {MAX_DIST} of "-showMatch" and "hideMatch" is returned in
    {*maxDistP}, and the options themselves in {*showMatchP} and {*hideMatchP}.
    If neither option is specified, {*maxDistP} is set to {-1}.
    
    The remaining parameters are returned in {*magnifyP} and {*perturbP}. */

#define dnae_vis_show_hide_match_option_HELP \
  "[ { -showMatch | -hideMatch } {MAX_DIST} ]"

#define dnae_vis_show_hide_match_option_INFO \
  "  -showMatch {MAX_DIST}\n" \
  "  -hideMatch {MAX_DIST}\n" \
  "    These mutually exclusive options have effect only if there are" \
  " two or more sequences and {MAX_DIST} is non-negative.  They specify" \
  " which parts of each sequence are to be shown in 'full' (opaque, full" \
  " color) or 'effaced' (transparent, ghostly whie).  For" \
  " \"-showMatch\", the program searches for the longest segment in" \
  " each sequence that matches a segment of some other sequence, and" \
  " shows that segment in full.  If there is no such segment, the entire" \
  " sequence is showsn effaced.   For \"-hideMatch\", the program finds" \
  " the longest prefix and the longest suffix in each sequence that matches" \
  " the same end segment of all the other sequences, and shows in full" \
  " only the datums that are between the prefix and suffix, or at the boundary.  If the" \
  " matched prefix and suffix of some sequence overlap by more than one" \
  " datum, the entire sequence is shown effaced.\n" \
  "    For both options, two sequence segments are assumed to match iff" \
  " they have the same length, and" \
  " the distance beyween any two corresponding datums in the two segments" \
  " is at most {MAX_DIST}.  The" \
  " distance here is the square root" \
  " of {dnae_datum_euc_distsq}.  If neither option is specified, or if {MAX_DIST} is" \
  " negative, all sequences are shown in full.  The datum distances" \
  " are computed after resampling, but before" \
  " applying the \"-magnify\" and \"-perturb\" options."

#define dna_vis_MAG_MIN (0.01)
#define dna_vis_MAG_MAX (10000.0)
  /* Allowed range for the magnification factor. */

#define dnae_vis_magnify_option_HELP \
  "[ -magnify {MAG_FACTOR} ]"
  
#define dnae_vis_magnify_option_INFO \
  "  -magnify {MAG_FACTOR}\n" \
  "    This optional argument specifies that each sequence model should" \
  " be magnified by {MAG_FACTOR} (keeping the barycenter of the tetrahedron fixed) for better" \
  " visualization.  If this argument omitted, the models are shown" \
  " with their actual size."
    
#define dnae_vis_perturb_DEFAULT (1.0e-5)
  /* Default perturbation magnitude, for {dnae_vis_parse_render_options}. */
 
#define dnae_vis_perturb_option_HELP \
  "[ -perturb {PERTURB_RAD} ]"
  
#define dnae_vis_perturb_option_INFO \
  "  -perturb {PERTURB_RAD}\n" \
  "    This optional argument specifies that each sequence model should" \
  " be displaced by {PERTURB_RAD} in a different direction, so that" \
  " coincindent sections can be distinguished. If this argument" \
  " omitted, {PERTURB_RAD} is set to " stringify(dnae_vis_perturb_DEFAULT) "."
 
#endif
