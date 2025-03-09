#ifndef pst_map_compare_H
#define pst_map_compare_H

/* pst_map_compare.h -- procedures for comparing maps (height, slopes, etc). */
/* Last edited on 2025-03-04 18:22:40 by stolfi */

#include <stdint.h>

#include <bool.h>
#include <float_image.h>

float_image_t *pst_map_compare
  ( float_image_t *A,
    float_image_t *B,
    int32_t NC,
    double avgA[],
    double devA[],
    double avgB[],
    double devB[],
    uint32_t undef[],
    double avgE[],
    double devE[],
    double devRelE[]
  );
  /* The two images {A,B} must have the same size, and the 
    channel count of each must be either {NC} or {NC-1}.
    
    Returns a map {E} of the same size with {NC} channels, where, for
    each pixel indices {X,Y} and each channel {c} other than {wch=NC-1},
    {E[c,X,Y]} is the difference {A[c,X,Y]-B[c,X,Y]}.
    
    If {wch} is a valid channel index for {A} and/or {B}, then {A[wch,X,Y]}
    and/or {B[wch,X,Y]} are assumed to be reliability weights for
    the other channels of {A} and {B}.  Otherwise the reliability weight 
    is assumed to be 1.0 at all pixels.  In any case, sample {E[wch,X,Y]} 
    is set to the harmonic sum of {A[wch,X,Y]} and {B[wch,X,Y]}.
    
    Also, for each channel {c} of {A} and {B}, INCLUDING channel {wch},
    returns in {avgA[c],avgB[c],avgE[c]} the weighted means of of the
    samples {A[c,X,Y]}, {B[c,X,Y]}, and {E[c,X,Y]}; in
    {devA[c],devB[c].devE[c]} the standard deviations of those samples;
    and in {devRelE[c]} the relative error {devE[c]/devAB[c]} where
    {devAB[c] = hypot(devA[c],devB[c])/sqrt(2)}. 
    
    In the above means and deviations, each sample {A[c,X,Y]},
    {B[c,X,Y]}, or {E[c,x,y]} gets the weight {E[wch,X,Y]}.
    
    The averages and deviations consider only the samples where where
    all three maps {A,B,E} are finite and have nonzero weight. Returns
    in {undef[c]} the number of samples that were excluded from the
    averages and deviations of each channel {c}. */
      
void pst_compare_check_sizes
  ( float_image_t *A,
    int32_t *NC_A,
    float_image_t *B,
    int32_t *NC_B,
    int32_t NC,
    int32_t *NX,
    int32_t *NY
  );
  /* Checks that maps {A} and {B} (which must not be {NULL}) have the
    same col and row counts, and that the number of channels in each map
    is either {NC} or {NC-1}. Returns in {*NC_A} and {*NC_B} the channel
    counts and in {*NX} and {*NY} the common col and row counts. */

/* REPORTING */

void pst_map_compare_analyze_and_write
  ( float_image_t *A,
    float_image_t *B, 
    int32_t NC,
    double change[],
    char *filePrefix,
    int32_t level,
    int32_t iter,
    char *tagA,
    char *tagB,
    char *tagE,
    bool_t verbose
  );
  /* A procedure for comparing two maps, writing out the differences
    and statistics.  Useful, for example, monitoring the state
    before, during, and after an iterative process, especially multiscale.
    
    The map {A} must be non-null and have either {NC} or {NC-1} channels.
    
    If {B} is non-null, it must have the
    same size as {A} and also either {NC} or {NC-1} channels. 
    In this case, the procedure
    calls {pst_map_compare} (q. v.) to computes an image {E=A-B}, and
    the statistical parameters
    {avgA,devA,avgB,devB,undefE,avgE,devE,devRelE} for each channel.
    
    If {wch=NC-1} is a valid channel index of {A} and/or {B}, it is assumed to
    be the index of the reliability weight channel (see
    {pst_map_compare}). 
    
    The maps {A}, {B}, and {E}, if they exist, are written to disk files
    "{filePrefix}-{level:02d}-{iter:09d}-{tag}.fni" where {tag} is
    {tagA}, {tagB}, or {tagE}, respectively. The {E} image is shifted so
    as to have zero mean before being written out.
    
    If {B} is not {NULL}, the procedure also writes a file called
    \"{filePrefix}-{level}-{iter}-{tagE}.txt\" See
    {pst_map_compare_error_summary_INFO} for the format.  
    
    Each of these output files is NOt created if the corresponding {tag}
    is {NULL}.
    
    When {iter} is zero, the procedure assumes that {A} and {B} are the
    intial state, before the first iteration. writing the state before
    the first iteration; in that case, it replaces the "{iter:09d}" part
    of the file names by "beg".
    
    If {iter} is {-1}, the pricedure assumes that {A} and {B} are the
    final result, after the last iteration; in that case, the {iter:09d}
    part of the file name is replaced by "end".
    
    If {level} is negative, the "-{level}" part is omitted from the file
    names. If {level} is non-negative, all messages to {stderr} are indented by
    {2*level+2} spaces. */
    
#define pst_map_compare_error_summary_INFO(A,B) \
  "The error summary file has a single line with the format\n" \
  "\n" \
  "    {level} {NX} {NY} {iter} {NC} {0} {STATS[0]} .. {NC-1} {STATS[NC-1]}\n" \
  "\n" \
  "  where  {level} is the recursion (hierarchy) level, {NX,NY} is the image size, {iter} is" \
  " the number of iterations performed at that level, and {NC} is the number of" \
  " channels (including the weight channel).  Each {STATS[c]} is a summary" \
  " for channel {c}, consisting of the fields\n" \
  "\n" \
  "   {change[c]} {avg" A "[c]} {dev" A "[c]} {avg" B "[c]} {dev" B "[c]} {avgE[c]} {devE[c]} {devRelE[c]}\n" \
  "\n" \
  "  where {change[c]} is the change in the last iteration, {avg" A ",avg" B ",avgE} are" \
  " the averages of {" A "}, {" B "}, and {E}, {dev" A ",dev" B ",devE}  are the" \
  " deviations, and {devRelE} is {devE} relative to {dev" A "} and {dev" B "}."

#endif
