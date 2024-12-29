/* Histogram of samples from a {float_image_t}.  */
/* Last edited on 2024-12-24 15:52:45 by stolfi */

#ifndef float_image_hist_H
#define float_image_hist_H

#include <stdio.h> 

#include <vec.h> 
#include <float_image.h>

void float_image_hist_build
  ( float_image_t *A,
    int32_t c,
    uint32_t N,
    double *hmin_P, 
    double *hmax_P, 
    double_vec_t *hist_P,
    double_vec_t *cumh_P,
    uint32_t *ngud_P,
    uint32_t *nbad_P
  );
  /* Generates a histogram of the samples in channel {c} of image {A}, excluding samples
    that are {±INF} or {NAN}.   The histogram will have {N} equal-width bins spanning an interval {[hmin_hmax]}
    that contains the true range {[smin_smax]} of the valid samples.  The chosen bounds {hmin} and {hmax}
    are returned in {*hmin_P} and {*hmax_P}.  The number {N} must be at least 2.
    
    The procedure also returns in {*ngud_P} and {*nbad_P} (if not {NULL})
    the counts {ngud} and {nbad} of valid and invalid samples found, respectively,.
    
    Specifically, for each sample value {s}, the procedure determines the bin {k0}
    of the histogram that contains {s}.  If {s} is the center of the bin, the count {hist.e[k0]} of that bin
    is incremented by 1.  Otherwise, if {s} lies {t} of the way from the center of nin k0 to the
    center of an adjacent bin {k1}, the procedure increment {hist[k0]} by {1-t} and {hist[k1]} by {t}.
    
    If {hist_P} is not {NULL}, the computed histogram {hist} is returned in {*hist_P}.
    Also, if {cumh_P} is not {NULL}, the procedure stores in {*cumh_P} a cumulative
    histogram {cumh} of the same size, such that {cumh.e[k]} is the sum of {hist.e[0..k]},
    for {k} in {0..N-1}.  Thus {cumh[N-1]} should be equal to {ngud}, apart from floating-point roundoff
    errors.
    
    If the actual range {[smin _ smax]} of valid samples is not trivial (that is, if {smin < smax}), the
    histogram range will be such that {smin} is the center of the first bin and {smax} is the center of
    the last one.  Thus no samples are lost in the count.  If {smin == smax}, the procedure 
    will choose {hmin} and {hmax} so that {smin} is the center of bin {N/2}. If there are
    no valid samples, the procedure sets {hmin} and {hmax} to {-1} and {+1}, respectively. */
  
void float_image_hist_write_file(FILE *wr, double hmin, double hmax, double_vec_t *hist);
  /* Writes to {wr} a file with five columns "{k} {slo[k]} {shi[k]} {hist[k]} {cum[k]}"
    where {k} varies in {0..N-1}, {N} is {hist.ne}, {slo[k]} and {shi[k]} are
    the lower and upper bounds of bin {k} of the histogram {hist.e[0..N-1]}, and {cumh[k]} is the sum of
    {hist[i]} for {i} in {0..k-1}. Assumes that the bins span the interval {[hmin _ hmax]}.
    
    Actually writes an extra line at the start with {k=-1}, {slo=shi=hmin} {hist[k]=cumh[k]=0},
    and an extra line at the end with {k=N}, lo-shi=hmax}, {hist[k]=0}, {cumh[k]=cumh[N-1]},
    to simplify the plotting of the histogram. */

void float_image_hist_write_named(char *fname, double hmin, double hmax, double_vec_t *hist);
  /* Same as {float_image_hist_write_file}, but writes to a new file with name "{fname}". */

#endif
