#ifndef msm_test_tools_H
#define msm_test_tools_H

/* Miscellaneous tools for test programs. */
/* Last edited on 2022-10-20 10:27:49 by stolfi */

#define msm_test_tools_H_COPYRIGHT \
  "Copyright © 2006  by the State University of Campinas (UNICAMP)"

#include <msm_basic.h>
#include <stdint.h>
#include <msm_pairing.h>
#include <msm_seq_desc.h>
#include <msm_cand.h>
#include <msm_dyn.h>

#include <vec.h>

void msm_debug_double_vec(double *v, int32_t nv, char *fmt);
void msm_debug_int32_vec(int32_t *v, int32_t nv, char *fmt);
  /* These procedures print {v[0..nv-1]} to {stderr}, each with format {fmt},
    separated by spaces and bracketed by '[' and ']'. */

msm_cand_t msm_test_tools_get_optimum_pairing
  ( msm_seq_desc_t *seq0, 
    msm_seq_desc_t *seq1, 
    msm_rung_t gini,
    msm_rung_t gfin, 
    int32_t delta, 
    msm_rung_step_score_proc_t *step_score, 
    bool_t verbose,
    msm_dyn_tableau_t *tb,
    int32_t maxIter
  );
  /* Finds an optimum pairing between the segments of 
    sequences {ap} and {bp} spanned by the two rungs {gini,gfin}. The
    procedure begins constructing the straightest pairing that covers
    both segments. The pairing is then refined {maxIter} times
    by incremental dynamic programming, where at each time the
    procedure considers a band of half-width {delta} around the
    current pairing. */

void msm_test_seq_write_and_plot_named
  ( msm_seq_desc_t *seq,
    double_vec_t *smp, 
    char *title,
    char *name,
    char *tag,
    double fontSize,
    bool_t orig
  );
  /* Writes a test sequence to a disk file called "{name}{tag}.txt".
    Also writes a plot of its sample values to a disk file called
    "{name}{tag}.eps". Assumes that the sequence samples are in the
    vector {*smp}, and that the sequence's descriptor is {*seq}. 
    The scale labels are in size {fontSize} (pt).
    
    If {orig} is false, the X coordinates of the plot will be the 
    indices of each sample in {seq}. If {orig} is true, the indices
    will be mapped to the corresponding indices in the original
    sequence. */

#endif

