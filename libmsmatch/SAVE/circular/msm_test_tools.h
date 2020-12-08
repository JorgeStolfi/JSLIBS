#ifndef msm_test_tools_H
#define msm_test_tools_H

/* Miscellaneous tools for test programs. */
/* Last edited on 2008-01-29 20:17:18 by hcgl */

#define msm_test_tools_H_COPYRIGHT \
  "Copyright © 2006  by the State University of Campinas (UNICAMP)"

#include <msm_basic.h>
#include <msm_pairing.h>
#include <msm_seq_desc.h>
#include <msm_cand.h>
#include <msm_dyn.h>

#include <vec.h>

void msm_debug_double_vec(double *v, int nv, char *fmt);
void msm_debug_int_vec(int *v, int nv, char *fmt);
  /* These procedures print {v[0..nv-1]} to {stderr}, each with format {fmt},
    separated by spaces and bracketed by '[' and ']'. */

msm_cand_t msm_test_tools_get_optimum_pairing
  ( msm_seq_desc_t *xp, 
    msm_seq_desc_t *yp, 
    int delta, 
    msm_rung_step_score_proc_t *step_score, 
    msm_dyn_tableau_t *tb,
    int maxIter
  );
  /* Finds an optimum pairing between the sequences {ap} and {bp}. The
    procedure begins constructing the straightest pairing that covers
    both sequences once. The pairing is then refined {maxIter} times
    by incremental dynamic programming, where at each time the
    procedure considers a band of half-width {delta} around the
    current pairing. */

void msm_test_seq_write_and_plot_named
  ( msm_seq_desc_t *seq,
    int den,
    double_vec_t *smp, 
    char *title,
    char *name,
    char *tag,
    double fontSize
  );
  /* Writes a test sequence to a disk file called "{name}{tag}.txt".
    Also writes a plot of its sample values to a disk file called
    "{name}{tag}.eps". Assumes that the sequence samples are in the
    vector {*smp}, and that the sequence's descriptor is {*seq}, with
    subsampling factor {den}. */

#endif

