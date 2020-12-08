/* See msm_cand.h */
/* Last edited on 2017-04-28 18:01:52 by stolfilocal */

#define msm_cand_C_COPYRIGHT \
  "Copyright © 2005  by the State University of Campinas (UNICAMP)"

#define _GNU_SOURCE
#include <assert.h>
#include <string.h>
#include <stdio.h>
#include <math.h>

#include <affirm.h>
#include <fget.h>
#include <nget.h>
#include <filefmt.h>
#include <jsrandom.h>
#include <jsmath.h>

#include <msm_basic.h>
#include <msm_rung.h>
#include <msm_pairing.h>
#include <msm_seq_desc.h>
#include <msm_cand.h>
#include <msm_cand_vec.h>

/* EXPORTED IMPLEMENTATIONS */

msm_cand_t msm_cand_from_pairing
  ( msm_seq_desc_t *seq0, 
    msm_seq_desc_t *seq1, 
    msm_pairing_t *pr,
    double score
  )
  { bool_t verbose = FALSE;
    (void)msm_cand_pairing_is_valid(pr, seq0->size, seq1->size, 0, FALSE, TRUE);
    msm_cand_t cd = (msm_cand_t)
      { .seq = { (*seq0), (*seq1) },
        .pr = pr,
        .score = score
      };
    if (verbose) { msm_cand_debug("created: ", &cd); }
    return cd;
  }

void msm_cand_free(msm_cand_t *cd)
  {
    msm_pairing_free(cd->pr);
  }

bool_t msm_cand_equivalent(msm_cand_t *cda, msm_cand_t *cdb, bool_t die)
  {
    int j;
    for (j = 0; j < 2; j++)
      { msm_seq_desc_t *saj = &(cda->seq[j]);
        msm_seq_desc_t *sbj = &(cdb->seq[j]);
        if(! msm_seq_desc_equal(saj, sbj, die)) { return FALSE; }
      }
    return msm_pairing_equal(cda->pr, cdb->pr, die);
  }

#define msm_cand_checking_level 1
  /* Define this as positive to get paranoid checking. */

#define CHKLEV(n) (msm_cand_checking_level >= (n))
  /* Use {if(CHKLEV(n)){...}} for checks of paranoia level {n}. */

#define msm_cand_debug_level 2
  /* Set this positive to get some debugging printouts. */

#define DEBLEV(n) (msm_cand_debug_level >= (n))
  /* Use {if(DEBLEV(n)){...}} for debug printouts of level {n}. */

int msm_cand_count_common_rungs(msm_cand_t *ca, msm_cand_t *cb)
  { int j;
    for (j = 0; j < 2; j++)
      { msm_seq_desc_t *saj = &(ca->seq[j]);
        msm_seq_desc_t *sbj = &(cb->seq[j]);
        if(! msm_seq_desc_equal(saj, sbj, FALSE)) { return 0; }
      }
    int n0 = ca->seq[0].size;
    int n1 = ca->seq[1].size;
    return msm_pairing_count_common_rungs(ca->pr, cb->pr, n0, n1);
  }
  
double msm_cand_compute_score
  ( msm_seq_desc_t *seq0, 
    msm_seq_desc_t *seq1,
    msm_pairing_t *pr,
    msm_rung_step_score_proc_t *step_score,
    bool_t ini_step,
    bool_t fin_step
  )
  { int ng = msm_pairing_num_rungs(pr);
    if (ng == 0) { /* Should not happen, but... */ return 0.0; }
    
    double sc = 0;  /* Total score. */
    int t_ini; /* Index of first step. */
    
    msm_rung_t g; /* Preceding rung. */
    if (ini_step)
      { /* Get the rung {g} that precedes rung 0: */
        g = msm_rung_none;
        t_ini = 0;
      }
    else
      { /* Start with rung 0: */
        g = msm_pairing_get_rung(pr, 0);
        t_ini = 1;
      }
    /* Add the scores of all steps that end in rungs {t_ini..ng-1}: */
    int t; 
    for (t = t_ini; t < ng; t++)
      { msm_rung_t h = msm_pairing_get_rung(pr, t);
        sc += step_score(seq0, seq1, &g, &h);
        g = h;
      }

    if (fin_step)
      { /* Add the score of the step from the last rung to nowhere: */ 
        sc += step_score(seq0, seq1, &g, NULL);
      }
    return sc;
  }

msm_cand_t msm_cand_throw
  ( msm_seq_desc_t *seq0, 
    msm_seq_desc_t *seq1, 
    int len,
    bool_t atomic,
    bool_t diagonal,
    double skipProb
  )
  { msm_cand_t cd;
    cd.score = dgaussrand();
    cd.seq[0] = (*seq0);
    cd.seq[1] = (*seq1);
    cd.pr = msm_pairing_throw
      ( seq0->size,
        seq1->size,
        len,
        atomic, diagonal, skipProb
      );
    return cd;
  }
  
msm_cand_t msm_cand_sub_throw(msm_cand_t *cd, int len)
  { int old_span = msm_pairing_span(cd->pr, -1); /* Total {s0} and {s1} span of {cd}: */
    int new_span = 2*len; /* Required total span. */
    double score = cd->score * ((double)new_span)/((double)old_span); /* Rough guess. */
    msm_cand_t cdn = (msm_cand_t)
      { .score = score,
        .seq = { cd->seq[0], cd->seq[1] },
        .pr = msm_pairing_sub_throw(cd->pr, len)
      };
    return cdn;
  }

void msm_cand_debug(char *tag, msm_cand_t *cd)
  {
    fprintf(stderr, "    %s ", tag);
    msm_cand_write(stderr, NULL, cd, 5, 8, 6, NULL,TRUE);
    fprintf(stderr, "\n");
  }

msm_cand_t msm_cand_map(msm_cand_t *cd, msm_seq_desc_t *s0_new, msm_seq_desc_t *s1_new)
  {
    bool_t debug = FALSE;

    if (debug)
      { fprintf(stderr, "  ----------------------------------------------------------\n");
        fprintf(stderr, "    mapping candidate\n");
        msm_cand_debug("old: ", cd);
      }

    msm_seq_desc_t *s0_old = &(cd->seq[0]);
    msm_seq_desc_t *s1_old = &(cd->seq[1]);

    msm_cand_t cdnew;
    cdnew.seq[0] = (*s0_new);
    cdnew.seq[1] = (*s1_new);
    cdnew.score = cd->score; 
    cdnew.pr = msm_pairing_map(cd->pr, s0_old, s1_old, s0_new, s1_new);

    if (debug) 
      { fprintf(stderr, "\n");
        msm_cand_debug("new: ", &cdnew);
        fprintf(stderr, "  ----------------------------------------------------------\n");
      }

    return cdnew; 
  }


msm_cand_t msm_cand_interpolate(msm_cand_t *cd)
  { 
    bool_t debug = FALSE;

    if (debug)
      { fprintf(stderr, "  ----------------------------------------------------------\n");
        fprintf(stderr, "    interpolating candidate\n");
        msm_cand_debug("old: ", cd);
      }
    
    msm_cand_t cdnew; 
    int j;
    cdnew.score = cd->score; 
    for (j = 0; j <= 1; j++) { cdnew.seq[j] = cd->seq[j]; }
    cdnew.pr = msm_pairing_interpolate(cd->pr);

    if (debug) 
      { fprintf(stderr, "\n");
        msm_cand_debug("new: ", &cdnew);
        fprintf(stderr, "  ----------------------------------------------------------\n");
      }

    return cdnew; 
  }

msm_cand_t msm_cand_make_increasing(msm_cand_t *cd, int minIncrEach, int minIncrSum)
  { 
    bool_t debug = FALSE;

    if (debug)
      { fprintf(stderr, "  ----------------------------------------------------------\n");
        fprintf(stderr, "    making candidate increasing\n");
        msm_cand_debug("old: ", cd);
      }
    
    msm_cand_t cdnew; 
    int j;
    cdnew.score = cd->score; 
    for (j = 0; j <= 1; j++) { cdnew.seq[j] = cd->seq[j]; }
    cdnew.pr = msm_pairing_make_increasing(cd->pr, minIncrEach, minIncrSum);

    if (debug) 
      { fprintf(stderr, "\n");
        msm_cand_debug("new: ", &cdnew);
        fprintf(stderr, "  ----------------------------------------------------------\n");
      }

    return cdnew; 
  }

void msm_cand_write
  ( FILE *wr, 
    char *pre, 
    msm_cand_t *cd, 
    int idSize, 
    int nameSize, 
    int ixSize, 
    char *suf,
    bool_t writePairing
  )
  { 
    if ((pre != NULL) && ((*pre) != 0)) { fputs(pre, wr); }
    int j;
    for (j = 0; j < 2; j++)
      { msm_seq_desc_t *sj = &(cd->seq[j]);
        msm_seq_desc_write(wr, "( ", sj, idSize, nameSize, ixSize, " )  ");
      }
    fprintf(wr, "%16.8f  ", cd->score); 
    msm_pairing_write(wr, cd->pr, ixSize, writePairing);
     
    if ((suf != NULL) && ((*suf) != 0)) { fputs(suf, wr); }
    fflush(wr);
  }

msm_cand_t msm_cand_read(FILE *rd, char *pre, char *suf)
  { if ((pre != NULL) && ((*pre) != 0)) { fget_skip_spaces(rd); fget_match(rd, pre); }
    msm_cand_t cd;
    int j;
    for (j = 0; j < 2; j++)
      { cd.seq[j] = msm_seq_desc_read(rd, "(", ")"); 
        
      }
    cd.score = fget_double(rd);
    cd.pr = msm_pairing_read(rd);
    if ((suf != NULL) && ((*suf) != 0)) { fget_skip_spaces(rd); fget_match(rd, suf); }
    return cd;
  }

bool_t msm_cand_pairing_is_valid(msm_pairing_t *pr, int n0, int n1, int minIncrEach, bool_t atomic, bool_t die)
  { if (! msm_pairing_is_valid(pr, die)) { return FALSE; }
    int ng = msm_pairing_num_rungs(pr);
    assert(ng >= 0);
    if (ng != 0) 
      { /* The pairing must be non-decreasing on both sides, non-stuttering: */
        if (minIncrEach < 0) { minIncrEach = 0; } /* Require non-decreasing steps in any case. */
        int minIncrSum = 1; /* Require non-stuttering steps in any case. */
        if (! msm_pairing_is_increasing(pr, minIncrEach, minIncrSum, atomic, die)) { return FALSE; }

        /* Check if first and last rungs are in the valid range, rest follows: */
        msm_rung_t gini = msm_pairing_get_rung(pr, 0);
        msm_rung_t gfin = msm_pairing_get_rung(pr, ng-1);
        if ((gini.c[0] < 0) || (gfin.c[0] >= n0)) { fail_test(die, "invalid 0-side indices in pairing"); }
        if ((gini.c[1] < 0) || (gfin.c[1] >= n1)) { fail_test(die, "invalid 1-side indices in pairing"); }
      }
    return TRUE;
  }

bool_t msm_cand_is_valid(msm_cand_t *cd, int minIncrEach, bool_t atomic, bool_t die)
  { 
    /* Check the descriptors: */
    if (! msm_seq_desc_is_valid(&(cd->seq[0]), FALSE)) { fail_test(die, "invalid sequence 0"); }
    if (! msm_seq_desc_is_valid(&(cd->seq[1]), FALSE)) { fail_test(die, "invalid sequence 1"); }
    /* Get and check the pairing: */
    msm_pairing_t *p = cd->pr;
    if (! msm_cand_pairing_is_valid(p, cd->seq[0].size, cd->seq[1].size, minIncrEach, atomic, die)) { return FALSE; }
    return TRUE;
  }      
