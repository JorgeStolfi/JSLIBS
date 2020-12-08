/* See msm_cand.h */
/* Last edited on 2008-04-19 07:35:36 by stolfi */

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

bool_t msm_cand_is_valid(msm_cand_t *cd, bool_t die)
  { /* Get and check the pairing: */
    msm_pairing_t *p = cd->pr;
    if (! msm_pairing_is_valid(p, die)) { return FALSE; }
    if (! msm_pairing_is_increasing(p, die)) { return FALSE; }
    int ng = msm_pairing_fund_rung_count(p);
    if (ng != 0) 
      { /* Get position counts of the two sequences: */
        int nx = cd->seq[0].npos;
        int ny = cd->seq[1].npos;

        /* The sequences must be non-empty: */
        if (nx <= 0) { fail_test(die, "invalid {nx} in cand"); }
        if (ny <= 0) { fail_test(die, "invalid {ny} in cand"); }

        /* The initial rung must lie in the initial cell: */
        msm_rung_t gini = msm_pairing_get_rung(p, 0);
        int ix = gini.c[0], iy = gini.c[1];
        if ((ix < 0) || (ix >= nx)) { fail_test(die, "invalid {ix} in cand"); }
        if ((iy < 0) || (iy >= ny)) { fail_test(die, "invalid {iy} in cand"); }

        if (msm_pairing_is_circular(p))
          { /* Circular pairing. */
            /* The loopback rung must be greater than the first one: */
            msm_rung_t gfin = msm_pairing_get_rung(p, ng);
            int fx = gfin.c[0], fy = gfin.c[1];
            if (fx <= ix) { fail_test(die, "invalid cand - {fx <= ix}"); }
            if (fy <= iy) { fail_test(die, "invalid cand - {fy <= iy}"); }
            /* Both sequences must be circular: */
            if ((! cd->seq[0].circ) || (! cd->seq[1].circ))  
              { fail_test(die, "circular cand, non circular seqs"); }
            /* Get the coordinate period vector {(DX,DY)}: */
            msm_rung_t rper = msm_pairing_period(p);
            int DX = rper.c[0], DY = rper.c[1];
            /* It must match the difference {(fx,fy) - (ix,iy)}. */
            if (DX != fx - ix) { fail_test(die, "loopback rung X inconsistent with period"); }
            if (DY != fy - iy) { fail_test(die, "loopback rung Y inconsistent with period"); }
            /* stationary or retrograde cands are not allowed. */
            if ((DX <= 0) || (DY <= 0)) 
              { fail_test(die, "stationary or retrograde circular cand"); }
            /* Each coordinate period must be one or more sweeps over the sequence: */
            if (DX % nx != 0) { fail_test(die, "circular cand does not close in X"); }
            if (DY % ny != 0) { fail_test(die, "circular cand does not close in Y"); }
          }
        else
          { /* Open pairing. */
            /* The final rung must be greater than or equal to the first one: */
            msm_rung_t gfin = msm_pairing_get_rung(p, ng-1);
            int fx = gfin.c[0], fy = gfin.c[1];
            if (fx < ix) { fail_test(die, "invalid cand - {fx < ix}"); }
            if (fy < iy) { fail_test(die, "invalid cand - {fy < iy}"); }
            /* If either sequence is open, the pairing must not fall off its end: */
            if ((! cd->seq[0].circ) && (fx >= nx)) { fail_test(die, "invalid open cand - {fx >= nx}"); }
            if ((! cd->seq[1].circ) && (fy >= ny)) { fail_test(die, "invalid open cand - {fy >= ny}"); }
          }
      }
    return TRUE;
  }      

msm_cand_t msm_cand_from_pairing
  ( msm_seq_desc_t *xp, 
    msm_seq_desc_t *yp, 
    msm_pairing_t *pr,
    double score
  )
  { demand(xp->level == yp->level, "cannot pair seqs of different levels");
    msm_cand_t cd;
    cd.seq[0] = (*xp);
    cd.seq[1] = (*yp);
    cd.pr = pr;
    cd.score = score;
    return cd;
  }

void msm_cand_free(msm_cand_t *cd)
  {
    msm_pairing_free(cd->pr);
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
        if(! msm_seq_desc_same_seq(saj, sbj, FALSE)) { return 0; }
        demand(saj->npos == sbj->npos, "inconsistent seq sizes"); 
      }
    int nx = ca->seq[0].npos, ny = ca->seq[1].npos;
    return msm_pairing_count_common_rungs(ca->pr, cb->pr, nx, ny);
  }
  
double msm_cand_compute_score
  ( msm_seq_desc_t *xp, 
    msm_seq_desc_t *yp,
    msm_pairing_t *pr,
    msm_rung_step_score_proc_t *step_score
  )
  { int ng = msm_pairing_fund_rung_count(pr);
    if (ng == 0) { /* Should not happen, but... */ return 0.0; }
    /* Get the rung {g} that precedes rung 0: */
    bool_t circp = msm_pairing_is_circular(pr);
    msm_rung_t g = (circp ?  msm_pairing_get_rung(pr, -1) : msm_rung_none);
    int t; 
    /* Add the scores of all steps that end in rungs {0..ng-1}: */
    double sc = 0;
    for (t = 0; t < ng; t++)
      { msm_rung_t h = msm_pairing_get_rung(pr, t);
        sc += step_score(xp, yp, g, h);
        g = h;
      }
    if (! circp) 
      { /* Add the score of the step from the last rung to nowhere: */ 
        sc += step_score(xp, yp, g, msm_rung_none);
      }
    return sc;
  }

msm_cand_t msm_cand_throw
  ( msm_seq_desc_t *xp, 
    msm_seq_desc_t *yp, 
    int len,
    bool_t circp,
    bool_t atomic,
    bool_t diagonal,
    double skipProb
  )
  { demand(xp->level == yp->level, "cannot pair seqs of different levels");
    msm_cand_t cd;
    cd.score = dgaussrand();
    cd.seq[0] = (*xp);
    cd.seq[1] = (*yp);
    cd.pr = msm_pairing_throw
      ( xp->npos, xp->circ,
        yp->npos, yp->circ,
        len, circp,
        atomic, diagonal, skipProb
      );
    return cd;
  }
  
msm_cand_t msm_cand_sub_throw(msm_cand_t *cd, int len)
  { int old_span = msm_pairing_span(cd->pr, -1); /* Total {x} and {y} span of {cd}: */
    int new_span = 2*len; /* Required total span. */
    msm_cand_t cdn;
    cdn.score = cd->score * ((double)new_span)/((double)old_span); /* Rough guess. */
    cdn.seq[0] = cd->seq[0];
    cdn.seq[1] = cd->seq[1];
    cdn.pr = msm_pairing_sub_throw(cd->pr, len);
    return cdn;
  }

void msm_cand_debug(char *tag, msm_cand_t *cd)
  {
    fprintf(stderr, "    %s ", tag);
    msm_cand_write(stderr, NULL, cd, 5, 8, 6, NULL);
    fprintf(stderr, "\n");
  }

msm_cand_t msm_cand_map_to_finer
  ( msm_cand_t *cd, 
    int nxnew, 
    int nynew,
    bool_t circx,
    bool_t circy, 
    int nwtb
  )
  { 
    bool_t debug = TRUE;

    if (debug)
      { fprintf(stderr, "  ----------------------------------------------------------\n");
        fprintf(stderr, "    mapping candidate\n");
        msm_cand_debug("old: ", cd);
      }
    
    int nold[2]; /* Counts of positions in old sequences. */
    msm_cand_t cdnew;
    cdnew.score = cd->score; 
    int j;
    for (j = 0; j <= 1; j++)
      { nold[j] = cd->seq[j].npos;
        cdnew.seq[j] = cd->seq[j];
        cdnew.seq[j].npos = (j == 0 ? nxnew : nynew);
        cdnew.seq[j].level --;
      }
    cdnew.pr = msm_pairing_map_to_finer
      ( cd->pr, nold[0], nold[1], nxnew, nynew, circx, circy, nwtb );

    if (debug) 
      { fprintf(stderr, "\n");
        msm_cand_debug("new: ", &cdnew);
        fprintf(stderr, "  ----------------------------------------------------------\n");
      }

    return cdnew; 
  }

msm_cand_t msm_cand_interpolate(msm_cand_t *cd)
  { 
    bool_t debug = TRUE;

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

void msm_cand_write
  ( FILE *wr, 
    char *pre, 
    msm_cand_t *cd, 
    int idSize, 
    int nameSize, 
    int ixSize, 
    char *suf
  )
  { 
    if ((pre != NULL) && ((*pre) != 0)) { fputs(pre, wr); }
    int j;
    for (j = 0; j < 2; j++)
      { msm_seq_desc_t *sj = &(cd->seq[j]);
        msm_seq_desc_write(wr, "( ", sj, idSize, nameSize, ixSize, " )  ");
      }
    fprintf(wr, "%16.8f  ", cd->score); 
    msm_pairing_write(wr, cd->pr, ixSize);
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
