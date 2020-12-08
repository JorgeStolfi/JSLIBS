/* See msm_seq_desc.h */
/* Last edited on 2008-02-01 12:43:09 by hcgl */

#define msm_seq_desc_C_COPYRIGHT \
  "Copyright © 2005  by the State University of Campinas (UNICAMP)" \
  " and the Federal Fluminense University (UFF)"

#define _GNU_SOURCE
#include <stdlib.h>
#include <math.h>
#include <stdio.h>
#include <string.h>
#include <ctype.h>
#include <assert.h>

#include <vec.h>
#include <nget.h>
#include <fget.h>
#include <filefmt.h>
#include <jsmath.h>
#include <jsrandom.h>

#include <msm_basic.h>

#include <msm_seq_desc.h>

msm_seq_desc_t msm_seq_desc_make
  ( msm_seq_id_t id,  /* Internal identifier of the sequence. */
    char *name,       /* External identifier of the sequence. */
    int level,        /* Filtering level. */
    int npos,         /* Number of fundamental indices in the sequence. */
    bool_t circ      /* Tells whether the sequence is circular. */
  )
  { msm_seq_desc_t seq;
    seq.id = id;
    seq.name = name;
    seq.level = level;
    seq.npos = npos;
    seq.circ = circ;
    return seq;
  }

bool_t msm_seq_desc_same_seq(msm_seq_desc_t *sa, msm_seq_desc_t *sb, bool_t die)
  { if (sa->id != sb->id)
      { fail_test(die, "sequence ids differ"); }
    if ((sa->name != sb->name) && (strcmp(sa->name,sb->name) != 0))
      { fail_test(die, "sequence names differ"); }
    if (sa->level != sb->level)
      { fail_test(die, "sequence levels differ"); }
    return TRUE;
  }

bool_t msm_seq_desc_equal(msm_seq_desc_t *sa, msm_seq_desc_t *sb, bool_t die)
  { if (sa->id != sb->id)
      { fail_test(die, "sequence ids differ"); }
    if ((sa->name != sb->name) && (strcmp(sa->name,sb->name) != 0))
      { fail_test(die, "sequence names differ"); }
    if (sa->level != sb->level)
      { fail_test(die, "sequence levels differ"); }
    if (sa->npos != sb->npos)
      { fail_test(die, "sequence lengths differ"); }
    if (sa->circ != sb->circ)
      { fail_test(die, "sequence circularities differ"); }
    return TRUE;
  }

void msm_seq_desc_throw_segment(int npos, bool_t circ, int nseg, int *ini, int *fin)
  { int ix = abrandom(0, npos - (circ ? 1 : nseg)); 
    int fx = ix + nseg;
    assert((ix >= 0) && (ix < npos));
    assert(ix <= fx); 
    if (! circ) { assert(fx < npos); }
    (*ini) = ix;
    (*fin) = fx;
  }
 
int msm_seq_desc_map_index_to_coarser(int ixFine, int nFine, int nCoarse, bool_t circ, int nwtb)
  { demand(nFine >= nwtb, "nFine is too small for filter kernel");
    demand(nCoarse >= 1, "cannot map to empty sequence");
    int ixCoarse;
    /* Apply the convolution shift: */
    ixFine -= (nwtb-1)/2;
    if (circ)
      { /* Reduce the index {ixFine} (which may be negative) modulo {nFine}: */
        int kxFine = imod(ixFine, nFine);
        /* Map {kxFine} from the fund. range of {sFine} to the fund. range of {sCoarse}: */
        int kxCoarse = (nCoarse * kxFine)/nFine;
        affirm((kxCoarse >= 0) && (kxCoarse < nCoarse), "bug");
        /* Add the original number {q} of full turns: */
        int q = (ixFine - kxFine)/nFine;
        ixCoarse = kxCoarse + q*nCoarse;
      }
    else
      { /* Compute max shifted position in {sFine}, accounting for filter kernel: */
        int ixFineMax = nFine - nwtb;
        /* Clip the index {ixFine} to the valid range: */
        if (ixFine < 0) { ixFine = 0; }
        if (ixFine > ixFineMax) { ixFine = ixFineMax; }
        /* Map {ixFine} from {0..ixFineMax} to {0..nCoarse-1}: */
        ixCoarse = (ixFineMax == 0 ? 0 : ((nCoarse - 1) * ixFine)/ixFineMax);
      }
    return ixCoarse;
  }

int msm_seq_desc_map_index_to_finer(int ixCoarse, int nCoarse, int nFine, bool_t circ, int nwtb)
  { demand(nFine >= nwtb, "nFine is too small for filter kernel");
    demand(nCoarse >= 1, "cannot map to empty sequence");
    int ixFine;
    if (circ)
      { /* Reduce the index {ixCoarse} (which may be negative) modulo {nCoarse}: */
        int kxCoarse = imod(ixCoarse, nCoarse);
        /* Map {kxCoarse} from the fund. range of {sCoarse} to the fund. range of {sFine}: */
        int kxFine = (nFine * kxCoarse)/nCoarse;
        affirm((kxFine >= 0) && (kxFine < nFine), "bug");
        /* Add the original number {q} of full turns: */
        int q = (ixCoarse - kxCoarse)/nCoarse;
        ixFine = kxFine + q*nFine;
      }
    else
      { /* Clip the index {ixCoarse} to the valid range: */
        if (ixCoarse < 0) { ixCoarse = 0; }
        if (ixCoarse >= nCoarse) { ixCoarse = nCoarse - 1; }
        /* Compute max shifted position in {sFine}, accounting for filter kernel: */
        int ixFineMax = nFine - nwtb;
        /* Map {ixCoarse} from {0..nCoarse-1} to {0..ixFineMax}: */
        ixFine = (nCoarse == 1 ? 0 : (ixFineMax * ixCoarse)/(nCoarse -1));
      }
    /* Undo the convolution shift: */
    ixFine += (nwtb-1)/2;
    return ixFine;
  }

void msm_seq_desc_write
  ( FILE *wr, 
    char *pre, 
    msm_seq_desc_t *s, 
    int idSize, 
    int nameSize, 
    int nposSize,
    char *suf
  )
  { if ((pre != NULL) & ((*pre) != 0)) { fputs(pre, wr); }
    fprintf(wr, "%*d", idSize, s->id);
    fprintf(wr, " %-*s", nameSize, s->name);
    fprintf(wr, " %02d", s->level);
    fprintf(wr, " %*d", nposSize, s->npos);
    fprintf(wr, " %c", (s->circ ? 'T' : 'F'));
    if ((suf != NULL) & ((*suf) != 0)) { fputs(suf, wr); }
  }

msm_seq_desc_t msm_seq_desc_read(FILE *rd, char *pre, char *suf)
  { if ((pre != NULL) & ((*pre) != 0)) { fget_skip_spaces(rd); fget_match(rd, pre); }
    msm_seq_desc_t s;
    s.id = fget_int(rd);
    s.name = fget_string(rd);
    s.level = fget_int(rd);
    s.npos = fget_int(rd);
    s.circ = fget_bool(rd);
    if ((suf != NULL) & ((*suf) != 0)) { fget_skip_spaces(rd); fget_match(rd, suf); }
    return s;
  }
