/* See msm_seq_desc.h */
/* Last edited on 2008-01-12 08:27:30 by stolfi */

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
    int nbas,         /* Number of fundamental indices in the sequence. */
    bool_t circ      /* Tells whether the sequence is circular. */
  )
  { msm_seq_desc_t seq;
    seq.id = id;
    seq.name = name;
    seq.level = level;
    seq.nbas = nbas;
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
    if (sa->nbas != sb->nbas)
      { fail_test(die, "sequence lengths differ"); }
    if (sa->circ != sb->circ)
      { fail_test(die, "sequence circularities differ"); }
    return TRUE;
  }

void msm_seq_desc_throw_segment(int ns, bool_t circ, int len, int *ini, int *fin)
  { int ix = abrandom(0, ns - (circ ? 1 : len)); 
    int fx = ix + len;
    assert((ix >= 0) && (ix < ns));
    assert(ix <= fx); 
    if (! circ) { assert(fx < ns); }
    (*ini) = ix;
    (*fin) = fx;
  }
 
int msm_seq_desc_filtered_size(int nFine, bool_t circ, int nwtb)
  { return (nFine - (circ ? 0 : nwtb - 1) + 1)/2; }

int msm_seq_desc_map_index_to_coarser(int ixFine, int nFine, int nCoarse, int nwtb)
  { /* Apply the convolution shift: */
    ixFine -= (nwtb-1)/2;
    /* Reduce the index {ixFine} (which may be negative) modulo {nFine}: */
    int kxFine = imod(ixFine, nFine);
    /* Apply the downsampling: */
    int kxCoarse = kxFine/2;
    affirm((kxCoarse >= 0) && (kxCoarse < nCoarse), "bug");
    /* Add the original number {q} of full turns: */
    int q = (ixFine - kxFine)/nFine;
    int ixCoarse = kxCoarse + q*nCoarse;
    return ixCoarse;
  }

int msm_seq_desc_map_index_to_finer(int ixCoarse, int nCoarse, int nFine, int nwtb)
  { /* Reduce the index {ixCoarse} (which may be negative) modulo {nCoarse}: */
    int kxCoarse = imod(ixCoarse, nCoarse);
    /* Reverse the downsampling: */
    int kxFine = 2*kxCoarse;
    affirm((kxFine >= 0) && (kxFine < nFine), "bug");
    /* Add the original number {q} of full turns: */
    int q = (ixCoarse - kxCoarse)/nCoarse;
    int ixFine = kxFine + q*nFine;
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
    int nbasSize,
    char *suf
  )
  { if ((pre != NULL) & ((*pre) != 0)) { fputs(pre, wr); }
    fprintf(wr, "%*d", idSize, s->id);
    fprintf(wr, " %-*s", nameSize, s->name);
    fprintf(wr, " %02d", s->level);
    fprintf(wr, " %*d", nbasSize, s->nbas);
    fprintf(wr, " %c", (s->circ ? 'T' : 'F'));
    if ((suf != NULL) & ((*suf) != 0)) { fputs(suf, wr); }
  }

msm_seq_desc_t msm_seq_desc_read(FILE *rd, char *pre, char *suf)
  { if ((pre != NULL) & ((*pre) != 0)) { fget_skip_spaces(rd); fget_match(rd, pre); }
    msm_seq_desc_t s;
    s.id = fget_int(rd);
    s.name = fget_string(rd);
    s.level = fget_int(rd);
    s.nbas = fget_int(rd);
    s.circ = fget_bool(rd);
    if ((suf != NULL) & ((*suf) != 0)) { fget_skip_spaces(rd); fget_match(rd, suf); }
    return s;
  }
