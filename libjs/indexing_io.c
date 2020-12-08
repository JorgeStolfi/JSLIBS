/* See indexing_io.h */
/* Last edited on 2011-12-21 21:36:47 by stolfi */
/* Copyright © 2007 by Jorge Stolfi, from University of Campinas, Brazil. */
/* See the rights and conditions notice at the end of this file. */

#include <indexing_io.h>
#include <indexing.h>
#include <stdio.h>
#include <stdint.h>
#include <assert.h>
#include <jswsize.h>

void ix_print_dim ( FILE *wr, ix_dim_t d )
  { fprintf(wr, "%u", d); }

void ix_print_pos ( FILE *wr, ix_pos_t p )
{ fprintf(wr, ("%" uint64_u_fmt), p); }

void ix_print_indices ( FILE *wr, char *lp, ix_dim_t d, ix_index_t ix[], int wd, char *sp, char *rp )
  { int i; 
    if (sp == NULL) { sp = " "; }
    char *epref = "";
    if (lp != NULL) { fprintf(wr, "%s", lp); }
    for (i = 0; i < d; i++) { fprintf(wr, ("%s%*" int64_d_fmt), epref, wd, ix[i]); epref = sp; }
    if (rp == NULL) { fprintf(wr, "%s", rp); }
  }

void ix_print_sizes ( FILE *wr, char *lp, ix_dim_t d, ix_size_t sz[], int wd, char *sp, char *rp )
  { int i; 
    if (sp == NULL) { sp = " "; }
    char *epref = "";
    if (lp != NULL) { fprintf(wr, "%s", lp); }
    for (i = 0; i < d; i++) { fprintf(wr, ("%s%*" uint64_u_fmt), epref, wd, sz[i]); epref = sp; }
    if (rp == NULL) { fprintf(wr, "%s", rp); }
  }

void ix_print_steps ( FILE *wr, char *lp, ix_dim_t d, ix_step_t st[], int wd, char *sp, char *rp )
  { int i; 
    if (sp == NULL) { sp = " "; }
    char *epref = "";
    if (lp != NULL) { fprintf(wr, "%s", lp); }
    for (i = 0; i < d; i++)
      { ix_step_t sti = st[i];
        fprintf(wr, "%s", epref); epref = sp; 
        if (sti == 0)
          { fprintf(wr, ("%*" int64_d_fmt), wd, st[i]); }
        else
          { fprintf(wr, ("%+*" int64_d_fmt), wd, st[i]); }
      }
    if (rp == NULL) { fprintf(wr, "%s", rp); }
  }

void ix_print_parms
 ( FILE *wr, 
   char *pre, 
   ix_dim_t d, 
   ix_pos_t *bp, 
   ix_size_t sz[], 
   ix_step_t st[], 
   int wd, 
   char *suf
 )
  { if (pre == NULL) { pre = ""; }
    if (suf == NULL) { suf = "\n"; }
    fprintf(wr, "%saxes = ", pre); ix_print_dim(wr, d); fprintf(wr, "%s", suf);
    if (bp != NULL) 
      { fprintf(wr, "%sbase =", pre); ix_print_pos(wr, (*bp)); fprintf(wr, "%s", suf); }
    if (st != NULL)
      { fprintf(wr, "%sstep =", pre);
        ix_print_sizes(wr, "", d, sz, wd, " ", ""); 
        fprintf(wr, "%s", suf);
      }
    if (sz != NULL)
      { fprintf(wr, "%ssize =", pre);
        ix_print_steps(wr, "", d, st, wd, " ", ""); 
        fprintf(wr, "%s", suf);
      }
  }
