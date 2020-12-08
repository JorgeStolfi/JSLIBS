/* See {nmsim_write.h} */
/* Last edited on 2019-07-04 17:46:47 by jstolfi */

#define _GNU_SOURCE
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <stdint.h>
#include <assert.h>
#include <math.h>

#include <affirm.h>
#include <fget.h>
#include <nget.h>
#include <jsmath.h>

#include <nmsim_basic.h>

#include <nmsim_write.h>
 
void nmsim_write_double_value(FILE *wr, double v, double prec, bool_t sgn, bool_t fudge_0, bool_t fudge_1)
  { bool_t debug = FALSE;
  
    /* Check for NAN and infinities: */
    char *vbad = NULL;
    if (isnan(v))
      { vbad = "nan"; }
    else if (v == -INF)
      { vbad = "-inf"; } 
    else if (v == +INF)
      { vbad = (sgn ? "+inf" : "+inf"); }
    if (vbad != NULL)
      { /* Special value: */
        fputs(vbad, wr);
        return;
      }
    else
      { /* Choose the number of decimal fraction digits {N}: */
        int32_t N = 0;
        demand(prec > 0.0, "invalid precision");
        while (prec < 0.5) { N++; prec *= 10; }
        /* If we are to fudge non-exact 0 and/or 1, we need at least one decimal fraction digit: */ 
        if ((N == 0) && (fudge_0 || fudge_1)){ N = 1; }
        /* Format the value: */
        char *vx = NULL; /* Full formatted parameter. */
        asprintf(&vx, " %.*f", N, v); 
        if (debug) { fprintf(stderr, "  fudge = %c%c  original = [%s]", "FT"[fudge_0], "FT"[fudge_1], vx); }
        assert(vx != NULL);
        /* Find first non-blank char {*vbeg}: */
        char *vbeg = vx; 
        while ((*vbeg) == ' ') { vbeg++; }
        /* If a sign, blank it for now: */
        if (*vbeg == '-') { (*vbeg) = ' '; vbeg++; }
        assert((*vbeg) != 0); /* Cannot be blank or just a sign... */
        if (strchr(vbeg, '.') != NULL)
          { /* Blank out trailing zeros from fraction: */
            char *vend = vbeg + strlen(vbeg) - 1; /* Last char.*/
            while ((*vend) == '0') { (*vend) = ' '; vend--; }
            if ((*vend) == '.') { (*vend) = ' '; vend--; }
          }
        /* If the number printed as ".00000", the above will blank out everything.  Fix: */
        if ((*vbeg) == ' ') { vbeg = "0"; }
        /* Decide whether printed value must be fudged: */
        bool_t fudge = FALSE; /* If true, must fudge the printed value. */
        char fudge_first = 0, fudge_fill = 0, fudge_last = 0;  /* Fudging params. */
        if ((fudge_0) && (v != 0.0) && ((*vbeg) == '0') && ((*(vbeg+1)) == ' '))
          { /* Print a nonzero value instead: */
            fudge = TRUE;
            fudge_first = '0'; fudge_fill = '0'; fudge_last = '1';
          }
        else if ((fudge_1) && (fabs(v) != 1.0) && ((*vbeg) == '1') && ((*(vbeg+1)) == ' '))
          { /* Print a non-unit value instead: */
            fudge = TRUE;
            if (v > 1.0)
              { fudge_first = '1'; fudge_fill = '0'; fudge_last = '1'; }
            else
              { fudge_first = '0'; fudge_fill = '9'; fudge_last = '9'; }
          }
        if (fudge)
          { /* Replace {*vbeg} ("0     ","1     ") by the fudged value ("0.0001", "1.0001", "0.9999"): */
            /* At this point {vbeg} points to the units digit. */
            char *vend = vbeg + strlen(vbeg) - 1; /* Last char.*/
            assert(vend >= (vbeg+2));
            assert(((*vend) == ' ') && ((*(vend+1)) == 0));
            (*vend) = fudge_last; vend--;
            while (vend > vbeg+1) { (*vend) = fudge_fill; vend--; }
            (*vend) = '.';
            (*vbeg) = fudge_first;
          }
        assert(((*vbeg) != 0) && ((*vbeg) != ' '));
        /* Check if output is just "0": */
        bool_t zout = ((*vbeg) == '0') && (((*(vbeg+1)) == ' ') || ((*(vbeg+1)) == 0));
        /* Print sign, if appropriate: */
        bool_t needs_sgn = ((v < 0) || sgn) && (! zout);
        if (debug) { fprintf(stderr, " cooked = ["); }
        if (needs_sgn) 
          { int32_t sgn = (v < 0 ? '-' : '+');
            fputc(sgn, wr);
            if (debug) { fprintf(stderr, "%c", sgn); }
          }
        /* Print the absolute value: */
        while (((*vbeg) != 0) && ((*vbeg) != ' ')) 
          { fputc((*vbeg), wr); 
            if (debug) { fprintf(stderr, "%c", (*vbeg)); }
            vbeg++;
          }
        if (debug) { fprintf(stderr, "]\n"); }
        /* Reclaim the temporary storage: */
        free(vx);
      }
  }

void nmsim_write_int64_param(FILE *wr, char *pref, char *name, int64_t v, char *fmt)
  { fprintf(wr, "%s%s = ", pref, name);
    fprintf(wr, fmt, v);
    fprintf(wr, "\n");
  }

void nmsim_write_double_param
  ( FILE *wr, 
    char *pref, 
    char *name, 
    double v, 
    double prec,
    bool_t sgn, 
    bool_t fudge_0, 
    bool_t fudge_1
  )
  { 
    fprintf(wr, "%s%s = ", pref, name);
    nmsim_write_double_value(wr, v, prec, sgn, fudge_0, fudge_1);
    fprintf(wr, "\n");
  }

