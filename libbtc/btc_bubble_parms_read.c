/* See {btc_bubble_parms_read.h} */
/* Last edited on 2024-12-05 10:23:11 by stolfi */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include <bool.h>
#include <fget.h>
#include <jsfile.h>
#include <affirm.h>

#include <btc_bubble_t.h>
#include <btc_date_read.h>
#include <btc_bubble_parms_validate.h>

#include <btc_bubble_parms_read.h>

void btc_bubble_parms_read(char* fName, int nd, char* dt[], int* nbP, btc_bubble_t** bpP)
  {    
    bool_t debug = TRUE;
    
    int nb_max = 100;
    
    FILE* rd = open_read(fName, TRUE);
    btc_bubble_t* bp = notnull(malloc(nb_max*sizeof(btc_bubble_t)), "no mem"); /* Bubble parameters, indexed {0..nb-1}. */
    
    int nb = 0;    /* Number of bubbles found. */
    int nlin = 0;  /* Line number, starting from 1. */
    bool_t ok = TRUE;
    while (! fget_skip_and_test_char(rd, EOF))
      { 
        /* One more line in file; should have skipped leading blanks: */
        nlin++;
        
        /* Check for comments and blank lines: */
        fget_skip_spaces(rd);
        if (fget_test_char(rd, '#')) { fget_skip_to_eol(rd); continue; }
        if (fget_test_char(rd, '\n')) { continue; }
        
        /* One more bubble: */
        demand(nb < nb_max, "too many bubbles");

        btc_bubble_t* bpj = &(bp[nb]);

        if (debug) { fprintf(stderr, "  %02d", nb); }

        bpj->coef = fget_double(rd); /* Coefficient in linear comination. */

        bpj->id_ini_sg = btc_date_read(rd, "dt_ini_sg", nd, dt, debug); /* Date of start of relevant period */

        bpj->rt_up = fget_double(rd); /* Rate of rally. */
        if (debug) { fprintf(stderr, "  %8.5f", bpj->rt_up); }

        bpj->id_fin_up = btc_date_read(rd, "dt_fin_up", nd, dt, debug); /* Date of end of rally */

        int id_ini_dn = btc_date_read(rd, "dt_ini_dn", nd, dt, debug); /* Date of start of decay */
        bpj->wd_plat = id_ini_dn - bpj->id_fin_up; /* Width of plateau. */

        bpj->rt_dn = fget_double(rd);  /* rate of decay. */
        if (debug) { fprintf(stderr, " %8.5f", bpj->rt_dn); }

        bpj->id_fin_sg = btc_date_read(rd, "dt_fin_sg", nd, dt, debug); /* Date of end of relevant period */

        bpj->tag = fget_string(rd); /* Tag of bubble. */
        if (debug) { fprintf(stderr, "  %-12s", bpj->tag); }

        bpj->color = fget_string(rd); /* Color for plotting. */
        if (debug) { fprintf(stderr, " %s", bpj->color); }

        if (debug) { fprintf(stderr, "\n"); }

        /* Validity checks: */
        ok &= btc_bubble_parms_validate(fName,nlin, nd,dt, bpj);

        nb++;
        fget_comment_or_eol(rd, '#', NULL);
      }
    fclose(rd);
    if (debug) { fprintf(stderr, "bubbles = %d\n", nb); }
      
    demand(ok, "aborted");
    
    (*nbP) = nb;
    (*bpP) = notnull(realloc(bp, nb*sizeof(btc_bubble_t)), "no mem");
  } 
  
