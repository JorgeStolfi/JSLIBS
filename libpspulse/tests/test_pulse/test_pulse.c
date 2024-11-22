/* Tests the unidimensional mother pulses */
/* Last edited on 2024-11-20 05:32:37 by stolfi */

#define _GNU_SOURCE
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include <jsfile.h>
#include <psp_pulse.h>
#include <psp_pulse_mother.h>

#define NSTEPS 30
  /* Number of samples per interval. */

#define MAXC 2
  /* Maximum continuity order. */
  
#define MAXDER 2
  /* Maximum derivative order to evaluate. */
  
#define BIGSZ 64
  /* A power of 2 greater than the support size of any spline. */

#define DIFF_STEP 0.0001
  /* Step for numerical derivation. */

#define DELTA ((MAXDER+1)*DIFF_STEP)
  /* Min distance to keep from knots. Must be greater than {MAXDER*DIFF_STEP}.  */
  
#define OUT_DIR "out"
  /* Directory for output files. */

/* INTERNAL PROTOS */

int main (int argc, char **argv);

void show_family ( psp_pulse_family_t *fam );

void show_pulse ( psp_pulse_family_t *fam, psp_pulse_t *p );

char *tsp_make_name
  ( psp_pulse_family_t *fam,
    psp_pulse_t *p
  );
    
void sample
  ( FILE *wr, 
    psp_pulse_family_t *fam,
    psp_pulse_t *p,
    double x, 
    psp_cont_t ord,
    double e[]
  );
  /* Ouputs to {wr} one sample point of the pulse {p} of family {fam}. 
    The argument {x} is the sampling point relative to the 
    root interval, assumed to be {[0 _ p->gsz]}. Plots derivatives up to {ord}
    and accumulates the analytical--numerical error in {e[0..ord-1]}. */

/* IMPLEMENTATIONS */

int main (int argc, char **argv)
  { int c;
    psp_pulse_kind_t pkind;
    for (pkind = psp_pulse_kind_MIN; pkind <= psp_pulse_kind_MAX; pkind++)
      { for (c = -1; c <= MAXC; c++)
          { int maxg = 2*c + 3;
            int g;
            for (g = 0; g <= maxg; g++)
              { psp_pulse_family_t fam = psp_pulse_family(pkind, c, g);
                show_family(&fam);
              }
          }
      }
    fprintf(stderr, "done.\n");
    return 0;
  }
                
void show_family ( psp_pulse_family_t *fam )
  {
    fprintf(stderr, "=== family ");
    psp_pulse_family_print(stderr, fam);
    fprintf(stderr, "  (%d mother pulses)", (int)fam->nmp);
    fprintf(stderr, " ===\n");
    
    psp_pulse_mother_index_t pix;
    for (pix = 0; pix < fam->nmp; pix++)
      { psp_grid_size_t msz = psp_pulse_mother_supp_count(fam, pix);
        psp_grid_size_t gsz = 1;
        do
          { psp_pulse_t p = psp_pulse(fam, pix, gsz, 0); 
            show_pulse(fam, &p);
            gsz = 2*gsz;
          }
        while (gsz <= 2*msz);
      }
  }
  
void show_pulse ( psp_pulse_family_t *fam, psp_pulse_t *p )
  {
    fprintf(stderr, "--- pulse ");
    psp_pulse_print(stderr, fam, p);
    fprintf(stderr, " ---\n");
    
    char *fname = tsp_make_name(fam, p);
    FILE *wr;
    /* Max derivation order to ask: */
    psp_cont_t ord = fam->g+1; 
    if (ord > MAXDER) { ord = MAXDER; }
    /* Max numerical-analytical discrepancy in each derivative: */
    double e[ord+1];
    { int k; for (k = 0; k <= ord; k++) { e[k] = 0; } }
    wr = open_write(fname, TRUE);
    /* Number of plot intervals in grid: */
    int NPLOT = (p->gsz)*NSTEPS; 
    int iz; /* Index of plot sampling point. */
    for (iz = -2; iz <= NPLOT+2; iz++)
      { /* Compute position of plot sampling point relative to mother's support: */
        double x = ((double)iz)/((double)NSTEPS);
        if ((iz % NSTEPS) == 0)
          { sample(wr, fam, p, x-DELTA, ord, e);
            sample(wr, fam, p, x+DELTA, ord, e);
          }
        else
          { sample(wr, fam, p, x, ord, e); }
      }
    if (wr == stdout) { fflush(wr); } else { fclose(wr); }
    fprintf(stderr, "max errors in derivatives: \n");
    { int k;
      for (k = 0; k <= MAXDER; k++) 
        { fprintf(stderr, " %24.16e", e[k]); }
    }
    fprintf(stderr, "\n");
  }

char *tsp_make_name
  ( psp_pulse_family_t *fam,
    psp_pulse_t *p
  )
  {
    char pkind_txt = psp_pulse_kind_to_char(fam->pkind);
    char *c_txt;
    if (fam->c >= 0) 
      { c_txt = jsprintf("%d", fam->c); } 
    else
      { c_txt = "-"; }
    char *gsz_txt;
    if (p->gsz < p->msz) 
      { gsz_txt = jsprintf("%02lld", p->gsz); }
    else 
      { gsz_txt = "--"; }
    char *fname = jsprintf(
        "%s/pu-z%s-%c-c%s-g%02d-m%02d.plt", 
        OUT_DIR, gsz_txt, pkind_txt, c_txt, fam->g, p->pix
      );
    if (fam->c >= 0) { free(c_txt); }
    return fname;
  }

void sample
  ( FILE *wr, 
    psp_pulse_family_t *fam, 
    psp_pulse_t *p, 
    double x, 
    psp_cont_t ord,
    double e[]
  )
  {
    int k, j;
    interval_t R = (interval_t){{ 0.0, (double)p->gsz }};
    
    /* Get pulse value and analyic derivatives {f[0..MAXDER]}: */
    double f[ord+1];
    for (k = 0; k <= ord; k++) { f[k] = 0; }
    psp_pulse_eval(fam, p, &R, x, ord, f);
    
    /* Get pulse value and numerical derivatives {d[0..MAXDER]}: */
    double d[ord+1];
    for (k = 0; k <= ord; k++)
      { psp_pulse_eval(fam, p, &R, x + k*DIFF_STEP, 0, &(d[k])); }
    for (k = 1; k <= ord; k++)
      { for (j = ord; j >= k; j--) { d[j] = (d[j]-d[j-1])/DIFF_STEP; } }

    /* Print values and derivs, and update the max error {e}. */
    fprintf(wr, "%8.5f", x);
    for (k = 0; k <= ord; k++) 
      { fprintf(wr, "  %10.5f %10.5f", f[k], d[k]);
        double err = fabs(f[k]-d[k]);
        if (err > e[k]) { e[k] = err; }
      }
    for (k = ord+1; k<= MAXDER; k++)
      { fprintf(wr, "  %10s %10s", "0", "0"); }
    fprintf(wr, "\n");
  }

