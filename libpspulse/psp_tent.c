/* See psp_tent.h */
/* Last edited on 2011-09-18 13:51:33 by stolfilocal */

#include <math.h>
#include <assert.h>

#include <affirm.h>
#include <vec.h>

#include <psp_basic.h>
#include <psp_pulse.h>
#include <psp_tent.h>

/* IMPLEMENTATIONS */

void psp_tent_eval
  ( psp_dim_t d, 
    psp_pulse_family_t fam[], 
    psp_pulse_t p[],
    interval_t R[],
    double x[],
    psp_degree_t ord[],
    double f[]
  )
  { double pf[psp_pulse_MAX_DEGREE + 1];
    f[0] = 1.0;
    int nf = 1;
    int i;
    for (i = 0; i < d; i++) 
      { interval_t *Ri = (R == NULL ? NULL : &(R[i]));
        psp_pulse_eval(&(fam[i]), &(p[i]), Ri, x[i], ord[i], pf);
        int j;
        for (j = 0; j < nf; j++)
          { int k; 
            for (k = 0; k <= ord[i]; k++)
              { f[k*nf + j] = f[j]*pf[k]; }
          }
        nf *= (ord[i] + 1);
      }
  }

void psp_tent_eval_total
  ( psp_dim_t d, 
    psp_pulse_family_t fam[], 
    psp_pulse_t p[],
    interval_t R[],
    double x[],
    psp_degree_t ord,
    double f[]
  )
  { double pf[psp_pulse_MAX_DEGREE + 1];
    f[0] = 1.0;
    int i;
    for (i = 0; i < d; i++) 
      { interval_t *Ri = (R == NULL ? NULL : &(R[i]));
        psp_pulse_eval(&(fam[i]), &(p[i]), Ri, x[i], ord, pf);
        affirm(FALSE, "not implemented");
      }
  }

/* PRINTOUT: */

void psp_tent_print
  ( FILE *wr, 
    psp_dim_t d,
    psp_pulse_family_t fam[], 
    psp_pulse_t p[]
  )
  { int i;
    /* Print pulse indices: */
    for (i = 0; i < d; i++) 
      { fprintf(wr, "[");
        psp_pulse_print(wr, &(fam[i]), &(p[i]));
        fprintf(wr, "]");
      }
  }

