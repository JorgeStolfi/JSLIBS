/* See {nmsim_test.h} */
/* Last edited on 2020-12-04 21:22:01 by jstolfi */

#define _GNU_SOURCE
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <assert.h>
#include <math.h>

#include <bool.h>
#include <affirm.h>
#include <jsfile.h>

#include <nmsim_firing_func.h>

#include <nmsim_test.h>

void nmsim_test_firing_func
  ( nmsim_firing_func_t *Phi,
    double r,
    int32_t ns
  )
  {
    char class = Phi->class;
    double V_M = Phi->V_M;
    double V_D = Phi->V_D;

    fprintf(stderr, "--- testing firing function class = %c ---\n", class);
    fprintf(stderr, "key potentials V_M = %6.2f V_D = %6.2f\n", V_M, V_D);
    fprintf(stderr, "relative plot radius = %.2f plotting steps = %d\n", r, ns);
    
    /* Filename: */
    char *fname = NULL;
    asprintf(&fname, "out/Phi_%c_M%+06.2f_D%06.2f.txt", class, V_M, V_D);
    FILE *wr = open_write(fname, TRUE);

    /* Run simulation: */
    double eps = 0.0001; /* Step for numeric derivative computation. */
    int32_t nbad = 0; /* Number of errors detected. */
    double Vmin = V_M - r * V_D;
    double Vmax = V_M + r * V_D;
    
    for (int32_t k = -ns; k <= ns; k++)
      { /* Compute the plot potential {V}: */
        double z = ((double)k)/((double)ns);
        double V = V_M + z * r * V_D;
        
        /* Evaluate the firing function: */
        double pr, dpr;
        nmsim_firing_func_eval(Phi, V, &pr, &dpr);
        fprintf(wr, "%+10.5f %17.15f %17.15f\n", V, pr, dpr);
        
        /* Check result ranges: */
        demand((pr >= 0.0) && (pr <= 1.0), "invalid probability");
        demand(dpr >= 0.0, "invalid prob density");
        
        /* Check derivative: */
        if ((pr <= 1.0e-12) || (pr >= 1.0-1.0e-12))
          { demand(fabs(dpr) < 1.0e-12, "nonzero derivative on flat region"); }
        else
          { /* Check derivative numerically: */
            double pra, prb;
            nmsim_firing_func_eval(Phi, V-eps, &pra, NULL);
            nmsim_firing_func_eval(Phi, V+eps, &prb, NULL);
            double dpr_num = (prb - pra)/(2*eps);
            double dpr_err = dpr - dpr_num;
            if (fabs(dpr_err/(dpr + 1.0e-4)) > 1.0e-6)
              { fprintf(stderr, "** incorrect derivative V = %10.5f", V);
                fprintf(stderr, "  dpr = %17.5f dpr_num = %17.5f", dpr, dpr_num);
                fprintf(stderr, "  dpr_err = %17.15f\n", dpr_err);
                nbad++;
                /* demand(FALSE, "aborted"); */
              }
          }
        /* Check inverse: */
        double V_check = nmsim_firing_func_eval_inv(Phi, pr);
        double pr_check;
        nmsim_firing_func_eval(Phi, V_check, &pr_check, NULL);
        double pr_err = pr - pr_check;
        if (fabs(pr_err > 1.0e-10))
          { fprintf(stderr, "** incorrect inverse V = %10.5f", V);
            fprintf(stderr, "  pr = %17.5f V_check = %20.15f", pr, V_check);
            fprintf(stderr, "  pr_check = %17.15f pr_err = %17.15f\n", pr_check, pr_err);
            nbad++;
            /* demand(FALSE, "aborted"); */
          }
      }
    demand(nbad == 0, "aborted");

    /* Plot reference lines: */
    for (int32_t k = -1; k <= +1; k++)
      { double V = V_M + k*V_D;
        double pr;
        nmsim_firing_func_eval(Phi, V, &pr, NULL);
        if ((pr > 0.001) && (pr < 0.999))
          { fprintf(wr, "\n");
            fprintf(wr, "0 0 0 %10.5f  %17.15f\n", Vmin+eps, pr);
            fprintf(wr, "0 0 0 %10.5f  %17.15f\n", Vmax-eps, pr);
          }
        fprintf(wr, "\n");
        fprintf(wr, "0 0 0 %10.5f  %17.15f\n", V, 0.0001);
        fprintf(wr, "0 0 0 %10.5f  %17.15f\n", V, 0.9999);
      }
    
    fclose(wr);
    free(fname);
  }
     
double *nmsim_test_NAN_vector(size_t n)
  { 
    double *v = notnull(malloc(n*sizeof(double)), "no mem");
    for (uint64_t i = 0; i < n; i++) { v[i] = NAN; }
    return v;
  }
