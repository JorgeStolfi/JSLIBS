/* See {nmsim_firing_func.h} */
/* Last edited on 2019-02-13 16:40:46 by jstolfi */

#define _GNU_SOURCE
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include <assert.h>
#include <math.h>

#include <affirm.h>
#include <bool.h>
#include <nmsim_firing_func_gauss.h>
#include <nmsim_firing_func_lramp.h>

#include <nmsim_firing_func.h>

nmsim_firing_func_t nmsim_firing_func_make
  ( nmsim_firing_func_class_t class,
    double V_M,
    double V_D
  )
  { nmsim_firing_func_t Phi = (nmsim_firing_func_t) 
      { .class = class, .V_M = V_M, .V_D = V_D };
    return Phi;
  }

void nmsim_firing_func_eval
  ( nmsim_firing_func_t *Phi,
    double V, 
    double *prP,
    double *dprP
  )
  {
    if (Phi->class == 'G') 
      { nmsim_firing_func_gauss_eval(Phi->V_M, Phi->V_D, V, prP, dprP); }
    else if (Phi->class == 'L') 
      { nmsim_firing_func_lramp_eval(Phi->V_M, Phi->V_D, V, prP, dprP); }
    else 
      { demand(FALSE, "invalid firing function class"); }
  }
 
double nmsim_firing_func_eval_inv(nmsim_firing_func_t *Phi, double pr)
  {
    demand((pr >= 0) && (pr <= 1), "invalid probability");
    if (Phi->class == 'G') 
      { return nmsim_firing_func_gauss_eval_inv(Phi->V_M, Phi->V_D, pr); }
    else if (Phi->class == 'L') 
      { return nmsim_firing_func_lramp_eval_inv(Phi->V_M, Phi->V_D, pr); }
    else 
      { demand(FALSE, "invalid firing function class"); }
  }

double nmsim_firing_func_eval_gauss
  ( nmsim_firing_func_t *Phi,
    double V_avg, 
    double V_dev
  )
  {
    double pr;
    if (Phi->class == 'G') 
      { double V_D = hypot(V_dev, Phi->V_D);
        nmsim_firing_func_gauss_eval(Phi->V_M, V_D, V_avg, &pr, NULL);
       }
    else 
      { demand(FALSE, "not implemented for this class");
        pr = NAN;
      }
    return pr;
  }
