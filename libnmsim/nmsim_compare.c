/* See {nmsim_basic.h} */
/* Last edited on 2019-07-04 17:37:42 by jstolfi */

#define _GNU_SOURCE
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <stdint.h>
#include <assert.h>
#include <math.h>

#include <affirm.h>
#include <jsmath.h>
#include <jsrandom.h>

#include <nmsim_basic.h>
#include <nmsim_write.h>
 
#include <nmsim_compare.h>

bool_t nmsim_compare_double_param_special
  ( double dvr, 
    double dve,
    double prec
  );
  /* Comparison of a read value {vr} and an expected value {ve} when
    there is with some special value {vz} (usually in {-1,0,+1}) that
    must not be rounded to. Expects to receive {dvr = vr-vz} and {dve
    = ve-vz}.
    
    Requires equality if eithher {ve} or {vr} are equal to {vz}.
    Otherwise, {fabs(dvr - dvz)} should be less than {prec}, unless
    {ve} is too close to {vz}, in which case {vr} may have been
    rounded away from {vz}, resulting in an error up to {2*prec}. */

void nmsim_compare_int64_param(char *name, int64_t vr, int64_t ve)
  { if (vr != ve)
      { fprintf(stderr, "** parameter %s mismatch", name);
        fprintf(stderr, " read = %ld  expected = %ld", vr, ve);
        demand(FALSE, "aborted");
     }
  }

void nmsim_compare_double_param
  ( char *name, 
    double vr, 
    double ve,
    double prec,
    bool_t special_0, 
    bool_t special_1
  )
  { 
    bool_t ok; /* Comparison succeeded. */
      if ((isnan(vr)) || (isnan(ve)))
      { /* Both must be {NAN}: */
        ok = (isnan(vr) & isnan(ve));
      }
    else if ((fabs(vr) == INF) || (fabs(ve) == INF))
      { /* Require identity: */
        ok = (vr == ve);
      }
    else if ((fabs(ve) < 0.5) && (special_0))
      { 
        /* May have rounded {vr} away from 0, if {ve} was too small: */
        ok = nmsim_compare_double_param_special(vr, ve, prec);
      }
    else if ((fabs(ve) > 0.5) && (special_1))
      { 
        /* May have rounded {vr} away from {+1} or {-1}, if {ve} was too close to that: */
        double one = (vr < 0.0 ? -1.0 : +1.0);
        ok = nmsim_compare_double_param_special(vr-one, ve-one, prec);
      }
    else
      { /* Require approximate equality: */
        ok = (fabs(vr - ve) <= prec);
      }
    if (! ok)
      { fprintf(stderr, "** parameter %s mismatch", name);
        fprintf(stderr, " read = %23.15e  expected = %23.15e", vr, ve);
        fprintf(stderr, " diff = %.15f\n", vr - ve);
        demand(FALSE, "aborted");
     }
  }

bool_t nmsim_compare_double_param_special
  ( double dvr, 
    double dve,
    double prec
  )
  {
    bool_t ok;
    if ((dvr == 0.0) || (dve == 0.0))
      { /* Require identity: */
        ok = (dvr == dve);
      }
    else if (fabs(dve) < prec)
      { /* May have rounded away from special value, on the same side of it: */
        ok = (dve*dvr > 0) && (fabs(dvr) < 2*prec);
      }
    else
      { /* Require approximate equality: */
        ok = (fabs(dvr - dve) <= prec);
      }
    return ok;
 }
 
