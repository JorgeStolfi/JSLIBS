/* See {nmsim_basic.h} */
/* Last edited on 2019-04-25 18:13:41 by jstolfi */

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
#include <jsrandom.h>

#include <nmsim_basic.h>

void nmsim_int64_range_clip
  ( int64_t alo,
    int64_t ahi, 
    int64_t blo, 
    int64_t bhi, 
    int64_t *cloP, 
    int64_t *chiP,
    char *itname
  )
  { /* Compute raw intersection: */
    int64_t clo = imax(alo, blo);
    int64_t chi = imin(ahi, bhi);
    /* Regularize empty range: */
    if (clo > chi) { clo = 0; chi = -1; }
    if (itname != NULL)
      { /* Print warnings, if any: */
        if (clo >chi)
          { fprintf(stderr, "!! %s range %ld..%ld is empty or fully outside %ld..%ld\n", itname, alo, ahi, blo, bhi); }
        else 
          { if (clo > alo)
              { fprintf(stderr, "!! %s subrange %ld..%ld is outside  %ld..%ld\n", itname, alo, clo-1, blo, bhi); }
            if (chi < ahi)
              { fprintf(stderr, "!! %s subrange %ld..%ld is outside  %ld..%ld\n", itname, chi+1, ahi, blo, bhi); }
          }
      }
    /* Return: */
    (*cloP) = clo; (*chiP) = chi;
  }
 
void nmsim_int32_range_clip
  ( int32_t alo,
    int32_t ahi, 
    int32_t blo, 
    int32_t bhi, 
    int32_t *cloP, 
    int32_t *chiP,
    char *itname
  )
  { int64_t clo, chi;
    nmsim_int64_range_clip(alo, ahi, blo, bhi, &clo, &chi, itname);
    assert((INT32_MIN <= clo) && (clo <= INT32_MAX));
    assert((INT32_MIN <= chi) && (chi <= INT32_MAX));
    (*cloP) = (int32_t)clo; (*chiP) = (int32_t)chi;
  }
 
void nmsim_int64_range_unite
  ( int64_t alo,
    int64_t ahi, 
    int64_t blo, 
    int64_t bhi, 
    int64_t *uloP, 
    int64_t *uhiP
  )
  { /* Compute raw union: */
    int64_t ulo, uhi;
    if (alo > ahi)
      { ulo = blo; uhi = bhi; }
    else if (blo > bhi)
      { ulo = alo; uhi = ahi; }
    else
      { ulo = imin(alo, blo); uhi = imax(ahi, bhi); }
    /* Regularize empty range: */
    if (ulo > uhi) { ulo = 0; uhi = -1; }
    /* Return: */
    (*uloP) = ulo; (*uhiP) = uhi;
  }

void nmsim_int32_range_unite
  ( int32_t alo,
    int32_t ahi, 
    int32_t blo, 
    int32_t bhi, 
    int32_t *uloP, 
    int32_t *uhiP
  )
  { int64_t ulo, uhi;
    nmsim_int64_range_unite(alo, ahi, blo, bhi, &ulo, &uhi);
    assert((INT32_MIN <= ulo) && (ulo <= INT32_MAX));
    assert((INT32_MIN <= uhi) && (uhi <= INT32_MAX));
    (*uloP) = (int32_t)ulo; (*uhiP) = (int32_t)uhi;
  }

double nmsim_throw_double(double vlo, double vhi)
  { double v = dabrandom(vlo, vhi);
    if (vlo != vhi)
      { /* Round value to nice number: */
        assert(vlo < vhi);
        double mod = nmsim_select_rounding_mod(vlo, vhi);
        double q = floor((v + mod/2)/mod);
        v = q*mod;
      }
    return v;
  }

double nmsim_select_rounding_mod(double vlo, double vhi)
  {
    double dif = fabs(vhi - vlo);
    /* Compute a power of 10 {mod} between {dif/100} and {dif/10): */
    double mod = 1.0;
    while (mod < dif/100) { mod = mod*10; }
    while (mod > dif/10) { mod = mod/10; }
    return mod;
  }

void nmsim_throw_time_range(nmsim_time_t tlo, nmsim_time_t thi, nmsim_time_t *tloP, nmsim_time_t *thiP)
  {
    if (tlo > thi) 
      { /* Empty interval: */
        (*tloP) = 0;
        (*thiP) = -1;
      }
    else
      { /* Choose min and max possible times: */
        nmsim_step_count_t dt = thi - tlo;  /* Num of steps in given range. */
        nmsim_step_count_t xt = dt / 2;     /* How much we should extrapolate. */
        if (xt < 10) { xt = 10; }
        nmsim_time_t tmin = (tlo >= nmsim_time_MIN + xt ? tlo - xt : nmsim_time_MIN);
        nmsim_time_t tmax = (thi <= nmsim_time_MAX - xt ? thi + xt : nmsim_time_MAX);
        /* Choose a random subrange of {tmin..tmax}: */
        nmsim_time_t tlo_pk = (nmsim_time_t)int64_abrandom(tmin-1, tmax);
        nmsim_time_t thi_pk = (nmsim_time_t)int64_abrandom(tmin, tmax);
        if (tlo_pk > thi_pk) { tlo_pk = tmax - (tlo_pk + 1 - tmin); thi_pk = tmax - (thi_pk - tmin); }
        assert((nmsim_time_MIN <= tlo_pk) && (tlo_pk <= thi_pk) && (thi_pk <= nmsim_time_MAX));
        (*tloP) = tlo_pk;
        (*thiP) = thi_pk;
      }
  }
        
