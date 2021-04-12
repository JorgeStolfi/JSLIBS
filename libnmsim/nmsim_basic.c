/* see {nmsim_basic.h} */
/* Last edited on 2020-12-24 22:23:19 by jstolfi */

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

double nmsim_basic_tau_from_mu(double mu, double timeStep)
  { if (mu == 0)
      { return 0; }
    else if (mu == 1.0)
      { return INF; }
    else
      { return -timeStep/log(mu); }
  }
  
double nmsim_basic_mu_from_tau(double tau, double timeStep)
  { if (tau < 0.02*timeStep)
      { return 0.0; }
    else if (tau == INF)
      { return 1.0; }
    else
      { return exp(-timeStep/tau); } 
  }

void nmsim_int64_range_clip
  ( int64_t aLo,
    int64_t aHi, 
    int64_t bLo, 
    int64_t bHi, 
    int64_t *cLoP, 
    int64_t *cHiP,
    char *itname
  )
  { /* Compute raw intersection: */
    int64_t cLo = imax(aLo, bLo);
    int64_t cHi = imin(aHi, bHi);
    /* Regularize empty range: */
    if (cLo > cHi) { cLo = 0; cHi = -1; }
    if (itname != NULL)
      { /* Print warnings, if any: */
        if (cLo >cHi)
          { fprintf(stderr, "!! %s range %ld..%ld is empty or fully outside %ld..%ld\n", itname, aLo, aHi, bLo, bHi); }
        else 
          { if (cLo > aLo)
              { fprintf(stderr, "!! %s subrange %ld..%ld is outside  %ld..%ld\n", itname, aLo, cLo-1, bLo, bHi); }
            if (cHi < aHi)
              { fprintf(stderr, "!! %s subrange %ld..%ld is outside  %ld..%ld\n", itname, cHi+1, aHi, bLo, bHi); }
          }
      }
    /* Return: */
    (*cLoP) = cLo; (*cHiP) = cHi;
  }
 
void nmsim_int32_range_clip
  ( int32_t aLo,
    int32_t aHi, 
    int32_t bLo, 
    int32_t bHi, 
    int32_t *cLoP, 
    int32_t *cHiP,
    char *itname
  )
  { int64_t cLo, cHi;
    nmsim_int64_range_clip(aLo, aHi, bLo, bHi, &cLo, &cHi, itname);
    assert((INT32_MIN <= cLo) && (cLo <= INT32_MAX));
    assert((INT32_MIN <= cHi) && (cHi <= INT32_MAX));
    (*cLoP) = (int32_t)cLo; (*cHiP) = (int32_t)cHi;
  }
 
void nmsim_int64_range_unite
  ( int64_t aLo,
    int64_t aHi, 
    int64_t bLo, 
    int64_t bHi, 
    int64_t *uLoP, 
    int64_t *uHiP
  )
  { /* Compute raw union: */
    int64_t uLo, uHi;
    if (aLo > aHi)
      { uLo = bLo; uHi = bHi; }
    else if (bLo > bHi)
      { uLo = aLo; uHi = aHi; }
    else
      { uLo = imin(aLo, bLo); uHi = imax(aHi, bHi); }
    /* Regularize empty range: */
    if (uLo > uHi) { uLo = 0; uHi = -1; }
    /* Return: */
    (*uLoP) = uLo; (*uHiP) = uHi;
  }

void nmsim_int32_range_unite
  ( int32_t aLo,
    int32_t aHi, 
    int32_t bLo, 
    int32_t bHi, 
    int32_t *uLoP, 
    int32_t *uHiP
  )
  { int64_t uLo, uHi;
    nmsim_int64_range_unite(aLo, aHi, bLo, bHi, &uLo, &uHi);
    assert((INT32_MIN <= uLo) && (uLo <= INT32_MAX));
    assert((INT32_MIN <= uHi) && (uHi <= INT32_MAX));
    (*uLoP) = (int32_t)uLo; (*uHiP) = (int32_t)uHi;
  }

double nmsim_throw_double(double vLo, double vHi)
  { double v = dabrandom(vLo, vHi);
    if (vLo != vHi)
      { /* Round value to nice number: */
        assert(vLo < vHi);
        double mod = nmsim_select_rounding_mod(vLo, vHi);
        double q = floor((v + mod/2)/mod);
        v = q*mod;
      }
    return v;
  }

double nmsim_select_rounding_mod(double vLo, double vHi)
  {
    double dif = fabs(vHi - vLo);
    /* Compute a power of 10 {mod} between {dif/100} and {dif/10): */
    double mod = 1.0;
    while (mod < dif/100) { mod = mod*10; }
    while (mod > dif/10) { mod = mod/10; }
    return mod;
  }

void nmsim_throw_time_range(nmsim_time_t tLo, nmsim_time_t tHi, nmsim_time_t *tLoP, nmsim_time_t *tHiP)
  {
    if (tLo > tHi) 
      { /* Empty interval: */
        (*tLoP) = 0;
        (*tHiP) = -1;
      }
    else
      { /* Choose min and max possible times: */
        nmsim_step_count_t dt = tHi - tLo;  /* Num of steps in given range. */
        nmsim_step_count_t xt = dt / 2;     /* How much we should extrapolate. */
        if (xt < 10) { xt = 10; }
        nmsim_time_t tmin = (tLo >= nmsim_time_MIN + xt ? tLo - xt : nmsim_time_MIN);
        nmsim_time_t tmax = (tHi <= nmsim_time_MAX - xt ? tHi + xt : nmsim_time_MAX);
        /* Choose a random subrange of {tmin..tmax}: */
        nmsim_time_t tLo_pk = (nmsim_time_t)int64_abrandom(tmin-1, tmax);
        nmsim_time_t tHi_pk = (nmsim_time_t)int64_abrandom(tmin, tmax);
        if (tLo_pk > tHi_pk) { tLo_pk = tmax - (tLo_pk + 1 - tmin); tHi_pk = tmax - (tHi_pk - tmin); }
        assert((nmsim_time_MIN <= tLo_pk) && (tLo_pk <= tHi_pk) && (tHi_pk <= nmsim_time_MAX));
        (*tLoP) = tLo_pk;
        (*tHiP) = tHi_pk;
      }
  }
