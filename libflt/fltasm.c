/* Setting rounding-mode direction */
/* Last edited on 2017-01-02 11:45:29 by jstolfi */

#define _GNU_SOURCE
#include <fpu_control.h>
#include <fenv.h>

#include <flt.h>

void flt_round_down (void)
  { 
    fesetround(FE_DOWNWARD);
  }
  
void flt_round_up   (void)
  { 
    fesetround(FE_UPWARD);
  }
  
void flt_round_near (void)
  { 
    fesetround(FE_TONEAREST);
  }
  
void flt_round_zero (void)
  { 
    fesetround(FE_TOWARDZERO);
  }

int flt_get_fsr(void)
  {
    fpu_control_t cw;
    _FPU_GETCW(cw);
    return cw;
  }

void flt_single_prec (void)
  {
    fpu_control_t cw;
    _FPU_GETCW(cw);
    cw = (fpu_control_t)((cw & ~_FPU_EXTENDED & ~_FPU_DOUBLE) | _FPU_SINGLE);
    _FPU_SETCW(cw);
  }

void flt_double_prec (void)
  {
    fpu_control_t cw;
    _FPU_GETCW(cw);
    cw = (fpu_control_t)((cw & ~_FPU_EXTENDED & ~_FPU_SINGLE) | _FPU_DOUBLE);
    _FPU_SETCW(cw);
  }

void flt_extended_prec (void)
  {
    fpu_control_t cw;
    _FPU_GETCW(cw);
    cw = (fpu_control_t)((cw & ~_FPU_DOUBLE & ~_FPU_SINGLE) | _FPU_EXTENDED);
    _FPU_SETCW(cw);
  }
