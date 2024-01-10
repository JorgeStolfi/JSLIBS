/* Setting rounding-mode direction */
/* Last edited on 2002-11-08 03:42:54 by stolfi */

#include <fptest.h>
#include <fenv.h>
#include <fpu_control.h>

void set_rounding (int dir)
{ 
  fesetround(dir);
}

int get_fsr(void)
{
  fpu_control_t cw;
  _FPU_GETCW(cw);
  return cw;
}
