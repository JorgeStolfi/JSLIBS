/* basic floating-point definitions - Sun SPARC4/Solaris 5 */
/* Last edited on 2007-11-01 20:18:13 by stolfi */

#ifndef fptest_H
#define fptest_H

#include <fpu_control.h>

void set_rounding (int dir);
int get_fsr (void); 

#endif
