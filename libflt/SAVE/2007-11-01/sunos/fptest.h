/* basic floating-point definitions - Sun SPARC4/Solaris 5 */
/* Last edited on 2007-11-01 20:18:34 by stolfi */

#ifndef fptest_H
#define fptest_H

#include <sys/ieeefp.h>

void set_rounding (enum fp_direction_type dir);
int get_fsr (void); 

#endif
