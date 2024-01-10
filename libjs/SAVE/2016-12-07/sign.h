#ifndef sign_H
#define sign_H

/* Sign data type. */
/* Last edited on 2007-10-28 19:19:38 by stolfi */

#ifndef HAVE_SIGN
typedef enum { NEG=-1, ZER=0, POS=+1 } sign_t;
#define NEG NEG
#define ZER ZER
#define POS POS
#define HAVE_SIGN 1
#endif

#endif
