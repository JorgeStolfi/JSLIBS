#ifndef bool_H
#define bool_H

/* Boolean data type. */
/* Last edited on 2005-02-14 10:30:48 by stolfi */

#ifndef HAVE_BOOL
typedef enum {FALSE = 0, TRUE = 1} bool_t;
#define bool_t bool_t
#define FALSE FALSE
#define TRUE TRUE
#define HAVE_BOOL 1
#endif

#endif
