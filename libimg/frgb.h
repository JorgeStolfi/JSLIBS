#ifndef frgb_H
#define frgb_H

/* frgb.h - floating-point RGB color data type */
/* Last edited on 2024-12-05 10:30:18 by stolfi */

#include <stdint.h>
#include <limits.h>

/* Channels in a color image: */
#define frgb_CHANNELS 3

typedef struct frgb_t { float c[frgb_CHANNELS]; } frgb_t;
  /* An RGB triplet. */

#ifndef INF
#define INF INFINITY 
#endif

#define frgb_Zeros   (frgb_t){{0.0, 0.0, 0.0}}
#define frgb_Black   (frgb_t){{0.0, 0.0, 0.0}}
#define frgb_White   (frgb_t){{1.0, 1.0, 1.0}}
#define frgb_Ones    (frgb_t){{1.0, 1.0, 1.0}}
#define frgb_NoColor (frgb_t){{NAN, NAN, NAN}}

#endif
