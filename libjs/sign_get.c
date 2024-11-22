/* See {sign_get.h}. */
/* Last edited on 2024-11-15 19:15:18 by stolfi */

#include <stdint.h>

#include <sign_get.h>
#include <sign.h>

sign_t sign_int(int x)
  { if (x < 0)
      { return -1; }
    else if (x > 0)
      { return +1; }
    else
      { return 0; }
  }

sign_t sign_long_int(long int x)
  { if (x < 0)
      { return -1; }
    else if (x > 0)
      { return +1; }
    else
      { return 0; }
  }

sign_t sign_int32(int32_t x)
  { if (x < 0)
      { return -1; }
    else if (x > 0)
      { return +1; }
    else
      { return 0; }
  }

sign_t sign_int64(int64_t x)
  { if (x < 0)
      { return -1; }
    else if (x > 0)
      { return +1; }
    else
      { return 0; }
  }

sign_t sign_float(float x)
  { if (x < 0)
      { return -1; }
    else if (x > 0)
      { return +1; }
    else
      { return 0; }
  }

sign_t sign_double(double x)
  { if (x < 0)
      { return -1; }
    else if (x > 0)
      { return +1; }
    else
      { return 0; }
  }


