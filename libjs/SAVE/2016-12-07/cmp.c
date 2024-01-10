/* See {cmp.h}. */
/* Last edited on 2007-10-28 21:33:29 by stolfi */

#include <cmp.h>
#include <sign.h>
#include <stdint.h>

sign_t cmp_int(const int* x, const int *y)
  { if ((*x) < (*y))
      { return -1; }
    else if ((*x) > (*y))
      { return +1; }
    else
      { return 0; }
  }

sign_t cmp_long_int(const long int *x, const long int *y)
  { if ((*x) < (*y))
      { return -1; }
    else if ((*x) > (*y))
      { return +1; }
    else
      { return 0; }
  }

sign_t cmp_unsigned_int(const unsigned int* x, const unsigned int *y)
  { if ((*x) < (*y))
      { return -1; }
    else if ((*x) > (*y))
      { return +1; }
    else
      { return 0; }
  }

sign_t cmp_unsigned_long_int(const unsigned long int *x, const unsigned long int *y)
  { if ((*x) < (*y))
      { return -1; }
    else if ((*x) > (*y))
      { return +1; }
    else
      { return 0; }
  }

sign_t cmp_int8(const int8_t *x, const int8_t *y)
  { if ((*x) < (*y))
      { return -1; }
    else if ((*x) > (*y))
      { return +1; }
    else
      { return 0; }
  }

sign_t cmp_int16(const int16_t *x, const int16_t *y)
  { if ((*x) < (*y))
      { return -1; }
    else if ((*x) > (*y))
      { return +1; }
    else
      { return 0; }
  }

sign_t cmp_int32(const int32_t *x, const int32_t *y)
  { if ((*x) < (*y))
      { return -1; }
    else if ((*x) > (*y))
      { return +1; }
    else
      { return 0; }
  }

sign_t cmp_int64(const int64_t *x, const int64_t *y)
  { if ((*x) < (*y))
      { return -1; }
    else if ((*x) > (*y))
      { return +1; }
    else
      { return 0; }
  }

sign_t cmp_uint8(const uint8_t *x, const uint8_t *y)
  { if ((*x) < (*y))
      { return -1; }
    else if ((*x) > (*y))
      { return +1; }
    else
      { return 0; }
  }

sign_t cmp_uint16(const uint16_t *x, const uint16_t *y)
  { if ((*x) < (*y))
      { return -1; }
    else if ((*x) > (*y))
      { return +1; }
    else
      { return 0; }
  }

sign_t cmp_uint32(const uint32_t *x, const uint32_t *y)
  { if ((*x) < (*y))
      { return -1; }
    else if ((*x) > (*y))
      { return +1; }
    else
      { return 0; }
  }

sign_t cmp_uint64(const uint64_t *x, const uint64_t *y)
  { if ((*x) < (*y))
      { return -1; }
    else if ((*x) > (*y))
      { return +1; }
    else
      { return 0; }
  }

sign_t cmp_float(const float *x, const float *y)
  { if ((*x) < (*y))
      { return -1; }
    else if ((*x) > (*y))
      { return +1; }
    else
      { return 0; }
  }

sign_t cmp_double(const double *x, const double *y)
  { if ((*x) < (*y))
      { return -1; }
    else if ((*x) > (*y))
      { return +1; }
    else
      { return 0; }
  }


