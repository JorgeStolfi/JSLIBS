/* See {btc_is_trivial_double_range.h} */
/* Last edited on 2024-12-05 10:23:37 by stolfi */

#include <stdio.h>

#include <bool.h>
#include <affirm.h>

#include <btc_is_trivial_double_range.h>

bool_t btc_is_trivial_double_range(double vlo, double v, double vhi, int ib, char* vname, bool_t die)
  { if ((vlo > v) || (v > vhi) || ((vlo != vhi) && die))
      { fprintf(stderr, "** bubble %d fields {%s}: [%25.16e _ (%25.16e) _ %25.16e]\n", ib, vname, vlo, v, vhi);
        demand((vlo <= vhi), "invalid parameter range");
        demand((vlo <= v) && (v <= vhi), "nominal parameter is out of range");
        demand((! die) || (vlo == vhi), "this parameter cannot be adjusted");
      }
    return (vlo == vhi);
  }
