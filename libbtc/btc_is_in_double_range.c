/* See {btc_is_in_double_range.h} */
/* Last edited on 2024-12-05 10:23:33 by stolfi */

#include <stdio.h>

#include <bool.h>
#include <affirm.h>

#include <btc_is_in_double_range.h>

    
bool_t btc_is_in_double_range(double vlo, double v, double vhi, int ib, char* vname, bool_t die)
  { if ((vlo > vhi) || (((v < vlo) || (v > vhi)) && die))
      { fprintf(stderr, "** bubble %d fields {%s}: %25.16e..(%25.16e)..%25.16e\n", ib, vname, vlo, v, vhi);
        demand((vlo <= vhi), "invalid parameter range");
        demand((! die) || ((vlo <= v) && (v <= vhi)), "nominal parameter is out of range");
      }
    return (vlo <= v) && (v <= vhi);
  }
