/* See {test_r2_opt_basic.h}. */
/* Last edited on 2025-03-19 14:33:39 by stolfi */

#include <stdint.h>
#include <stdio.h>
#include <math.h>

#include <r2.h>
#include <bool.h>

#include <test_r2_opt_basic.h>

void tr2o_debug_points
  ( uint32_t indent,
    char *title,
    uint32_t NI, 
    char *pname,
    r2_t p[], 
    r2_t pini[], 
    r2_t popt[], 
    r2_t arad[],
    r2_t astp[],
    double f2p, 
    double b2p
  )
  { 
    auto void show_disp(char *qtag, r2_t *q, r2_t *o, r2_t *r, r2_t *s);
      /* Shows the difference {q} to {o}, absolute and relative to the 
        search radius {r} (if {r} is not {NULL}) and/or
        the search step {s} (if {s} is not {NULL}). */

    auto void show_rel_disp(char *qtag, r2_t *d, char *rtag, r2_t *r);
      /* If {r} is not {NULL} and not {(0,0)}, Shows the displacement
        {d} from the point named {qtag} as a fraction of the vector {r}
        named {rtag}. */
       
    fprintf(stderr, "%*s%s\n", indent, "", title);
    for (uint32_t i = 0;  i < NI; i++)
      { fprintf(stderr, "%*s  %s[%02d] =       ( %9.4f %9.4f )\n", indent, "", pname, i, p[i].c[0], p[i].c[1]);
        r2_t *r = (arad != NULL ? &(arad[i]) : NULL);
        r2_t *s = (astp != NULL ? &(astp[i]) : NULL);
        if (pini != NULL) { show_disp("ini", &(p[i]), &(pini[i]), r, s); }
        if (popt != NULL) { show_disp("opt", &(p[i]), &(popt[i]), r, s); }
        fprintf(stderr, "\n");
      }
    if ((! isnan(f2p)) || (! isnan(b2p)))
      { fprintf(stderr, "%*s", indent, "");
        if (! isnan(f2p)) { fprintf(stderr, "  f2(p) = %22.10f", f2p); }
        if (! isnan(b2p)) { fprintf(stderr, "  b2(p) = %22.10f", b2p); }
        if ((! isnan(f2p)) && (! isnan(b2p))) { fprintf(stderr, "  f2_full(p) = %22.10f\n", f2p+b2p); }
      }
    fprintf(stderr, "\n");
    return;

    void show_disp(char *qtag, r2_t *q, r2_t *o, r2_t *r, r2_t *s) 
      { r2_t d; r2_sub(q, o, &d);
        fprintf(stderr, "%*s           = %s + ( %+9.4f %+9.4f )\n", indent, "", qtag, d.c[0], d.c[1]);
        show_rel_disp(qtag, &d, "rad", r);
        show_rel_disp(qtag, &d, "stp", s);
      }
      
    void show_rel_disp(char *qtag, r2_t *d, char *rtag, r2_t *r)
      { if (r != NULL)
          { double rx = r->c[0];
            double ry = r->c[1];
            if ((rx != 0.0) || (ry != 0.0))
              { double dx = (rx == 0 ? 0.0 : d->c[0]/rx);
                double dy = (ry == 0 ? 0.0 : d->c[1]/ry);
                fprintf(stderr, "%*s           = %s + ( %+9.4f %+9.4f )*%s\n", indent, "", qtag, dx, dy, rtag);
              }
          }
      }
   }

void tr2o_debug_r2
  ( uint32_t indent,
    char *pname,
    uint32_t i,
    r2_t *p
  )
  { 
    fprintf(stderr, "%*s %s[%02d] = ( %20.15f %20.15f )\n", indent, "", pname, i, p->c[0], p->c[1]); 
  }

void tr2o_debug_params
  ( uint32_t indent,
    char *pname,
    uint32_t NI, 
    r2_t p[]
  )
  { for (uint32_t i = 0;  i < NI; i++) { tr2o_debug_r2(indent, pname, i, &(p[i])); }
    fprintf(stderr, "\n");
    return;
  }
  
