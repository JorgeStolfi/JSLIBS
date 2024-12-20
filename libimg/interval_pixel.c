/* See {interval_pixel.h}. */
/* Last edited on 2024-12-04 23:33:37 by stolfi */

#include <assert.h>
#include <limits.h>
#include <string.h>
#include <math.h>
 
#include <bool.h>
#include <affirm.h>
#include <interval.h>
#include <interval_pixel.h>

#define INF INFINITY

/* INTERNAL PROTOTYPES */

/* IMPLEMENTATIONS */

void ivpix_accum_pixel(int32_t chns, interval_t vs[], double wt, interval_t v[], double *wtotP)
  {
    if (wt != 0.0)
      { int32_t ich;
        for (ich = 0; ich < chns; ich++)
          { interval_t *ek = &(vs[ich]); 
            interval_t *vk = &(v[ich]);
            if (wt > 0.0)
              { LO(*vk) += wt*LO(*ek); HI(*vk) += wt*HI(*ek); }
            else 
              { LO(*vk) += wt*HI(*ek); HI(*vk) += wt*LO(*ek); }
          }
        (*wtotP) += wt;
      }
  }

void ivpix_make_pixel_undef(int32_t chns, interval_t v[])
  {
    int32_t ich;
    for (ich = 0; ich < chns; ich++) 
      { v[ich] = (interval_t){{ -INF, +INF }}; }
  }

float ivpix_floatize_interval(interval_t *v)
  {
    if (LO(*v) == +INF)
      { return +INF; }
    else if (HI(*v) == -INF)
      { return -INF; }
    else
      { double fval = 0.5*LO(*v) + 0.5*HI(*v);
        if (isnan(fval))
          { return 0.5; }
        else 
          { return (float)fval; }
      }
  }

void ivpix_scale_pixel(int32_t chns, double s, interval_t v[])
  {
    int32_t ich;
    for (ich = 0; ich < chns; ich++) 
      { interval_t *vk = &(v[ich]);
        if (s >= 0) 
          { LO(*vk) *= s; HI(*vk) *= s; }
        else
          { double t = LO(*vk); HI(*vk) = s*LO(*vk); LO(*vk) = s*t; }
      }
  }

void ivpix_debug_itv_pixel(char *label, double x, double y, int32_t chns, interval_t v[], char *tail)
  { 
    int32_t ich;
    fprintf(stderr, "    %s(%9.4f,%9.4f) = (", label, x, y);
    for (ich = 0; ich < chns; ich++) 
      { fprintf(stderr, " ");
        ivpix_print_interval(stderr, &(v[ich]), 7, 4);
      }
    fprintf(stderr, " )%s", tail);
  }
  
void ivpix_print_interval(FILE *wr, interval_t *v, int32_t width, int32_t prec)
  { fprintf(stderr, "[");
    ivpix_print_bound(stderr, LO(*v), width, prec);
    fprintf(stderr, " _ ");
    ivpix_print_bound(stderr, HI(*v), width, prec);
    fprintf(stderr, "]");
  }
  
void ivpix_print_bound(FILE *wr, double v, int32_t width, int32_t prec)
  { if (v == +INF) 
      { fprintf(stderr, "%*s", width, "+oo"); }
    else if (v == -INF) 
      { fprintf(stderr, "%*s", width, "-oo"); }
    else
      { fprintf(stderr, "%*.*f", width, prec, v); }
  }

