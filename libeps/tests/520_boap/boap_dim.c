/* See {boap_dim.h} */

#define _GNU_SOURCE
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <stdint.h>
#include <math.h>
#include <assert.h>

#include <bool.h>
#include <sign.h>
#include <affirm.h>
#include <vec.h>
#include <r3.h>
#include <epswr.h>
#include <epswr_dim.h>
#include <jsfile.h>

#include <boap_dim.h>

boap_dim_t *boap_dim_new
  ( double Xa, double Ya, double Za,
    double Xb, double Yb, double Zb,
    double agap,
    double bgap, 
    double elen,
    bool_t inner,
    int8_t decimals
  )
  {
    boap_dim_t *dim = notnull(malloc(sizeof(boap_dim_t)), "no mem");
    
    dim->a = (r3_t){{ Xa, Ya, Za }};
    dim->b = (r3_t){{ Xb, Yb, Zb }};
    
    dim->agap = agap;
    dim->bgap = bgap;
    dim->elen = elen;
    dim->inner = inner;
    dim->decimals = decimals;

    return dim;
  }
  
void boap_dim_draw (epswr_figure_t *epsf,  boap_dim_t *dim, int8_t uax)
  {
    /* Project points {A,B} on axes perpendicular to {vax}: */
    int8_t vax = (int8_t)((uax + 1) % 3);
    int8_t wax = (int8_t)((vax + 1) % 3);
    double xa = dim->a.c[vax], ya = dim->a.c[wax];
    double xb = dim->b.c[vax], yb = dim->b.c[wax];
    if (hypot(xb-xa, yb-ya) < 0.5) { return; }
    
    sign_t dir = ( dim->elen < 0 ? -1 : +1); /* Direction of extension lines rel {a-->b}. */
    double dpos = 1.0;
    double dlen = (dim->inner ? 100000 : 6.0);
    double hoff = 0.0;
    double voff = dir*(dim->inner ? 2.5 : 3.0);

    double dab, xr, yr, rot;
    epswr_set_fill_color(epsf, 0.000, 0.000, 0.000);
    epswr_set_pen(epsf, 0.000,0.000,0.000, 0.25, 0.0,0.0);
    epswr_dim_linear
      ( epsf, xa, ya, xb, yb, &dab,
        dim->agap, dim->bgap, dim->elen, 
        dim->inner, dpos, dlen, 
        hoff, voff, &xr, &yr, &rot
      );
    char *label = NULL;
    asprintf(&label, "%.*f", dim->decimals, dab);
    epswr_label(epsf, label, xr, yr, rot, FALSE, 0.5, 0.5, TRUE, FALSE);
    free(label);
  }

void boap_dim_print(FILE *wr, boap_dim_t *dim)
  {
    fprintf(wr, "DIM");
    fprintf(wr, " dist = %7.2f", r3_dist(&(dim->a), &(dim->b)));
    fprintf(wr, " a = (%7.1f,%7.1f,%7.1f) gap %7.1f", dim->a.c[0], dim->a.c[1],dim->a.c[2], dim->agap);
    fprintf(wr, " b = (%7.1f,%7.1f,%7.1f) gap %7.1f", dim->b.c[0], dim->b.c[1],dim->b.c[2], dim->bgap);
    fprintf(wr, " elen %7.1f", dim->elen);
    fprintf(wr, " inner = %c", "FT"[dim->inner]);
    fprintf(wr, " dec %d", dim->decimals);
    fprintf(wr, "\n");
  }

vec_typeimpl(boap_dim_ref_vec_t, boap_dim_ref_vec, boap_dim_ref_t);
