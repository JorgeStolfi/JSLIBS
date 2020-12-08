/* See hi2.h */
/* Last edited on 2013-11-21 02:55:04 by stolfilocal */ 

#define _GNU_SOURCE
#include <math.h>
#include <assert.h>

#include <i2.h>
#include <i3.h>
/* #include <i3x3.h> */
#include <affirm.h>
#include <jsmath.h>
#include <sign.h>
#include <sign_get.h>

#include <hi2.h>

#define NH 3
  /* Number of homogeneous coordinates in a point. */

#define NC 2
  /* Number of Cartesian coordinates in a point. */

hi2_point_t hi2_from_i2(i2_t *c)
  {
    return (hi2_point_t){{{1, c->c[0], c->c[1]}}};
  }

sign_t hi2_side(hi2_point_t *p, hi2_line_t *L)
  {
    int64_t dd = i3_dot(&(p->c), &(L->f));
    return sign_int64(dd);
  }

sign_t hi2_orient(hi2_point_t *p, hi2_point_t *q, hi2_point_t *r)
  {
    int64_t dd = i3_det(&(p->c), &(q->c), &(r->c));
    return sign_int64(dd);
  }

hi2_line_t hi2_join(hi2_point_t *p, hi2_point_t *q)
  {
    i3_t f;
    i3_cross(&(p->c), &(q->c), &f);
    return (hi2_line_t){f};
  }

hi2_point_t hi2_meet(hi2_line_t *K, hi2_line_t *L)
  {
    i3_t c;
    i3_cross(&(K->f), &(L->f), &c);
    return (hi2_point_t){c};
  }

i2_t hi2_line_dir(hi2_line_t *L)
  {
    int32_t dx = L->f.c[2];
    int32_t dy = -L->f.c[1];
    return (i2_t){{dx, dy}};
  }

i2_t hi2_line_normal(hi2_line_t *L)
  {
    int32_t nx = L->f.c[1];
    int32_t ny = L->f.c[2];
    return (i2_t){{nx, ny}};
  }

i2_t hi2_point_point_dir(hi2_point_t *frm, hi2_point_t *tto)
  {
    int32_t fw = frm->c.c[0];
    int32_t fx = frm->c.c[1];
    int32_t fy = frm->c.c[2];

    int32_t tw = tto->c.c[0];
    int32_t tx = tto->c.c[1];
    int32_t ty = tto->c.c[2];

    int32_t dx = (int32_t)(fw*(int64_t)tx - tw*(int64_t)fx);
    int32_t dy = (int32_t)(fw*(int64_t)ty - tw*(int64_t)fy);
    
    return (i2_t){{dx, dy}};
  }

urat64_t hi2_dist_sqr(hi2_point_t *a, hi2_point_t *b)
  {
    int32_t M0 = 38967; /* Max abs inputs to avoid overflow. */
  
    /* Range: {-M .. +M} */
    int32_t wa = a->c.c[0]; demand((-M0 <= wa) && (wa <= +M0), "input coord wa too big");
    int32_t xa = a->c.c[1]; demand((-M0 <= xa) && (xa <= +M0), "input coord xa too big");
    int32_t ya = a->c.c[2]; demand((-M0 <= ya) && (ya <= +M0), "input coord ya too big");
                                                                          
    int32_t wb = b->c.c[0]; demand((-M0 <= wb) && (wb <= +M0), "input coord wb too big");
    int32_t xb = b->c.c[1]; demand((-M0 <= xb) && (xb <= +M0), "input coord xb too big");
    int32_t yb = b->c.c[2]; demand((-M0 <= yb) && (yb <= +M0), "input coord yb too big");
                              
    int64_t dx = xa*(int64_t)wb - xb*(int64_t)wa;   /* {-2*M^2 .. +2*M^2} */
    int64_t dy = ya*(int64_t)wb - yb*(int64_t)wa;   /* {-2*M^2 .. +2*M^2} */
    int64_t dw = wa*(int64_t)wb;                    /* {-M^2 .. +M^2} */
    
    uint64_t num = dx*dx + dy*dy; /* {0 .. +8*M^4} */
    uint64_t den = dw*dw;         /* {0 .. +M^4} */
    
    urat64_t d = (urat64_t){ num, den };
    return d;
  }

sign_t hi2_in_circle(hi2_point_t *a, hi2_point_t *b, hi2_point_t *c, hi2_point_t *d)
  {
    int32_t M0 = 32767; /* Max abs inputs to avoid overflow. */

    /* Range: {-M .. +M} */
    int32_t wa = a->c.c[0]; demand((-M0 <= wa) && (wa <= +M0), "input coord wa too big");
    int32_t xa = a->c.c[1]; demand((-M0 <= xa) && (xa <= +M0), "input coord xa too big");
    int32_t ya = a->c.c[2]; demand((-M0 <= ya) && (ya <= +M0), "input coord ya too big");
                                                                          
    int32_t wb = b->c.c[0]; demand((-M0 <= wb) && (wb <= +M0), "input coord wb too big");
    int32_t xb = b->c.c[1]; demand((-M0 <= xb) && (xb <= +M0), "input coord xb too big");
    int32_t yb = b->c.c[2]; demand((-M0 <= yb) && (yb <= +M0), "input coord yb too big");
                                                                          
    int32_t wc = c->c.c[0]; demand((-M0 <= wc) && (wc <= +M0), "input coord wc too big");
    int32_t xc = c->c.c[1]; demand((-M0 <= xc) && (xc <= +M0), "input coord xc too big");
    int32_t yc = c->c.c[2]; demand((-M0 <= yc) && (yc <= +M0), "input coord yc too big");
                                                                          
    int32_t wd = d->c.c[0]; demand((-M0 <= wd) && (wd <= +M0), "input coord wd too big");
    int32_t xd = d->c.c[1]; demand((-M0 <= xd) && (xd <= +M0), "input coord xd too big");
    int32_t yd = d->c.c[2]; demand((-M0 <= yd) && (yd <= +M0), "input coord yd too big");
    
    /* Range: {-2*M^2 .. +2*M^2} */
    int32_t M1 = 2147483647;

    int32_t yda = yd*wa - ya*wd; assert((-M1 <= yda) && (yda <= +M1));
    int32_t xda = xd*wa - xa*wd; assert((-M1 <= xda) && (xda <= +M1));
    int32_t xbc = xb*wc - xc*wb; assert((-M1 <= xbc) && (xbc <= +M1));
    int32_t ybc = yb*wc - yc*wb; assert((-M1 <= ybc) && (ybc <= +M1));
                                                             
    int32_t ydc = yd*wc - yc*wd; assert((-M1 <= ydc) && (ydc <= +M1));
    int32_t xdc = xd*wc - xc*wd; assert((-M1 <= xdc) && (xdc <= +M1));
    int32_t xba = xb*wa - xa*wb; assert((-M1 <= xba) && (xba <= +M1));
    int32_t yba = yb*wa - ya*wb; assert((-M1 <= yba) && (yba <= +M1));
    
    /* Range: {-8*M^4 .. +8*M^4} */
    int64_t uda_bc = yda*(int64_t)xbc + xda*(int64_t)ybc;
    int64_t udc_ba = ydc*(int64_t)xba + xdc*(int64_t)yba;
    int64_t vda_bc = xda*(int64_t)xbc - yda*(int64_t)ybc;
    int64_t vdc_ba = xdc*(int64_t)xba - ydc*(int64_t)yba;
    
    /* Range: {-64*M^8 .. +64*M^8} */
    int64_t A_1; uint64_t A_0; int64_mul(uda_bc, vdc_ba, &A_1, &A_0); /* Will not overflow. */
    int64_t B_1; uint64_t B_0; int64_mul(udc_ba, vda_bc, &B_1, &B_0); /* Will not overflow. */
    
    if (A_1 > B_1)
      { return +1; }
    else if (A_1 < B_1)
      { return -1; }
    else if (A_0 > B_0)
      { return +1; }
    else if (A_0 < B_0)
      { return -1; }
    else
      { return 0; }
  }

// bool_t hi2_pmap_is_identity(hi2_pmap_t *m)
//   {
//     return i3x3_is_unif_scaling(&(m->dir), m->dir.c[0][0]); 
//   }
// 
// hi2_point_t hi2_map_point(hi2_point_t *p, hi2_pmap_t *m)
//   {
//     i3_t q; 
//     i3x3_map_row(&(p->c), &(m->dir), &q);
//     return (hi2_point_t){q};
//   }
// 
// hi2_line_t hi2_map_line(hi2_line_t *L, hi2_pmap_t *m)
//   {
//     i3_t f;
//     i3x3_map_col(&(m->inv), &(L->f), &f);
//     return (hi2_line_t){f};
//   }
// 
// hi2_pmap_t hi2_comp_map(hi2_pmap_t *m, hi2_pmap_t *n)
//   {
//     i3x3_t dir, inv;
//     i3x3_mul(&(m->dir), &(n->dir), &dir);
//     i3x3_mul(&(n->inv), &(m->inv), &inv);
//     return (hi2_pmap_t){dir, inv};
//   }
// 
// hi2_pmap_t hi2_inv_map(hi2_pmap_t *m)
//   {
//     return (hi2_pmap_t){m->inv, m->dir};
//   }

// hi2_pmap_t hi2_pmap_from_points(hi2_point_t *p, hi2_point_t *q, hi2_point_t *r, hi2_point_t *u)
//   {
//     int i;
//     hi2_pmap_t m; /* The resulting map. */
//     
//     /* Compute weights {(a,b,c)=w.c[0..2]} for rows of {Q} that map {[1,1,1]} to {u}: */
//     i3_t w;
//     { /* Compute a matrix {Q} that maps the cardinal points to {p,q,r} as given: */
//       i3x3_t Q;
//       Q.c[0][0] = p->c.c[0]; Q.c[0][1] = p->c.c[1]; Q.c[0][2] = p->c.c[2];
//       Q.c[1][0] = q->c.c[0]; Q.c[1][1] = q->c.c[1]; Q.c[1][2] = q->c.c[2];
//       Q.c[2][0] = r->c.c[0]; Q.c[2][1] = r->c.c[1]; Q.c[2][2] = r->c.c[2];
//       i3x3_inv(&Q, &Q);
//       w.c[0] = u->c.c[0];  w.c[1] = u->c.c[1];  w.c[2] = u->c.c[2];
//       i3x3_map_row(&w, &Q, &w);
//     }
//     
//     /* Make the weights positive, so that {p,q,r} are strictly honored: */
//     for (i = 0; i < NH; i++) { w.c[i] = fabs(w.c[i]); }
// 
//     /* Ensure that {m.dir} maps the cardinal points to {p,q,r} and some unit point to {u}: */
//     i3x3_t *P = &(m.dir);
//     P->c[0][0] = w.c[0]*p->c.c[0]; P->c[0][1] = w.c[0]*p->c.c[1]; P->c[0][2] = w.c[0]*p->c.c[2];
//     P->c[1][0] = w.c[1]*q->c.c[0]; P->c[1][1] = w.c[1]*q->c.c[1]; P->c[1][2] = w.c[1]*q->c.c[2];
//     P->c[2][0] = w.c[2]*r->c.c[0]; P->c[2][1] = w.c[2]*r->c.c[1]; P->c[2][2] = w.c[2]*r->c.c[2];
// 
//     /* Compute the inverse map: */
//     i3x3_inv(&(m.dir), &(m.inv));
// 
//     return m;
// 
//   }

