/* See hr3.h */
/* Last edited on 2024-11-20 12:03:09 by stolfi */ 

#include <stdint.h>
#include <math.h>
#include <assert.h>
#include <stdio.h>

#include <sign_get.h>
#include <sign.h>
#include <affirm.h>

#include <r4x4.h>
#include <r3x3.h>
#include <r6.h>
#include <r4.h>
#include <r3.h>

#include <hr3.h>

/* Based on HR3.m3, created 1993-05-18 by Marcos C. Carrard. */
/* Based on H3.pas by J. Stolfi. */

#define NH 4
  /* Number of homogeneous coordinates in a point or coefficients in a plane. */

#define NG 6
  /* Number of Grassman coordinates of a line. */

#define NC 3
  /* Number of Cartesian coordinates in a point. */

r3_t r3_from_hr3_nan(hr3_point_t *p)
  {
    double w = p->c.c[0];
    double m = fmax(fmax(fabs(p->c.c[1]), fabs(p->c.c[2])), fabs(p->c.c[3]));
    r3_t c;
    if (fabs(w) <= m*1e-200) 
      { c = (r3_t){{ NAN, NAN, NAN }}; }
    else
      { c = (r3_t){{ p->c.c[1]/w, p->c.c[2]/w, p->c.c[3]/w }}; } 
    return c;
  }

hr3_point_t hr3_from_r3(r3_t *c)
  {
    return (hr3_point_t){{{1.0, c->c[0], c->c[1], c->c[2]}}};
  }
    
r3_t r3_from_hr3(hr3_point_t *p)
  {
    double w = p->c.c[0];
    demand(w != 0.0, "point at infinity");
    return (r3_t){{ p->c.c[1]/w, p->c.c[2]/w, p->c.c[3]/w }};
  }

hr3_point_t hr3_point_at_infinity(r3_t *dir)
  { hr3_point_t p;
    p.c.c[0] = 0.0;
    p.c.c[1] = dir->c[0];
    p.c.c[2] = dir->c[1];
    p.c.c[3] = dir->c[2];
    return p;
  }

double hr3_pt_pt_diff(hr3_point_t *p, hr3_point_t *q)
  {
    return r4_angle(&(p->c), &(q->c));
  }

sign_t hr3_side(hr3_point_t *p, hr3_plane_t *Q)
  {
    double dd = r4_dot(&(p->c), &(Q->f));
    return sign_double(dd);
  }       

sign_t hr3_orient(hr3_point_t *p, hr3_point_t *q, hr3_point_t *r, hr3_point_t *s)
  {
    double dd = r4_det(&(p->c), &(q->c), &(r->c), &(s->c));
    return sign_double(dd);
  }
  
hr3_line_t hr3_line_from_two_points(hr3_point_t *p, hr3_point_t *q)
  {
    double a0 = p->c.c[0];
    double a1 = p->c.c[1];
    double a2 = p->c.c[2];
    double a3 = p->c.c[3];
                      
    double b0 = q->c.c[0];
    double b1 = q->c.c[1];
    double b2 = q->c.c[2];
    double b3 = q->c.c[3];
    return 
      (hr3_line_t)
        {{{ a0*b1 - a1*b0,
            a0*b2 - a2*b0,
            a1*b2 - a2*b1,
            a0*b3 - a3*b0,
            a1*b3 - a3*b1,
            a2*b3 - a3*b2
        }}};
  }
  
hr3_plane_t hr3_plane_from_three_points(hr3_point_t *p, hr3_point_t *q, hr3_point_t *r)
  {
    r4_t f;
    r4_cross(&(p->c), &(q->c), &(r->c), &f);
    return (hr3_plane_t){f};
  }   

hr3_plane_t hr3_plane_from_line_and_point(hr3_line_t *L, hr3_point_t *r)
  {
    double a01 = L->k.c[0];
    double a02 = L->k.c[1];
    double a12 = L->k.c[2];
    double a03 = L->k.c[3];
    double a13 = L->k.c[4];
    double a23 = L->k.c[5];
    
    double b0 = r->c.c[0];
    double b1 = r->c.c[1];
    double b2 = r->c.c[2];
    double b3 = r->c.c[3];
    return 
      (hr3_plane_t)
        {{{
            - a12*b3 + a13*b2 - a23*b1,
            + a02*b3 - a03*b2 + a23*b0,
            - a01*b3 + a03*b1 - a13*b0,
            + a01*b2 - a02*b1 + a12*b0
        }}};
  }

hr3_line_t hr3_line_from_two_planes(hr3_plane_t *P, hr3_plane_t *Q)
  {
    double a0 = P->f.c[0];
    double a1 = P->f.c[1];
    double a2 = P->f.c[2];
    double a3 = P->f.c[3];
                      
    double b0 = Q->f.c[0];
    double b1 = Q->f.c[1];
    double b2 = Q->f.c[2];
    double b3 = Q->f.c[3];
    return 
      (hr3_line_t)
        {{{
            + a2*b3 - a3*b2,
            - a1*b3 + a3*b1,
            + a0*b3 - a3*b0,
            + a1*b2 - a2*b1,
            - a0*b2 + a2*b0,
            + a0*b1 - a1*b0
        }}};
  }
  
hr3_point_t hr3_point_from_three_planes(hr3_plane_t *P, hr3_plane_t *Q, hr3_plane_t *R)
  {
    r4_t c;
    r4_cross(&(R->f), &(Q->f), &(P->f), &c);
    return (hr3_point_t){c};
  } 
  
hr3_point_t hr3_point_from_line_and_plane(hr3_line_t *L, hr3_plane_t *R)
  {
    double a01 = L->k.c[0];
    double a02 = L->k.c[1];
    double a12 = L->k.c[2];
    double a03 = L->k.c[3];
    double a13 = L->k.c[4];
    double a23 = L->k.c[5];
     
    double b123 = R->f.c[0];
    double b023 = R->f.c[1];
    double b013 = R->f.c[2];
    double b012 = R->f.c[3];
    return 
      (hr3_point_t)
        {{{
            + a01*b023 + a02*b013 + a03*b012,
            - a01*b123 + a12*b013 + a13*b012,
            - a02*b123 - a12*b023 + a23*b012,
            - a03*b123 - a13*b023 - a23*b013
        }}};
  }


r3_t hr3_line_dir(hr3_line_t *L)
  {
    double a01 = L->k.c[0];
    double a02 = L->k.c[1];
    double a03 = L->k.c[3];
    r3_t d = (r3_t){{a01, a02, a03}};
    r3_dir(&d, &d);
    return d;
  }

r3_t hr3_plane_normal(hr3_plane_t *P)
  {
    double nx = P->f.c[1];
    double ny = P->f.c[2];
    double nz = P->f.c[3];
    double length = hypot(hypot(nx, ny), nz);
    return (r3_t){{nx/length, ny/length, nz/length}};
  }    

r3_t hr3_point_point_dir(hr3_point_t *frm, hr3_point_t *tto)
  {
    double fw = frm->c.c[0];
    double tw = tto->c.c[0];
    
    double fx = frm->c.c[1];
    double tx = tto->c.c[1];
    double dx = fw * tx - tw * fx;
    
    double fy = frm->c.c[2];
    double ty = tto->c.c[2];
    double dy = fw * ty - tw * fy;
    
    double fz = frm->c.c[3];
    double tz = tto->c.c[3];
    double dz = fw * tz - tw * fz;
    
    double length = hypot(hypot(dx, dy), dz);
    
    return (r3_t){{dx/length, dy/length, dz/length}};
  }  
  
double hr3_dist(hr3_point_t *a, hr3_point_t *b)
  {
    double aw = 1.0/a->c.c[0];
    double bw = 1.0/b->c.c[0];
    
    double ax = a->c.c[1];
    double bx = b->c.c[1];
    double dx = ax*aw - bx*bw;
    
    double ay = a->c.c[2];
    double by = b->c.c[2];
    double dy = ay*aw - by*bw;
    
    double az = a->c.c[3];
    double bz = b->c.c[3];
    double dz = az*aw - bz*bw;
    return hypot(hypot(dx, dy), dz);
  }
    
double hr3_dist_sqr(hr3_point_t *a, hr3_point_t *b)
  {
    double aw = 1.0/a->c.c[0];
    double bw = 1.0/b->c.c[0];
    
    double ax = a->c.c[1];
    double bx = b->c.c[1];
    double dx = ax*aw - bx*bw;
    
    double ay = a->c.c[2];
    double by = b->c.c[2];
    double dy = ay*aw - by*bw;
    
    double az = a->c.c[3];
    double bz = b->c.c[3];
    double dz = az*aw - bz*bw;
    return dx*dx + dy*dy + dz*dz;
  }

hr3_point_t hr3_point_mix(double pt, hr3_point_t *p, double qt, hr3_point_t *q)
  { hr3_point_t r;
    r.c.c[0] = pt * p->c.c[0] + qt * q->c.c[0];
    r.c.c[1] = pt * p->c.c[1] + qt * q->c.c[1];
    r.c.c[2] = pt * p->c.c[2] + qt * q->c.c[2];
    r.c.c[3] = pt * p->c.c[3] + qt * q->c.c[3];
    return r;
  }

void hr3_L_inf_normalize_point(hr3_point_t *p)
  { double m = r4_L_inf_norm(&(p->c));
    if (m != 1.0)
      {	p->c.c[0] /= m;
	p->c.c[1] /= m;
	p->c.c[2] /= m;
	p->c.c[3] /= m;
      }
  }
  
void hr3_L_inf_normalize_plane(hr3_plane_t *P)
  { double m = r4_L_inf_norm(&(P->f));
    if (m != 1.0)
      {	P->f.c[0] /= m;
	P->f.c[1] /= m;
	P->f.c[2] /= m;
	P->f.c[3] /= m;
      }
  }

void hr3_L_inf_normalize_line(hr3_line_t *L)
  { double m = r6_L_inf_norm(&(L->k));
    if (m != 1.0)
      {	L->k.c[0] /= m;
	L->k.c[1] /= m;
	L->k.c[2] /= m;
	L->k.c[3] /= m;
	L->k.c[4] /= m;
	L->k.c[5] /= m;
      }
  }
  
hr3_point_t hr3_point_throw(void)
  { hr3_point_t p;
    r4_throw_dir(&(p.c));
    return p;
  }
   
hr3_plane_t hr3_plane_throw(void)
  { hr3_plane_t A;
    r4_throw_dir(&(A.f));
    return A;
  }

void hr3_point_print(FILE *wr, char *pre, hr3_point_t *p, char *fmt, char *suf)
  { if ((pre != NULL) && ((*pre) != 0)) { fputs(pre, wr); }
    fputs("[ ", wr);
    if (fmt == NULL) { fmt = "24.26e"; }
    for (int32_t i = 0; i < NH; i++)
      { if (i != 0) { fputc(' ', wr); }
        fprintf(wr, fmt, p->c.c[i]);
      }
    fputs(" ]", wr);
    if ((suf != NULL) && ((*suf) != 0)) { fputs(suf, wr); }
  }

void hr3_plane_print(FILE *wr, char *pre, hr3_plane_t *P, char *fmt, char *suf)
  { if ((pre != NULL) && ((*pre) != 0)) { fputs(pre, wr); }
    fputs("< ", wr);
    if (fmt == NULL) { fmt = "24.26e"; }
    for (int32_t i = 0; i < NH; i++)
      { if (i != 0) { fputc(' ', wr); }
        fprintf(wr, fmt, P->f.c[i]);
      }
    fputs(" >", wr);
    if ((suf != NULL) && ((*suf) != 0)) { fputs(suf, wr); }
  }

void hr3_line_print(FILE *wr, char *pre, hr3_line_t *L, char *fmt, char *suf)
  { if ((pre != NULL) && ((*pre) != 0)) { fputs(pre, wr); }
    fputs("[[ ", wr);
    if (fmt == NULL) { fmt = "24.26e"; }
    for (int32_t i = 0; i < NG; i++)
      { if (i != 0) { fputc(' ', wr); }
        fprintf(wr, fmt, L->k.c[i]);
      }
    fputs(" ]]", wr);
    if ((suf != NULL) && ((*suf) != 0)) { fputs(suf, wr); }
  }
