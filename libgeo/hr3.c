/* See hr3.h */
/* Last edited on 2021-06-09 20:15:55 by jstolfi */ 

#define _GNU_SOURCE
#include <stdint.h>
#include <hr3.h>

#include <r3.h>
#include <r4.h>
#include <r6.h>
#include <r3x3.h>
#include <r4x4.h>
#include <affirm.h>
#include <sign.h>
#include <sign_get.h>

#include <math.h>

/* Based on HR3.m3, created 1993-05-18 by Marcos C. Carrard. */
/* Based on H3.pas by J. Stolfi. */

hr3_point_t hr3_from_r3(r3_t *c)
  {
    return (hr3_point_t){{{1.0, c->c[0], c->c[1], c->c[2]}}};
  }
    
r3_t r3_from_hr3(hr3_point_t *p)
  {
    double w = p->c.c[0];
    demand(w != 0.0, "point at infinity");
    return (r3_t){{p->c.c[1]/w, p->c.c[2]/w, p->c.c[3]/w}};
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

hr3_plane_t hr3_plane_from_line_and_point(hr3_line_t *n, hr3_point_t *r)
  {
    double a01 = n->k.c[0];
    double a02 = n->k.c[1];
    double a12 = n->k.c[2];
    double a03 = n->k.c[3];
    double a13 = n->k.c[4];
    double a23 = n->k.c[5];
    
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
  
hr3_point_t hr3_point_from_line_and_plane(hr3_line_t *n, hr3_plane_t *R)
  {
    double a01 = n->k.c[0];
    double a02 = n->k.c[1];
    double a12 = n->k.c[2];
    double a03 = n->k.c[3];
    double a13 = n->k.c[4];
    double a23 = n->k.c[5];
     
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

bool_t hr3_pmap_is_identity(hr3_pmap_t *m)
  {
    return r4x4_is_unif_scaling(&(m->dir), m->dir.c[0][0]); 
  }

hr3_point_t hr3_pmap_point(hr3_point_t *p, hr3_pmap_t *m)
  {
    r4_t c;
    r4x4_map_row(&(p->c), &(m->dir), &c);
    return (hr3_point_t){c};
  }
        
hr3_point_t hr3_pmap_inv_point(hr3_point_t *p, hr3_pmap_t *m)
  {
    r4_t c;
    r4x4_map_row(&(p->c), &(m->inv), &c);
    return (hr3_point_t){c};
  }
        
hr3_plane_t hr3_pmap_plane(hr3_plane_t *P, hr3_pmap_t *m)
  {
    r4_t f;
    r4x4_map_col(&(m->inv), &(P->f), &f);
    return (hr3_plane_t){f};
  }
    
hr3_plane_t hr3_pmap_inv_plane(hr3_plane_t *P, hr3_pmap_t *m)
  {
    r4_t f;
    r4x4_map_col(&(m->dir), &(P->f), &f);
    return (hr3_plane_t){f};
  }
    

hr3_pmap_t hr3_pmap_translation(r3_t *v)
  {
    hr3_pmap_t m; 
    r4x4_ident(&(m.dir));
    m.dir.c[0][1] = +v->c[0];
    m.dir.c[0][2] = +v->c[1];
    m.dir.c[0][3] = +v->c[2];
    r4x4_ident(&(m.inv));
    m.inv.c[0][1] = -v->c[0];
    m.inv.c[0][2] = -v->c[1];
    m.inv.c[0][3] = -v->c[2];
    return m;
  }

hr3_pmap_t hr3_pmap_u_v_rotation(r3_t *u, r3_t *v)
{ /* Compute the matrix {r} for {R^3} to {R^3}: */
    r3x3_t r;
    r3x3_u_v_rotation(u, v, &r);
    
    /* Convert to a projective map (note that inverse is just transpose): */
    hr3_pmap_t m;
    int32_t i, j;
    for (i = 0; i < 4; i++)
      { for (j = 0; j < 4; j++)
          { if ((i == 0) && (j == 0))
              { m.dir.c[i][j] = m.inv.c[i][j] = 1.0; }
            else if ((i == 0) || (j == 0))
              { m.dir.c[i][j] = m.inv.c[i][j] = 0.0; }
            else
              { m.dir.c[i][j] = r.c[i-1][j-1]; 
                m.inv.c[i][j] = r.c[j-1][i-1]; 
              }
          }
      }
    return m;
  }

hr3_pmap_t hr3_pmap_comp(hr3_pmap_t *m, hr3_pmap_t *n)
  {
    r4x4_t dir, inv; 
    r4x4_mul(&(m->dir), &(n->dir), &dir);
    r4x4_mul(&(n->inv), &(m->inv), &inv);
    return (hr3_pmap_t){dir, inv};
  }    
  
hr3_pmap_t hr3_pmap_inv(hr3_pmap_t *m)
  {
    return (hr3_pmap_t){m->inv, m->dir};
  }

r3_t hr3_line_dir(hr3_line_t *n)
  {
    double a01 = n->k.c[0];
    double a02 = n->k.c[1];
    double a03 = n->k.c[3];
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
    return dx * dx + dy * dy + dz * dz;
  }

hr3_point_t hr3_point_at_infinity(r3_t *dir)
  { hr3_point_t p;
    p.c.c[0] = 0.0;
    p.c.c[1] = dir->c[0];
    p.c.c[2] = dir->c[1];
    p.c.c[3] = dir->c[2];
    return p;
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

void hr3_L_inf_normalize_line(hr3_line_t *n)
  { double m = r6_L_inf_norm(&(n->k));
    if (m != 1.0)
      {	n->k.c[0] /= m;
	n->k.c[1] /= m;
	n->k.c[2] /= m;
	n->k.c[3] /= m;
	n->k.c[4] /= m;
	n->k.c[5] /= m;
      }
  }

hr3_pmap_t hr3_pmap_persp(hr3_point_t *obs, hr3_point_t *foc, double rad, hr3_point_t *upp)
  { hr3_pmap_t m;
    r4x4_t mt, mr, mc, mrt;
    int32_t i, j;
    
    demand(foc->c.c[0] > 0.0, "focus must be finite and hither");
    demand(obs->c.c[0] >= 0.0, "observer must be hither or infinite");
    demand(upp->c.c[0] >= 0.0, "zenith must be hither or infinite"); 

    /* Start with a translation from {foc} to the origin: */
    for (i = 0; i < 4; i++){
      for (j = 0; j < 4; j++){
        if (i == j)
          { mt.c[i][j] = foc->c.c[0]; }
        else if (i == 0)
          { mt.c[i][j] = -(foc->c.c[j]); }
        else
          { mt.c[i][j] = 0.0; }
      }
    }
    
    /* Compute the image frame vectors {r,s,t}: */
    r3_t t = hr3_point_point_dir(foc, obs); /* Vector {t} points out of image towards {obs}. */
    r3_t u = hr3_point_point_dir(foc, upp); /* Direction from {foc} towards zenith. */
    r3_t v, s;
    r3_decomp(&u, &t, &v, &s);
    /* Zenith reference point {upp} must not be on imagesys Z axis: */
    demand(r3_norm_sqr(&s) >= 1.0e-12, "bad zenith"); 
    r3_dir(&s, &s); /* Vector {s} is the image's vertical axis. */
    r3_t r;
    r3_cross(&s, &t, &r);
    r3_dir(&r, &r); /* Vector {r} is the image's horizontal axis. */

    /* Append the rotation matrix that moves {r,s,t} to X,Y,Z: */
    mr.c[0][0] = 1.0;
    for (i = 1; i < 4; i++)
      { mr.c[0][i] = 0.0;
        mr.c[i][0] = 0.0;
        mr.c[i][1] = r.c[i-1];
        mr.c[i][2] = s.c[i-1];
        mr.c[i][3] = t.c[i-1];
      }
    r4x4_mul(&mt, &mr, &mrt);

    /* Do we need a conical projection step? */
    if (obs->c.c[0] == 0.0)
      { /* Observer is at infinity; cilindrical projection. */
        m.dir = mrt;
      }
    else
      { /* Observer is finite; add conical projection step. */
        double d = hr3_dist(foc, obs);
        double uno = -1.0;
        if (fabs(d) > 1.0) { uno = -1.0/d; d = 1.0; }
        for (i = 0; i < 4; i++)
          for (j = 0; j < 4; j++)
            { if (i == j)
                { mc.c[i][j] = d; }
              else
                { mc.c[i][j] = 0.0; }
            }
        mc.c[3][0] = uno;
        r4x4_mul(&mrt, &mc, &(m.dir));
      }

    /* Asked for scaling? */
    if (rad > 0.0)
      { /* Combine {map} with a uniform scale of {1/rad}: */
        for (i = 0; i < 4; i++)
          for (j = 1; j < 4; j++)
            { m.dir.c[i][j] /= rad; }
      }

    /* Compute inverse matrix: */
    r4x4_inv(&m.dir, &m.inv);
    return m;
  }

void hr3_point_print (FILE *f, char *pre, hr3_point_t *p, char *fmt, char *suf)
  { if ((pre != NULL) && ((*pre) != 0)) { fputs(pre, f); }
    fputc('[', f);
    if (fmt == NULL) { fmt = "24.26e"; }
    int32_t i;
    for (i = 0; i < 4; i++)
      { if (i != 0) { fputc(' ', f); }
        fprintf(f, fmt, p->c.c[i]);
      }
    fputc(']', f);
    if ((suf != NULL) && ((*suf) != 0)) { fputs(suf, f); }
  }

void hr3_plane_print (FILE *f, char *pre, hr3_plane_t *P, char *fmt, char *suf)
  { if ((pre != NULL) && ((*pre) != 0)) { fputs(pre, f); }
    fputc('[', f);
    if (fmt == NULL) { fmt = "24.26e"; }
    int32_t i;
    for (i = 0; i < 4; i++)
      { if (i != 0) { fputc(' ', f); }
        fprintf(f, fmt, P->f.c[i]);
      }
    fputc(']', f);
    if ((suf != NULL) && ((*suf) != 0)) { fputs(suf, f); }
  }

void hr3_line_print (FILE *f, char *pre, hr3_line_t *n, char *fmt, char *suf)
  { if ((pre != NULL) && ((*pre) != 0)) { fputs(pre, f); }
    fputc('[', f);
    if (fmt == NULL) { fmt = "24.26e"; }
    int32_t i;
    for (i = 0; i < 4; i++)
      { if (i != 0) { fputc(' ', f); }
        fprintf(f, fmt, n->k.c[i]);
      }
    fputc(']', f);
    if ((suf != NULL) && ((*suf) != 0)) { fputs(suf, f); }
  }
