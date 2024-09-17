/* See hr2_pmap.h */
/* Last edited on 2024-09-17 16:24:20 by stolfi */ 

#define _GNU_SOURCE
#include <stdint.h>
#include <math.h>
#include <assert.h>
#include <stdio.h>
#include <string.h>

#include <sign_get.h>
#include <sign.h>
#include <affirm.h>

#include <r2x2.h>
#include <r3x3.h>
#include <rn.h>
#include <r3.h>
#include <r2.h>
#include <hr2.h>


#include <hr2_pmap.h>

#define NH 3
  /* Number of homogeneous coordinates in a point. */

#define NC 2
  /* Number of Cartesian coordinates in a point. */


hr2_point_t hr2_pmap_point(hr2_point_t *p, hr2_pmap_t *M)
  {
    r3_t q; 
    r3x3_map_row(&(p->c), &(M->dir), &q);
    return (hr2_point_t){q};
  }

hr2_point_t hr2_pmap_inv_point(hr2_point_t *p, hr2_pmap_t *M)
  {
    r3_t q; 
    r3x3_map_row(&(p->c), &(M->inv), &q);
    return (hr2_point_t){q};
  }

hr2_line_t hr2_pmap_line(hr2_line_t *L, hr2_pmap_t *M)
  {
    r3_t f;
    r3x3_map_col(&(M->inv), &(L->f), &f);
    return (hr2_line_t){f};
  }

hr2_line_t hr2_pmap_inv_line(hr2_line_t *L, hr2_pmap_t *M)
  {
    r3_t f;
    r3x3_map_col(&(M->dir), &(L->f), &f);
    return (hr2_line_t){f};
  }

r2_t hr2_pmap_r2_point(r2_t *p, hr2_pmap_t *M)
  {
    hr2_point_t ph = (hr2_point_t){{{ 1.0, p->c[0], p->c[1] }}};
    hr2_point_t qh;
    r3x3_map_row(&(ph.c), &(M->dir), &(qh.c));
    r2_t qc; hr2_point_to_r2_nan(&qh, &qc);
    return qc;
  }

r2_t hr2_pmap_inv_r2_point(r2_t *p, hr2_pmap_t *M)
  {
    hr2_point_t ph = (hr2_point_t){{{ 1, p->c[0], p->c[1] }}};
    hr2_point_t qh;
    r3x3_map_row(&(ph.c), &(M->inv), &(qh.c));
    r2_t qc; hr2_point_to_r2_nan(&qh, &qc);
    return qc;
  }

hr2_pmap_t hr2_pmap_identity(void)
  { hr2_pmap_t M;
    r3x3_ident(&(M.dir));
    M.inv = M.dir;
    return M;
  }

hr2_pmap_t hr2_pmap_mirror(sign_t sgnx, sign_t sgny)
  { demand(sgnx != 0, "{sgnx} must be {-1} or {+1}");
    demand(sgny != 0, "{sgny} must be {-1} or {+1}");
    hr2_pmap_t M;
    r3x3_ident(&(M.dir));
    M.dir.c[1][1] = sgnx;
    M.dir.c[2][2] = sgny;
    M.inv = M.dir;
    return M;
  }

hr2_pmap_t hr2_pmap_xy_swap(void)
  { hr2_pmap_t M;
    r3x3_ident(&(M.dir));
    M.dir.c[1][1] = 0; M.dir.c[1][2] = 1;
    M.dir.c[2][1] = 1; M.dir.c[2][2] = 0;
    M.inv = M.dir;
    return M;
  }

hr2_pmap_t hr2_pmap_translation(r2_t *vec)
  { hr2_pmap_t M;
    r3x3_ident(&(M.dir));

    M.dir.c[0][1] = +vec->c[0];
    M.dir.c[0][2] = +vec->c[1];

    r3x3_ident(&(M.inv));
    M.inv.c[0][1] = -vec->c[0];
    M.inv.c[0][2] = -vec->c[1];
    return M;
  }

hr2_pmap_t hr2_pmap_rotation(double ang)
  {
    double c = cos(ang);
    double s = sin(ang);

    hr2_pmap_t M;
    r3x3_ident(&(M.dir));
    r3x3_ident(&(M.inv));

    M.dir.c[1][1] = +c;
    M.dir.c[1][2] = +s;
    M.dir.c[2][1] = -s;
    M.dir.c[2][2] = +c;

    M.inv.c[1][1] = +c;
    M.inv.c[1][2] = -s;
    M.inv.c[2][1] = +s;
    M.inv.c[2][2] = +c;
    return M;
  }

hr2_pmap_t hr2_pmap_scaling(r2_t *scale)
  {
    for (int32_t j = 0; j < 2; j++)
      { demand(scale->c[j] != 0.0, "cannot scale by zero"); }
    
    hr2_pmap_t M;
    r3x3_ident(&(M.dir));
    r3x3_ident(&(M.inv));

    M.dir.c[1][1] = scale->c[0];
    M.dir.c[2][2] = scale->c[1];

    M.inv.c[1][1] = 1/scale->c[0];
    M.inv.c[2][2] = 1/scale->c[1];
    return M;
  }

hr2_pmap_t hr2_pmap_rotation_and_scaling(double ang, double scale)
  {
    double c = cos(ang);
    double s = sin(ang);

    hr2_pmap_t M;
    M.dir.c[0][0] = 1.0;       M.inv.c[0][0] = 1.0;     
    M.dir.c[0][1] = 0.0;       M.inv.c[0][1] = 0.0;     
    M.dir.c[0][2] = 0.0;       M.inv.c[0][2] = 0.0;     

    M.dir.c[1][0] = 0.0;       M.inv.c[1][0] = 0.0;     
    M.dir.c[1][1] = +c*scale;  M.inv.c[1][1] = +c/scale;
    M.dir.c[1][2] = +s*scale;  M.inv.c[1][2] = -s/scale;

    M.dir.c[2][0] = 0.0;       M.inv.c[2][0] = 0.0;     
    M.dir.c[2][1] = -s*scale;  M.inv.c[2][1] = +s/scale;
    M.dir.c[2][2] = +c*scale;  M.inv.c[2][2] = +c/scale;

    return M;
  }

hr2_pmap_t hr2_pmap_compose(hr2_pmap_t *M, hr2_pmap_t *N)
  {
    r3x3_t dir, inv;
    r3x3_mul(&(M->dir), &(N->dir), &dir);
    r3x3_mul(&(N->inv), &(M->inv), &inv);
    return (hr2_pmap_t){dir, inv};
  }

hr2_pmap_t hr2_pmap_inv(hr2_pmap_t *M)
  {
    return (hr2_pmap_t){M->inv, M->dir};
  }

hr2_pmap_t hr2_pmap_inv_compose(hr2_pmap_t *M, hr2_pmap_t *N)
  {
    r3x3_t dir, inv;
    r3x3_mul(&(M->inv), &(N->dir), &dir);
    r3x3_mul(&(N->inv), &(M->dir), &inv);
    return (hr2_pmap_t){dir, inv};
  }

hr2_pmap_t hr2_pmap_congruence_from_point_and_dir(r2_t *p, r2_t *u, sign_t sgn)
  { 
    /* Normalize {u} to unit length to get the image of the {(1,0)} vector: */
    r2_t du;  
    double mu = r2_dir(u, &du);
    demand(mu != 0, "invalid direction {u}");

    /* Get the unit vector {dv} that is to be the image of the {(0,1)} vector: */
    r2_t dv;
    if (sgn < 0)
      { dv = (r2_t){{ +du.c[1], -du.c[0] }}; }
    else
      { dv = (r2_t){{ -du.c[1], +du.c[0] }}; }
      
    /* Assemble the matrix: */
    hr2_pmap_t M; /* The resulting map. */
    M.dir.c[0][0] = 1.0;     
    M.dir.c[0][1] = p->c[0]; 
    M.dir.c[0][2] = p->c[1]; 

    M.dir.c[1][0] = 0.0;     
    M.dir.c[1][1] = du.c[0]; 
    M.dir.c[1][2] = du.c[1]; 

    M.dir.c[2][0] = 0.0;     
    M.dir.c[2][1] = dv.c[0]; 
    M.dir.c[2][2] = dv.c[1]; 
    
    r3x3_inv(&(M.dir), &(M.inv));

    return M;
  }

hr2_pmap_t hr2_pmap_similarity_from_two_points(r2_t *p, r2_t *q, sign_t sgn)
  {
    /* Compute the vector {u} that is the image of the {(1,0)} vector: */
    r2_t u; r2_sub(q, p, &u);

    /* Get the squared length {d2} of {u}: */
    double d2 = r2_norm_sqr(&u);
    demand(d2 != 0, "points {p,q} coincide");

    /* Get the vector {v} that is to be the image of the {(0,1)} vector: */
    r2_t v;
    if (sgn < 0)
      { v = (r2_t){{ +u.c[1], -u.c[0] }}; }
    else
      { v = (r2_t){{ -u.c[1], +u.c[0] }}; }
      
    /* Assemble the matrix: */
    hr2_pmap_t M; /* The resulting map. */
    M.dir.c[0][0] = 1.0;        
    M.dir.c[0][1] = p->c[0];    
    M.dir.c[0][2] = p->c[1];    

    M.dir.c[1][0] = 0.0;         
    M.dir.c[1][1] = +u.c[0];    
    M.dir.c[1][2] = +u.c[1];    

    M.dir.c[2][0] = 0.0;         
    M.dir.c[2][1] = +v.c[0];    
    M.dir.c[2][2] = +v.c[1];    
    
    r3x3_inv(&(M.dir), &(M.inv));

    return M;
  }

hr2_pmap_t hr2_pmap_from_four_points(hr2_point_t *p, hr2_point_t *q, hr2_point_t *r, hr2_point_t *u)
  {
    hr2_pmap_t M; /* The resulting map. */
    
    /* Compute weights {(a,b,c)=w.c[0..2]} for rows of {Q} that map {[1,1,1]} to {u}: */
    r3_t w;
    { /* Compute a matrix {Q} that maps the cardinal points to {p,q,r} as given: */
      r3x3_t Q;
      Q.c[0][0] = p->c.c[0]; Q.c[0][1] = p->c.c[1]; Q.c[0][2] = p->c.c[2];
      Q.c[1][0] = q->c.c[0]; Q.c[1][1] = q->c.c[1]; Q.c[1][2] = q->c.c[2];
      Q.c[2][0] = r->c.c[0]; Q.c[2][1] = r->c.c[1]; Q.c[2][2] = r->c.c[2];
      /* Map {u} by the inverse of {Q}: */
      r3x3_inv(&Q, &Q);
      r3x3_map_row(&(u->c), &Q, &w);
    }
    
    /* Make the weights positive, so that {p,q,r} are strictly honored: */
    for (int32_t i = 0; i < NH; i++) { w.c[i] = fabs(w.c[i]); }

    /* Ensure that {M.dir} maps the cardinal points to {p,q,r} and some unit point to {u}: */
    r3x3_t *P = &(M.dir);
    P->c[0][0] = w.c[0]*p->c.c[0]; P->c[0][1] = w.c[0]*p->c.c[1]; P->c[0][2] = w.c[0]*p->c.c[2];
    P->c[1][0] = w.c[1]*q->c.c[0]; P->c[1][1] = w.c[1]*q->c.c[1]; P->c[1][2] = w.c[1]*q->c.c[2];
    P->c[2][0] = w.c[2]*r->c.c[0]; P->c[2][1] = w.c[2]*r->c.c[1]; P->c[2][2] = w.c[2]*r->c.c[2];

    /* Compute the inverse map: */
    r3x3_inv(&(M.dir), &(M.inv));

    return M;
  }

hr2_pmap_t hr2_pmap_from_four_r2_points(r2_t *p, r2_t *q, r2_t *r, r2_t *u)
  {
    hr2_point_t hp = hr2_from_r2(p);
    hr2_point_t hq = hr2_from_r2(q);
    hr2_point_t hr = hr2_from_r2(r);
    hr2_point_t hu = hr2_from_r2(u);
    hr2_pmap_t M = hr2_pmap_from_four_points(&hp, &hq, &hr, &hu);
    for (int32_t j = 0; j < NH; j++)
      { M.dir.c[1][j] -= M.dir.c[0][j];  
        M.dir.c[2][j] -= M.dir.c[0][j];
        M.inv.c[j][0] += M.inv.c[j][1] + M.inv.c[j][2];
      }
    return M;
  }

hr2_pmap_t hr2_pmap_aff_from_mat_and_disp(r2x2_t *E, r2_t *d)
  {
    hr2_pmap_t M;

    M.dir.c[0][0] = 1.0;
    M.dir.c[0][1] = d->c[0];
    M.dir.c[0][2] = d->c[1];
    
    M.dir.c[1][0] = 0.0;
    M.dir.c[1][1] = E->c[0][0];
    M.dir.c[1][2] = E->c[0][1];
    
    M.dir.c[2][0] = 0.0;
    M.dir.c[2][1] = E->c[1][0];
    M.dir.c[2][2] = E->c[1][1];
    
    r2x2_t F;
    r2x2_inv(E, &(F));
    r2_t psid;
    r2x2_map_row(d, &(F), &psid);
    
    M.inv.c[0][0] = 1.0;
    M.inv.c[0][1] = - psid.c[0];
    M.inv.c[0][2] = - psid.c[1];
    
    M.inv.c[1][0] = 0.0;
    M.inv.c[1][1] = F.c[0][0];
    M.inv.c[1][2] = F.c[0][1];
    
    M.inv.c[2][0] = 0.0;
    M.inv.c[2][1] = F.c[1][0];
    M.inv.c[2][2] = F.c[1][1];
    
    return M;
  }

hr2_pmap_t hr2_pmap_aff_from_three_points(r2_t *o, r2_t *p, r2_t *q)
  {
    hr2_pmap_t M; 
    
    M.dir.c[0][0] = 1.0;
    M.dir.c[0][1] = o->c[0];
    M.dir.c[0][2] = o->c[1];

    M.dir.c[1][0] = 0.0;
    M.dir.c[1][1] = p->c[0] - o->c[0];
    M.dir.c[1][2] = p->c[1] - o->c[1];

    M.dir.c[2][0] = 0.0;
    M.dir.c[2][1] = q->c[0] - o->c[0];
    M.dir.c[2][2] = q->c[1] - o->c[1];

    r3x3_inv(&(M.dir), &(M.inv));

    /* Just in case: */
    M.inv.c[1][0] = 0.0;
    M.inv.c[2][0] = 0.0;
    
    return M;
  }

hr2_pmap_t hr2_pmap_r2_from_class(int32_t class)
  {
    hr2_pmap_t M; 
    
    M.dir.c[0][0] = 1;
    M.dir.c[0][1] = 0;
    M.dir.c[0][2] = 0;

    int32_t b = ((class & 2) == 0 ? +1 : -1);
    M.dir.c[1][0] = b - 1;
    M.dir.c[1][1] = b;
    M.dir.c[1][2] = 0;

    int32_t c = ((class & 1) == 0 ? +1 : -1);
    M.dir.c[2][0] = c - 1;
    M.dir.c[2][1] = 0;
    M.dir.c[2][2] = c;

    r3x3_inv(&(M.dir), &(M.inv));
    
    return M;
  }

double hr2_pmap_diff_sqr(hr2_pmap_t *M, hr2_pmap_t *N)
  { double sum_d2 = 0;
    for (int32_t sense = 0; sense < 2; sense++)
      { r3x3_t *A = (sense == 0 ? &(M->dir) : &(M->inv));
        double Am = r3x3_norm(A) + 1.0e-200;
        r3x3_t *B = (sense == 0 ? &(N->dir) : &(N->inv));
        double Bm = r3x3_norm(B) + 1.0e-200;
        for (int32_t i = 0; i < NH; i++)
          { for (int32_t j = 0; j < NH; j++)
             { double Aij = A->c[i][j]/Am;
               double Bij = B->c[i][j]/Bm;
               double dij = Aij - Bij;
               sum_d2 += dij*dij;
             }
          }
      }
    return sum_d2;
  }

double hr2_pmap_mismatch_sqr(hr2_pmap_t *M, int32_t np, r2_t p1[], r2_t p2[], double w[])
  {
    double sum_wD2 = 0.0;
    double sum_w = 1.0e-200; /* In case there are no points with positive weight. */
    for (int32_t k = 0; k < np; k++)
      { r2_t *p1k = &(p1[k]);
        r2_t *p2k = &(p2[k]);
        double wk = (w == NULL ? 1.0 : w[k]);
        r2_t q1k = hr2_pmap_r2_point(p1k, M);
        r2_t q2k = hr2_pmap_inv_r2_point(p2k, M);
        double D2k_dir = r2_dist_sqr(&q1k, p2k);
        double D2k_inv = r2_dist_sqr(p1k, &q2k);
        sum_wD2 += wk*(D2k_dir + D2k_inv);
        sum_w += wk;
      }
    return 0.5*sum_wD2/sum_w;
  }

double hr2_pmap_deform_sqr(r2_t ph[], hr2_pmap_t *M)
  {
    int32_t nk = (1 << NC); /* Number of corners of the quadrilateral. */
    r2_t qh[nk];
    for (int32_t k = 0; k < nk; k++)
      { qh[k] = hr2_pmap_r2_point(&(ph[k]), M); }
    
    int32_t nd = 4 + 2; /* Number of distances to probe. */
    double logr[nd];
    int32_t kd = 0;
    for (int32_t ik = 1; ik < nk; ik++)
      { for (int32_t jk = 0; jk < ik; jk++)
          { double dp2 = r2_dist_sqr(&(ph[ik]), &(ph[jk]));
            double dq2 = r2_dist_sqr(&(qh[ik]), &(qh[jk]));
            assert(kd < nd);
            logr[kd] = log(dq2/dp2)/2;
            kd++;
          }
      }
    assert(kd == nd);
    
    /* Compute the variance of the logs: */
    double sum = 0;
    for (int32_t kd = 0; kd < nd; kd++) { sum += logr[kd]; }
    double avg = sum/nd;
    double sum2 = 0;
    for (int32_t kd = 0; kd < nd; kd++) { double dk = logr[kd] - avg; sum2 += dk*dk; }
    double var = sum2/(nd-1);
    return var;
  }

double hr2_pmap_aff_discr_sqr(hr2_pmap_t *M, hr2_pmap_t *N)
  {
    demand((M->dir.c[1][0] == 0) && (M->dir.c[2][0] == 0), "{M} is not affine");
    demand(M->dir.c[0][0] > 0, "map {M} does not preserve side");
    r3x3_t A; double wA = M->dir.c[0][0]; r3x3_scale(1/wA, &(M->dir), &A);
   
    demand((N->dir.c[1][0] == 0) && (N->dir.c[2][0] == 0), "{N} is not affine");
    demand(M->dir.c[0][0] > 0, "map {N} does not preserve side");
    r3x3_t B; double wB = N->dir.c[0][0]; r3x3_scale(1/wB, &(N->dir), &B);
   
    /* Hope the math is right: */
    r3x3_t E, H;
    r3x3_sub(&A, &B, &E);
    r3x3_mul_tr(&E, &E, &H);
    double h2 = (H.c[1][1] + H.c[2][2])/2;
    double d2 = H.c[0][0];
    return h2 + d2;
  }

bool_t hr2_pmap_is_valid(hr2_pmap_t *M, double tol)
  { bool_t debug = FALSE;
    
    if (debug) { fprintf(stderr, "   . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .\n"); }
    
    r3x3_t P; r3x3_mul(&(M->dir), &(M->inv), &P);

    if (debug) 
      { fprintf(stderr, "  M = \n");
        hr2_pmap_gen_print(stderr, M, "%+24.16e", NULL,  NULL,NULL,NULL, NULL,NULL,NULL, "\n");
        fprintf(stderr, "  P = \n");
        r3x3_gen_print(stderr, &P, "%+24.16e", "  ","\n  ","\n", NULL,NULL,NULL);
      } 
    double mv = r3x3_L_inf_normalize(&P);
    if (mv == 0) { return FALSE; }
    for (int32_t i = 0; i < NH; i++)
      for (int32_t j = 0; j < NH; j++)
        { double vexp = (i == j ? 1.0 : 0.0);
          double err = P.c[i][j] - vexp;
          if (fabs(err) > tol) { return FALSE; }
        }
    return TRUE;
  }

sign_t hr2_pmap_sign(hr2_pmap_t *M)
  { double det = r3x3_det(&(M->dir));
    if (det < 0)
      { return -1; }
    else if (det > 0)
      { return +1; }
    else
      { return 0; }
  }

bool_t hr2_pmap_is_identity(hr2_pmap_t *M, double tol)
  { 
    for (int32_t dir = 0; dir <= 1; dir++)
      { r3x3_t *A = (dir == 0 ? &(M->dir) : &(M->inv));
        double A00 = A->c[0][0];
        if (A00 < 1.0e-200) { return FALSE; }
        double Atol = A00*tol;
        
        if (fabs(A->c[0][1]) > Atol) { return FALSE; }
        if (fabs(A->c[0][2]) > Atol) { return FALSE; }
        
        if (fabs(A->c[1][0]) > Atol) { return FALSE; }
        if (fabs(A->c[1][1] - A00) > Atol) { return FALSE; }
        if (fabs(A->c[1][2]) > Atol) { return FALSE; }
        
        if (fabs(A->c[2][0]) > Atol) { return FALSE; }
        if (fabs(A->c[2][1]) > Atol) { return FALSE; }
        if (fabs(A->c[2][2] - A00) > Atol) { return FALSE; }
      }
    return TRUE;
  }

bool_t hr2_pmap_is_translation(hr2_pmap_t *M, double tol)
  {
    if (! hr2_pmap_is_affine(M, tol)) { return FALSE; }
    for (int32_t dir = 0; dir <= 1; dir++)
      { r3x3_t *A = (dir == 0 ? &(M->dir) : &(M->inv));
        double A00 = A->c[0][0];
        if (A00 < 1.0e-200) { return FALSE; }
        double Atol = A00*tol;
        
        if (fabs(A->c[1][1] - A00) > Atol) { return FALSE; }
        if (fabs(A->c[1][2]) > Atol) { return FALSE; }
        
        if (fabs(A->c[2][1]) > Atol) { return FALSE; }
        if (fabs(A->c[2][2] - A00) > Atol) { return FALSE; }
      }
    return TRUE;
  }
  
bool_t hr2_pmap_is_congruence(hr2_pmap_t *M, double tol)
  {
    return hr2_pmap_is_similarity(M, 1.0, tol);
  }
    
bool_t hr2_pmap_is_similarity(hr2_pmap_t *M, double scale, double tol)
  {
    bool_t debug = FALSE;
    
    if (debug) { fprintf(stderr, "   . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .\n"); }
    
    if (debug) 
      { fprintf(stderr, "  M = \n");
        hr2_pmap_gen_print(stderr, M, "%+24.16e", NULL,  NULL,NULL,NULL, NULL,NULL,NULL, "\n");
      }
    
    if (! hr2_pmap_is_affine(M, tol)) { return FALSE; }
    for (int32_t dir = 0; dir <= 1; dir++)
      { r3x3_t *A = (dir == 0 ? &(M->dir) : &(M->inv));
        double A00 = A->c[0][0];
        double A2tol = A00*A00*tol;
        
        double uu = A->c[1][1]*A->c[1][1] + A->c[1][2]*A->c[1][2];
        double uv = A->c[1][1]*A->c[2][1] + A->c[1][2]*A->c[2][2];
        double vv = A->c[2][1]*A->c[2][1] + A->c[2][2]*A->c[2][2];
        
        if (debug) 
          { fprintf(stderr, "    uu/AA = %24.16e", uu/(A00*A00)); 
            fprintf(stderr, " vv/AA = %24.16e",    vv/(A00*A00)); 
            fprintf(stderr, " uv/AA = %24.16e\n",  uv/(A00*A00));
          } 
        
        if (fabs(uv) > A2tol) { return FALSE; }
        if (fabs((uu - vv)) > A2tol) { return FALSE; }
        
        if (! isnan(scale))
          { demand(scale > 1.0e-100, "invalid {scale}");
            double sc = (dir == 0 ? scale : 1.0/scale);
            double det = A->c[1][1]*A->c[2][2] - A->c[1][2]*A->c[2][1];
            /* The determinant should be the same as {uu} and {vv}, so: */
            double err = fabs(fabs(det) - sc*sc*A00*A00);
            if (debug) 
              { fprintf(stderr, "    det/AA = %24.16e sc2 = %24.16e",  fabs(det)/(A00*A00), sc*sc); 
                fprintf(stderr, "    err/AA = %24.16e\n",  err/(A00*A00)); 
              } 
            assert(fabs((2*fabs(det) - uu - vv)) < 10*A2tol);
            if (err > 4*A2tol) { return FALSE; }
          }
      }
    return TRUE;
  }

bool_t hr2_pmap_is_affine(hr2_pmap_t *M, double tol)
  {
    for (int32_t dir = 0; dir <= 1; dir++)
      { r3x3_t *A = (dir == 0 ? &(M->dir) : &(M->inv));
        double A00 = A->c[0][0];
        if (A00 < 1.0e-200) { return FALSE; }
        double Atol = A00*tol;
        
        if (fabs(A->c[1][0]) > Atol) { return FALSE; }
        if (fabs(A->c[2][0]) > Atol) { return FALSE; }
        
        double det = A->c[1][1]*A->c[2][2] - A->c[1][2]*A->c[2][1];
        if (fabs(det) < A00*A00*tol) { return FALSE; }
      }
    return TRUE;
  }
    
bool_t hr2_pmap_is_generic(hr2_pmap_t *M, double tol)
  {
    for (int32_t dir = 0; dir <= 1; dir++)
      { r3x3_t *A = (dir == 0 ? &(M->dir) : &(M->inv));
        double Amag = r3x3_norm(A);
        if (Amag < 1.0e-200) { return FALSE; }
        
        double det = r3x3_det(A);
        if (fabs(det) < Amag*Amag*Amag*tol) { return FALSE; }
      }
    return TRUE;
  }

void hr2_pmap_print (FILE *wr, hr2_pmap_t *M, char *pref, char *suff)
  { 
    hr2_pmap_gen_print(wr, M, "%12.6f", pref, "  ", "  ", "\n", "[ ", " ", " ]", suff);
    fflush(wr);
  }

void hr2_pmap_gen_print
  ( FILE *wr, hr2_pmap_t *M,
    char *fmt, char *pref,                /* Overall prefix. */
    char *rpref, char *rsep, char *rsuff, /* Row prefix, matrix separator, and suffix. */
    char *elp, char *esep, char *erp,     /* Delimiters for each row. */
    char *suff                            /* Overall sufffix. */
  )
  {
    /* Defaults: */
    if (fmt == NULL) { fmt = "%12.6f"; }
    if (elp == NULL) { elp = "[ "; }
    if (esep == NULL) { esep = " "; }
    if (erp == NULL) { erp = " ]"; }
    if (rpref == NULL) { rpref = "  "; }
    if (rsep == NULL) { rsep = " "; }
    if (rsuff == NULL) { rsuff = "\n"; }
    
    if (pref != NULL) { fputs(pref, wr); }
    for (int32_t i = 0; i < NH; i++)
      { fputs(rpref, wr);
        for (int32_t k = 0; k < 2; k++)
          { if (k != 0) { fputs(rsep, wr); }
            fputs(elp, wr);
            for (int32_t j = 0; j < NH; j++)
              { if (j != 0) { fputs(esep, wr); }
                fprintf(wr, fmt, (k == 0 ? M->dir : M->inv).c[i][j]);
              }
            fputs(erp, wr);
          }
        fputs(rsuff, wr);
      }
    if (suff != NULL) { fputs(suff, wr); }
  }

bool_t hr2_pmap_type_from_string(char *tname, hr2_pmap_type_t *type_P)
  { if ((strcmp(tname, "IDENTITY") == 0) || (strcmp(tname, "identity") == 0))
      { (*type_P) = hr2_pmap_type_IDENTITY; }
    else if ((strcmp(tname, "IDENTITY") == 0) || (strcmp(tname, "identity") == 0))
      { (*type_P) = hr2_pmap_type_IDENTITY; }
    else if ((strcmp(tname, "TRANSLATION") == 0) || (strcmp(tname, "translation") == 0))
      { (*type_P) = hr2_pmap_type_TRANSLATION; }
    else if ((strcmp(tname, "CONGRUENCE") == 0) || (strcmp(tname, "congruence") == 0))
      { (*type_P) = hr2_pmap_type_CONGRUENCE; }
    else if ((strcmp(tname, "SIMILARITY") == 0) || (strcmp(tname, "similarity") == 0))
      { (*type_P) = hr2_pmap_type_SIMILARITY; }
    else if ((strcmp(tname, "AFFINE") == 0) || (strcmp(tname, "affine") == 0))
      { (*type_P) = hr2_pmap_type_AFFINE; }
    else if ((strcmp(tname, "GENERIC") == 0) || (strcmp(tname, "projective") == 0))
      { (*type_P) = hr2_pmap_type_GENERIC; }
    else if ((strcmp(tname, "NONE") == 0) || (strcmp(tname, "none") == 0))
      { (*type_P) = hr2_pmap_type_NONE; }
    else 
      { return FALSE; }
    return TRUE;
  }

char *hr2_pmap_type_to_string(hr2_pmap_type_t type)
  { 
    switch(type)
      { case hr2_pmap_type_IDENTITY:
          return "IDENTITY";
        case hr2_pmap_type_TRANSLATION:
          return "TRANSLATION";
        case hr2_pmap_type_CONGRUENCE:
          return "CONGRUENCE";
        case hr2_pmap_type_SIMILARITY:
          return "SIMILARITY";
        case hr2_pmap_type_AFFINE:
          return "AFFINE";
        case hr2_pmap_type_NONE:
          return "NONE";
        case hr2_pmap_type_GENERIC:
          return "GENERIC";
      }
    assert(FALSE);
  }


void hr2_pmap_set_sign(hr2_pmap_t *M, sign_t sgn)
  { demand((sgn == -1) || (sgn == +1), "invalid {sgn}"); 
    double det = r3x3_det(&(M->dir));
    demand(fabs(det) > 1.0e-250, "invalid map - determinant is zero");
    if (sgn*det < 0)
      { hr2_pmap_t S = hr2_pmap_xy_swap();
        (*M) = hr2_pmap_compose(&S, M);
      }
  }

bool_t hr2_pmap_is_type(hr2_pmap_t *M, hr2_pmap_type_t type, sign_t sgn, double tol)
  { 
    double det = r3x3_det(&(M->dir));
    if (sgn*det < 0) { return FALSE; }
    
    hr2_pmap_t N = (*M);
    hr2_pmap_set_sign(&N, +1);

    switch (type)
      { case hr2_pmap_type_IDENTITY:
          return hr2_pmap_is_identity(&N, tol);
        case hr2_pmap_type_TRANSLATION:
          return hr2_pmap_is_translation(&N, tol);
        case hr2_pmap_type_CONGRUENCE:
          return hr2_pmap_is_congruence(&N, tol);
        case hr2_pmap_type_SIMILARITY:
          return hr2_pmap_is_similarity(&N, NAN, tol);
        case hr2_pmap_type_AFFINE:
          return hr2_pmap_is_affine(&N, tol);
        case hr2_pmap_type_GENERIC:
          return hr2_pmap_is_generic(&N, tol);
        default:
          demand(FALSE, "unrecognized or invalid map type");
      }
  }

