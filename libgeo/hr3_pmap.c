/* See hr3_pmap.h */
/* Last edited on 2025-01-05 00:35:21 by stolfi */ 

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

#include <hr3_pmap.h>

#define NH 4
  /* Number of homogeneous coordinates in a point or coefficients in a plane. */

#define NG 6
  /* Number of Grassman coordinates of a line. */

#define NC 3
  /* Number of Cartesian coordinates in a point. */

bool_t hr3_pmap_is_identity(hr3_pmap_t *M)
  {
    return r4x4_is_unif_scaling(&(M->dir), M->dir.c[0][0]); 
  }

hr3_point_t hr3_pmap_point(hr3_point_t *p, hr3_pmap_t *M)
  {
    r4_t c;
    r4x4_map_row(&(p->c), &(M->dir), &c);
    return (hr3_point_t){c};
  }
        
hr3_point_t hr3_pmap_inv_point(hr3_point_t *p, hr3_pmap_t *M)
  {
    r4_t c;
    r4x4_map_row(&(p->c), &(M->inv), &c);
    return (hr3_point_t){c};
  }
        
hr3_plane_t hr3_pmap_plane(hr3_plane_t *P, hr3_pmap_t *M)
  {
    r4_t f;
    r4x4_map_col(&(M->inv), &(P->f), &f);
    return (hr3_plane_t){f};
  }
    
hr3_plane_t hr3_pmap_inv_plane(hr3_plane_t *P, hr3_pmap_t *M)
  {
    r4_t f;
    r4x4_map_col(&(M->dir), &(P->f), &f);
    return (hr3_plane_t){f};
  }

r3_t hr3_pmap_r3_point(r3_t *p, hr3_pmap_t *M)
  {
    hr3_point_t ph = (hr3_point_t){{{ 1.0, p->c[0], p->c[1], p->c[2] }}};
    hr3_point_t qh;
    r4x4_map_row(&(ph.c), &(M->dir), &(qh.c));
    r3_t qc = r3_from_hr3_nan(&qh);
    return qc;
  }

r3_t hr3_pmap_inv_r3_point(r3_t *p, hr3_pmap_t *M)
  {
    hr3_point_t ph = (hr3_point_t){{{ 1, p->c[0], p->c[1], p->c[2] }}};
    hr3_point_t qh;
    r4x4_map_row(&(ph.c), &(M->inv), &(qh.c));
    r3_t qc = r3_from_hr3_nan(&qh);
    return qc;
  }

hr3_pmap_t hr3_pmap_identity(void)
  { hr3_pmap_t M;
    r4x4_ident(&(M.dir));
    M.inv = M.dir;
    return M;
  }

hr3_pmap_t hr3_pmap_translation(r3_t *v)
  {
    hr3_pmap_t M; 
    r4x4_ident(&(M.dir));
    M.dir.c[0][1] = +v->c[0];
    M.dir.c[0][2] = +v->c[1];
    M.dir.c[0][3] = +v->c[2];
    r4x4_ident(&(M.inv));
    M.inv.c[0][1] = -v->c[0];
    M.inv.c[0][2] = -v->c[1];
    M.inv.c[0][3] = -v->c[2];
    return M;
  }

hr3_pmap_t hr3_pmap_scaling(r3_t *scale)
  {
    hr3_pmap_t M;
    for (uint32_t i = 0;  i < NH; i++)
      { for (uint32_t j = 0;  j < NH; j++)
          { if ((i == 0) && (j == 0))
              { M.dir.c[i][j] = M.inv.c[i][j] = 1.0; }
            else if (i == j)
              { M.dir.c[i][j] = scale->c[j-1];
                M.inv.c[i][j] = 1/scale->c[j-1];
              }
            else
              { M.dir.c[i][j] = M.inv.c[i][j] = 0.0; }
          }
      }
    return M;
  
  }

hr3_pmap_t hr3_pmap_u_to_v_rotation(r3_t *u, r3_t *v)
  { /* Compute the matrix {r} for {R^3} to {R^3}: */
    r3x3_t r;
    r3x3_u_to_v_rotation(u, v, &r);
    
    /* Convert to a projective map (note that inverse is just transpose): */
    hr3_pmap_t M;
    for (uint32_t i = 0;  i < NH; i++)
      { for (uint32_t j = 0;  j < NH; j++)
          { if ((i == 0) && (j == 0))
              { M.dir.c[i][j] = M.inv.c[i][j] = 1.0; }
            else if ((i == 0) || (j == 0))
              { M.dir.c[i][j] = M.inv.c[i][j] = 0.0; }
            else
              { M.dir.c[i][j] = r.c[i-1][j-1]; 
                M.inv.c[i][j] = r.c[j-1][i-1]; 
              }
          }
      }
    return M;
  }

hr3_pmap_t hr3_pmap_compose(hr3_pmap_t *M, hr3_pmap_t *N)
  {
    r4x4_t dir, inv; 
    r4x4_mul(&(M->dir), &(N->dir), &dir);
    r4x4_mul(&(N->inv), &(M->inv), &inv);
    return (hr3_pmap_t){dir, inv};
  }    
  
hr3_pmap_t hr3_pmap_inv(hr3_pmap_t *M)
  {
    return (hr3_pmap_t){M->inv, M->dir};
  }

hr3_pmap_t hr3_pmap_inv_compose(hr3_pmap_t *M, hr3_pmap_t *N)
  {
    r4x4_t dir, inv;
    r4x4_mul(&(M->inv), &(N->dir), &dir);
    r4x4_mul(&(N->inv), &(M->dir), &inv);
    return (hr3_pmap_t){dir, inv};
  }

hr3_pmap_t hr3_pmap_aff_from_mat_and_disp(r3x3_t *E, r3_t *d)

  {
    hr3_pmap_t M;

    r3x3_t F; r3x3_inv(E, &(F));
    r3_t p;
    r3x3_map_row(d, &(F), &p);

    for (uint32_t i = 0;  i < NH; i++)
      { for (uint32_t j = 0;  j < NH; j++)
          { double Pij, Qij; /* Elements of {M.dir} and {M.inv}. */
            if (i == 0)
              { if (j == 0)
                  { Pij = 1.0; Qij = 1.0; }
                else
                  { Pij = d->c[j-1]; Qij = p.c[j-1]; }
              }
            else
              { if (j == 0)
                  { Pij = 0.0; Qij = 0.0; }
                else
                  { Pij = E->c[i-1][j-1]; Qij = F.c[i-1][j-1]; }
              }
            M.dir.c[i][j] = Pij; 
            M.inv.c[i][j] = Qij;
          }
      }
    
    return M;
  }
  
hr3_pmap_t hr3_pmap_aff_from_four_points(r3_t *o, r3_t *p, r3_t *q, r3_t *r)
  {
    hr3_pmap_t M; 

    for (uint32_t j = 0;  j < NH; j++)
      { M.dir.c[0][j] = (j == 0 ? 1.0 : o->c[j-1]);
        M.dir.c[1][j] = (j == 0 ? 0.0 : p->c[j-1] - o->c[j-1]);
        M.dir.c[2][j] = (j == 0 ? 0.0 : q->c[j-1] - o->c[j-1]);
        M.dir.c[3][j] = (j == 0 ? 0.0 : r->c[j-1] - o->c[j-1]);
      }

    r4x4_inv(&(M.dir), &(M.inv));
    
    return M;
  }
  
hr3_pmap_t hr3_pmap_from_five_points(hr3_point_t *p, hr3_point_t *q, hr3_point_t *r, hr3_point_t *s, hr3_point_t *u)
  {
    hr3_pmap_t M; /* The resulting map. */
    
    /* Compute weights {(a,b,c)=w.c[0..2]} for rows of {Q} that map {[1,1,1]} to {u}: */
    r4_t w;
    { /* Compute a matrix {Q} that maps the cardinal points to {p,q,r,s} as given: */
      r4x4_t Q;
      for (uint32_t j = 0;  j < NH; j++)
        { M.dir.c[0][j] = p->c.c[j];
          M.dir.c[1][j] = q->c.c[j];
          M.dir.c[2][j] = r->c.c[j];
          M.dir.c[3][j] = s->c.c[j];
        }
      /* Map {u} by the inverse of {Q}: */
      r4x4_inv(&(M.dir), &Q);
      w = u->c;
      r4x4_map_row(&w, &Q, &w);
    }
    
    /* Make the weights positive, so that {p,q,r,s} are strictly honored: */
    for (uint32_t i = 0;  i < NH; i++) { w.c[i] = fabs(w.c[i]); }

    /* Ensure that {M.dir} maps the cardinal points to {p,q,r,s} and some unit point to {u}: */
    for (uint32_t i = 0;  i < NH; i++)
      { for (uint32_t j = 0;  j < NH; j++)
          { M.dir.c[i][j] *= w.c[i];  }
      }

    /* Compute the inverse map: */
    r4x4_inv(&(M.dir), &(M.inv));

    return M;
  }

hr3_pmap_t hr3_pmap_persp(hr3_point_t *obs, hr3_point_t *foc, double rad, hr3_point_t *upp)
  { hr3_pmap_t M;
    r4x4_t Mt, Mr, Mc, Mrt;
    
    demand(foc->c.c[0] > 0.0, "focus must be finite and hither");
    demand(obs->c.c[0] >= 0.0, "observer must be hither or infinite");
    demand(upp->c.c[0] >= 0.0, "zenith must be hither or infinite"); 

    /* Start with a translation from {foc} to the origin: */
    for (uint32_t i = 0;  i < NH; i++){
      for (uint32_t j = 0;  j < NH; j++){
        if (i == j)
          { Mt.c[i][j] = foc->c.c[0]; }
        else if (i == 0)
          { Mt.c[i][j] = -(foc->c.c[j]); }
        else
          { Mt.c[i][j] = 0.0; }
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
    Mr.c[0][0] = 1.0;
    for (uint32_t i = 1;  i < NH; i++)
      { Mr.c[0][i] = 0.0;
        Mr.c[i][0] = 0.0;
        Mr.c[i][1] = r.c[i-1];
        Mr.c[i][2] = s.c[i-1];
        Mr.c[i][3] = t.c[i-1];
      }
    r4x4_mul(&Mt, &Mr, &Mrt);

    /* Do we need a conical projection step? */
    if (obs->c.c[0] == 0.0)
      { /* Observer is at infinity; cilindrical projection. */
        M.dir = Mrt;
      }
    else
      { /* Observer is finite; add conical projection step. */
        double d = hr3_dist(foc, obs);
        double uno = -1.0;
        if (fabs(d) > 1.0) { uno = -1.0/d; d = 1.0; }
        for (uint32_t i = 0;  i < NH; i++)
          for (uint32_t j = 0;  j < NH; j++)
            { if (i == j)
                { Mc.c[i][j] = d; }
              else
                { Mc.c[i][j] = 0.0; }
            }
        Mc.c[3][0] = uno;
        r4x4_mul(&Mrt, &Mc, &(M.dir));
      }

    /* Asked for scaling? */
    if (rad > 0.0)
      { /* Combine {M} with a uniform scale of {1/rad}: */
        for (uint32_t i = 0;  i < NH; i++)
          for (uint32_t j = 1;  j < NH; j++)
            { M.dir.c[i][j] /= rad; }
      }

    /* Compute inverse matrix: */
    r4x4_inv(&M.dir, &M.inv);
    return M;
  }

bool_t hr3_pmap_is_affine(hr3_pmap_t *M)
  { r4x4_t *A = &(M->dir);
    if (A->c[0][0] <= 0) { return FALSE; }
    if ((A->c[1][0] != 0.0) || (A->c[2][0] != 0.0) || (A->c[3][0] != 0.0))
      { return FALSE; }
    return TRUE;
  }
  
double hr3_pmap_diff_sqr(hr3_pmap_t *M, hr3_pmap_t *N)
  { double sum_d2 = 0;
    for (uint32_t sense = 0;  sense < 2; sense++)
      { r4x4_t *A = (sense == 0 ? &(M->dir) : &(M->inv));
        double Am = r4x4_norm(A) + 1.0e-200;
        r4x4_t *B = (sense == 0 ? &(N->dir) : &(N->inv));
        double Bm = r4x4_norm(B) + 1.0e-200;
        for (uint32_t i = 0;  i < NH; i++)
          { for (uint32_t j = 0;  j < NH; j++)
             { double Aij = A->c[i][j]/Am;
               double Bij = B->c[i][j]/Bm;
               double dij = Aij - Bij;
               sum_d2 += dij*dij;
             }
          }
      }
    return sum_d2;
  }

double hr3_pmap_mismatch_sqr(hr3_pmap_t *M, uint32_t np, r3_t p1[], r3_t p2[])
  {
    bool_t debug = FALSE;
    
    double sum2 = 0.0;
    for (uint32_t k = 0;  k < np; k++)
      { r3_t *p1k = &(p1[k]);
        r3_t *p2k = &(p2[k]);
        r3_t q1k = hr3_pmap_r3_point(p1k, M);
        r3_t q2k = hr3_pmap_inv_r3_point(p2k, M);
        double d2 = r3_dist_sqr(&q1k, &q2k);
        if (debug)
          { r3_gen_print(stderr, p1k, "%+8.5f", "  @ ( ", " ", " )");
            r3_gen_print(stderr, &q1k, "%+12.8f", " -> ( ", " ", " )");
            fprintf(stderr, " |%.6f| ", d2);
            r3_gen_print(stderr, &q2k, "%+12.8f", "( ", " ", " ) <- ");
            r3_gen_print(stderr, p2k, "%+8.5f", "( ", " ", " )\n");
          }
        sum2 += d2;
      }
    return sum2/np;
  }

double hr3_pmap_deform_sqr(r3_t ph[], hr3_pmap_t *M)
  {
    uint32_t nk = (1 << NC); /* Number of corners of the cuboid. */
    r3_t qh[nk];
    for (uint32_t k = 0;  k < nk; k++)
      { qh[k] = hr3_pmap_r3_point(&(ph[k]), M); }
    
    uint32_t nd = 12 + 4; /* Number of distances to probe. */
    double logr[nd];
    uint32_t kd = 0;
    for (uint32_t ik = 1;  ik < nk; ik++)
      { for (uint32_t jk = 0;  jk < ik; jk++)
          { uint32_t eij = (uint32_t)(ik ^ jk); /* Exclusive OR of indices. */
            uint32_t hij = (eij & 1) + (eij & 2) + (eij & 4); /* Hamming dist. */
            if ((hij == 1) || (hij == 3))
              { /* Side or main diagonal: */
                double dp2 = r3_dist_sqr(&(ph[ik]), &(ph[jk]));
                double dq2 = r3_dist_sqr(&(qh[ik]), &(qh[jk]));
                assert(kd < nd);
                logr[kd] = log(dq2/dp2)/2;
                kd++;
              }
          }
      }
    assert(kd == nd);
    
    /* Compute the variance of the logs: */
    double sum = 0;
    for (uint32_t kd = 0;  kd < nd; kd++) { sum += logr[kd]; }
    double avg = sum/nd;
    double sum2 = 0;
    for (uint32_t kd = 0;  kd < nd; kd++) { double dk = logr[kd] - avg; sum2 += dk*dk; }
    double var = sum2/(nd-1);
    return var;
  }

double hr3_pmap_aff_discr_sqr(hr3_pmap_t *M, hr3_pmap_t *N)
  {
    demand((M->dir.c[1][0] == 0) && (M->dir.c[2][0] == 0), "{M} is not affine");
    demand(M->dir.c[0][0] > 0, "map {M} does not preserve side");
    r4x4_t A; double wA = M->dir.c[0][0]; r4x4_scale(1/wA, &(M->dir), &A);
   
    demand((N->dir.c[1][0] == 0) && (N->dir.c[2][0] == 0), "{N} is not affine");
    demand(M->dir.c[0][0] > 0, "map {N} does not preserve side");
    r4x4_t B; double wB = N->dir.c[0][0]; r4x4_scale(1/wB, &(N->dir), &B);
   
    /* Hope the math is right: */
    r4x4_t E, H;
    r4x4_sub(&A, &B, &E);
    r4x4_mul_tr(&E, &E, &H);
    double h2 = (H.c[1][1] + H.c[2][2])/2;
    double d2 = H.c[0][0];
    return h2 + d2;
  }
  
void hr3_pmap_print(FILE *wr, hr3_pmap_t *M, char *pref, char *suff)
  { 
    hr3_pmap_gen_print(wr, M, "%12.6f", pref, "  ", "  ", "\n", "[ ", " ", " ]", suff);
    fflush(wr);
  }

void hr3_pmap_gen_print
  ( FILE *wr, hr3_pmap_t *M,
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
    for (uint32_t i = 0;  i < NH; i++)
      { fputs(rpref, wr);
        for (uint32_t k = 0;  k < 2; k++)
          { if (k != 0) { fputs(rsep, wr); }
            fputs(elp, wr);
            for (uint32_t j = 0;  j < NH; j++)
              { if (j != 0) { fputs(esep, wr); }
                fprintf(wr, fmt, (k == 0 ? M->dir : M->inv).c[i][j]);
              }
            fputs(erp, wr);
          }
        fputs(rsuff, wr);
      }
    if (suff != NULL) { fputs(suff, wr); }
  }
