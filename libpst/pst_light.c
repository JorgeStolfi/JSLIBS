/* See pst_lamp.h */
/* Last edited on 2016-03-16 16:08:48 by stolfilocal */

#define _GNU_SOURCE
#include <stdio.h>
#include <math.h>
#include <values.h>
#include <stdlib.h>
#include <assert.h>

#include <bool.h>
#include <r3.h> 
#include <r3x3.h>
#include <affirm.h>
#include <argparser.h>

#include <pst_basic.h>
#include <pst_lamp.h>
#include <pst_light.h>
#include <pst_argparser.h>

/* INTERNAL PROTOTYPES */

void pst_light_adjust_lamp_radii(pst_lamp_vec_t *lmpv, int IS, int NS);
  /* Adjust the angular radii of each lamp {lmpv[IS..NS-1]}
    to be approximately half the average distance to its nearest
    neighbors among those lamps. */

/* IMPLEMENTATIONS */

pst_light_t pst_light_new(int NS, int NC)
  { pst_light_t lht;
    lht.lmpv = pst_lamp_vec_new(NS);
    r3_t zer = (r3_t) {{ 0, 0, 0 }};
    int i;
    for (i = 0; i < NS; i++)
      { lht.lmpv.e[i] = pst_lamp_new(NC, &zer, 0.0, +1.0); }
    pst_lamp_vec_trim(&(lht.lmpv), NS);
    return lht;
  }

pst_light_t pst_light_from_lamps(pst_lamp_vec_t lmpv)
  { pst_light_t lht;
    lht.lmpv = lmpv;
    return lht;
  }

pst_light_t pst_light_copy(pst_light_t *lht)
  { 
    /* Get lamp vectro of original light field: */
    pst_lamp_vec_t *olmpv = &(lht->lmpv);
    int NS = olmpv->ne;
    
    /* Create new light field: */
    pst_light_t new;
    pst_lamp_vec_t *nlmpv = &(new.lmpv);
    (*nlmpv) = pst_lamp_vec_new(NS);
    
    /* Copy lamps: */
    int i, c;
    for (i = 0; i < NS; i++)
      { /* Get original lamp: */
        pst_lamp_t *osrc = olmpv->e[i];
        /* Get lamp parameters from original descriptor: */
        int NC = osrc->pwr.ne;
        r3_t *dir = &(osrc->dir);
        double crad = osrc->crad;
        /* Create new lamp descriptor and set its parameters: */
        pst_lamp_t *nsrc = pst_lamp_new(NC, dir, 0.0, crad);
        for (c = 0; c < NC; c++) { nsrc->pwr.e[c] = osrc->pwr.e[c]; }
        /* Save in new light model: */
        nlmpv->e[i] = nsrc;
      }

    /* Just to be sure: */
    pst_lamp_vec_trim(nlmpv, NS);
    return new;
  }

void pst_light_ensure_one_lamp(pst_light_t *lht, double cmin, pst_lamp_t **src)
  { pst_lamp_vec_t *lmpv = &(lht->lmpv);
    int NS = lmpv->ne;

    /* Look for one bounded lamp: */
    int isrc = -1;
    int i;
    for (i = 0; i < NS; i++)
      { pst_lamp_t *src = lmpv->e[i];
        if (src->crad > cmin)
          { isrc = i; 
            cmin = src->crad;
          }
      }
      
    /* Append a lamp if missing: */
    if (isrc < 0)
      { pst_lamp_vec_expand(lmpv, NS);
        r3_t def = (r3_t){{ 0,0,1 }};
        lmpv->e[NS] = pst_lamp_new(1, &def, 0.0, cmin);
        isrc = NS; NS++;
        pst_lamp_vec_trim(lmpv, NS);
      }

    (*src) = lmpv->e[isrc];
  }

void pst_light_add_single(pst_lamp_vec_t *lmpv, int *NSP, double crad)
  { pst_lamp_vec_expand(lmpv, (*NSP));
    r3_t one = (r3_t) {{ +1, 0, 0 }};
    r3_t *dir = (crad <= -1.0 ? NULL : &one );
    lmpv->e[(*NSP)] = pst_lamp_new(1, dir, 0.0, crad);
    (*NSP)++;
  }

void pst_light_add_pair(pst_lamp_vec_t *lmpv, int *NSP, double crad)
  { /* A wall towards +X: */
    { r3_t dir = (r3_t){{ +1, 0, 0 }};
      pst_lamp_vec_expand(lmpv, (*NSP));
      lmpv->e[(*NSP)] = pst_lamp_new(1, &dir, 0.0, crad);
      (*NSP)++;
    }
    /* A wall towards -X: */
    { r3_t dir = (r3_t){{ -1, 0, 0 }};
      pst_lamp_vec_expand(lmpv, (*NSP));
      lmpv->e[(*NSP)] = pst_lamp_new(1, &dir, 0.0, crad);
      (*NSP)++;
    }
  }

void pst_light_add_tetra(pst_lamp_vec_t *lmpv, int *NSP, bool_t dual)
  { /* Even or odd vertices of a cube. */
    int i0, i1, i2;
    for (i0 = 0; i0 < 2; i0++)
      { for (i1 = 0; i1 < 2; i1++)
          { for (i2 = 0; i2 < 2; i2++)
              { int parity = (i0 + i1 + i2) % 2;
                if ((parity == 0) == dual)
                  { /* A corner of the cube: */
                    r3_t dir = (r3_t){{ 2*i0-1, 2*i1-1, 2*i2-1 }};
                    (void)r3_dir(&dir, &dir);
                    pst_lamp_vec_expand(lmpv, (*NSP));
                    lmpv->e[(*NSP)] = pst_lamp_new(1, &dir, 0.0, 1.0);
                    (*NSP)++;
                  }
              }
          }
      }
  }

void pst_light_add_octa(pst_lamp_vec_t *lmpv, int *NSP)
  { /* Vertices on the coordinates axes. */
    int i, s;
    for (i = 0; i < 3; i++)
      { for (s = 0; s < 2; s++)
          { /* A vertex of the octahedron: */
            r3_t dir = (r3_t){{ 0, 0, 0 }};
            dir.c[i] = 2*s-1;
            (void)r3_dir(&dir, &dir);
            pst_lamp_vec_expand(lmpv, (*NSP));
            lmpv->e[(*NSP)] = pst_lamp_new(1, &dir, 0.0, 1.0);
            (*NSP)++;
          }
      }
  }

void pst_light_add_icosa(pst_lamp_vec_t *lmpv, int *NSP)
  { /* Stabbed by the axes across the midpoints of opposite edges. */
    double phi = (sqrt(5)-1)/2;
    int i, s, t;
    for (i = 0; i < 3; i++)
      { for (s = 0; s < 2; s++)
          { for (t = 0; t < 2; t++)
              { /* A vertex of the icosahedron: */
                r3_t dir = (r3_t){{ 0, 0, 0 }};
                dir.c[i] = 2*s-1;
                dir.c[(i+1)%3] = (2*t-1)*phi;
                (void)r3_dir(&dir, &dir);
                pst_lamp_vec_expand(lmpv, (*NSP));
                lmpv->e[(*NSP)] = pst_lamp_new(1, &dir, 0.0, 1.0);
                (*NSP)++;
              }
          }
      }
  }

void pst_light_add_dodeca(pst_lamp_vec_t *lmpv, int *NSP)
  { /* Stabbed by the axes across the midpoints of opposite edges. */
    /* Two vertices on each face of a cube: */
    double phi = (sqrt(5)-1)/2;
    int i, s, t;
    for (i = 0; i < 3; i++)
      { for (s = 0; s < 2; s++)
          { for (t = 0; t < 2; t++)
              { /* A vertex of the dodecahedron: */
                r3_t dir = (r3_t){{ 0, 0, 0 }};
                dir.c[i] = 2*s-1;
                dir.c[(i+2)%3] = (2*t-1)*(2+phi);
                (void)r3_dir(&dir, &dir);
                pst_lamp_vec_expand(lmpv, (*NSP));
                lmpv->e[(*NSP)] = pst_lamp_new(1, &dir, 0.0, 1.0);
                (*NSP)++;
              }
          }
      }
    /* The corners of a cube: */
    int i0, i1, i2;
    for (i0 = 0; i0 < 2; i0++)
      { for (i1 = 0; i1 < 2; i1++)
          { for (i2 = 0; i2 < 2; i2++)
              { /* A corner of the cube: */
                r3_t dir = (r3_t){{ 2*i0-1, 2*i1-1, 2*i2-1 }};
                (void)r3_dir(&dir, &dir);
                pst_lamp_vec_expand(lmpv, (*NSP));
                lmpv->e[(*NSP)] = pst_lamp_new(1, &dir, 0.0, 1.0);
                (*NSP)++;
              }
          }
      }
  }

void pst_light_add_uniform_array
  ( pst_lamp_vec_t *lmpv, 
    int *NSP, 
    int NA,
    r3_t *dir0, 
    r3_t *dir1, 
    double_vec_t *pwr,
    double crad
  )
  { 
    int IS = (*NSP); /* Save initial number of lamps. */
    switch(NA)
      {
      case 1: 
        /* Single ambient lamp. */
        pst_light_add_single(lmpv, NSP, -1.0); 
        break;
      case 2: 
        /* Pair of opposite wall lamps. */
        pst_light_add_pair(lmpv, NSP, 0.0); 
        break;
      case 3: 
        /* Pair of opposite point lamps and an ambient light. */
        pst_light_add_single(lmpv, NSP, -1.0); 
        pst_light_add_pair(lmpv, NSP, 0.0); 
        break;
      case 4: 
        /* Vertices of a tetrahedron. */
        pst_light_add_tetra(lmpv, NSP, FALSE); 
        break;
      case 6: 
        /* Vertices of an octahedron. */
        pst_light_add_octa(lmpv, NSP); 
        break;
      case 8: 
        /* Vertices of a cube. */
        pst_light_add_tetra(lmpv, NSP, FALSE); 
        pst_light_add_tetra(lmpv, NSP, TRUE); 
        break;
      case 12: 
        /* Vertices of an icosahedron. */
        pst_light_add_icosa(lmpv, NSP); 
        break;
      case 14: 
        /* Octahedron and cube. */
        pst_light_add_octa(lmpv, NSP); 
        pst_light_add_tetra(lmpv, NSP, FALSE); 
        pst_light_add_tetra(lmpv, NSP, TRUE); 
        break;
      case 20: 
        /* Dodecahedron. */
        pst_light_add_dodeca(lmpv, NSP); 
        break;
      case 32: 
        /* Dodecahedron and icosahedron. */
        pst_light_add_icosa(lmpv, NSP); 
        pst_light_add_dodeca(lmpv, NSP); 
        break;
      default:
        /* Not implemented yet */
        demand(FALSE, "cannot generate mode with that many lamps");
      }
    assert(IS + NA == (*NSP));

    /* Adjust lamp parameters as appropriate: */
    int i;
    if ((dir0 != NULL) && (r3_L_inf_norm(dir0) != 0.0))
      { /* Adjust directions of all lamps in array. */
        /* Get directions of first two lamps: */
        r3_t *sdir0 = (IS < (*NSP) ? &(lmpv->e[IS]->dir) : NULL); 
        r3_t *sdir1 = (IS+1 < (*NSP) ? &(lmpv->e[IS+1]->dir) : NULL);
        /* Compute the rotation matrices: */
        r3x3_t S, D;
        pst_light_alignment_matrix(sdir0, sdir1, &S);
        pst_light_alignment_matrix(dir0, dir1, &D);
        /* Apply it to all lamps: */
        for (i = IS; i < (*NSP); i++)
          { pst_lamp_t *src = lmpv->e[i];
            r3_t tdir;
            /* Map {src.dir} by {S}, then by the inverse of {D}: */
            r3x3_map_row(&(src->dir), &S, &(tdir));
            r3x3_map_col(&D, &(tdir), &(src->dir));
          }
      }
    
    if ((pwr != NULL) && (pwr->ne != 0))
      { /* Set all power vectors to copies of {pwr}. */
        int KC = pwr->ne;
        int c;
        for (i = IS; i < (*NSP); i++)
          { pst_lamp_t *src = lmpv->e[i];
            src->pwr = double_vec_new(KC);
            for (c = 0; c < KC; c++) { src->pwr.e[c] = pwr->e[c]; }
          }
      }
    
    if (fabs(crad) != INF)
      { /* Set all cosine-radius fields to {crad}. */
        for (i = IS; i < (*NSP); i++)
          { pst_lamp_t *src = lmpv->e[i];
            src->crad = crad;
          }
      }
    else if (NA >= 4)
      { /* Set the radii according to neighbor's distances: */
        pst_light_adjust_lamp_radii(lmpv, IS, (*NSP));
      }
  }
      
void pst_light_alignment_matrix(r3_t *u, r3_t *v, r3x3_t *M)
  {
    r3x3_ident(M); /* By default. */
    if ((u == NULL) || (r3_L_inf_norm(u) == 0.0))
      { /* If {u} is indeterminate, there's nothing to do. */ }
    else
      { /* Set {M} to the shortest rotation matrix that sends {u} to {x}. */
        /* Get the unit +X vector {x}: */
        r3_t x = (r3_t) {{ 1, 0, 0 }};
        /* Make sure {u} has unit length: */
        r3_t uu; (void)r3_dir(u, &uu);
        /* We need the the direction {p} of {u + x}: */
        r3_t p; r3_add(&uu, &x, &p);
        if (r3_L_inf_norm(&p) < 1.0e-8)
          { /* The vectors are practically opposite, rotate 180 deg around Z: */
            M->c[0][0] = M->c[1][1] = -1;
          }
        else
          { /* Scale {p} to unit norm: */
            (void)r3_dir(&p, &p);
            /* Set {M} to a mirror in direction {p}, followed by negation of {x}: */
            int i, j;
            for (i = 0; i < 3; i++)
              { for (j = 0; j < 3; j++) { M->c[i][j] += - 2*p.c[i]*p.c[j]; }
                M->c[i][0] = - M->c[i][0];
              }
          }

        /* Now fixt the effect on {v}: */
        if ((v == NULL) && (r3_L_inf_norm(v) == 0.0))
          { /* Nothing else to do. */ }
        else
          { /* We need a vector perpendicular to {u,v}: */
            r3_t w; r3_cross(u, v, &w);
            if (r3_L_inf_norm(v) == 0.0)
              { /* {u} and {v} are practically collinear, nothing to do. */ }
            else
              { /* Map {w} by {M}, result is {q}: */
                r3_t q; r3x3_map_row(&w, M, &q);
                /* Scale {q} to unit norm: */
                (void)r3_dir(&q, &q);
                /* The vector {q} should be perpendicular to {x}: */
                assert(fabs(q.c[0]) < 1.0e-7);
                /* Make a rotation matrix around {x} that takes {q} to {(0,0,1)}: */
                r3x3_t V; r3x3_ident(&V);
                V.c[1][1] = V.c[2][2] = q.c[2];
                V.c[1][2] = +q.c[1]; 
                V.c[2][1] = -q.c[1]; 
                /* Compose with {M} (first {M} then {V}): */
                r3x3_mul(M, &V, M);
              }
          }
      }
  }

void pst_light_regularize_channels(pst_light_t *lht, int NC)
  { pst_lamp_vec_t *lmpv = &(lht->lmpv);
    int NS = lmpv->ne;
    if (NS > 0)
      { double defpwr = 1.0/NS;
        int i;
        for (i = 0; i < NS; i++)
          { pst_lamp_t *src = lmpv->e[i];
            pst_double_vec_regularize(&(src->pwr), NC, defpwr);
            /* The regularization must have succeded: */
            assert(src->pwr.ne == NC);
          }
      }
  }

void pst_light_adjust_lamp_radii(pst_lamp_vec_t *lmpv, int IS, int NS)
  { int i, j;
    for (i = IS; i < NS; i++)
      { pst_lamp_t *srci = lmpv->e[i];
        /* Find mean distance to nearest neighbors: */
        double twsum = 0, wsum = 0;
        for (j = IS; j < NS; j++)
          { /* Check all other lamps: */
            if (j == i) { continue; }
            pst_lamp_t *srcj = lmpv->e[j];
            /* Compute angle {ang} between the two: */
            double cang = r3_dot(&(srci->dir), &(srcj->dir));
            demand (cang < 1.0, "two lights with same direction"); 
            double ang = acos(fmin(+1.0, fmax(-1.0, cang)));
            if (cang > - 1.0 + 1.0e-5) 
              { /* Consider the lamp's disk as extending to the midpoint. */
                double thaf = tan(ang/2);
                /* Assign to this radius a weight that decays quickly with distance: */
                double wt = exp(-thaf)/thaf;
                twsum += thaf*wt;
                wsum += wt;
              }
          }
        demand(wsum > 0, "lights are too far apart");
        /* Radius is the weighted average: */
        double trad = twsum/wsum;
        double rad = atan(trad);
        assert(rad > 0);
        srci->crad = cos(rad);
      }
  }

pst_light_t pst_light_spec_parse(argparser_t *pp, bool_t next, int *NCP)
  { 
    int NS = 0; /* Number of light sources. */
    pst_lamp_vec_t lmpv = pst_lamp_vec_new(0);
    while(TRUE)
      { if (pst_keyword_present(pp, "-lamp", next))
          { pst_lamp_t *src = pst_lamp_new(0, NULL, 0.0, +INF);
            (*src) = pst_lamp_spec_params_next_parse(pp, NCP);
            /* pst_lamp_spec_write(stderr, src); */
            pst_lamp_vec_expand(&lmpv, NS);
            lmpv.e[NS] = src;
            NS++;
          }
        else if (pst_keyword_present(pp, "-array", next))
          { int NA = (int)argparser_get_next_int(pp, 1, MAXINT);
            pst_lamp_t def = pst_lamp_spec_params_next_parse(pp, NCP);
            pst_light_add_uniform_array
              ( &lmpv, &NS, NA, &(def.dir), NULL, &(def.pwr), def.crad );
          }
        else
          {
            /* Assume it is the end of the light field spec: */
            break;
          }
        /* The {next} argument applies only to the first lamp/array: */
        next = FALSE;
      }
    pst_lamp_vec_trim(&lmpv, NS);
    return pst_light_from_lamps(lmpv);
  }

void pst_light_spec_write(FILE *wr, pst_light_t *lht)
  { pst_lamp_vec_t *lmpv = &(lht->lmpv);
    int NS = lmpv->ne;
    int i;
    for(i = 0; i < NS; i++)
      { pst_lamp_spec_write(wr, lmpv->e[i]); }
  }

vec_typeimpl(pst_light_vec_t,pst_light_vec,pst_light_t);
