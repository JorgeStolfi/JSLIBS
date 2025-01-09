/* See pst_lamp.h */
/* Last edited on 2025-01-05 10:28:31 by stolfi */

#include <stdio.h>
#include <stdint.h>
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

void pst_light_adjust_lamp_radii(pst_lamp_t *lmpv[], uint32_t NL);
  /* Adjust the angular radii of each lamp {lmpv[0..NL-1]}
    to be approximately half the average distance to its nearest
    neighbors among those lamps. */

void pst_light_set_default_power(pst_lamp_t *lmpv[], uint32_t NL);
  /* Sets the power of each lamp that has {NAN} power to {1/NL}. */


/* IMPLEMENTATIONS */

pst_light_t *pst_light_new(uint32_t NL)
  { pst_light_t *lht = talloc(1, pst_light_t);
    lht->lmpv = pst_lamp_vec_new(NL);
    frgb_t zer = (frgb_t) {{ 0, 0, 0 }};
    for (uint32_t i = 0; i < NL; i++)
      { lht->lmpv.e[i] = pst_lamp_new(NULL, &zer, +1.0); }
    lht->NL = NL;
    return lht;
  }

pst_light_t *pst_light_from_lamps(uint32_t NL, pst_lamp_vec_t lmpv)
  { pst_light_t *lht = talloc(1, pst_light_t);
    demand(NL <= lmpv.ne, "invalid source count {NL}");
    lht->lmpv = lmpv;
    lht->NL = NL;
    return lht;
  }

pst_light_t *pst_light_copy(pst_light_t *lht)
  { 
    
    /* Get lamp vectro of original light field: */
    pst_lamp_vec_t *olmpv = &(lht->lmpv);
    uint32_t NL = olmpv->ne;
    
    /* Create new light field: */
    pst_light_t *new = talloc(1, pst_light_t);
    new->lmpv = pst_lamp_vec_new(NL);
     
    /* Copy lamps: */
    for (uint32_t i = 0; i < NL; i++)
      { /* Get original lamp: */
        pst_lamp_t *osrc = lht->lmpv.e[i];
        /* Get lamp parameters from original descriptor: */
        r3_t *dir = &(osrc->dir);
        frgb_t *pwr = &(osrc->pwr);
        double crad = osrc->crad;
        pst_lamp_t *nsrc = pst_lamp_new(dir, pwr, crad);
        new->lmpv.e[i] = nsrc;
      }

    /* Just to be sure: */
    pst_lamp_vec_trim(&(new->lmpv), NL);
    return new;
  }

void pst_light_ensure_one_lamp(pst_light_t *lht, double cmin, pst_lamp_t **src)
  { pst_lamp_vec_t *lmpv = &(lht->lmpv);
    uint32_t NL = lmpv->ne;

    /* Look for one bounded lamp: */
    int32_t isrc = -1;
    for (uint32_t i = 0; i < NL; i++)
      { pst_lamp_t *src = lmpv->e[i];
        if (src->crad > cmin)
          { isrc = (int32_t)i; 
            cmin = src->crad;
          }
      }
      
    /* Append a lamp if missing: */
    if (isrc < 0)
      { pst_lamp_vec_expand(lmpv, (vec_index_t)NL);
        r3_t zdir = (r3_t){{ 0, 0, 1 }};
        frgb_t pwr = (frgb_t){{ 1, 1, 1 }};
        lmpv->e[NL] = pst_lamp_new(&zdir, &pwr, cmin);
        isrc = (int32_t)NL; NL++;
      }

    (*src) = lmpv->e[isrc];
  }

void pst_light_add_single(pst_light_t *lht, double crad)
  { pst_lamp_vec_t *lmpv = &(lht->lmpv);
    pst_lamp_vec_expand(lmpv, (vec_index_t)lht->NL);
    r3_t one = (r3_t) {{ +1, 0, 0 }};
    r3_t *dir = (crad <= -1.0 ? NULL : &one );
    lmpv->e[lht->NL] = pst_lamp_new(dir, NULL, crad);
    lht->NL++;
  }

void pst_light_add_pair(pst_light_t *lht, double crad)
  { /* A wall towards +X: */
    pst_lamp_vec_t *lmpv = &(lht->lmpv);
    { r3_t dir = (r3_t){{ +1, 0, 0 }};
      pst_lamp_vec_expand(lmpv, (vec_index_t)lht->NL);
      lmpv->e[lht->NL] = pst_lamp_new(&dir, NULL, crad);
      lht->NL++;
    }
    /* A wall towards -X: */
    { r3_t dir = (r3_t){{ -1, 0, 0 }};
      pst_lamp_vec_expand(lmpv, (vec_index_t)lht->NL);
      lmpv->e[lht->NL] = pst_lamp_new(&dir, NULL, crad);
      lht->NL++;
    }
  }

void pst_light_add_tetra(pst_light_t *lht, bool_t dual)
  { /* Even or odd vertices of a cube. */
    pst_lamp_vec_t *lmpv = &(lht->lmpv);
    for (int32_t i0 = 0; i0 < 2; i0++)
      { for (int32_t i1 = 0; i1 < 2; i1++)
          { for (int32_t i2 = 0; i2 < 2; i2++)
              { int32_t parity = (i0 + i1 + i2) % 2;
                if ((parity == 0) == dual)
                  { /* A corner of the cube: */
                    r3_t dir = (r3_t){{ 2*i0-1, 2*i1-1, 2*i2-1 }};
                    (void)r3_dir(&dir, &dir);
                    pst_lamp_vec_expand(lmpv, (vec_index_t)lht->NL);
                    lmpv->e[lht->NL] = pst_lamp_new(&dir, NULL, 1.0);
                    lht->NL++;
                  }
              }
          }
      }
  }

void pst_light_add_octa(pst_light_t *lht)
  { /* Vertices on the coordinates axes. */
    pst_lamp_vec_t *lmpv = &(lht->lmpv);
    for (int32_t i = 0; i < 3; i++)
      { for (int32_t s = 0; s < 2; s++)
          { /* A vertex of the octahedron: */
            r3_t dir = (r3_t){{ 0, 0, 0 }};
            dir.c[i] = 2*s-1;
            (void)r3_dir(&dir, &dir);
            pst_lamp_vec_expand(lmpv, (vec_index_t)lht->NL);
            lmpv->e[lht->NL] = pst_lamp_new(&dir, NULL, 1.0);
            lht->NL++;
          }
      }
  }

void pst_light_add_icosa(pst_light_t *lht)
  { /* Stabbed by the axes across the midpoints of opposite edges. */
    pst_lamp_vec_t *lmpv = &(lht->lmpv);
    double phi = (sqrt(5)-1)/2;
    for (int32_t i = 0; i < 3; i++)
      { for (int32_t s = 0; s < 2; s++)
          { for (int32_t t = 0; t < 2; t++)
              { /* A vertex of the icosahedron: */
                r3_t dir = (r3_t){{ 0, 0, 0 }};
                dir.c[i] = 2*s-1;
                dir.c[(i+1)%3] = (2*t-1)*phi;
                (void)r3_dir(&dir, &dir);
                pst_lamp_vec_expand(lmpv, (vec_index_t)lht->NL);
                lmpv->e[lht->NL] = pst_lamp_new(&dir, NULL, 1.0);
                lht->NL++;
              }
          }
      }
  }

void pst_light_add_dodeca(pst_light_t *lht)
  { /* Stabbed by the axes across the midpoints of opposite edges. */
    /* Icosahedron vertices, two on each face of a cube: */
    pst_lamp_vec_t *lmpv = &(lht->lmpv);
    double phi = (sqrt(5)-1)/2;
    for (int32_t i = 0; i < 3; i++)
      { for (int32_t s = 0; s < 2; s++)
          { for (int32_t t = 0; t < 2; t++)
              { /* A vertex of the dodecahedron: */
                r3_t dir = (r3_t){{ 0, 0, 0 }};
                dir.c[i] = 2*s-1;
                dir.c[(i+2)%3] = (2*t-1)*(2+phi);
                (void)r3_dir(&dir, &dir);
                pst_lamp_vec_expand(lmpv, (vec_index_t)lht->NL);
                lmpv->e[lht->NL] = pst_lamp_new(&dir, NULL, 1.0);
                lht->NL++;
              }
          }
      }
    /* Cuber vertices: */
    for (int32_t i0 = 0; i0 < 2; i0++)
      { for (int32_t i1 = 0; i1 < 2; i1++)
          { for (int32_t i2 = 0; i2 < 2; i2++)
              { /* A corner of the cube: */
                r3_t dir = (r3_t){{ 2*i0-1, 2*i1-1, 2*i2-1 }};
                (void)r3_dir(&dir, &dir);
                pst_lamp_vec_expand(lmpv, (vec_index_t)lht->NL);
                lmpv->e[lht->NL] = pst_lamp_new(&dir, NULL, 1.0);
                lht->NL++;
              }
          }
      }
  }

void pst_light_add_uniform_array
  ( pst_light_t *lht, 
    uint32_t NA,
    r3_t *dir0, 
    r3_t *dir1, 
    frgb_t *pwr,
    double crad
  )
  { 
    uint32_t IS = lht->NL; /* Initial number of lights. */
    switch(NA)
      {
      case 1: 
        /* Single ambient lamp. */
        pst_light_add_single(lht, -1.0); 
        break;
      case 2: 
        /* Pair of opposite wall lamps. */
        pst_light_add_pair(lht, 0.0); 
        break;
      case 3: 
        /* Pair of opposite point lamps and an ambient light. */
        pst_light_add_single(lht, -1.0); 
        pst_light_add_pair(lht, 0.0); 
        break;
      case 4: 
        /* Vertices of a tetrahedron. */
        pst_light_add_tetra(lht, FALSE); 
        break;
      case 6: 
        /* Vertices of an octahedron. */
        pst_light_add_octa(lht); 
        break;
      case 8: 
        /* Vertices of a cube. */
        pst_light_add_tetra(lht, FALSE); 
        pst_light_add_tetra(lht, TRUE); 
        break;
      case 12: 
        /* Vertices of an icosahedron. */
        pst_light_add_icosa(lht); 
        break;
      case 14: 
        /* Octahedron and cube. */
        pst_light_add_octa(lht); 
        pst_light_add_tetra(lht, FALSE); 
        pst_light_add_tetra(lht, TRUE); 
        break;
      case 20: 
        /* Dodecahedron. */
        pst_light_add_dodeca(lht); 
        break;
      case 32: 
        /* Dodecahedron and icosahedron. */
        pst_light_add_icosa(lht); 
        pst_light_add_dodeca(lht); 
        break;
      default:
        /* Not implemented yet */
        demand(FALSE, "cannot generate mode with that many lamps");
      }
    assert(lht->NL == IS + NA);

    /* Adjust lamp parameters as appropriate: */
    pst_lamp_vec_t *lmpv = &(lht->lmpv);
    if ((dir0 != NULL) && (r3_L_inf_norm(dir0) != 0.0))
      { /* Adjust directions of all lamps in array. */
        /* Get directions of first two lamps: */
        r3_t *sdir0 = (IS < lht->NL ? &(lmpv->e[IS]->dir) : NULL); 
        r3_t *sdir1 = (IS+1 < lht->NL ? &(lmpv->e[IS+1]->dir) : NULL);
        /* Compute the rotation matrices: */
        r3x3_t S, D;
        r3x3_u_v_to_x_y_rotation(sdir0, sdir1, &S);
        r3x3_u_v_to_x_y_rotation(dir0, dir1, &D);
        /* Apply it to all lamps: */
        for (uint32_t i = IS; i < lht->NL; i++)
          { pst_lamp_t *src = lmpv->e[i];
            r3_t tdir;
            /* Map {src.dir} by {S}, then by the inverse of {D}: */
            r3x3_map_row(&(src->dir), &S, &(tdir));
            r3x3_map_col(&D, &(tdir), &(src->dir));
          }
      }
    
    if (pwr != NULL)
      { /* Set power fields of new lamps to {pwr}. */
        for (uint32_t i = IS; i < lht->NL; i++)
          { pst_lamp_t *src = lmpv->e[i];
            src->pwr = (*pwr);
          }
      }
    
    if (isfinite(crad))
      { /* Set all cosine-radius fields to {crad}. */
        for (uint32_t i = IS; i < lht->NL; i++)
          { pst_lamp_t *src = lmpv->e[i];
            src->crad = crad;
          }
      }
    else if (NA >= 4)
      { /* Set the radii according to neighbor's distances: */
        pst_light_adjust_lamp_radii(&(lht->lmpv.e[IS]), NA);
      }
  }

void pst_light_adjust_lamp_radii(pst_lamp_t *lmpv[], uint32_t NL)
  { bool_t debug = TRUE;
    for (uint32_t i = 0; i < NL; i++)
      { pst_lamp_t *srci = lmpv[i];
        /* Find mean distance to nearest neighbors: */
        if (debug) { r3_gen_print(stderr, &(srci->dir), "%+7.4f", "  srci = ( ", " ", " )\n"); }
        double twsum = 0, wsum = 0;
        for (uint32_t j = 0; j < NL; j++)
          { if (j == i) { /* Skip self: */ continue; }
            pst_lamp_t *srcj = lmpv[j];
            /* Compute angle {ang} between the two: */
            double cang = r3_dot(&(srci->dir), &(srcj->dir));
            if (debug) { r3_gen_print(stderr, &(srcj->dir), "%+7.4f", "    srcj = ( ", " ", " )\n"); }
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
  
void pst_light_set_default_power(pst_lamp_t *lmpv[], uint32_t NL)
  { float defpwr = (float)(1.0/((double)NL));
    for (int32_t i = 0; i < NL; i++)
      { pst_lamp_t *lmp = lmpv[i];
        for (int32_t c = 0; c < 3; c++)
          { if (isnan(lmp->pwr.c[c])) { lmp->pwr.c[c] = defpwr; } }
      }
  }
         

frgb_t pst_light_shading(pst_light_t *lht, r3_t *nrm)
  { uint32_t NL = lht->NL;
    double val[3] = { 0.0, 0.0, 0.0 };
    for (int32_t i = 0; i < NL; i++)
      { pst_lamp_t *src = lht->lmpv.e[i];
        double gf = pst_lamp_geom_factor(nrm, &(src->dir), src->crad);
        for (int32_t c = 0; c < 3; c++)
          { val[c] += gf*src->pwr.c[c]; }
      }
    return (frgb_t){{ (float)(val[0]), (float)(val[1]), (float)(val[2]) }}; 
  }

pst_light_t *pst_light_spec_parse(argparser_t *pp, bool_t next)
  { 
    pst_light_t *lht = pst_light_new(0);
    pst_lamp_vec_t *lmpv = &(lht->lmpv);
    while(TRUE)
      { uint32_t N;
        pst_lamp_t *src = pst_lamp_spec_parse(pp, &N);
        if (src == NULL) { break; }
        if (N == 1)
          { /* Single lamp: */
            pst_lamp_vec_expand(lmpv, (vec_index_t)lht->NL);
            lmpv->e[lht->NL] = src;
            lht->NL++;
          }
        else 
          { /* Lamp array: */
            pst_light_add_uniform_array
              ( lht, N, &(src->dir), NULL, &(src->pwr), src->crad );
            free(src);
          }
      }
    pst_lamp_vec_trim(lmpv, lht->NL);
    pst_light_set_default_power(lmpv->e, lmpv->ne);
    return lht;
  }

void pst_light_spec_write(FILE *wr, pst_light_t *lht)
  { pst_lamp_vec_t *lmpv = &(lht->lmpv);
    uint32_t NL = lmpv->ne;
    for(uint32_t i = 0; i < NL; i++)
      { pst_lamp_spec_write(wr, lmpv->e[i]); }
  }
