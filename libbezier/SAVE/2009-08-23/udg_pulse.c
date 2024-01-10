/* See udg_pulse.h */
/* Last edited on 2009-08-22 19:11:03 by stolfi */

#include <math.h>

#include <bz_basic.h>
#include <affirm.h>

#include <udg_pulse.h>
#include <udg_pulse_B.h>
#include <udg_pulse_H.h>
#include <udg_pulse_N.h>

/* IMPLS */

char udg_pulse_kind_to_char(udg_pulse_kind_t pkind)
  { switch(pkind)
      { case udg_PK_B: 
          return 'B';
        case udg_PK_H: 
          return 'H';
        case udg_PK_N: 
          return 'N';
        default:
          affirm(FALSE, "unknown pulse kind");
          return '?'; /* To pacify the compiler. */
      }
  }  

udg_pulse_family_t udg_pulse_family(udg_pulse_kind_t pkind, udg_cont_t c, udg_degree_t g)
  { unsigned int nmp;
    switch(pkind)
      { case udg_PK_B: 
          nmp = udg_pulse_B_num_mothers(c, g); break;
        case udg_PK_H: 
          nmp = udg_pulse_H_num_mothers(c, g); break;
        case udg_PK_N: 
          nmp = udg_pulse_N_num_mothers(c, g); break;
        default:
          affirm(FALSE, "unknown pulse kind");
          nmp = 0;
      }
    return (udg_pulse_family_t){ .pkind = pkind, .c = c, .g = g, .nmp = nmp }; 
  }

udg_pulse_t udg_pulse
  ( udg_pulse_family_t *fam, 
    udg_pulse_mother_index_t pix, 
    udg_grid_size_t gsz, 
    udg_grid_pos_t pos
  )
  { udg_grid_size_t msz = udg_pulse_mother_supp_count(fam, pix);
    return (udg_pulse_t){ .pix = pix, .gsz = gsz, .pos = pos, .msz = msz };
  }

void udg_pulse_eval
  ( udg_pulse_family_t *fam,
    udg_pulse_t *p,
    interval_t *R,
    double x,
    udg_degree_t ord,
    double *f
  )
  { /* Scale argument from {R} to {[0 _ p.gsz)} and shift it from {p.pos} to 0: */
    double scale = (double)p->gsz, y = x;
    if (R != NULL) 
      { double xlo = (R == NULL ? 0 : LO(*R));
        double xhi = (R == NULL ? 1 : HI(*R));
        scale /= xhi - xlo;  y -= xlo;
      }
    y = y*scale - p->pos;
    
    /* Evaluate the pulse and derivatives w.r.t. {y}: */
    udg_pulse_eval_cell_relative(fam, p->pix, p->gsz, y, ord, f);

    if ((ord > 0) && (scale != 1))
      { /* Adjust derivatives to compensate for scaling: */
        double fac = scale;
        int k;
        for (k = 1; k <= ord; k++) { f[k] *= fac; fac *= scale; }
      }
  }

void udg_pulse_eval_cell_relative
  ( udg_pulse_family_t *fam,
    udg_pulse_mother_index_t pix, 
    udg_grid_size_t gsz,
    double x,
    udg_cont_t ord,
    double *f
  )
  { affirm(fam->g >= 0, "negative degree");
    affirm(fam->g <= udg_pulse_MAX_DEGREE, "excessive degree");
    affirm(ord <= udg_pulse_MAX_DIFF_ORD, "diff order too high");
    
    /* Reduce {x} to {[0 _ gsz)} modulo {gsz}: */
    double z = x - floor(x/gsz)*gsz;
    while (z >= gsz) { z -= gsz; }
    while (z < 0) { z += gsz; }
    
    /* Clear out the result vector: */
    int k;
    for (k = 0; k <= ord; k++) { f[k] = 0; }
    
    /* Discover the high endpoint {hi} of mother pulse's support: */
    udg_grid_pos_t hi = udg_pulse_mother_supp_count(fam, pix);

    /* Add mother's values for all {z} in {(0 _ hi)} congruent to {x} mod {gsz}: */
    if (z == 0) { z += gsz; } 
    double w[udg_pulse_MAX_DEGREE + 1];
    while (z < hi) 
      { udg_pulse_mother_eval(fam, pix, z, ord, w);
        for (k = 0; k <= ord; k++) { f[k] += w[k]; }
        z += gsz;
      }
  }

udg_grid_size_t udg_pulse_max_supp_count ( udg_pulse_family_t *fam )
  {
    int pix, nsmax = 0;
    for (pix = 0; pix < fam->nmp; pix++)
      { int nsm = udg_pulse_mother_supp_count(fam, pix);
        if (nsm > nsmax) { nsmax = nsm; }
      }
    return nsmax;
  }

void udg_pulse_to_bezier
  ( udg_pulse_family_t *fam,
    udg_pulse_t *p, 
    udg_grid_pos_t ix,
    double *bz
  )
  { affirm(fam->g >= 0, "negative degree");
    affirm(fam->g <= udg_pulse_MAX_DEGREE, "excessive degree");
    
    /* Shift pulse from {pos} to 0: */
    ix -= p->pos;
    
    /* Compute {z} = {ix} reduced to {0..p->gsz-1} modulo {p->gsz}: */
    int z = ix % p->gsz;
    while (z < 0) { z += p->gsz; }
    while (z >= p->gsz) { z -= p->gsz; };
    
    /* Clear out the Bézier coeff vector: */
    int k;
    for (k = 0; k <= fam->g; k++) { bz[k] = 0; }

    /* Add mother's coeffs for all {z} in {0..p->msz-1} congruent to {ix} mod {gsz}: */
    double wbz[udg_pulse_MAX_DEGREE + 1];
    while (z < p->msz) 
      { udg_pulse_mother_to_bezier(fam, p->pix, z, wbz);
        for (k = 0; k <= fam->g; k++) { bz[k] += wbz[k]; }
        z += p->gsz;
      }
  }

/* STANDARD MOTHER PULSES */

udg_grid_size_t udg_pulse_mother_supp_count
  ( udg_pulse_family_t *fam, 
    udg_pulse_mother_index_t pix
  )
  { switch(fam->pkind)
      { case udg_PK_B: 
          return udg_pulse_B_mother_supp_count(fam->c, fam->g, pix);
        case udg_PK_H: 
          return udg_pulse_H_mother_supp_count(fam->c, fam->g, pix);
        case udg_PK_N: 
          return udg_pulse_N_mother_supp_count(fam->c, fam->g, pix);
        default:
          affirm(FALSE, "unknown pulse kind");
          return 0.0; /* To keep the compiler happy. */
      }
  }

void udg_pulse_mother_eval
  ( udg_pulse_family_t *fam, 
    udg_pulse_mother_index_t pix, 
    double z,
    udg_cont_t ord,
    double *w
  )
  { affirm(fam->g >= 0, "negative degree");
    affirm(fam->g <= udg_pulse_MAX_DEGREE, "excessive degree");
    affirm(ord <= udg_pulse_MAX_DIFF_ORD, "diff order too high");
    /* Get Bézier coeffs for interval containing {z}: */
    double bz[udg_pulse_MAX_DEGREE + 1];
    int k = (int) floor(z); 
    udg_pulse_mother_to_bezier(fam, pix, k, bz);
    /* Apply DeCasteljau's algorithm: */
    bz_eval(fam->g, bz, 0.0, 1.0, z - k, ord, w);
  }

void udg_pulse_mother_to_bezier
  ( udg_pulse_family_t *fam, 
    udg_pulse_mother_index_t pix,
    udg_grid_pos_t ix,
    double *bz
  )
  { affirm(fam->g >= 0, "negative degree");
    affirm(fam->g <= udg_pulse_MAX_DEGREE, "excessive degree");
    switch(fam->pkind)
      { case udg_PK_B: 
          udg_pulse_B_mother_to_bezier(fam->c, fam->g, pix, ix, bz);
          return;
        case udg_PK_H: 
          udg_pulse_H_mother_to_bezier(fam->c, fam->g, pix, ix, bz);
          return;
        case udg_PK_N: 
          udg_pulse_N_mother_to_bezier(fam->c, fam->g, pix, ix, bz);
          return;
        default:
          affirm(FALSE, "unknown pulse kind");
      }
  }

void udg_pulse_family_print(FILE *wr, udg_pulse_family_t *fam)
  { /* Pulse kind: */
    fprintf(wr, "%c", udg_pulse_kind_to_char(fam->pkind));
    /* Continuity and degree: */
    fprintf(wr, "_%d^%d", fam->c, fam->g);
  }

void udg_pulse_print
  ( FILE *wr, 
    udg_pulse_family_t *fam, 
    udg_pulse_t *p
  )
  { udg_pulse_family_print(wr, fam); /* Pulse family: */
    /* Mother pulse index: */
    if (fam->nmp != 1)
      { fprintf(wr, "[");
        fprintf(wr, "%d", p->pix);
        if (p->pix < 0) { fprintf(wr, "!?"); }
        if (p->pix >= fam->nmp) { fprintf(wr, "!?(%d)", fam->nmp); }
        fprintf(wr, "]");
      }
    /* Support endpoints: */
    fprintf(wr, ":%lld:%lld/%lld", p->pos, p->msz, p->gsz);
  }
