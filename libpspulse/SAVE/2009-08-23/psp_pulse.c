/* See psp_pulse.h */
/* Last edited on 2009-08-23 19:03:20 by stolfi */

#include <math.h>

#include <bz_basic.h>
#include <affirm.h>

#include <psp_pulse_mother.h>

#include <psp_pulse.h>

/* IMPLS */

char psp_pulse_kind_to_char(psp_pulse_kind_t pkind)
  { switch(pkind)
      { case psp_PK_B: 
          return 'B';
        case psp_PK_H: 
          return 'H';
        case psp_PK_N: 
          return 'N';
        default:
          affirm(FALSE, "unknown pulse kind");
          return '?'; /* To pacify the compiler. */
      }
  }  

psp_pulse_family_t psp_pulse_family(psp_pulse_kind_t pkind, psp_cont_t c, psp_degree_t g)
  { psp_pulse_mother_count_t nmp = psp_pulse_mother_count(pkind,c,g);
    return (psp_pulse_family_t){ .pkind = pkind, .c = c, .g = g, .nmp = nmp }; 
  }

psp_pulse_t psp_pulse
  ( psp_pulse_family_t *fam, 
    psp_pulse_mother_index_t pix, 
    psp_grid_size_t gsz, 
    psp_grid_pos_t pos
  )
  { psp_grid_size_t msz = psp_pulse_mother_supp_count(fam, pix);
    return (psp_pulse_t){ .pix = pix, .gsz = gsz, .pos = pos, .msz = msz };
  }

void psp_pulse_eval
  ( psp_pulse_family_t *fam,
    psp_pulse_t *p,
    interval_t *R,
    double x,
    psp_degree_t ord,
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
    psp_pulse_eval_cell_relative(fam, p->pix, p->gsz, y, ord, f);

    if ((ord > 0) && (scale != 1))
      { /* Adjust derivatives to compensate for scaling: */
        double fac = scale;
        int k;
        for (k = 1; k <= ord; k++) { f[k] *= fac; fac *= scale; }
      }
  }

void psp_pulse_eval_cell_relative
  ( psp_pulse_family_t *fam,
    psp_pulse_mother_index_t pix, 
    psp_grid_size_t gsz,
    double x,
    psp_cont_t ord,
    double *f
  )
  { affirm(fam->g >= 0, "negative degree");
    affirm(fam->g <= psp_pulse_MAX_DEGREE, "excessive degree");
    affirm(ord <= psp_pulse_MAX_DIFF_ORD, "diff order too high");
    
    /* Reduce {x} to {[0 _ gsz)} modulo {gsz}: */
    double z = x - floor(x/gsz)*gsz;
    while (z >= gsz) { z -= gsz; }
    while (z < 0) { z += gsz; }
    
    /* Clear out the result vector: */
    int k;
    for (k = 0; k <= ord; k++) { f[k] = 0; }
    
    /* Discover the high endpoint {hi} of mother pulse's support: */
    psp_grid_pos_t hi = psp_pulse_mother_supp_count(fam, pix);

    /* Add mother's values for all {z} in {(0 _ hi)} congruent to {x} mod {gsz}: */
    if (z == 0) { z += gsz; } 
    double w[psp_pulse_MAX_DEGREE + 1];
    while (z < hi) 
      { psp_pulse_mother_eval(fam, pix, z, ord, w);
        for (k = 0; k <= ord; k++) { f[k] += w[k]; }
        z += gsz;
      }
  }

psp_grid_size_t psp_pulse_max_supp_count ( psp_pulse_family_t *fam )
  {
    int pix, nsmax = 0;
    for (pix = 0; pix < fam->nmp; pix++)
      { int nsm = psp_pulse_mother_supp_count(fam, pix);
        if (nsm > nsmax) { nsmax = nsm; }
      }
    return nsmax;
  }

void psp_pulse_to_bezier
  ( psp_pulse_family_t *fam,
    psp_pulse_t *p, 
    psp_grid_pos_t ix,
    double *bz
  )
  { affirm(fam->g >= 0, "negative degree");
    affirm(fam->g <= psp_pulse_MAX_DEGREE, "excessive degree");
    
    /* Shift pulse from {pos} to 0: */
    ix -= p->pos;
    
    /* Compute {z} = {ix} reduced to {0..p->gsz-1} modulo {p->gsz}: */
    int z = ix % p->gsz;
    while (z < 0) { z += p->gsz; }
    while (z >= p->gsz) { z -= p->gsz; };
    
    /* Clear out the B�zier coeff vector: */
    int k;
    for (k = 0; k <= fam->g; k++) { bz[k] = 0; }

    /* Add mother's coeffs for all {z} in {0..p->msz-1} congruent to {ix} mod {gsz}: */
    double wbz[psp_pulse_MAX_DEGREE + 1];
    while (z < p->msz) 
      { psp_pulse_mother_to_bezier(fam, p->pix, z, wbz);
        for (k = 0; k <= fam->g; k++) { bz[k] += wbz[k]; }
        z += p->gsz;
      }
  }

/* STANDARD MOTHER PULSES */

void psp_pulse_family_print(FILE *wr, psp_pulse_family_t *fam)
  { /* Pulse kind: */
    fprintf(wr, "%c", psp_pulse_kind_to_char(fam->pkind));
    /* Continuity and degree: */
    fprintf(wr, "_%d^%d", fam->c, fam->g);
  }

void psp_pulse_print
  ( FILE *wr, 
    psp_pulse_family_t *fam, 
    psp_pulse_t *p
  )
  { psp_pulse_family_print(wr, fam); /* Pulse family: */
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
