/* See dg_tent.h */
/* Last edited on 2009-05-18 14:20:54 by stolfi */

#include <math.h>
#include <assert.h>

#include <affirm.h>
#include <vec.h>

#include <bz_basic.h>
/* #include <bz_patch.h> */

#include <mdg_grid.h>
#include <udg_pulse.h>
#include <dg_tent.h>

/* IMPLEMENTATIONS */

void dg_tent_eval
  ( mdg_dim_t d, 
    udg_pulse_family_t fam[], 
    udg_pulse_t p[],
    interval_t R[],
    double x[],
    dg_degree_t ord[],
    double f[]
  )
  { double pf[udg_pulse_MAX_DEGREE + 1];
    f[0] = 1.0;
    int nf = 1;
    int i;
    for (i = 0; i < d; i++) 
      { interval_t *Ri = (R == NULL ? NULL : &(R[i]));
        udg_pulse_eval(&(fam[i]), &(p[i]), Ri, x[i], ord[i], pf);
        int j;
        for (j = 0; j < nf; j++)
          { int k; 
            for (k = 0; k <= ord[i]; k++)
              { f[k*nf + j] = f[j]*pf[k]; }
          }
        nf *= (ord[i] + 1);
      }
  }

void dg_tent_eval_total
  ( mdg_dim_t d, 
    udg_pulse_family_t fam[], 
    udg_pulse_t p[],
    interval_t R[],
    double x[],
    dg_degree_t ord,
    double f[]
  )
  { double pf[udg_pulse_MAX_DEGREE + 1];
    f[0] = 1.0;
    int i;
    for (i = 0; i < d; i++) 
      { interval_t *Ri = (R == NULL ? NULL : &(R[i]));
        udg_pulse_eval(&(fam[i]), &(p[i]), Ri, x[i], ord, pf);
        affirm(FALSE, "not implemented");
      }
  }

dg_tent_vec_t dg_tents_from_cells
  ( mdg_dim_t d, 
    dg_tent_mother_index_t tix,
    mdg_cell_index_vec_t cell
  )
  { dg_tent_vec_t tv = dg_tent_vec_new(cell.ne);
    int nt = 0;
    int i; 
    for (i = 0; i < cell.ne; i++)
      { mdg_cell_index_t ck = cell.e[i];
        tv.e[nt] = (dg_tent_t){tix,ck};
        nt++;
      }
    assert(nt == tv.ne);
    return tv;
  }

dg_tent_mother_index_t dg_tent_mother_index_max(mdg_dim_t d, udg_pulse_family_t fam[])
  { dg_tent_mother_index_t tix = 0;
    int i;
    for (i = d-1; i >= 0; i--) { tix = tix * fam[i].nmp + (fam[i].nmp - 1); }
    return tix; 
  }

dg_tent_mother_index_t dg_tent_mother_index_pack
  ( mdg_dim_t d,
    udg_pulse_family_t fam[],
    udg_pulse_mother_index_t pix[]
  )
  { dg_tent_mother_index_t tix = 0;
    int i;
    for (i = d-1; i >= 0; i--) { tix = tix * fam[i].nmp + pix[i]; }
    return tix;
  }

void dg_tent_mother_index_unpack
  ( mdg_dim_t d, 
    udg_pulse_family_t fam[],
    dg_tent_mother_index_t tix,
    udg_pulse_mother_index_t pix[]
  )
  { int i;
    for (i = 0; i < d; i++) { pix[i] = tix % fam[i].nmp; tix /= fam[i].nmp; }
  }

/* PRINTOUT: */

void dg_tent_print
  ( FILE *wr, 
    mdg_dim_t d,
    udg_pulse_family_t fam[], 
    udg_pulse_t p[]
  )
  { int i;
    /* Print pulse indices: */
    for (i = 0; i < d; i++) 
      { fprintf(wr, "[");
        udg_pulse_print(wr, &(fam[i]), &(p[i]));
        fprintf(wr, "]");
      }
  }

void dg_tent_packed_print
  ( FILE *wr, 
    mdg_dim_t d,
    udg_pulse_family_t fam[], 
    dg_tent_t *t
  )
  { /* Unpack puse indices from tent: */
    udg_pulse_t p[d];
    dg_tent_unpack(d, fam, t, p);
    dg_tent_print(wr, d, fam, p);
  }

/* MANIPULATION OF TENT VECTORS */

vec_typeimpl(dg_tent_vec_t,dg_tent_vec,dg_tent_t);

