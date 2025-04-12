/* See {pst_gr.h} */
/* Last edited on 2025-03-11 14:54:56 by stolfi */
/* Created by Rafael F. V. Saracchini */

#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include <math.h>
#include <assert.h>

#include <jsfile.h>
#include <jsprintf.h>
#include <affirm.h>
#include <rn.h>

#include <pst_imgsys.h>
#include <pst_gr.h>

#include <pst_gr_integrate.h>

#define DELETED pst_gr_mark_DELETED
#define MAX_COEFFS pst_imgsys_MAX_COEFFS
#define FLUFF (1.0e-140)

pst_imgsys_t *pst_gr_integrate_build_system
  ( pst_gr_t *gr,
    int32_t NX_Z,
    int32_t NY_Z,
    int32_t kz_from_kv[],
    bool_t verbose
  )
  { 
    uint32_t N_gone = 0; /* Vertex records marked as {DELETED}. */
    uint32_t N_null = 0; /* Equations with zero weight. */
    uint32_t N_hors = 0; /* Equations without pixel indices but nonzero weight. */
    uint32_t N_good = 0; /* Equations with pixel indices and nonzero weight. */
  
    /* Build the table {kz_from_kv} that maps vertex index to equation/height index, or {-1}: */
    uint32_t NZ = 0; /* Number of equations and variables. */
    for (uint32_t kv = 0; kv < gr->NV ; kv++)
      { pst_gr_vertex_data_t *vdk = &(gr->vdata[kv]);
        if (vdk->vmark == DELETED) 
          { kz_from_kv[kv] = -1; N_gone++; }
        else
          { kz_from_kv[kv] = (int32_t)NZ; NZ++; }
      }
    
    /* Now build the equations: */
    
    pst_imgsys_equation_t *eq = talloc(NZ, pst_imgsys_equation_t);
    for (uint32_t kv = 0; kv < gr->NV ; kv++)
      { int32_t kz = kz_from_kv[kv];
        if (kz != -1)
          { assert((kz >= 0) && (kz < NZ));
            pst_imgsys_equation_t *eqk = &(eq[kz]);
            eqk->nt = 0; /* Number of terms in equation. */
            eqk->rhs = 0.0;
            eqk->uid[eqk->nt] = (uint32_t)kz;  
            eqk->cf[eqk->nt] = 0.00; 
            eqk->nt++;
            pst_gr_vertex_data_t *vdk = &(gr->vdata[kv]);
            assert((vdk->x == -1) == (vdk->y == -1));
            pst_gr_arc_t a0 = vdk->aout;
            pst_gr_arc_t aj = a0;
            do
              { double dj = pst_gr_arc_delta(gr, aj);
                double wj = pst_gr_arc_weight(gr, aj);
                uint32_t dstj = pst_gr_arc_org(gr,pst_gr_arc_sym(aj)); /* Index of {aj} dst vertex. */
                assert(dstj < gr->NV);
                int32_t jz = kz_from_kv[dstj];
                assert((jz >= 0) && (jz < NZ));
                assert(isfinite(wj) && (fabs(wj) >= FLUFF));
                assert(isfinite(dj));
                uint32_t j = eqk->nt;
                demand(j < MAX_COEFFS, "too many edges out of a vertex");
                eqk->uid[j] = (uint32_t)jz;
                eqk->cf[j] = -wj;
                eqk->cf[0] += wj;
                eqk->rhs += -wj*dj;
                eqk->wtot += wj;
                eqk->nt++;
                aj = pst_gr_arc_onext(aj);
              } while(aj != a0);
            assert(eqk->nt <= MAX_COEFFS);
            if (fabs(eqk->wtot) < FLUFF)
              { assert(eqk->nt == 1);
                N_null++;
              }
            else 
              { assert(eqk->nt >= 2);
                N_good++;
              }
          }
      }

    /* Build the pixel index tables {col[0..N-1],row[0..N-1]} and the inverse {zid[0..NX*NY-1]}: */
    int32_t *col = talloc(NZ, int32_t);
    int32_t *row = talloc(NZ, int32_t);
    int32_t NXY_Z = NX_Z*NY_Z;
    int32_t *zid_from_xy = talloc(NXY_Z, int32_t);
    for (int32_t xy = 0; xy < NXY_Z; xy++) { zid_from_xy[xy] = -1; }
    for (int32_t kv = 0; kv < gr->NV; kv++)
      { int32_t kz = kz_from_kv[kv]; 
        if (kz_from_kv[kv] != -1)
          { assert((kz >= 0) && (kz < NZ));
            pst_gr_vertex_data_t *vdk = &(gr->vdata[kv]);
            int32_t x = vdk->x;
            int32_t y = vdk->y;
            assert((x == -1) == (y == -1));
            if (x == -1)
              { N_hors++; }
            else
              { demand((x >= 0) && (x < NX_Z), "invalid {x} pixel index in vertex");
                demand((y >= 0) && (y < NY_Z), "invalid {y} pixel index in vertex");
                int32_t xy = x + y*NX_Z;
                zid_from_xy[xy] = kz;
              }
            col[kv] = x;
            row[kv] = y;
          }
      }
    assert(N_gone + N_null + N_good == gr->NV);
    assert(N_null + N_good == NZ);
    if (verbose)
      { fprintf(stderr, "%5d deleted vertices, ignored\n", N_gone);
        fprintf(stderr, "%5d height values with null equations\n", N_null);
        fprintf(stderr, "%5d height values with non-null equations\n", N_good);
        fprintf(stderr, "%5d height values are not height map pixels\n", N_hors);
      }
   
    /* Now package the equations as a system: */
    pst_imgsys_t *S = pst_imgsys_from_eqs(NZ, eq, col, row);  
    return S;
  }
