/* See pst_imgsys.h  */

/* Created on 2005-10-01 by Jorge Stolfi, unicamp, <stolfi@dcc.unicamp.br> */
/* Based on the work of Rafael Saracchini, U.F.Fluminense. */
/* Last edited on 2025-03-03 14:25:05 by stolfi */
/* See the copyright and authorship notice at the end of this file.  */

#include <stdint.h>
#include <stdio.h>
#include <assert.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include <bool.h>
#include <affirm.h>
#include <jsfile.h>
#include <float_image.h>
#include <float_image_mscale.h>
#include <filefmt.h>

#include <pst_imgsys.h>

#define MAX_COEFFS pst_imgsys_MAX_COEFFS

#define FLUFF (1.0e-140)
  /* Practically zero. */

pst_imgsys_t *pst_imgsys_new_grid(int32_t NX, int32_t NY)
  {
    uint32_t N = (uint32_t)(NX*NY);
    pst_imgsys_t *S = talloc(1, pst_imgsys_t);
    S->N = N;
    S->eq = talloc(N, pst_imgsys_equation_t);
    S->NX = NX;
    S->NY = NY;
    S->col = talloc(N, int32_t);
    S->row = talloc(N, int32_t);
    for (int32_t uid = 0; uid < N; uid++)
      { S->eq[uid].nt = 1;
        S->eq[uid].wtot = 0;
        S->eq[uid].rhs = 0;
        S->eq[uid].cf[0] = 0;
        S->eq[uid].uid[0] = (uint32_t)uid;
        int32_t xy = uid;
        S->col[uid] = xy % NX;
        S->row[uid] = xy / NX;
      }
    return S;
  }

pst_imgsys_t *pst_imgsys_from_eqs
  ( uint32_t N,
    pst_imgsys_equation_t *eq,
    int32_t *col,
    int32_t *row
  )
  {
    pst_imgsys_t *S = talloc(1, pst_imgsys_t);
    S->N = N;
    S->eq = eq;
    S->col = col;
    S->row = row;
    pst_imgsys_check_valid(S, N, -1, -1);
    return S;
  }

void pst_imgsys_check_valid(pst_imgsys_t *S, uint32_t N, int32_t NX, int32_t NY)
  { 
    demand(S->N == N, "mismatch {N, S.N}");
    /* Check consitency of {S.uid,S.col,S.row} */
    for (int32_t uid = 0; uid < N; uid++)
      { int32_t x = S->col[uid];
        int32_t y = S->row[uid];
        demand((x == -1) == (y == -1), "inconsistent {S.row[k],S.col[k]}");
        if (x != -1)
          { demand((x >= 0) && ((NX == -1 ) || (x < NX)), "invalid {S.col[k]}");
            demand((y >= 0) && ((NY == -1 ) || (y < NY)), "invalid {S.row[k]}");
          }
      }
    /* Check consistency of {S.eq}: */
    for (int32_t eid = 0; eid < N; eid++)
      { pst_imgsys_equation_t *eqi = &(S->eq[eid]); 
        demand((eqi->nt >= 1) && (eqi->nt <= MAX_COEFFS), "invalid term count in equation");
        int32_t uid = eid;
        demand(eqi->uid[0] == uid, "equation's first term is not its main variable");
      }
  }

void pst_imgsys_free(pst_imgsys_t *S)
  {
    free(S->eq);
    free(S->col);
    free(S->row);
    free(S);
  }

bool_t pst_imgsys_equation_is_null(uint32_t uid, pst_imgsys_equation_t *eq, uint32_t N)
  { assert((uid >= 0) &&(uid < N));
    assert(eq->uid[0] == uid);
    assert((eq->nt >= 1) && (eq->nt <= MAX_COEFFS));
    if (fabs(eq->wtot) < FLUFF)
      { /* The equation must be indeterminate: */
        assert(eq->nt == 1);
        assert(eq->rhs == 0.0);
        assert(eq->cf[0] == 0);
        return TRUE;
      }
    else
      { /* The equation must be determinate: */
        assert(eq->wtot > 0);
        assert(eq->nt >= 2);
        assert(fabs(eq->cf[0]) >= FLUFF);
        double sum = 0;
        for (int32_t j = 1; j < eq->nt; j++) { sum += fabs(eq->cf[j]); }
        assert(sum >= FLUFF);
        return FALSE;
      }
  }
    
void pst_imgsys_write(FILE *wr, pst_imgsys_t *S, char *fmt)
  {
    filefmt_write_header(wr, pst_imgsys_FILE_TYPE, pst_imgsys_FILE_VERSION);
    
    /* Provide a default format: */
    if ((fmt == NULL) || (strlen(fmt) == 0)) { fmt = "%+24.16e"; }
    
    fprintf(wr, "N = %d\n", S->N);
    fprintf(wr, "NX = %d\n", S->NX);
    fprintf(wr, "NY = %d\n", S->NY);

    for(int32_t k = 0; k < S->N; k++)
      {
        /* Index of equation's main variable: */
        int32_t uidk = k;
        
        /* Get indices of pixel corresponding to the variable {Z[k]}: */
        int32_t xk = S->col[uidk]; int32_t yk = S->row[uidk];
         
        /* Print indices of main pixel: */
        fprintf(wr, "  %8d = %4d %4d", k, xk, yk);
        assert((xk == -1) == (yk == -1));
        if (xk != -1)
          { assert((xk >= 0) && (xk < S->NX));
            assert((yk >= 0) && (yk < S->NY));
          }
            
        /* Get equation {k}: */
        pst_imgsys_equation_t *eqk = &(S->eq[k]);
        
        /* Print the total weight: */
        fprintf(wr, " w = %8.6f ", eqk->wtot);

        /* Print the right-hand side and weight: */
        fprintf(wr, " rhs = "); fprintf(wr, fmt, eqk->rhs);
        
        /* Print the number of terms: */
        fprintf(wr, " nt = %2d", eqk->nt);
        
        assert((eqk->nt >= 1) && (eqk->nt <= MAX_COEFFS));
        assert(eqk->uid[0] == k);
        
        /* Print the equation coefficients: */
        for(int32_t j = 0; j < eqk->nt; j++)
          { /* Get hold of a variable {Z[uidj]} in equation {k}, and its coeff: */
            uint32_t uidj = eqk->uid[j];
            assert((uidj >= 0) && (uidj < S->N));
            double cfj = eqk->cf[j];
            fprintf(wr, " "); fprintf(wr, fmt, cfj); 
            fprintf(wr, "*Z[%d", uidj);
            /* Get indices of the corresponding pixel: */
            int32_t xj = S->col[uidj]; int32_t yj = S->row[uidj];
            assert((xj == -1) == (yj == -1));
            if (xj != -1)
              { assert((xj >= 0) && (xj < S->NX));
                assert((yj >= 0) && (yj < S->NY));
                fprintf(wr, "=%d,%d", xj, yj);
              }
            fprintf(wr, "]");
          }
        fprintf(wr, "\n");
      }
    filefmt_write_footer(wr, pst_imgsys_FILE_TYPE);
    fflush(wr);
  }

float_image_t* pst_imgsys_make_weight_image(pst_imgsys_t *S)
  {
    float_image_t *SW = float_image_new(1, S->NX, S->NY);
    
    int32_t nout = 0; /* Equations without valid pixel indices. */
    for (int32_t k = 0; k < S->N; k++)
      {
        /* Index of equation's main variable: */
        int32_t uidk = k;
        
        /* Get equation {k}: */
        pst_imgsys_equation_t *eqk = &(S->eq[k]);
        
        /* Get indices of pixel corresponding to the variable {Z[k]}: */
        int32_t xk = S->col[uidk]; int32_t yk = S->row[uidk];
         
        /* Print the total weight: */
        if ((xk >= 0) && (xk < S->NX) && (yk >= 0) && (yk < S->NY))
          { float_image_set_sample(SW, 0, xk, yk, (float)(eqk->wtot)); }
        else
          { nout++; }
      }
    if (nout > 0)
      { fprintf(stderr, "{pst_imgsys_make_weight_image}: there were %d equations without pixel indices\n", nout); }

    return SW;
  }

void pst_imgsys_write_named(char *fileName, pst_imgsys_t *S, char *fmt, int32_t indent)
  {
    if (S == NULL) { return; }
    fprintf(stderr, "%*swriting %s ...\n", indent, "", fileName);
    FILE* wr = open_write(fileName, FALSE);
    assert(wr != NULL); 
    pst_imgsys_write(wr, S, fmt);
    if (wr == stdout) { fflush(wr); } else { fclose(wr); }
  }

void pst_imgsys_copy_image_to_sol_vec(pst_imgsys_t *S, float_image_t *Z, double z[], double vdef)
  {
    uint32_t N = S->N;
    int32_t NX = S->NX, NY = S->NY;
    float_image_check_size(Z, 2, NX, NY, "bad height map");
    for (uint32_t k = 0; k < N; k++)
      { int32_t x = S->col[k];
        int32_t y = S->row[k];
        assert((x == -1) == (y == -1));
        if (x != -1)
          { assert((x >= 0) && (x < NX));
            assert((y >= 0) && (y < NY));
            z[k] = float_image_get_sample(Z, 0, x, y);
          }
        else
          { /* Height Value is not associated with any height map pixel: */
            z[k] = vdef;
          }
      }
  }

void pst_imgsys_copy_sol_vec_to_image(pst_imgsys_t *S, double z[], float_image_t *Z, float vdef)
  {
    int32_t NC, NX, NY;
    float_image_get_size(Z, &NC, &NX, &NY);
    float_image_fill_channel(Z, 0, vdef);
    if (NC >= 2) { float_image_fill_channel(Z, 1, 0.0); }
    for (uint32_t uid = 0; uid < S->N; uid++)
      { int32_t x = S->col[uid];
        int32_t y = S->row[uid];
        assert((x == -1) == (y == -1));
        if (x != -1)
          { demand((x >= 0) && (x < NX), "invalid image col in system");
            demand((y >= 0) && (y < NY), "invalid image row in system");
            float_image_set_sample(Z, 0, x, y, (float)z[uid]);
            if (NC >= 2) { float_image_set_sample(Z, 1, x, y, (float)S->eq[uid].wtot); }
          }
      }
  }


/*
**
** Copyright (C) Jorge Stolfi, Unicamp.
**
** Permission to use, copy, modify, and distribute this software and its
** documentation for any purpose and without fee is hereby granted, provided
** that the above copyright notice appear in all copies and that both that
** copyright notice and this permission notice appear in supporting
** documentation.  This software is provided "as is" without express or
** implied warranty. Neither the author nor its employers are liable to
** any damages which may result from its use.
*/
