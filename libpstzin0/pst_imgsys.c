/* See pst_imgsys.h  */

/* Created on 2005-10-01 by Jorge Stolfi, unicamp, <stolfi@dcc.unicamp.br> */
/* Based on the work of Rafael Saracchini, U.F.Fluminense. */
/* Last edited on 2025-01-18 12:36:19 by stolfi */
/* See the copyright and authorship notice at the end of this file.  */

#include <stdio.h>
#include <assert.h>
#include <stdint.h>
#include <stdlib.h>
#include <math.h>

#include <bool.h>
#include <affirm.h>
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
  
void pst_imgsys_remove_holes(pst_imgsys_t *S, int32_t *nuid_from_ouid)
  { bool_t debug = TRUE;
  
    uint32_t N_in = S->N; 
      
    /* Provide a local table if the client did not give one: */
    int32_t *nfo_old = nuid_from_ouid;
    if (nuid_from_ouid == NULL) { nuid_from_ouid = talloc(N_in, int32_t); }
    for (int32_t k = 0; k < N_in; k++) { nuid_from_ouid[k] = -1; }
    
    /* Get number {N-ot} of non-null equations and set {nuid_from_ouid}: */
    uint32_t N_ot = 0; /* Non-null equations. */
    uint32_t N_ex = 0; /* Number of equations excluded. */
    for (uint32_t k = 0; k < N_in; k++)
      { uint32_t uidk_old = k;
        pst_imgsys_equation_t *eqk_old = &(S->eq[k]);
        if (pst_imgsys_equation_is_null(uidk_old, eqk_old, N_in))
          { /* Will exclude it: */
            if (debug) { fprintf(stderr, "    excluding null height/equation %d\n", uidk_old); }
            N_ex++;
          }
        else
          { /* Will keep it: */
            uint32_t uidk_new = N_ot;
            nuid_from_ouid[uidk_old] = (int32_t)uidk_new;
            N_ot++; 
          }
      }          
    assert(N_ot + N_ex == N_in); /* Paranoia. */
    
    if (N_ex > 0)
      { /* Now create the equation, colum, row tables for a new system and copy the non-null equations into them: */
        pst_imgsys_equation_t *eq_ot = talloc(N_ot, pst_imgsys_equation_t);
        int32_t *col_ot = talloc(N_ot, int32_t);
        int32_t *row_ot = talloc(N_ot, int32_t);
        uint32_t k_new = 0;
        for (uint32_t k = 0; k < N_in; k++)
          { uint32_t uidk_old = k;
            pst_imgsys_equation_t *eqk = &(S->eq[k]);
            if (! pst_imgsys_equation_is_null(uidk_old, eqk, N_in))
              { /* Will keep it: */
                uint32_t uidk_new = k_new;
                assert(nuid_from_ouid[uidk_old] == uidk_new);
                eq_ot[k_new] = (*eqk);
                /* Remap terms; for good measure, discard null ones: */
                uint32_t nt = eqk->nt;
                assert((nt >= 2) && (nt <= MAX_COEFFS));
                uint32_t nt_new = 0;
                for (int32_t j = 0; j < nt; j++)
                  { if ((j == 0) || (fabs(eqk->cf[j]) >= FLUFF))
                      { uint32_t uidj_old = eqk->uid[j];
                        int32_t j_new = (int32_t)nt_new;
                        int32_t uidj_new = nuid_from_ouid[uidj_old];
                        assert(uidj_new != -1);
                        assert((uidj_new >= 0) && (uidj_new < N_ot));
                        eqk->uid[j_new] = (uint32_t)uidj_new;
                        eqk->cf[j_new] = eqk->cf[j];
                        nt_new++;
                      }
                  }
                assert(nt_new >= 2);
                eqk->nt = nt_new;
                assert(eqk->uid[0] == uidk_new);
                /* Remap {col,row,uid}: */
                int32_t x = S->col[uidk_old];
                int32_t y = S->row[uidk_old];
                col_ot[uidk_new] = x;
                row_ot[uidk_new] = y;
                /* One more: */
                k_new++;
              }
          }
        assert(k_new == N_ot);
    
        /* Replace the tables: */
        S->N = N_ot;
        free(S->eq); S->eq = eq_ot;
        free(S->col); S->col = col_ot;
        free(S->row); S->row = row_ot;
        pst_imgsys_check_valid(S, N_ot, -1, -1);
      }
    
    if (nuid_from_ouid != nfo_old) { free(nuid_from_ouid); }
    return;
  }
    
void pst_imgsys_write(FILE *wr, pst_imgsys_t *S)
  {
    filefmt_write_header(wr, pst_imgsys_FILE_TYPE, pst_imgsys_FILE_VERSION);

    for(int32_t k = 0; k < S->N; k++)
      {
        /* Index of equation's main variable: */
        int32_t uidk = k;
        
        /* Compute indices of pixel corresponding to the variable {Z[k]}: */
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

        /* Print the number of terms: */
        fprintf(wr, " nt = %2d", eqk->nt);
        
        /* Print the right-hand side and weight: */
        fprintf(wr, " eq = %+12.7f", eqk->rhs);
        
        assert((eqk->nt >= 1) && (eqk->nt <= MAX_COEFFS));
        assert(eqk->uid[0] == k);
        
        /* Print the equation coefficients: */
        for(int32_t j = 0; j < eqk->nt; j++)
          { /* Get hold of a variable {Z[uidj]} in equation {k}, and its coeff: */
            uint32_t uidj = eqk->uid[j];
            assert((uidj >= 0) && (uidj < S->N));
            double cfj = eqk->cf[j];
            fprintf(wr, " %+12.7f*Z[%d", cfj, uidj);
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

void pst_imgsys_write_report(pst_imgsys_t *S, char *filePrefix, int32_t level, char *tag)
  {
    if (S == NULL) { return; }
    char *fileName = float_image_mscale_file_name(filePrefix, level, -1, tag, "sys");
    int32_t indent = (level < -1 ? 0 : 2*level+2);
    fprintf(stderr, "%*swriting %s ...", indent, "", fileName);
    FILE* wr = fopen(fileName, "wt");
    assert(wr != NULL); 
    pst_imgsys_write(wr, S);
    if (wr == stdout) { fflush(wr); } else { fclose(wr); }
    fprintf(stderr, "\n");
    free(fileName);
  }

void pst_imgsys_copy_image_to_sol_vec(pst_imgsys_t *S, float_image_t *Z, double h[], double vdef)
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
            h[k] = float_image_get_sample(Z, 0, x, y);
          }
        else
          { /* Height Value is not associated with any height map pixel: */
            h[k] = vdef;
            ??? say that channel 1 is ignored
          }
      }
  }

void pst_imgsys_copy_sol_vec_to_image(pst_imgsys_t *S, double h[], float_image_t *Z, float vdef)
  {
    int32_t NC, NX, NY;
    float_image_get_size(Z, &NC, &NX, &NY);
    demand(NC == 1, "image should have only one channel");
    float_image_fill_channel(Z, 0, vdef);
    for (uint32_t uid = 0; uid < S->N; uid++)
      { int32_t x = S->col[uid];
        int32_t y = S->row[uid];
        assert((x == -1) == (y == -1));
        if (x != -1)
          { demand((x >= 0) && (x < NX), "invalid image col in system");
            demand((y >= 0) && (y < NY), "invalid image row in system");
            float_image_set_sample(Z, 0, x, y, (float)h[uid]);
          }
      }
  }

void pst_imgsys_extract_system_eq_tot_weight_image(pst_imgsys_t *S, float_image_t *U, float vdef)
  { int32_t NC, NX, NY;
    float_image_get_size(U, &NC, &NX, &NY);
    demand(NC == 1, "image should have only one channel");
    float_image_fill_channel(U, 0, vdef);
    
    for (uint32_t uid = 0; uid < S->N; uid++)
      { int32_t x = S->col[uid];
        int32_t y = S->row[uid];
        assert((x == -1) == (y == -1));
        if (x != -1)
          { demand((x >= 0) && (x < NX), "invalid image col in system");
            demand((y >= 0) && (y < NY), "invalid image row in system");
            double wtot = S->eq[uid].wtot;
            float_image_set_sample(U, 0, x,y, (float)wtot);
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
