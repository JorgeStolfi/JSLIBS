/* See cpk_lopt.h */
/* Last edited on 2024-12-31 15:57:42 by stolfi */ 

#include <stdint.h>
#include <assert.h>
#include <math.h>
#include <time.h>
#include <stdlib.h>

#include <vec.h>
#include <bool.h>

#include <cpk_io.h>
#include <cpk_debug.h>
#include <cpk_valid.h>
#include <cpk_mis.h>
#include <cpk_graph.h>
#include <cpk_lopt_multi.h>
#include <cpk_lopt_simple.h>
#include <cpk_basic.h>

#include <cpk_lopt.h>

/* PROTOTIPOS INTERNOS */

/* IMPLEMENTA��ES */

void cpk_local_opt
  ( cpk_mis_state_t *SQ, /* Solu��o corrente {S} e estruturas auxiliares. */
    uint32_t maxRem,        /* N�mero m�ximo de v�rtice a trocar em cada tentativa. */
    uint32_t *nSM,          /* Tamanho da melhor solu��o j� vista. */
    uint32_t SM[],          /* {SM[0..nSM-1]} � a melhor solu��o j� vista. */
    double *WSM,       /* Peso total da melhor solu��o {SM}. */
    bool_t verbose     /* TRUE para imprimir mensagens de diagn�stico. */
  )
  { if (maxRem == 1) 
      { cpk_local_opt_simple(SQ, nSM, SM, WSM, verbose); }
    else
      { cpk_local_opt_multiple(SQ, maxRem, nSM, SM, WSM, verbose); }
  }

void cpk_lopt_update_best_solution
  ( cpk_mis_state_t *SQ, /* Solu��o corrente {S} e estruturas auxiliares. */
    uint32_t *nSM,          /* Tamanho da melhor solu��o j� vista. */
    uint32_t SM[],          /* {SM[0..nSM-1]} � a melhor solu��o j� vista. */
    double *WSM,       /* Peso total da melhor solu��o {SM}. */
    bool_t verbose     /* TRUE para imprimir mensagens de diagn�stico. */
  )
  {
    if (cpk_weight_cmp(SQ->WS, *WSM) > 0)
      { /* O resultado � melhor que {SM}; copia {S} para {SM}: */
        for (int32_t i = 0; i < SQ->nS; i++) { SM[i] = SQ->S[i]; }
        (*nSM) = SQ->nS;
        (*WSM) = SQ->WS; 
        if (verbose)
          { TRACE_SOL("  Melhorou solu��o global", (*nSM), (*WSM));
            cpk_print_vertex_set(stderr, "", (*nSM), SM, NULL, 30);
            if (verbose) { fprintf(stderr, "\n");  }
          }
      }
  }

void cpk_lopt_restore_solution
  ( cpk_mis_state_t *SQ, /* Solu��o corrente {S} e estruturas auxiliares. */
    uint32_t nDel,        /* N�mero de v�rtices a retirar. */
    uint32_t nAdd,        /* N�mero de v�rtices a recolocar. */
    uint32_t Add[],       /* {Add[0..nAdd-1]} s�o os v�rtices a recolocar. */
    bool_t verbose   /* TRUE para imprimir mensagens de diagn�stico. */
  )
  {
    if (nDel > 0)
      { if (verbose) { fprintf(stderr, "      Retirando os �ltimos %d v�rtices", nDel); }
        uint32_t nSred = SQ->nS - nDel;  /* N�mero de v�rtices a deixar: */
        while(SQ->nS > nSred) { uint32_t v = SQ->S[SQ->nS-1];  cpk_mis_delete_indep(SQ, v, verbose); }
        if (verbose) { fprintf(stderr, "\n");  }
      }
    if (nAdd > 0)
      { if (verbose) { fprintf(stderr, "      Recolocando %d v�rtices", nAdd); }
        while (nAdd > 0) { nAdd--; uint32_t v = Add[nAdd]; cpk_mis_add_indep(SQ, v, verbose); }
        if (verbose) { fprintf(stderr, "\n");  }
      }
  }

bool_t cpk_lopt_same_set(uint32_t nX, uint32_t X[], uint32_t nY, uint32_t Y[])
  {
    if (nX != nY) { return FALSE; }
    for (int32_t iX = 0; iX < nX; iX++)
      { int32_t iY = 0;
        while ((iY < nY) && (X[iX] != Y[iY])) { iY++; }
        if (iY >= nY) { return FALSE; }
      }
    return TRUE;
  }
