/* See cpk_lopt_simple.h */
/* Last edited on 2024-12-31 16:27:03 by stolfi */ 

#include <stdint.h>
#include <assert.h>
#include <math.h>
#include <time.h>
#include <stdlib.h>

#include <vec.h>
#include <bool.h>

#include <cpk_basic.h>
#include <cpk_io.h>
#include <cpk_debug.h>
#include <cpk_valid.h>
#include <cpk_weight.h>
#include <cpk_mis.h>
#include <cpk_graph.h>
#include <cpk_greedy.h>
#include <cpk_lopt.h>

#include <cpk_lopt_simple.h>

/* PROTOTIPOS INTERNOS */

bool_t cpk_local_opt_twiddle_one
  ( cpk_mis_state_t *SQ, /* Solução corrente {S} e estruturas auxiliares. */
    uint32_t u,             /* Vértice a retirar de {S}. */
    bool_t verbose,    /* TRUE para imprimir mensagens de diagnóstico. */
    bool_t temp[]      /* Área de trabalho */
  );
  /* Tenta melhorar a solução {S} eliminando o vértice {u} de {S} e
    cobrindo o buraco de forma a deixar {S} maximal novamente. na
    saída, o peso total {WS} não deve ser menor do que na entrada.
    Devolve {TRUE} se e somente se {S} foi alterado.
    
    A área {temp} deve ter pelo menos {nV} elementos; deve ser toda FALSE
    na entrada e será toda FALSE na saída. */

/* IMPLEMENTAÇÕES */

void cpk_local_opt_simple
  ( cpk_mis_state_t *SQ, /* Solução corrente {S} e estruturas auxiliares. */
    uint32_t *nSM,          /* Tamanho da melhor solução já vista. */
    uint32_t SM[],          /* {SM[0..nSM-1]} é a melhor solução já vista. */
    double *WSM,       /* Peso total da melhor solução {SM}. */
    bool_t verbose     /* TRUE para imprimir mensagens de diagnóstico. */
  )
  {
    uint32_t nV = SQ->G->nV;
    /* Areas de trabalho: */
    bool_t temp[nV];
    for (int32_t v = 0; v < nV; v++) { temp[v] = FALSE; }
    
    /* Repete o passo base até tentar todos os vértices de {S} em seguida sem sucesso: */
    uint32_t iu = 0; /* Índice em {S} do último vértice experimentado: */
    uint32_t nTentou = 0; /* Número de tentativas fracassadas consecutivas. */
    uint32_t nPass = 0; /* Número de tentativas que deram certo. */
    do
      { if (nTentou == 0) { if (verbose) { fprintf(stderr, "\n  Passada %d\n\n", nPass); } }
        /* Ajusta {iu}: */
        iu--; if (iu < 0) { iu = SQ->nS-1; }
        /* Pega o próximo vértice: */
        uint32_t u = SQ->S[iu];
        if (verbose) { fprintf(stderr, "  Tentando %d\n", nPass); }
        /* Lembra do peso atual: */
        double WSold = SQ->WS;
        uint32_t nSold = SQ->nS;
        /* Tente melhorar {S} retirando o vértice {u}: */
        if (cpk_local_opt_twiddle_one(SQ, u, verbose, temp))
          { /* Mudou; deve ser para melhor: */
            assert(cpk_weight_cmp(WSold, SQ->WS) < 0);
            /* Tente atualizar a melhor solução. */
            cpk_lopt_update_best_solution(SQ, nSM, SM, WSM, verbose);
            /* comece tudo de novo: */ 
            nTentou = 0; nPass++;
          }
        else
          { /* Não mudou; deve ter mantido peso e tamanho: */
            assert(nSold == SQ->nS);
            assert(cpk_weight_cmp(WSold, SQ->WS) == 0);
            /* Conte mais um fracasso: */
            nTentou++;
          }
      }
    while (nTentou < SQ->nS);
  }

bool_t cpk_local_opt_twiddle_one
  ( cpk_mis_state_t *SQ, /* Solução corrente {S} e estruturas auxiliares. */
    uint32_t u,             /* Vértice a retirar de {S}. */
    bool_t verbose,    /* TRUE para imprimir mensagens de diagnóstico. */
    bool_t temp[]      /* Área de trabalho */
  )
  { pqueue_t *Q = SQ->Q; 
    assert(SQ->locS[u] != locNONE);
    /* Lembra tamanho e pesos originais: */
    uint32_t nSold = SQ->nS;
    double WSold = SQ->WS;
    /* Retira o vértice {u} do conjunto {S}, atualizando {Q} e pesos: */
    if (verbose) { fprintf(stderr, "  Retirando vértice %d:", u); }
    cpk_mis_delete_indep(SQ, u, verbose);
    if (verbose) { fprintf(stderr, "\n"); }
    assert(pqueue_count(Q) > 0);
    assert(SQ->nS == nSold-1);
    /* Expande {S} para maximal, gulosamente mas com lokahead de 2: */
    if (verbose) { fprintf(stderr, "  Expansão gulosa..."); }
    cpk_greedy_maximize_lookahead(SQ, verbose, temp);
    if (verbose) { fprintf(stderr, " recolocou %d vértices\n", SQ->nS - (nSold-1)); }
    /* Consistências: */
    assert(pqueue_count(Q) == 0);
    if (cpk_weight_cmp(SQ->WS, WSold) > 0)
      { if (verbose) { TRACE_SOL("  Solução corrente melhorou", SQ->nS, SQ->WS); }
        return TRUE;
      }
    else if ((SQ->nS == nSold) && (SQ->S[SQ->nS-1] == u))
      { /* Pelas propriedades de {cpk_maximize_greedy_lookahead}, {S} não mudou: */
        if (verbose) { TRACE_SOL("  Voltou para mesma solução", SQ->nS, SQ->WS); }
        return FALSE;
      }
    else 
      { /* Solução mudou mas não melhorou. Restaure mesmo que igual, para evitar loops: */
        uint32_t v = SQ->S[SQ->nS-1]; /* Vértice acrescentado. */
        if (verbose) 
          { if (cpk_weight_cmp(SQ->WS, WSold) < 0)
              { TRACE_SOL("  Solução corrente piorou", SQ->nS, SQ->WS);
                fprintf(stderr, "   "); 
                fprintf(stderr, " W[%d] = " WT_FFMT, u, SQ->W[u]);
                fprintf(stderr, " W[%d] = " WT_FFMT, v, SQ->W[v]);
              }
            else 
              { TRACE_SOL("  Solução mudou mas não melhorou", SQ->nS, SQ->WS); }
          }
        if (verbose) { fprintf(stderr, "    Retirando vértice %d:", v); }
        cpk_mis_delete_indep(SQ, v, verbose);
        if (verbose) { fprintf(stderr, "\n"); }
        if (verbose) { fprintf(stderr, "    Recolocando vértice %d:", u); }
        cpk_mis_add_indep(SQ, u, verbose);
        if (verbose) { fprintf(stderr, "\n"); }
        if (verbose) { TRACE_SOL("  Solução restaurada", SQ->nS, SQ->WS); }
        assert(cpk_weight_cmp(SQ->WS, WSold) == 0);
        return FALSE;
      }   
  }
