/* See cpk_lopt_multi.h */
/* Last edited on 2024-12-31 16:26:54 by stolfi */ 

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
#include <cpk_weight.h>
#include <cpk_mis.h>
#include <cpk_graph.h>
#include <cpk_greedy.h>
#include <cpk_lopt.h>
#include <cpk_basic.h>
#include <cpk_lopt_multi.h>

/* PROTOTIPOS INTERNOS */

bool_t cpk_local_opt_at_pivot_1
  ( cpk_mis_state_t *SQ,     /* Solução corrente {S} e estruturas auxiliares. */
    uint32_t nRem,        /* Número de vértices a trocar em cada tentativa. */
    uint32_t x,           /* Vértice pivô. */
    bool_t verbose,  /* TRUE para imprimir mensagens de diagnóstico. */
    bool_t piv[],    /* {piv[v] = TRUE} se vale a pena tentar {v} como pivô. */
    uint32_t T[],
    bool_t temp[]
  );
  /* Se {x} está em {S} ou não tem exatamente {nRem} vizinhos em {S},
    não faz nada e retorna FALSE. Senão, tenta melhorar {S} retirando todos
    os vértices {T} de {S} que são adjacentes a {x} (que não deve esta
    em {S}) e recobrindo o buraco assim aberto com um algoritmo guloso
    com lookahead. Se isso melhorar a solução, devolve TRUE; senão
    restaura {S} como antes e devolve FALSE.
    
    O procedimento também faz {piv[v]= FALSE} para {x} e todo vértice
    {v} que ficou descoberto com a retirada de {T}, ou seja, para {x}
    e para todo {v} em {G-S} cujos vizinhos em {S} estão todos em {T}.
    
    Os parâmetros {T} e {temp} são áreas de trabalho que devem ser
    alocadas pelo cliente, com pelo menos {nV} elementos. O vetor
    {temp} deve ser todo FALSE na entrada, e será todo FALSE na
    saída. */

bool_t cpk_local_opt_at_pivot_2
  ( cpk_mis_state_t *SQ,     /* Solução corrente {S} e estruturas auxiliares. */
    uint32_t nRem,        /* Número de vértices a trocar em cada tentativa. */
    uint32_t x,           /* Vértice pivô. */
    bool_t verbose,  /* TRUE para imprimir mensagens de diagnóstico. */
    bool_t piv[],    /* {piv[v] = TRUE} se vale a pena tentar {v} como pivô. */
    uint32_t T[],
    bool_t temp[]
  );
  /* Semelhante a {}, exceto que o pivô {x} é incluído na nova solução como vértice
    inicial da busca gulosa, e ele é o único vértice que é marcado com {piv = FALSE}. */

/* IMPLEMENTAÇÕES */

#define VALIDATE TRUE

void cpk_local_opt_multiple
  ( cpk_mis_state_t *SQ, /* Solução corrente {S} e estruturas auxiliares. */
    uint32_t maxRem,        /* Número máximo de vértice a trocar em cada tentativa. */
    uint32_t *nSM,          /* Tamanho da melhor solução já vista. */
    uint32_t SM[],          /* {SM[0..nSM-1]} é a melhor solução já vista. */
    double *WSM,       /* Peso total da melhor solução {SM}. */
    bool_t verbose     /* TRUE para imprimir mensagens de diagnóstico. */
  )
  { if (verbose) { cpk_trace_entry(); }
    
    /* A estratégia é retirar de {S} grupos de {nRem} vértices
    próximos, com {nRem} em {1..maxRem}, e depois completar {S}
    para maximal, gulosamente.

    A cada retirada, escolhe-se um vértice /pivô/ {x} de {G-S} com
    exatamente {nRem} vizinhos em {S}. Seja {T} o conjunto desses
    vizinhos. Retira-se {T} de {S}. Seja {Q} o conjunto dos vértices
    de {G} que, depois dessa retirada, não estão em {S} nem são
    adjacentes a {S}; este conjunto não é vazio pois contém pelo
    menos {x}. O conjunto {S} é então re-expandido gulosamente às
    custas de {Q}, até se tornar maximal.

    Para expandir {S-T}, usamos uma busca gulosa com lookahead (vide
    abaixo). Note que é possível que o novo conjunto maximal seja pior
    que o anterior. Nesse caso a mexida é desfeita, restaurando-se o
    conjunto {S} original.

    Sempre que um vértice {w} é incluído em {Q}, ele é marcado
    fazendo-se {piv[w] = FALSE}, e excluído de
    futuras escolhas de pivô para este {nRem}. Quando todos os 
    vértices foram excluídos, recomeça-se tudo de novo com 
    {nRem = nRem+1} e {piv[0..nV-1] = TRUE}.
    Quando a solução corrente é melhorada,
    recomeça-se tudo de novo com {nRem = 1} e {piv[0..nV-1] = TRUE}. */

    cpk_graph_t *G = SQ->G;
    uint32_t nV = G->nV;
    
    /* Areas de trabalho: */
    bool_t piv[nV];
    uint32_t T[nV];
    bool_t temp[nV];
    for (int32_t v = 0; v < nV; v++) { temp[v] = FALSE; }
    
    /* Inicializa o vetor de vértices pivotáveis: */
    for (int32_t v = 0; v < nV; v++) { piv[v] = TRUE; }

    /* Repete busca na vizinhança enquanto houver melhoria. */
    uint32_t x = 0; /* Último vértice de {G} usado como pivô. */
    uint32_t nRem = 1; /* Número de {S}-vizinhos do pivô a remover: */
    uint32_t nPass = 0; /* Número de tentativas que deram certo. */
    uint32_t nTentou = 0;  /* Número de tentativas fracassadas consecutivas p/ {nRem}. */
    do 
      { if (nTentou == 0) 
          { if (verbose) { fprintf(stderr, "\n  Passada %d nRem = %d\n\n", nPass, nRem); } }
        /* Pega próximo pivô: */
        x--; if (x < 0) { x = nV-1; }

        /* Lembra do peso atual: */
        double WSold = SQ->WS;
        uint32_t nSold = SQ->nS;
        /* Tenta troca com {nRem} vértices na vizinhança de {x}: */
        if (piv[x] && cpk_local_opt_at_pivot_1(SQ, nRem, x, verbose, piv, T, temp))
          { /* Mudou; deve ser para melhor: */
            assert(cpk_weight_cmp(WSold, SQ->WS) < 0);
            /* Tente atualizar a melhor solução. */
            cpk_lopt_update_best_solution(SQ, nSM, SM, WSM, verbose);
            /* Comece tudo de novo: */ 
            nTentou = 0; nRem = 1; nPass++;
            /* Reinicializa o vetor {piv}. */
            for (int32_t v = 0; v < nV; v++) { piv[v] = TRUE; }
          }
        else
          { /* Não mudou; deve ter mantido peso e tamanho: */
            assert(nSold == SQ->nS);
            assert(cpk_weight_cmp(WSold, SQ->WS) == 0);
            /* Por via das dúvidas: */
            piv[x] = FALSE;
            /* Mais uma tentativa fracassada: */
            nTentou++;
            /* Esgotamos este {nRem}? */
            if (nTentou >= nV)
              { /* Vamos para o seguinte: */
                nRem++; nTentou = 0; 
                /* Reinicializa {piv}: */
                for (int32_t v = 0; v < nV; v++) { piv[v] = TRUE; }
              }
          } 
      } 
    while (nRem <= maxRem);
    if (VALIDATE) { cpk_check_indep_set(G, SQ->nS, SQ->S, NULL); }
    if (verbose) { cpk_trace_exit(); }
  }

bool_t cpk_local_opt_at_pivot_1
  ( cpk_mis_state_t *SQ,     /* Solução corrente {S} e estruturas auxiliares. */
    uint32_t nRem,        /* Número de vértices a trocar em cada tentativa. */
    uint32_t x,           /* Vértice pivô. */
    bool_t verbose,  /* TRUE para imprimir mensagens de diagnóstico. */
    bool_t piv[],    /* {piv[v] = TRUE} se vale a pena tentar {v} como pivô. */
    uint32_t T[],
    bool_t temp[]
  )
  {
    cpk_graph_t *G = SQ->G;
    pqueue_t *Q = SQ->Q;
    
    assert(piv[x]); /* Só deve chamar se {x} ainda está marcado como possível pivô. */
    assert(pqueue_count(Q) == 0); /* Neste ponto {S} deveria ser maximal. */

    /* O vértice {x} realmente satisfaz as condições para pivô? */
    if ((SQ->locS[x] != locNONE) || (SQ->degS[x] != nRem)) { return FALSE; } 
    
    if (verbose) { fprintf(stderr, "    Mexendo nas redondezas do vértice %d\n", x); }

    /* Lembra do peso da solução corrente: */
    double WSold = SQ->WS;

    /* Retire de {S} todos os vizinhos de {x}, salvando-os em {T[0..nT-1]}. */
    /* Coloque os vértices que ficaram admissíveis em {Q}. */
    if (verbose) { fprintf(stderr, "      Retirando %d vértices:", SQ->degS[x]); }
    uint32_t nT = 0;
    { uint32_t dx = G->deg[x]; uint32_t *nbx = &(G->nbr[G->fnb[x]]);
      for (int32_t i = 0; i < dx; i++)
        { uint32_t v = nbx[i]; 
          if (SQ->locS[v] != locNONE) 
            { cpk_mis_delete_indep(SQ, v, verbose); T[nT] = v; nT++; }
        }
    }
    assert(SQ->degS[x] == 0); /* Objetivo da retirada. */
    if (verbose) { fprintf(stderr, "\n"); }

    /* Marca todos os vértices de {Q} como não mais pivotáveis: */
    uint32_t nQ = pqueue_count(Q);
    if (verbose) { fprintf(stderr, "      Criados %d vértices admissíveis:", nQ);}
    assert(nQ > 0); /* Pois pelo menos {x} deve ter ficado admissível. */
    { for (uint32_t i = 0; i < nQ; i++)
        { pqueue_item_t q = pqueue_item(Q, i);
          piv[q] = FALSE; 
          if (verbose) { fprintf(stderr, " %d", q);}
        }
    }
    if (verbose) { fprintf(stderr, "\n");}
    /* TRACE_Q("     Q =", Q); */

    /* Lembra posição corrente em {S}: */
    uint32_t nSred = SQ->nS; 

    /* Roda algoritmo greedy sobre {G[Q]}: */
    if (verbose) { fprintf(stderr, "      Expansão gulosa..."); }
    cpk_greedy_maximize_lookahead(SQ, verbose, temp);
    if (verbose) { fprintf(stderr, " recolocou %d vértices\n", SQ->nS - nSred);  }

    /* Houve melhora? */
    if (cpk_weight_cmp(SQ->WS, WSold) > 0)
      { /* Houve melhora! */
        if (verbose) { TRACE_SOL("    Solução corrente melhorou", SQ->nS, SQ->WS); }
        return TRUE;
      }
    else if (cpk_lopt_same_set(nT, T, SQ->nS - nSred, &(SQ->S[nSred])))
      { /* Voltamos à solução de partida. */
        if (verbose) { TRACE_SOL("    Ficou na mesma solução", SQ->nS, SQ->WS); }
        return FALSE;
      }
    else
      { /* A emenda não foi melhor que o soneto, o négócio é desfazer a troca. */
        uint32_t nDel = SQ->nS - nSred;
        if (verbose)
          {  if (cpk_weight_cmp(SQ->WS, WSold) < 0)
              { TRACE_SOL("      Solução corrente piorou", SQ->nS, SQ->WS);
                fprintf(stderr, "     ");
                for (uint32_t i = 0; i < nDel; i++) 
                  { uint32_t w = SQ->S[nSred+i]; fprintf(stderr, " W[%d]=" WT_FFMT, w, SQ->W[w]); }
                fprintf(stderr, "\n");
                fprintf(stderr, "     ");
                for (int32_t i = 0; i < nT; i++) 
                  { uint32_t w = T[i]; fprintf(stderr, " W[%d]=" WT_FFMT, w, SQ->W[w]); }
                fprintf(stderr, "\n");
              }
            else 
              { TRACE_SOL("      Solução mudou mas não melhorou", SQ->nS, SQ->WS); }
          }
        cpk_lopt_restore_solution(SQ, nDel, nT, T, verbose);
        if (verbose) { TRACE_SOL("    Solução restaurada", SQ->nS, SQ->WS); }
        assert(cpk_weight_cmp(SQ->WS, WSold) == 0);
        return FALSE;
      }
  }
  
bool_t cpk_local_opt_at_pivot_2
  ( cpk_mis_state_t *SQ,     /* Solução corrente {S} e estruturas auxiliares. */
    uint32_t nRem,        /* Número de vértices a trocar em cada tentativa. */
    uint32_t x,           /* Vértice pivô. */
    bool_t verbose,  /* TRUE para imprimir mensagens de diagnóstico. */
    bool_t piv[],    /* {piv[v] = TRUE} se vale a pena tentar {v} como pivô. */
    uint32_t T[],
    bool_t temp[]
  )
  {
    cpk_graph_t *G = SQ->G;
    pqueue_t *Q = SQ->Q;
    
    assert(piv[x]); /* Só deve chamar se {x} ainda está marcado como possível pivô. */
    assert(pqueue_count(Q) == 0); /* Neste ponto {S} deveria ser maximal. */

    /* O vértice {x} realmente satisfaz as condições para pivô? */
    if ((SQ->locS[x] != locNONE) || (SQ->degS[x] != nRem)) { return FALSE; } 
    
    if (verbose) { fprintf(stderr, "    Mexendo nas redondezas do vértice %d\n", x); }

    /* Lembra do peso da solução corrente: */
    double WSold = SQ->WS;

    /* Retire de {S} todos os vizinhos de {x}, salvando-os em {T[0..nT-1]}. */
    /* Coloque os vértices que ficaram admissíveis em {Q}. */
    if (verbose) { fprintf(stderr, "      Retirando %d vértices:", SQ->degS[x]); }
    uint32_t nT = 0;
    { uint32_t dx = G->deg[x]; uint32_t *nbx = &(G->nbr[G->fnb[x]]);
      for (int32_t i = 0; i < dx; i++)
        { uint32_t v = nbx[i]; 
          if (SQ->locS[v] != locNONE) 
            { cpk_mis_delete_indep(SQ, v, verbose); T[nT] = v; nT++; }
        }
    }
    assert(SQ->degS[x] == 0); /* Objetivo da retirada. */
    if (verbose) { fprintf(stderr, "\n"); }

    /* Lembra posição corrente em {S}: */
    uint32_t nSred = SQ->nS; 

    /* Coloca o vértice {x} no conjunto {S} na marra: */
    if (verbose) { fprintf(stderr, "      Colocando %d na solução...", x); }
    cpk_mis_add_indep(SQ, x, verbose);
    if (verbose) { fprintf(stderr, "\n");  }

    /* Roda algoritmo greedy sobre {G[Q]}: */
    if (verbose) { fprintf(stderr, "      Expansão gulosa..."); }
    cpk_greedy_maximize_lookahead(SQ, verbose, temp);
    if (verbose) { fprintf(stderr, " recolocou mais %d vértices\n", SQ->nS - nSred);  }
    if (verbose) { TRACE_SOL("      Solução corrente", SQ->nS, SQ->WS); }

    /* Marca os vértices acrescentados como não mais pivotáveis: */
    if (verbose) { fprintf(stderr, "      Eliminando futuros pivôs:"); }
    { for (uint32_t i = nSred; i < SQ->nS; i++)
        { uint32_t s = SQ->S[i];
          piv[s] = FALSE; 
          if (verbose) { fprintf(stderr, " %d", s);}
        }
    }
    if (verbose) { fprintf(stderr, "\n");}

    /* Houve melhora? */
    if (cpk_weight_cmp(SQ->WS, WSold) > 0)
      { /* Houve melhora! */
        if (verbose) { TRACE_SOL("    Solução corrente melhorou", SQ->nS, SQ->WS); }
        return TRUE;
      }
    else if (cpk_lopt_same_set(nT, T, SQ->nS - nSred, &(SQ->S[nSred])))
      { /* Voltamos à solução de partida. */
        if (verbose) { TRACE_SOL("    Ficou na mesma solução", SQ->nS, SQ->WS); }
        return FALSE;
      }
    else
      { /* A emenda não foi melhor que o soneto, o négócio é desfazer a troca. */
        uint32_t nDel = SQ->nS - nSred;
        if (verbose)
          {  if (cpk_weight_cmp(SQ->WS, WSold) < 0)
              { TRACE_SOL("      Solução corrente piorou", SQ->nS, SQ->WS);
                fprintf(stderr, "     ");
                for (uint32_t i = 0; i < nDel; i++) 
                  { uint32_t w = SQ->S[nSred+i]; fprintf(stderr, " W[%d]=" WT_FFMT, w, SQ->W[w]); }
                fprintf(stderr, "\n");
                fprintf(stderr, "     ");
                for (int32_t i = 0; i < nT; i++) 
                  { uint32_t w = T[i]; fprintf(stderr, " W[%d]=" WT_FFMT, w, SQ->W[w]); }
                fprintf(stderr, "\n");
              }
            else 
              { TRACE_SOL("      Solução mudou mas não melhorou", SQ->nS, SQ->WS); }
          }
        cpk_lopt_restore_solution(SQ, nDel, nT, T, verbose);
        if (verbose) { TRACE_SOL("    Solução restaurada", SQ->nS, SQ->WS); }
        assert(cpk_weight_cmp(SQ->WS, WSold) == 0);
        return FALSE;
      }
  }
  
    
