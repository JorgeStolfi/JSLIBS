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
  ( cpk_mis_state_t *SQ,     /* Solu��o corrente {S} e estruturas auxiliares. */
    uint32_t nRem,        /* N�mero de v�rtices a trocar em cada tentativa. */
    uint32_t x,           /* V�rtice piv�. */
    bool_t verbose,  /* TRUE para imprimir mensagens de diagn�stico. */
    bool_t piv[],    /* {piv[v] = TRUE} se vale a pena tentar {v} como piv�. */
    uint32_t T[],
    bool_t temp[]
  );
  /* Se {x} est� em {S} ou n�o tem exatamente {nRem} vizinhos em {S},
    n�o faz nada e retorna FALSE. Sen�o, tenta melhorar {S} retirando todos
    os v�rtices {T} de {S} que s�o adjacentes a {x} (que n�o deve esta
    em {S}) e recobrindo o buraco assim aberto com um algoritmo guloso
    com lookahead. Se isso melhorar a solu��o, devolve TRUE; sen�o
    restaura {S} como antes e devolve FALSE.
    
    O procedimento tamb�m faz {piv[v]= FALSE} para {x} e todo v�rtice
    {v} que ficou descoberto com a retirada de {T}, ou seja, para {x}
    e para todo {v} em {G-S} cujos vizinhos em {S} est�o todos em {T}.
    
    Os par�metros {T} e {temp} s�o �reas de trabalho que devem ser
    alocadas pelo cliente, com pelo menos {nV} elementos. O vetor
    {temp} deve ser todo FALSE na entrada, e ser� todo FALSE na
    sa�da. */

bool_t cpk_local_opt_at_pivot_2
  ( cpk_mis_state_t *SQ,     /* Solu��o corrente {S} e estruturas auxiliares. */
    uint32_t nRem,        /* N�mero de v�rtices a trocar em cada tentativa. */
    uint32_t x,           /* V�rtice piv�. */
    bool_t verbose,  /* TRUE para imprimir mensagens de diagn�stico. */
    bool_t piv[],    /* {piv[v] = TRUE} se vale a pena tentar {v} como piv�. */
    uint32_t T[],
    bool_t temp[]
  );
  /* Semelhante a {}, exceto que o piv� {x} � inclu�do na nova solu��o como v�rtice
    inicial da busca gulosa, e ele � o �nico v�rtice que � marcado com {piv = FALSE}. */

/* IMPLEMENTA��ES */

#define VALIDATE TRUE

void cpk_local_opt_multiple
  ( cpk_mis_state_t *SQ, /* Solu��o corrente {S} e estruturas auxiliares. */
    uint32_t maxRem,        /* N�mero m�ximo de v�rtice a trocar em cada tentativa. */
    uint32_t *nSM,          /* Tamanho da melhor solu��o j� vista. */
    uint32_t SM[],          /* {SM[0..nSM-1]} � a melhor solu��o j� vista. */
    double *WSM,       /* Peso total da melhor solu��o {SM}. */
    bool_t verbose     /* TRUE para imprimir mensagens de diagn�stico. */
  )
  { if (verbose) { cpk_trace_entry(); }
    
    /* A estrat�gia � retirar de {S} grupos de {nRem} v�rtices
    pr�ximos, com {nRem} em {1..maxRem}, e depois completar {S}
    para maximal, gulosamente.

    A cada retirada, escolhe-se um v�rtice /piv�/ {x} de {G-S} com
    exatamente {nRem} vizinhos em {S}. Seja {T} o conjunto desses
    vizinhos. Retira-se {T} de {S}. Seja {Q} o conjunto dos v�rtices
    de {G} que, depois dessa retirada, n�o est�o em {S} nem s�o
    adjacentes a {S}; este conjunto n�o � vazio pois cont�m pelo
    menos {x}. O conjunto {S} � ent�o re-expandido gulosamente �s
    custas de {Q}, at� se tornar maximal.

    Para expandir {S-T}, usamos uma busca gulosa com lookahead (vide
    abaixo). Note que � poss�vel que o novo conjunto maximal seja pior
    que o anterior. Nesse caso a mexida � desfeita, restaurando-se o
    conjunto {S} original.

    Sempre que um v�rtice {w} � inclu�do em {Q}, ele � marcado
    fazendo-se {piv[w] = FALSE}, e exclu�do de
    futuras escolhas de piv� para este {nRem}. Quando todos os 
    v�rtices foram exclu�dos, recome�a-se tudo de novo com 
    {nRem = nRem+1} e {piv[0..nV-1] = TRUE}.
    Quando a solu��o corrente � melhorada,
    recome�a-se tudo de novo com {nRem = 1} e {piv[0..nV-1] = TRUE}. */

    cpk_graph_t *G = SQ->G;
    uint32_t nV = G->nV;
    
    /* Areas de trabalho: */
    bool_t piv[nV];
    uint32_t T[nV];
    bool_t temp[nV];
    for (int32_t v = 0; v < nV; v++) { temp[v] = FALSE; }
    
    /* Inicializa o vetor de v�rtices pivot�veis: */
    for (int32_t v = 0; v < nV; v++) { piv[v] = TRUE; }

    /* Repete busca na vizinhan�a enquanto houver melhoria. */
    uint32_t x = 0; /* �ltimo v�rtice de {G} usado como piv�. */
    uint32_t nRem = 1; /* N�mero de {S}-vizinhos do piv� a remover: */
    uint32_t nPass = 0; /* N�mero de tentativas que deram certo. */
    uint32_t nTentou = 0;  /* N�mero de tentativas fracassadas consecutivas p/ {nRem}. */
    do 
      { if (nTentou == 0) 
          { if (verbose) { fprintf(stderr, "\n  Passada %d nRem = %d\n\n", nPass, nRem); } }
        /* Pega pr�ximo piv�: */
        x--; if (x < 0) { x = nV-1; }

        /* Lembra do peso atual: */
        double WSold = SQ->WS;
        uint32_t nSold = SQ->nS;
        /* Tenta troca com {nRem} v�rtices na vizinhan�a de {x}: */
        if (piv[x] && cpk_local_opt_at_pivot_1(SQ, nRem, x, verbose, piv, T, temp))
          { /* Mudou; deve ser para melhor: */
            assert(cpk_weight_cmp(WSold, SQ->WS) < 0);
            /* Tente atualizar a melhor solu��o. */
            cpk_lopt_update_best_solution(SQ, nSM, SM, WSM, verbose);
            /* Comece tudo de novo: */ 
            nTentou = 0; nRem = 1; nPass++;
            /* Reinicializa o vetor {piv}. */
            for (int32_t v = 0; v < nV; v++) { piv[v] = TRUE; }
          }
        else
          { /* N�o mudou; deve ter mantido peso e tamanho: */
            assert(nSold == SQ->nS);
            assert(cpk_weight_cmp(WSold, SQ->WS) == 0);
            /* Por via das d�vidas: */
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
  ( cpk_mis_state_t *SQ,     /* Solu��o corrente {S} e estruturas auxiliares. */
    uint32_t nRem,        /* N�mero de v�rtices a trocar em cada tentativa. */
    uint32_t x,           /* V�rtice piv�. */
    bool_t verbose,  /* TRUE para imprimir mensagens de diagn�stico. */
    bool_t piv[],    /* {piv[v] = TRUE} se vale a pena tentar {v} como piv�. */
    uint32_t T[],
    bool_t temp[]
  )
  {
    cpk_graph_t *G = SQ->G;
    pqueue_t *Q = SQ->Q;
    
    assert(piv[x]); /* S� deve chamar se {x} ainda est� marcado como poss�vel piv�. */
    assert(pqueue_count(Q) == 0); /* Neste ponto {S} deveria ser maximal. */

    /* O v�rtice {x} realmente satisfaz as condi��es para piv�? */
    if ((SQ->locS[x] != locNONE) || (SQ->degS[x] != nRem)) { return FALSE; } 
    
    if (verbose) { fprintf(stderr, "    Mexendo nas redondezas do v�rtice %d\n", x); }

    /* Lembra do peso da solu��o corrente: */
    double WSold = SQ->WS;

    /* Retire de {S} todos os vizinhos de {x}, salvando-os em {T[0..nT-1]}. */
    /* Coloque os v�rtices que ficaram admiss�veis em {Q}. */
    if (verbose) { fprintf(stderr, "      Retirando %d v�rtices:", SQ->degS[x]); }
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

    /* Marca todos os v�rtices de {Q} como n�o mais pivot�veis: */
    uint32_t nQ = pqueue_count(Q);
    if (verbose) { fprintf(stderr, "      Criados %d v�rtices admiss�veis:", nQ);}
    assert(nQ > 0); /* Pois pelo menos {x} deve ter ficado admiss�vel. */
    { for (uint32_t i = 0; i < nQ; i++)
        { pqueue_item_t q = pqueue_item(Q, i);
          piv[q] = FALSE; 
          if (verbose) { fprintf(stderr, " %d", q);}
        }
    }
    if (verbose) { fprintf(stderr, "\n");}
    /* TRACE_Q("     Q =", Q); */

    /* Lembra posi��o corrente em {S}: */
    uint32_t nSred = SQ->nS; 

    /* Roda algoritmo greedy sobre {G[Q]}: */
    if (verbose) { fprintf(stderr, "      Expans�o gulosa..."); }
    cpk_greedy_maximize_lookahead(SQ, verbose, temp);
    if (verbose) { fprintf(stderr, " recolocou %d v�rtices\n", SQ->nS - nSred);  }

    /* Houve melhora? */
    if (cpk_weight_cmp(SQ->WS, WSold) > 0)
      { /* Houve melhora! */
        if (verbose) { TRACE_SOL("    Solu��o corrente melhorou", SQ->nS, SQ->WS); }
        return TRUE;
      }
    else if (cpk_lopt_same_set(nT, T, SQ->nS - nSred, &(SQ->S[nSred])))
      { /* Voltamos � solu��o de partida. */
        if (verbose) { TRACE_SOL("    Ficou na mesma solu��o", SQ->nS, SQ->WS); }
        return FALSE;
      }
    else
      { /* A emenda n�o foi melhor que o soneto, o n�g�cio � desfazer a troca. */
        uint32_t nDel = SQ->nS - nSred;
        if (verbose)
          {  if (cpk_weight_cmp(SQ->WS, WSold) < 0)
              { TRACE_SOL("      Solu��o corrente piorou", SQ->nS, SQ->WS);
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
              { TRACE_SOL("      Solu��o mudou mas n�o melhorou", SQ->nS, SQ->WS); }
          }
        cpk_lopt_restore_solution(SQ, nDel, nT, T, verbose);
        if (verbose) { TRACE_SOL("    Solu��o restaurada", SQ->nS, SQ->WS); }
        assert(cpk_weight_cmp(SQ->WS, WSold) == 0);
        return FALSE;
      }
  }
  
bool_t cpk_local_opt_at_pivot_2
  ( cpk_mis_state_t *SQ,     /* Solu��o corrente {S} e estruturas auxiliares. */
    uint32_t nRem,        /* N�mero de v�rtices a trocar em cada tentativa. */
    uint32_t x,           /* V�rtice piv�. */
    bool_t verbose,  /* TRUE para imprimir mensagens de diagn�stico. */
    bool_t piv[],    /* {piv[v] = TRUE} se vale a pena tentar {v} como piv�. */
    uint32_t T[],
    bool_t temp[]
  )
  {
    cpk_graph_t *G = SQ->G;
    pqueue_t *Q = SQ->Q;
    
    assert(piv[x]); /* S� deve chamar se {x} ainda est� marcado como poss�vel piv�. */
    assert(pqueue_count(Q) == 0); /* Neste ponto {S} deveria ser maximal. */

    /* O v�rtice {x} realmente satisfaz as condi��es para piv�? */
    if ((SQ->locS[x] != locNONE) || (SQ->degS[x] != nRem)) { return FALSE; } 
    
    if (verbose) { fprintf(stderr, "    Mexendo nas redondezas do v�rtice %d\n", x); }

    /* Lembra do peso da solu��o corrente: */
    double WSold = SQ->WS;

    /* Retire de {S} todos os vizinhos de {x}, salvando-os em {T[0..nT-1]}. */
    /* Coloque os v�rtices que ficaram admiss�veis em {Q}. */
    if (verbose) { fprintf(stderr, "      Retirando %d v�rtices:", SQ->degS[x]); }
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

    /* Lembra posi��o corrente em {S}: */
    uint32_t nSred = SQ->nS; 

    /* Coloca o v�rtice {x} no conjunto {S} na marra: */
    if (verbose) { fprintf(stderr, "      Colocando %d na solu��o...", x); }
    cpk_mis_add_indep(SQ, x, verbose);
    if (verbose) { fprintf(stderr, "\n");  }

    /* Roda algoritmo greedy sobre {G[Q]}: */
    if (verbose) { fprintf(stderr, "      Expans�o gulosa..."); }
    cpk_greedy_maximize_lookahead(SQ, verbose, temp);
    if (verbose) { fprintf(stderr, " recolocou mais %d v�rtices\n", SQ->nS - nSred);  }
    if (verbose) { TRACE_SOL("      Solu��o corrente", SQ->nS, SQ->WS); }

    /* Marca os v�rtices acrescentados como n�o mais pivot�veis: */
    if (verbose) { fprintf(stderr, "      Eliminando futuros piv�s:"); }
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
        if (verbose) { TRACE_SOL("    Solu��o corrente melhorou", SQ->nS, SQ->WS); }
        return TRUE;
      }
    else if (cpk_lopt_same_set(nT, T, SQ->nS - nSred, &(SQ->S[nSred])))
      { /* Voltamos � solu��o de partida. */
        if (verbose) { TRACE_SOL("    Ficou na mesma solu��o", SQ->nS, SQ->WS); }
        return FALSE;
      }
    else
      { /* A emenda n�o foi melhor que o soneto, o n�g�cio � desfazer a troca. */
        uint32_t nDel = SQ->nS - nSred;
        if (verbose)
          {  if (cpk_weight_cmp(SQ->WS, WSold) < 0)
              { TRACE_SOL("      Solu��o corrente piorou", SQ->nS, SQ->WS);
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
              { TRACE_SOL("      Solu��o mudou mas n�o melhorou", SQ->nS, SQ->WS); }
          }
        cpk_lopt_restore_solution(SQ, nDel, nT, T, verbose);
        if (verbose) { TRACE_SOL("    Solu��o restaurada", SQ->nS, SQ->WS); }
        assert(cpk_weight_cmp(SQ->WS, WSold) == 0);
        return FALSE;
      }
  }
  
    
