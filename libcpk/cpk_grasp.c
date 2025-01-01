/* See cpk_grasp.h */
/* Last edited on 2024-12-31 16:09:57 by stolfi */

#include <stdint.h>
#include <math.h>
#include <stdlib.h>
#include <time.h>
#include <assert.h>

#include <cpk_grasp.h>
#include <cpk_mis.h>
#include <cpk_graph.h>
#include <cpk_greedy.h>
#include <cpk_lopt.h>
#include <cpk_valid.h>
#include <cpk_debug.h>
#include <cpk_io.h>
#include <cpk_basic.h>

#include <bool.h>
#include <vec.h>
#include <affirm.h>
#include <jsmath.h>

typedef struct cpk_seed_pair_t 
  { uint32_t v1;
    uint32_t v2;
    double escore;
  } cpk_seed_pair_t;
  /* Estrutura usada para encontrar pares de v�rtices independentes de
    alto peso e baixo consumo. */

/* PROT�TIPOS */

ui2_vec_t cpk_grasp_get_seed_pairs(cpk_graph_t *G, double W[], uint32_t maxPairs, bool_t verbose);
  /* Devolve uma lista {K} de pares de v�rtices independentes de {G},
     para serem usados como solu��es de partida da fase de constru��o
     do GRASP. A lista ter� no m�ximo {maxPairs} pares. */

double cpk_grasp_pair_score
  ( cpk_graph_t *G, 
    double W[], 
    uint32_t u, 
    uint32_t v, 
    bool_t temp[]
  );
  /* Calcula o escore de um par de v�rtices {u,v} considerado como
    conjunto inicial para a fase de constru��o do GRASP. O escore �
    -INF se {u} e {v} s�o vizinhos, caso contr�rio � a soma dos pesos
    de {u} e {v} dividida pelo tamanho da uni�o de suas vizinhan�as
    (incluindo {u} e {v}).
    
    O vetor de booleanos {temp} � usado internamente; ele deve ser
    todo FALSE na entrada e ser� todo FALSE na sa�da. */

/* IMPLEMENTA��ES */

#define VALIDATE TRUE

uint32_vec_t cpk_grasp_find_indep_set
  ( cpk_graph_t *G,
    double W[], 
    double maxClock,  /* Tempo limite para retornar a resposta. */
    uint32_t seed,    /* Semente para escolhas aleat�rias. */
    bool_t verbose
  )
  {
    /* As fases de constru��o e busca sao repetidas para varios pares
    de v�rtices independentes de partida, ciclicamente. */

    if (verbose) { cpk_trace_entry(); }

    uint32_t nV = G->nV;

    /* Melhor solu��o encontrada at� o momento (inicialmente vazia): */
    uint32_t *SM = talloc(nV, uint32_t);
    uint32_t nSM = 0;     /* A melhor solu��o � {SM[0..nSM-1]}. */
    double WSM = 0;  /* Peso total de {SM}. */

    /* Inicializa gerador de n�meros aleat�rios: */
    srand(seed);
    
    if (nV > 0)
      { /* Decide quantos pares-semente usar. Melhor ficar com poucos e bons, acho: */
        uint32_t maxSeedPairs = (uint32_t)ceil(sqrt((double)nV));
        /* Constr�i um conjunto de no m�ximo {maxSeedPairs} pares de v�rtices iniciais: */
        ui2_vec_t seedPair = cpk_grasp_get_seed_pairs(G, W, maxSeedPairs, verbose);
        if (verbose) { cpk_trace(); }
        
        /* Representa��o do conjunto independente {S} e admiss�vel {Q}: */
        cpk_mis_state_t *SQ = cpk_mis_state_new(G, W);

        /* Roda v�rias tentativas GRASP: */
        uint32_t nIter = 0;
        do
          { if (verbose)  { fprintf(stderr, "Itera��o %d\n", nIter); }

            /* FASE DE CONSTRU��O */

            /* Esvazia {S}, coloca todos os n�s em {Q}, ordenados com base em WR: */
            cpk_mis_init(SQ, verbose);

            if (nIter == 0)
              { /* Na itera��o 0, come�a com a pilha vazia: */
              }
            else
              { /* Incializa {S} com uma semente variada: */
                if (verbose) { fprintf(stderr, "\n  ------------------------\n"); }
                if (seedPair.ne > 0)
                  { /* Inicializa {S} com o pr�ximo par-semente: */
                    if (verbose) { fprintf(stderr, "  Inserindo um par inicial"); }
                    { uint32_t t = nIter % seedPair.ne;
                      uint32_t u = (uint32_t)seedPair.e[t].c[0];
                      uint32_t v = (uint32_t)seedPair.e[t].c[1];
                      cpk_mis_add_indep(SQ, u, verbose);
                      cpk_mis_add_indep(SQ, v, verbose);
                    }
                  }
                else 
                  { /* Inicializa {S} com o pr�ximo v�rtice: */
                    uint32_t u = nIter % nV;
                    cpk_mis_add_indep(SQ, u, verbose);
                  }
                if (verbose) { fprintf(stderr, "\n"); }
                if (verbose) { TRACE_SOL("  Solu��o semente", SQ->nS, SQ->WS);  }
              }
            /* Expande a solu��o aleat�rio-gulosamente at� ficar maximal: */
            if (verbose) { fprintf(stderr, "\n  ------------------------\n"); }
            if (verbose) { fprintf(stderr, "  Expans�o gulosa randomizada:"); }
            cpk_greedy_maximize_random(SQ, verbose);
            if (verbose) { fprintf(stderr, "\n");  }
            if (verbose) { TRACE_SOL("  Solu��o gulosa randomizada", SQ->nS, SQ->WS);  }

            if (VALIDATE) { cpk_check_indep_set(G, SQ->nS, SQ->S, NULL); }

            /* Atualiza a melhor solu��o: */
            cpk_lopt_update_best_solution(SQ, &nSM, SM, &WSM, verbose);

            /* FASE DE BUSCA */

            /* At� quantos v�rtices vale a pena tentar relocar por vez? */
            uint32_t maxRem = 3;

            /* Faz otimiza��o local com retirada de m�ltiplos v�rtices: */
            assert(pqueue_count(SQ->Q) == 0);
            if (verbose) { fprintf(stderr, "\n  ------------------------\n"); }
            if (verbose) { fprintf(stderr, "  Otimiza��o local:"); }
            cpk_local_opt(SQ, maxRem,  &nSM, SM, &WSM, verbose);
            if (verbose) { fprintf(stderr, "\n");  }
            if (verbose) { TRACE_SOL("  Solu��o localmente otimizada", SQ->nS, SQ->WS);  }

            if (VALIDATE)  
              { cpk_check_indep_set(G, SQ->nS, SQ->S, NULL);
                cpk_check_indep_set(G, nSM, SM, NULL);
              }

            /* Mais uma itera��o: */
            nIter++;
          }
        while(cpk_cpu_time_1() < maxClock);

        /* Copia melhor solu��o para solu��o corrente: */
        cpk_mis_copy_indep(SQ, nSM, SM, verbose);
        assert(cpk_weight_cmp(SQ->WS, WSM) == 0);
        if (verbose) { fprintf(stderr, "\n  ------------------------\n"); }
        if (verbose) { TRACE_SOL("  Solu��o �tima bruta", SQ->nS, SQ->WS);  }

        /* Faz busca local com max 1 remo��o, para garantir: */
        if (verbose) { fprintf(stderr, "\n  ------------------------\n"); }
        if (verbose) { fprintf(stderr, "  Re-otimiza��o local:");  }
        cpk_local_opt(SQ, 1, &nSM, SM, &WSM, verbose);
        if (verbose) { fprintf(stderr, "\n");  }

        /* Copia melhor solu��o para solu��o corrente: */
        cpk_mis_copy_indep(SQ, nSM, SM, verbose);
        assert(cpk_weight_cmp(SQ->WS, WSM) == 0);
        if (verbose) { TRACE_SOL("  Solu��o �tima otimizada", SQ->nS, SQ->WS);  }

        /* Limpeza: */
        cpk_mis_state_free(SQ);
      }

    if (VALIDATE) { cpk_check_indep_set(G, nSM, SM, NULL); }

    if (verbose)
      { TRACE_SOL("Solu��o GRASP final", nSM, WSM);
        cpk_print_vertex_set(stderr, "", nSM, SM, NULL, 30);
      }

    /* prepara a saida da rotina */
    uint32_vec_t J = uint32_vec_new(nSM); 
    for (int32_t i = 0; i < nSM; i++) { J.e[i] = SM[i]; }

    if (verbose) { cpk_trace_exit(); }
    return J;
  }

ui2_vec_t cpk_grasp_get_seed_pairs(cpk_graph_t *G, double W[], uint32_t maxPairs, bool_t verbose) 
  {
    if (verbose) { cpk_trace_entry(); }
    uint32_t nV = G->nV;
    uint32_t totPares = (nV-1)*nV/2;
    if (verbose) { fprintf(stderr, "  Grafo com %d v�rtices (%d pares)\n", nV, totPares); }

    /* Obt�m lista de v�rtices {X[0..nV-1]}} em ordem decrescente de escore: */
    uint32_t X[nV];
    for(uint32_t v = 0; v < nV; v++) { X[v] = v; }
    
    /* Ordena lista por escore decrescente: */
    auto int32_t compara_vertices(const void *a, const void *b);
    qsort(X, nV, sizeof(uint32_t), &compara_vertices);

    int32_t compara_vertices(const void *a, const void *b)
      { uint32_t u = *((uint32_t*)a); double eu = W[u]/(G->deg[u]+1);
        uint32_t v = *((uint32_t*)b); double ev = W[v]/(G->deg[v]+1);
        return -dblcmp(eu,ev); /* Como s�o {doubles} n�o podemos retornar {eu-ev}: */
      }

    /* Quantos pares vale a pena olhar no m�ximo? (n^2 � demais!) */
    uint32_t maxPairsExm = (uint32_t)imin(10*nV, totPares);
    if (verbose) { fprintf(stderr, "  Coletando %d pares-semente\n", maxPairsExm); }
    cpk_seed_pair_t *uvs;
    uvs = (cpk_seed_pair_t *)notnull(malloc(maxPairsExm*sizeof(cpk_seed_pair_t)), "no mem");
    uint32_t nPairs = 0;
    { bool_t temp[nV]; /* Tempor�rio para c�lculo dos escores. */
      for (int32_t v = 0; v < nV; v++) { temp[v] = FALSE; }
      /* Enumera pares : */
      for(int32_t s = 0; (s < nV) && (nPairs < maxPairsExm); s++)
        { for(int32_t j = 0; (j <= s) && (nPairs < maxPairsExm); j++)
            { uint32_t u = X[j];
              uint32_t v = X[s - j];
              double escore = cpk_grasp_pair_score(G, W, u, v, temp);
              if (escore > -INF)
                { uvs[nPairs] = (cpk_seed_pair_t){ u, v, escore};
                  nPairs++;
                }
            }
        }
      if (verbose) { fprintf(stderr, "  Encontrou %d pares-semente\n", nPairs); }
       
      /* Ordena o vetor {uvs[0..nPairs-1]} em ordem decrescente de escore: */
      auto int32_t compara_escores(const void *a, const void *b);
      qsort(uvs, nPairs, sizeof(cpk_seed_pair_t), &compara_escores);

      int32_t compara_escores(const void *a, const void *b)
        { double ea = ((cpk_seed_pair_t *)a)->escore;
          double eb = ((cpk_seed_pair_t *)b)->escore;
          /* Como escores s�o {double}s, n�o podemos simplesmente subtrair: */
          if (ea > eb)
            { return -1; }
          else if (ea < eb)
            { return +1; }
          else
            { return 0; }
        }
    }

    /* Devolve os melhores {maxPairs} pares: */
    nPairs = (uint32_t)imin(nPairs, maxPairs);
    ui2_vec_t par = ui2_vec_new(nPairs);
    { for (int32_t i = 0; i < nPairs; i++)
        { cpk_seed_pair_t *pi = &(uvs[i]);
          par.e[i] = (ui2_t){{ pi->v1, pi->v2 }};
        }
    }
    
    /* Limpeza: */
    free(uvs);
    
    if (verbose) { cpk_trace_exit(); }
    return par;
  }

double cpk_grasp_pair_score
  ( cpk_graph_t *G, 
    double W[], 
    uint32_t u, 
    uint32_t v, 
    bool_t temp[]
  )
  {
    /* Paran�ia: */
    if (u == v) { return -INF; }
    
    /* Pega graus e vizinhan�as de {u} e {v}: */
    uint32_t du = G->deg[u]; uint32_t *nbu = &(G->nbr[G->fnb[u]]);
    uint32_t dv = G->deg[v]; uint32_t *nbv = &(G->nbr[G->fnb[v]]);
    
    /* Sup�e que {temp} est� toda {FALSE} incialmente. */
    /* Liga {temp[w]} sse {w} � vizinho de {u}: */
    { int32_t j;  for (j = 0; j < du; j++) { assert(! temp[nbu[j]]); temp[nbu[j]] = TRUE; } }
    
    /* Calcula escore: */
    double escore;
    if (temp[v])
      { /* {v} � vizinho de {u}. */
        escore = -INF;
      }
    else
      { /* Conta vizinhos de {v} que s�o vizinhos de {u}: */
        uint32_t duv = 0;
        { int32_t j; for (j = 0; j < dv; j++) { if (temp[nbv[j]]) { duv++; } } }
        if (duv > 0)
          { /* Calcula tamanho {tamViz} da uni�o das vizinhan�as: */
            uint32_t tamViz = du + dv + 2 - duv; /* Come�a contando {u}, {v}, e os vizinhos de {u}. */
            /* Escore � raz�o peso/consumo: */
            escore = (W[u]+W[v]) / ((double)tamViz);
          }
        else
          { /* Dist�ncia de {u} a {v} � maior que 2: */
            escore = -INF;
          }
      }
      
    /* Limpa {temp} para futuras chamadas: */
    for (int32_t j = 0; j < du; j++) { temp[nbu[j]] = FALSE; }
    
    return escore; 
  }
