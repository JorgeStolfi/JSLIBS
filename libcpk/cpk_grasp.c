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
  /* Estrutura usada para encontrar pares de vértices independentes de
    alto peso e baixo consumo. */

/* PROTÓTIPOS */

ui2_vec_t cpk_grasp_get_seed_pairs(cpk_graph_t *G, double W[], uint32_t maxPairs, bool_t verbose);
  /* Devolve uma lista {K} de pares de vértices independentes de {G},
     para serem usados como soluções de partida da fase de construção
     do GRASP. A lista terá no máximo {maxPairs} pares. */

double cpk_grasp_pair_score
  ( cpk_graph_t *G, 
    double W[], 
    uint32_t u, 
    uint32_t v, 
    bool_t temp[]
  );
  /* Calcula o escore de um par de vértices {u,v} considerado como
    conjunto inicial para a fase de construção do GRASP. O escore é
    -INF se {u} e {v} são vizinhos, caso contrário é a soma dos pesos
    de {u} e {v} dividida pelo tamanho da união de suas vizinhanças
    (incluindo {u} e {v}).
    
    O vetor de booleanos {temp} é usado internamente; ele deve ser
    todo FALSE na entrada e será todo FALSE na saída. */

/* IMPLEMENTAÇÕES */

#define VALIDATE TRUE

uint32_vec_t cpk_grasp_find_indep_set
  ( cpk_graph_t *G,
    double W[], 
    double maxClock,  /* Tempo limite para retornar a resposta. */
    uint32_t seed,    /* Semente para escolhas aleatórias. */
    bool_t verbose
  )
  {
    /* As fases de construção e busca sao repetidas para varios pares
    de vértices independentes de partida, ciclicamente. */

    if (verbose) { cpk_trace_entry(); }

    uint32_t nV = G->nV;

    /* Melhor solução encontrada até o momento (inicialmente vazia): */
    uint32_t *SM = talloc(nV, uint32_t);
    uint32_t nSM = 0;     /* A melhor solução é {SM[0..nSM-1]}. */
    double WSM = 0;  /* Peso total de {SM}. */

    /* Inicializa gerador de números aleatórios: */
    srand(seed);
    
    if (nV > 0)
      { /* Decide quantos pares-semente usar. Melhor ficar com poucos e bons, acho: */
        uint32_t maxSeedPairs = (uint32_t)ceil(sqrt((double)nV));
        /* Constrói um conjunto de no máximo {maxSeedPairs} pares de vértices iniciais: */
        ui2_vec_t seedPair = cpk_grasp_get_seed_pairs(G, W, maxSeedPairs, verbose);
        if (verbose) { cpk_trace(); }
        
        /* Representação do conjunto independente {S} e admissível {Q}: */
        cpk_mis_state_t *SQ = cpk_mis_state_new(G, W);

        /* Roda várias tentativas GRASP: */
        uint32_t nIter = 0;
        do
          { if (verbose)  { fprintf(stderr, "Iteração %d\n", nIter); }

            /* FASE DE CONSTRUÇÃO */

            /* Esvazia {S}, coloca todos os nós em {Q}, ordenados com base em WR: */
            cpk_mis_init(SQ, verbose);

            if (nIter == 0)
              { /* Na iteração 0, começa com a pilha vazia: */
              }
            else
              { /* Incializa {S} com uma semente variada: */
                if (verbose) { fprintf(stderr, "\n  ------------------------\n"); }
                if (seedPair.ne > 0)
                  { /* Inicializa {S} com o próximo par-semente: */
                    if (verbose) { fprintf(stderr, "  Inserindo um par inicial"); }
                    { uint32_t t = nIter % seedPair.ne;
                      uint32_t u = (uint32_t)seedPair.e[t].c[0];
                      uint32_t v = (uint32_t)seedPair.e[t].c[1];
                      cpk_mis_add_indep(SQ, u, verbose);
                      cpk_mis_add_indep(SQ, v, verbose);
                    }
                  }
                else 
                  { /* Inicializa {S} com o próximo vértice: */
                    uint32_t u = nIter % nV;
                    cpk_mis_add_indep(SQ, u, verbose);
                  }
                if (verbose) { fprintf(stderr, "\n"); }
                if (verbose) { TRACE_SOL("  Solução semente", SQ->nS, SQ->WS);  }
              }
            /* Expande a solução aleatório-gulosamente até ficar maximal: */
            if (verbose) { fprintf(stderr, "\n  ------------------------\n"); }
            if (verbose) { fprintf(stderr, "  Expansão gulosa randomizada:"); }
            cpk_greedy_maximize_random(SQ, verbose);
            if (verbose) { fprintf(stderr, "\n");  }
            if (verbose) { TRACE_SOL("  Solução gulosa randomizada", SQ->nS, SQ->WS);  }

            if (VALIDATE) { cpk_check_indep_set(G, SQ->nS, SQ->S, NULL); }

            /* Atualiza a melhor solução: */
            cpk_lopt_update_best_solution(SQ, &nSM, SM, &WSM, verbose);

            /* FASE DE BUSCA */

            /* Até quantos vértices vale a pena tentar relocar por vez? */
            uint32_t maxRem = 3;

            /* Faz otimização local com retirada de múltiplos vértices: */
            assert(pqueue_count(SQ->Q) == 0);
            if (verbose) { fprintf(stderr, "\n  ------------------------\n"); }
            if (verbose) { fprintf(stderr, "  Otimização local:"); }
            cpk_local_opt(SQ, maxRem,  &nSM, SM, &WSM, verbose);
            if (verbose) { fprintf(stderr, "\n");  }
            if (verbose) { TRACE_SOL("  Solução localmente otimizada", SQ->nS, SQ->WS);  }

            if (VALIDATE)  
              { cpk_check_indep_set(G, SQ->nS, SQ->S, NULL);
                cpk_check_indep_set(G, nSM, SM, NULL);
              }

            /* Mais uma iteração: */
            nIter++;
          }
        while(cpk_cpu_time_1() < maxClock);

        /* Copia melhor solução para solução corrente: */
        cpk_mis_copy_indep(SQ, nSM, SM, verbose);
        assert(cpk_weight_cmp(SQ->WS, WSM) == 0);
        if (verbose) { fprintf(stderr, "\n  ------------------------\n"); }
        if (verbose) { TRACE_SOL("  Solução ótima bruta", SQ->nS, SQ->WS);  }

        /* Faz busca local com max 1 remoção, para garantir: */
        if (verbose) { fprintf(stderr, "\n  ------------------------\n"); }
        if (verbose) { fprintf(stderr, "  Re-otimização local:");  }
        cpk_local_opt(SQ, 1, &nSM, SM, &WSM, verbose);
        if (verbose) { fprintf(stderr, "\n");  }

        /* Copia melhor solução para solução corrente: */
        cpk_mis_copy_indep(SQ, nSM, SM, verbose);
        assert(cpk_weight_cmp(SQ->WS, WSM) == 0);
        if (verbose) { TRACE_SOL("  Solução ótima otimizada", SQ->nS, SQ->WS);  }

        /* Limpeza: */
        cpk_mis_state_free(SQ);
      }

    if (VALIDATE) { cpk_check_indep_set(G, nSM, SM, NULL); }

    if (verbose)
      { TRACE_SOL("Solução GRASP final", nSM, WSM);
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
    if (verbose) { fprintf(stderr, "  Grafo com %d vértices (%d pares)\n", nV, totPares); }

    /* Obtém lista de vértices {X[0..nV-1]}} em ordem decrescente de escore: */
    uint32_t X[nV];
    for(uint32_t v = 0; v < nV; v++) { X[v] = v; }
    
    /* Ordena lista por escore decrescente: */
    auto int32_t compara_vertices(const void *a, const void *b);
    qsort(X, nV, sizeof(uint32_t), &compara_vertices);

    int32_t compara_vertices(const void *a, const void *b)
      { uint32_t u = *((uint32_t*)a); double eu = W[u]/(G->deg[u]+1);
        uint32_t v = *((uint32_t*)b); double ev = W[v]/(G->deg[v]+1);
        return -dblcmp(eu,ev); /* Como são {doubles} não podemos retornar {eu-ev}: */
      }

    /* Quantos pares vale a pena olhar no máximo? (n^2 é demais!) */
    uint32_t maxPairsExm = (uint32_t)imin(10*nV, totPares);
    if (verbose) { fprintf(stderr, "  Coletando %d pares-semente\n", maxPairsExm); }
    cpk_seed_pair_t *uvs;
    uvs = (cpk_seed_pair_t *)notnull(malloc(maxPairsExm*sizeof(cpk_seed_pair_t)), "no mem");
    uint32_t nPairs = 0;
    { bool_t temp[nV]; /* Temporário para cálculo dos escores. */
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
          /* Como escores são {double}s, não podemos simplesmente subtrair: */
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
    /* Paranóia: */
    if (u == v) { return -INF; }
    
    /* Pega graus e vizinhanças de {u} e {v}: */
    uint32_t du = G->deg[u]; uint32_t *nbu = &(G->nbr[G->fnb[u]]);
    uint32_t dv = G->deg[v]; uint32_t *nbv = &(G->nbr[G->fnb[v]]);
    
    /* Supõe que {temp} está toda {FALSE} incialmente. */
    /* Liga {temp[w]} sse {w} é vizinho de {u}: */
    { int32_t j;  for (j = 0; j < du; j++) { assert(! temp[nbu[j]]); temp[nbu[j]] = TRUE; } }
    
    /* Calcula escore: */
    double escore;
    if (temp[v])
      { /* {v} é vizinho de {u}. */
        escore = -INF;
      }
    else
      { /* Conta vizinhos de {v} que são vizinhos de {u}: */
        uint32_t duv = 0;
        { int32_t j; for (j = 0; j < dv; j++) { if (temp[nbv[j]]) { duv++; } } }
        if (duv > 0)
          { /* Calcula tamanho {tamViz} da união das vizinhanças: */
            uint32_t tamViz = du + dv + 2 - duv; /* Começa contando {u}, {v}, e os vizinhos de {u}. */
            /* Escore é razão peso/consumo: */
            escore = (W[u]+W[v]) / ((double)tamViz);
          }
        else
          { /* Distância de {u} a {v} é maior que 2: */
            escore = -INF;
          }
      }
      
    /* Limpa {temp} para futuras chamadas: */
    for (int32_t j = 0; j < du; j++) { temp[nbu[j]] = FALSE; }
    
    return escore; 
  }
