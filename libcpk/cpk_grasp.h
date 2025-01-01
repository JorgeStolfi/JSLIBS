#ifndef cpk_grasp_H
#define cpk_grasp_H

/* Definições gerais para solução de MIS com GRASP. */
/* Last edited on 2024-12-31 13:23:31 by stolfi */ 

#include <stdint.h>
#include <vec.h>
#include <r2.h>
#include <bool.h>

#include <cpk_graph.h>
#include <cpk_basic.h>

/* 4. Constantes usadas no GRASP */
#define cpk_grasp_REL_WT_NOISE 0.10 
  /* Perturbação relativa máxima no peso para randomização da construção. */

uint32_vec_t cpk_grasp_find_indep_set
  ( cpk_graph_t *G,   /* The incompatibilty graph. */
    double *W,        /* Weight of each vertex. */
    double maxClock,  /* Deadline for result. */
    uint32_t seed,         /* Seed for random choices. */
    bool_t verbose    /* TRUE prints debugging diagnostics. */
  );
  /* Computa um conjunto independente (IS) de vértices de {G} com peso
    máximo (espera-se), usando uma adaptação do algoritmo GRASP
    proposto por Resende (``A Greedy Adaptative Search Procedure for
    the Maximum Independent Set'', Operations Research, ...).
    
    O procedimento repetidamente constrói uma solução maximal {S},
    (/fase de construção/), e aplica a cada uma delas um procedimento
    de otimização local (/fase de busca/). Este processo termina
    quando o relógio de CPU {now()} passar de {maxClock}.
    
    A fase de construção inicializa {S} com dois vértices "bons" e
    acrescenta mais vértices, um a um, até {S} ficar maximal. Os
    vértices são escolhidos de forma gulosa (pelo critério da razão
    benefício/custo), mas com uma certa perturbação aleatória.
    
    A fase de busca repetidamente retira até {maxRem} vértices da
    solução {S}, e tenta cobrir o "buraco" assim criado, de maneira
    gulosa ou exaustiva. */

#endif
