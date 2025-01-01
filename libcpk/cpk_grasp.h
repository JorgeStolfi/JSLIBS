#ifndef cpk_grasp_H
#define cpk_grasp_H

/* Defini��es gerais para solu��o de MIS com GRASP. */
/* Last edited on 2024-12-31 13:23:31 by stolfi */ 

#include <stdint.h>
#include <vec.h>
#include <r2.h>
#include <bool.h>

#include <cpk_graph.h>
#include <cpk_basic.h>

/* 4. Constantes usadas no GRASP */
#define cpk_grasp_REL_WT_NOISE 0.10 
  /* Perturba��o relativa m�xima no peso para randomiza��o da constru��o. */

uint32_vec_t cpk_grasp_find_indep_set
  ( cpk_graph_t *G,   /* The incompatibilty graph. */
    double *W,        /* Weight of each vertex. */
    double maxClock,  /* Deadline for result. */
    uint32_t seed,         /* Seed for random choices. */
    bool_t verbose    /* TRUE prints debugging diagnostics. */
  );
  /* Computa um conjunto independente (IS) de v�rtices de {G} com peso
    m�ximo (espera-se), usando uma adapta��o do algoritmo GRASP
    proposto por Resende (``A Greedy Adaptative Search Procedure for
    the Maximum Independent Set'', Operations Research, ...).
    
    O procedimento repetidamente constr�i uma solu��o maximal {S},
    (/fase de constru��o/), e aplica a cada uma delas um procedimento
    de otimiza��o local (/fase de busca/). Este processo termina
    quando o rel�gio de CPU {now()} passar de {maxClock}.
    
    A fase de constru��o inicializa {S} com dois v�rtices "bons" e
    acrescenta mais v�rtices, um a um, at� {S} ficar maximal. Os
    v�rtices s�o escolhidos de forma gulosa (pelo crit�rio da raz�o
    benef�cio/custo), mas com uma certa perturba��o aleat�ria.
    
    A fase de busca repetidamente retira at� {maxRem} v�rtices da
    solu��o {S}, e tenta cobrir o "buraco" assim criado, de maneira
    gulosa ou exaustiva. */

#endif
