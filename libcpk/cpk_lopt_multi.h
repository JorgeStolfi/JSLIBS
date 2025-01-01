#ifndef cpk_lopt_multi_H
#define cpk_lopt_multi_H

/* Otimização local de conjunto indep maximal - remocoes multiplas. */
/* Last edited on 2024-12-31 15:42:03 by stolfi */ 

#include <stdint.h>

#include <bool.h>

#include <cpk_basic.h>
#include <cpk_mis.h>
#include <cpk_graph.h>

void cpk_local_opt_multiple
  ( cpk_mis_state_t *SQ, /* Solução corrente {S} e estruturas auxiliares. */
    uint32_t maxRem,        /* Número máximo de vértice a trocar em cada tentativa. */
    uint32_t *nSM,          /* Tamanho da melhor solução já vista. */
    uint32_t SM[],          /* {SM[0..nSM-1]} é a melhor solução já vista. */
    double *WSM,       /* Peso total da melhor solução {SM}. */
    bool_t verbose     /* TRUE para imprimir mensagens de diagnóstico. */
  );
  /* Supõe que o conjunto {S} armazenado em {SQ} é um conjunto
    independente maximal. Examina alguns outros conjuntos maximais
    {S'} derivados de {S} pela retirada de até {maxRem} vértices e
    re-espansão gulosa para conjunto maximal; faz no máximo {maxTry}
    dessas tentativas. 
    
    As variáveis {SM[0..nSM-1]} e {WSM} devem conter respectivamente 
    os vértices e o peso total da melhor solução já examinada até o momento.
    Elas são atualizadas sempre que o procedimento encontra uma 
    solução {S} melhor que {SM}.
    
    O procedimento faz várias passadas sobre o grafo, tentando todos
    os subconjuntos de {S} de tamanho {nRem = 1}, que encontra, depois
    {nRem=2}, {nRem = 3}, etc.; retorna quando uma passada completa
    sobre todos os vértice e todos os tamanhos {nRem=1..maxRem} não
    consegue melhorar a solução.  */

#endif
