#ifndef cpk_lopt_simple_H
#define cpk_lopt_simple_H

/* Otimização local de conjunto indep maximal - remocoes multiplas. */
/* Last edited on 2024-12-31 15:40:59 by stolfi */ 

#include <stdint.h>

#include <bool.h>

#include <cpk_basic.h>
#include <cpk_mis.h>

void cpk_local_opt_simple
  ( cpk_mis_state_t *SQ, /* Solução corrente {S} e estruturas auxiliares. */
    uint32_t *nSM,          /* Tamanho da melhor solução já vista. */
    uint32_t SM[],          /* {SM[0..nSM-1]} é a melhor solução já vista. */
    double *WSM,       /* Peso total da melhor solução {SM}. */
    bool_t verbose     /* TRUE para imprimir mensagens de diagnóstico. */
  );
  /* Supõe que o conjunto {S} armazenado em {SQ} é um conjunto
    independente maximal. Examina alguns outros conjuntos maximais
    {S'} derivados de {S} pela retirada de um único vértice e
    re-expansão gulosa para conjunto maximal. 
    
    As variáveis {SM[0..nSM-1]} e {WSM} devem conter respectivamente 
    os vértices e o peso total da melhor solução já examinada até o momento.
    Elas são atualizadas sempre que o procedimento encontra uma 
    solução {S} melhor que {SM}.

    O procedimento faz várias passadas sobre os vértices de {S}, e
    retorna quando uma passada sobre todos os vértices de {S} não
    consegue melhorar a solução. */

#endif
