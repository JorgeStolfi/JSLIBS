#ifndef cpk_lopt_multi_H
#define cpk_lopt_multi_H

/* Otimiza��o local de conjunto indep maximal - remocoes multiplas. */
/* Last edited on 2024-12-31 15:42:03 by stolfi */ 

#include <stdint.h>

#include <bool.h>

#include <cpk_basic.h>
#include <cpk_mis.h>
#include <cpk_graph.h>

void cpk_local_opt_multiple
  ( cpk_mis_state_t *SQ, /* Solu��o corrente {S} e estruturas auxiliares. */
    uint32_t maxRem,        /* N�mero m�ximo de v�rtice a trocar em cada tentativa. */
    uint32_t *nSM,          /* Tamanho da melhor solu��o j� vista. */
    uint32_t SM[],          /* {SM[0..nSM-1]} � a melhor solu��o j� vista. */
    double *WSM,       /* Peso total da melhor solu��o {SM}. */
    bool_t verbose     /* TRUE para imprimir mensagens de diagn�stico. */
  );
  /* Sup�e que o conjunto {S} armazenado em {SQ} � um conjunto
    independente maximal. Examina alguns outros conjuntos maximais
    {S'} derivados de {S} pela retirada de at� {maxRem} v�rtices e
    re-espans�o gulosa para conjunto maximal; faz no m�ximo {maxTry}
    dessas tentativas. 
    
    As vari�veis {SM[0..nSM-1]} e {WSM} devem conter respectivamente 
    os v�rtices e o peso total da melhor solu��o j� examinada at� o momento.
    Elas s�o atualizadas sempre que o procedimento encontra uma 
    solu��o {S} melhor que {SM}.
    
    O procedimento faz v�rias passadas sobre o grafo, tentando todos
    os subconjuntos de {S} de tamanho {nRem = 1}, que encontra, depois
    {nRem=2}, {nRem = 3}, etc.; retorna quando uma passada completa
    sobre todos os v�rtice e todos os tamanhos {nRem=1..maxRem} n�o
    consegue melhorar a solu��o.  */

#endif
