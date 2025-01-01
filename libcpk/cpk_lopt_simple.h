#ifndef cpk_lopt_simple_H
#define cpk_lopt_simple_H

/* Otimiza��o local de conjunto indep maximal - remocoes multiplas. */
/* Last edited on 2024-12-31 15:40:59 by stolfi */ 

#include <stdint.h>

#include <bool.h>

#include <cpk_basic.h>
#include <cpk_mis.h>

void cpk_local_opt_simple
  ( cpk_mis_state_t *SQ, /* Solu��o corrente {S} e estruturas auxiliares. */
    uint32_t *nSM,          /* Tamanho da melhor solu��o j� vista. */
    uint32_t SM[],          /* {SM[0..nSM-1]} � a melhor solu��o j� vista. */
    double *WSM,       /* Peso total da melhor solu��o {SM}. */
    bool_t verbose     /* TRUE para imprimir mensagens de diagn�stico. */
  );
  /* Sup�e que o conjunto {S} armazenado em {SQ} � um conjunto
    independente maximal. Examina alguns outros conjuntos maximais
    {S'} derivados de {S} pela retirada de um �nico v�rtice e
    re-expans�o gulosa para conjunto maximal. 
    
    As vari�veis {SM[0..nSM-1]} e {WSM} devem conter respectivamente 
    os v�rtices e o peso total da melhor solu��o j� examinada at� o momento.
    Elas s�o atualizadas sempre que o procedimento encontra uma 
    solu��o {S} melhor que {SM}.

    O procedimento faz v�rias passadas sobre os v�rtices de {S}, e
    retorna quando uma passada sobre todos os v�rtices de {S} n�o
    consegue melhorar a solu��o. */

#endif
