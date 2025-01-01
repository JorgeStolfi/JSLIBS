#ifndef cpk_lopt_H
#define cpk_lopt_H

/* Otimiza��o local gulosa de um conjunto independente maximal. */
/* Last edited on 2024-12-31 13:19:19 by stolfi */ 

#include <stdint.h>
#include <bool.h>

#include <cpk_basic.h>
#include <cpk_mis.h>
#include <cpk_graph.h>

void cpk_local_opt
  ( cpk_mis_state_t *SQ, /* Solu��o corrente {S} e estruturas auxiliares. */
    uint32_t maxRem,        /* N�mero m�ximo de v�rtice a trocar em cada tentativa. */
    uint32_t *nSM,          /* Tamanho da melhor solu��o j� vista. */
    uint32_t SM[],          /* {SM[0..nSM-1]} � a melhor solu��o j� vista. */
    double *WSM,       /* Peso total da melhor solu��o {SM}. */
    bool_t verbose     /* TRUE para imprimir mensagens de diagn�stico. */
  );
  /* Sup�e que o conjunto {S} armazenado em {SQ} � um conjunto
    independente maximal. Examina alguns outros conjuntos maximais
    {S'} derivados de {S} pela retirada de um certo n�mero {nRem}
    de v�rtices e re-expans�o gulosa para conjunto maximal. 
    
    As vari�veis {SM[0..nSM-1]} e {WSM} devem conter respectivamente
    os v�rtices e o peso total da melhor solu��o j� examinada at� o
    momento. Elas s�o atualizadas sempre que o procedimento encontra
    uma solu��o {S} melhor que {SM}.
    
    O procedimento faz v�rias passadas sobre o grafo, come�ando com
    {nRem=1}, depois com {nRem=2}, etc; e voltando para {nRem=1} cada
    vez que consegue uma melhora. O procedimento retorna quando uma
    passada sobre todos os v�rtices e todos os {nRem} de 1 a {maxRem}
    n�o consegue melhorar a solu��o. */

bool_t cpk_lopt_same_set(uint32_t nX, uint32_t X[], uint32_t nY, uint32_t Y[]);
  /* Verdadeiro se e somente se o conjunto {X[0..nX-1]} � igual a {Y[0..nY-1]},
    ignorando a ordem. �til para verificar se v�rtices recolocados
    na solu��o s�o os que foram tirados. Tempo: O(nX*nY). */

void cpk_lopt_update_best_solution
  ( cpk_mis_state_t *SQ, /* Solu��o corrente {S} e estruturas auxiliares. */
    uint32_t *nSM,          /* Tamanho da melhor solu��o j� vista. */
    uint32_t SM[],          /* {SM[0..nSM-1]} � a melhor solu��o j� vista. */
    double *WSM,       /* Peso total da melhor solu��o {SM}. */
    bool_t verbose     /* TRUE para imprimir mensagens de diagn�stico. */
  );
  /* Se o peso {WS} da solu��o corrente {S} guardada em {SQ} � maior
    que {WSM}, copia {nS,S} para {nSM,SM}, e atualiza {WSM}. */
    
void cpk_lopt_restore_solution
  ( cpk_mis_state_t *SQ,  /* Solu��o corrente {S} e estruturas auxiliares. */
    uint32_t nDel,        /* N�mero de v�rtices a retirar. */
    uint32_t nAdd,        /* N�mero de v�rtices a recolocar. */
    uint32_t Add[],       /* {Add[0..nAdd-1]} s�o os v�rtices a recolocar. */
    bool_t verbose        /* TRUE para imprimir mensagens de diagn�stico. */
  );
  /* Modifica a solu��o {S} retirando dela os �ltimos {nDel} v�rtices,
    e acrescentando os v�rtices {Add[0..nAdd-1]}, que ficar�o no final
    de {S}. Atualiza o peso total {WS} e a fila {Q} adequadamente. */

#endif
