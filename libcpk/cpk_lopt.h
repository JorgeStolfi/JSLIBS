#ifndef cpk_lopt_H
#define cpk_lopt_H

/* Otimização local gulosa de um conjunto independente maximal. */
/* Last edited on 2024-12-31 13:19:19 by stolfi */ 

#include <stdint.h>
#include <bool.h>

#include <cpk_basic.h>
#include <cpk_mis.h>
#include <cpk_graph.h>

void cpk_local_opt
  ( cpk_mis_state_t *SQ, /* Solução corrente {S} e estruturas auxiliares. */
    uint32_t maxRem,        /* Número máximo de vértice a trocar em cada tentativa. */
    uint32_t *nSM,          /* Tamanho da melhor solução já vista. */
    uint32_t SM[],          /* {SM[0..nSM-1]} é a melhor solução já vista. */
    double *WSM,       /* Peso total da melhor solução {SM}. */
    bool_t verbose     /* TRUE para imprimir mensagens de diagnóstico. */
  );
  /* Supõe que o conjunto {S} armazenado em {SQ} é um conjunto
    independente maximal. Examina alguns outros conjuntos maximais
    {S'} derivados de {S} pela retirada de um certo número {nRem}
    de vértices e re-expansão gulosa para conjunto maximal. 
    
    As variáveis {SM[0..nSM-1]} e {WSM} devem conter respectivamente
    os vértices e o peso total da melhor solução já examinada até o
    momento. Elas são atualizadas sempre que o procedimento encontra
    uma solução {S} melhor que {SM}.
    
    O procedimento faz várias passadas sobre o grafo, começando com
    {nRem=1}, depois com {nRem=2}, etc; e voltando para {nRem=1} cada
    vez que consegue uma melhora. O procedimento retorna quando uma
    passada sobre todos os vértices e todos os {nRem} de 1 a {maxRem}
    não consegue melhorar a solução. */

bool_t cpk_lopt_same_set(uint32_t nX, uint32_t X[], uint32_t nY, uint32_t Y[]);
  /* Verdadeiro se e somente se o conjunto {X[0..nX-1]} é igual a {Y[0..nY-1]},
    ignorando a ordem. Útil para verificar se vértices recolocados
    na solução são os que foram tirados. Tempo: O(nX*nY). */

void cpk_lopt_update_best_solution
  ( cpk_mis_state_t *SQ, /* Solução corrente {S} e estruturas auxiliares. */
    uint32_t *nSM,          /* Tamanho da melhor solução já vista. */
    uint32_t SM[],          /* {SM[0..nSM-1]} é a melhor solução já vista. */
    double *WSM,       /* Peso total da melhor solução {SM}. */
    bool_t verbose     /* TRUE para imprimir mensagens de diagnóstico. */
  );
  /* Se o peso {WS} da solução corrente {S} guardada em {SQ} é maior
    que {WSM}, copia {nS,S} para {nSM,SM}, e atualiza {WSM}. */
    
void cpk_lopt_restore_solution
  ( cpk_mis_state_t *SQ,  /* Solução corrente {S} e estruturas auxiliares. */
    uint32_t nDel,        /* Número de vértices a retirar. */
    uint32_t nAdd,        /* Número de vértices a recolocar. */
    uint32_t Add[],       /* {Add[0..nAdd-1]} são os vértices a recolocar. */
    bool_t verbose        /* TRUE para imprimir mensagens de diagnóstico. */
  );
  /* Modifica a solução {S} retirando dela os últimos {nDel} vértices,
    e acrescentando os vértices {Add[0..nAdd-1]}, que ficarão no final
    de {S}. Atualiza o peso total {WS} e a fila {Q} adequadamente. */

#endif
