#ifndef pst_signature_cube_H
#define pst_signature_cube_H

/* pst_signature_cube.h -- bucket grid search of light signatures. */
/* Last edited on 2006-04-03 15:27:11 by stolfi */

#include <pst_basic.h>

typedef struct ltn_cube_t /* A bucket grid for light signatures. */
  {
    int NS;           /* Number of entries in each light signature. */
  } ltn_cube_t;

typedef struct Lista lista;

/*
ltn_cube_t* gera_cubo(int resolucao);
void processa_tabela(ltn_table_t* tab, ltn_cube_t* cubo);
long int busca_cubo(ltn_table_t* tab,pixel* cores,int canal,ltn_cube_t* cubo);
*/

#endif
