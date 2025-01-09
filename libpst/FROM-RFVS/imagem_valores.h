/* Imagens com elementos {float}. */
/* Last edited on 2010-03-17 17:35:38 by stolfi */ 

#ifndef imagem_valores_H
#define imagem_valores_H

typedef struct imagem_valores_t
  { int tx; /* Numero de colunas. */
    int ty; /* Numero de linhas. */
    float** pixel;  /* {pixel[y][x]} est� na linha {y} e coluna {x}. */
  } imagem_valores_t;
  
imagem_valores_t *cria_imagem_valores(int tx, int ty);
  /* Aloca uma imagem de normais com {tx} colunas e {ty} linhas.
    os pixels n�o s�o inicializados. */

imagem_valores_t *le_imagem_valores(char* nome);
  /* Le uma imagem de valores normais do arquivo "{nome}".
    As duas primeiras linhas devem ser "tx = {tx}" e 
    "ty = {ty}" onde {tx} � o n�mero de colunas e {ty}
    � o n�mero de linhas. A seguir devem ocorrer {tx*ty}
    linhas contendo os campos "{x} {y} {ux} {uy} {uz}",
    onde {(ux,uy,uz)} � o valor normal no pixel da 
    coluna {x} e linha {y}, linha por linha. */

void escreve_imagem_valores(char *nome, imagem_valores_t *I);
  /* Escreve a imagem {I} no arquivo de nome "{arq}",
    num formato compat�vel com {le_imagem_valores}. */

#endif

