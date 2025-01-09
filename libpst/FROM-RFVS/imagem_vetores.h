/* Imagens com elementos vetoriais. */
/* Last edited on 2010-03-17 17:35:57 by stolfi */ 

#ifndef imagem_vetores_H
#define imagem_vetores_H

typedef struct vetor_t { float x, y, z; } vetor_t;

typedef struct imagem_vetores_t
  { int tx; /* Numero de colunas. */
    int ty; /* Numero de linhas. */
    vetor_t** pixel;  /* {pixel[y][x]} está na linha {y} e coluna {x}. */
  } imagem_vetores_t;
  
imagem_vetores_t *cria_imagem_vetores(int tx, int ty);
  /* Aloca uma imagem de normais com {tx} colunas e {ty} linhas.
    os pixels não são inicializados. */

imagem_vetores_t *le_imagem_vetores(char* nome);
  /* Le uma imagem de vetores normais do arquivo "{nome}".
    As duas primeiras linhas devem ser "tx = {tx}" e 
    "ty = {ty}" onde {tx} é o número de colunas e {ty}
    é o número de linhas. A seguir devem ocorrer {tx*ty}
    linhas contendo os campos "{x} {y} {ux} {uy} {uz}",
    onde {(ux,uy,uz)} é o vetor normal no pixel da 
    coluna {x} e linha {y}, linha por linha. */

void escreve_imagem_vetores(char *nome, imagem_vetores_t *I);
  /* Escreve a imagem {I} no arquivo de nome "{arq}",
    num formato compatível com {le_imagem_vetores}. */

#endif

