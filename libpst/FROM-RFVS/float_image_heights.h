/* Last edited on 2010-03-17 17:40:30 by stolfi */

#ifndef float_image_heights_H
#define float_image_heights_H

#define _GNU_SOURCE
#include <sistema.h>
#include <float_image.h>
#include <r3.h>
/* 
  IN comments, one assumes que as linhas s� numeradas DE BAIXO PARA
  CIMA e as colunas s� numeradas da esquerda para a direita.
  Portanto o pixel {[0][0]} �o canto INFERIOR esquerdo da imagem.
  */

//imagem_valores_t *calcula_derivadas(imagem_vetores_t *IN, char eixo, imagem_valores_t *IW);
float_image_t *calcula_derivadas(float_image_t *IN, char eixo, float_image_t* IW);
  /* Dada a imagem de vetores normais {IN}, devolve uma imagem de
    valores {ID}, com mesmo tamanho, tal que {ID[y][x]} �a derivada
    da altura no centro do pixel {[y][x]}, isto � no ponto
    {(x+0.5,y+0.5)}. A derivada �{dz/dx} ou {dz/dy}, conforme {eixo}
    for 'x' ou 'y', respectivamente; ela �calculada a partir da
    normal {IN[y][x]} no mesmo ponto. 

    Se a normal for {0,0,0} a derivada é NAN

    Se a normal for inválida (horizontal ou apontando para baixo), 
    supõe que a superficie é horizontal, mas atribui peso 0 ao pixel. */

/*imagem_valores_t *calcula_alturas_recursiva
  ( imagem_valores_t *IDX, 
    imagem_valores_t *IDY, 
    imagem_valores_t *IW,
    int prof,
    char *prefDebug
  );*/
float_image_t *calcula_alturas_recursiva
  ( float_image_t *IDX, 
    float_image_t *IDY, 
    float_image_t *IW,
    int prof,
    char *prefDebug
  );

  /* Calcula uma imagem de alturas {IZ} dadas imagens de derivadas {IDX}
    (horizontal) e {IDY} (vertical), a imagem de pesos {IW}. A imagem {IZ} deve
    ter�uma coluna e uma linha a mais do que {IDX,IDY,IW} (que devem ter o mesmo tamanho).

    A imagem �calculada resolvendo-se um sistema linear (vide abaixo)
    pelo m�odo de Gauss-Jordan. Este procedimento usa a t�nica
    (recursiva) de multiplas escalas para acelerar a converg�cia do
    mesmo. O par�etro {prof} indica a profundidade da recurs�. Se
    {prefDebug} n� �NULL, escreve as imagens de derivadas e alturas
    usadas em cada escala, bem como o sistema linear correspondente. */


void interpola_derivadas_recursivo(float_image_t *IDX, 
    float_image_t *IDY, 
    float_image_t *IW,
    int prof,
    char *prefDebug
  );


//imagem_valores_t *reduz_derivadas(imagem_valores_t *ID);
float_image_t *reduz_derivadas(float_image_t *ID);
  /* Dada uma imagem de derivadas {ID} de uma imagem de alturas {IZ},
    retorna outra imagem {JD} com as derivadas para outra imagem {JZ}
    de alturas, com metade do tamanho e com as alturas divididas pela metade. 
    
    Se o tamanho de {ID} for �par, a ltima linha e/ou a ltima coluna 
    s� duplicadas implicitamente antes da redu�o. */


//imagem_valores_t *reduz_pesos(imagem_valores_t *ID);
float_image_t *reduz_pesos(float_image_t *ID);
/* Dada uma imagem de pesos {ID}
    retorna outra imagem de pesos{JD}, com metade do tamanho.
    
    Se o tamanho de {ID} for �par, a ltima linha e/ou a ltima coluna 
    s� duplicadas implicitamente antes da redu�o. */

void escreve_sistema(sistema_t *S, char *nome, int tx, int ty);
  /* Escreve o sistema {S} no arquivo com o {nome} dado,
    explicitando os �dices dos pixels a que se refere cada equa�o 
    e cada inc�nita. */

//void preenche_sistema(imagem_valores_t *IDX, imagem_valores_t *IDY, imagem_valores_t *IW ,sistema_t* S);
void preenche_sistema(float_image_t *IDX, float_image_t *IDY, float_image_t *IW ,sistema_t* S);
  /*
    Dadas duas matrizes {IDX} e {IDY} com as derivadas experimentais de
    uma imagem de alturas,e uma imagem de pesos IW preenche os coeficientes do sistema linear
    {S} cuja solu�o s� as tais alturas.
    
    As imagens {IDX},{IDY}e {IW} devem ter o mesmo tamanho {tx} colunas e
    {ty} linhas. O valor de {IDX[y][x]} deve ser a derivada horizontal
    {Dz/Dx} da altura, medida no centro do pixel da linha {y} e coluna
    {x}. O valor de {IDY[y][x]} deve ser a derivada vertical {Dz/Dy},
    no mesmo ponto.
    
    O sistema tem a forma {M Z = B} onde {M} �uma matriz conhecida,
    {Z} �um vetor de inc�nitas, e {B} um vetor conhecido. As
    inc�nitas {Z[j]} s� as alturas nos v�tices (cantos dos pixels)
    de {IDX} e {IDY}; e portanto o tamanho do sistema �{N =
    (tx+1)+(ty+1)}.
    
    Mais precisamente, a inc�nita {Z[j]}, para todo {j = 0..N-1}, �a
    altura no ponto {(x,y)} tal que {j = x + (tx+1)*y}, com {x} em
    {0..tx} e {y} em {0..ty}. Este ponto �o canto inferior esquerdo
    do pixel {IDX[y][x]}. */
      

void preenche_sistema_interpolacao(float_image_t* ID,float_image_t*IW,sistema_t* S);

//void extrai_valores_da_imagem(imagem_valores_t *IZ, double *Z);
void extrai_valores_da_imagem(float_image_t *IZ, double *Z);            
//void coloca_valores_na_imagem(double *Z, imagem_valores_t *IZ);
void coloca_valores_na_imagem(double *Z, float_image_t *IZ);
  /* Copia os valores das inc�nitas {Z[i]} de ou para os 
    elementos da imagem de alturas {IZ->elem[y][x]}. */

//imagem_valores_t *expande_imagem(imagem_valores_t *JZ, int tx, int ty);
float_image_t *expande_alturas(float_image_t *JZ, int tx, int ty);
  /* Dada uma imagem de alturas {JZ}, retorna outra imagem {IZ} com
    {tx} colunas e {ty} linhas.
    
    A imagem original deve ter tamanho {tx/2 + 1} por {ty/2 + 1},
    arredondado para baixo. Na nova imagem, pixels de �dices pares
    s� copiados do original, e pixels �pares s� interpolados. */

float_image_t *expande_derivadas(float_image_t *JD, int tx, int ty);
  /* Dada uma imagem de alturas {JD}, retorna outra imagem {IZ} com
    {tx} colunas e {ty} linhas.
    
    A imagem original deve ter tamanho {(tx+1)/2} por {(ty+1)/2},
    arredondado para baixo. Na nova imagem, cada pixel é duplicado por 4 pixels. */


//void ajusta_pesos(imagem_valores_t* IW, double expPeso);
void ajusta_pesos_e_normais(float_image_t* IW, float_image_t* IN, double expPeso);
/*
	Eleva imagem de pesos a um expoente expPeso
*/

void elimina_termo_linear(float_image_t* IZ, float_image_t* IDX, float_image_t*IDY, float_image_t*IW);
/*
	Soma ao mapa de alturas IZ, um termo linear {Cx*x +Cy*y} de modo a
	minimizar a diferença entre as normais dadas {IDX} e {IDY} e as normais calculadas por
	diferença finita do mapa de alturas {IZ}
*/
#endif
