/* Veja imagem_valores.h */
/* Last edited on 2010-03-17 17:35:31 by stolfi */ 

#include <imagem_valores.h>

#include <assert.h>
#include <stdio.h>
#include <stdlib.h>

imagem_valores_t *cria_imagem_valores(int tx, int ty)
  {
    imagem_valores_t *IM = (imagem_valores_t*)malloc(sizeof(imagem_valores_t));
    int y;

    IM-> pixel = (float**)malloc(sizeof(float*)*ty);
    for(y = 0; y < ty; y++)
      { IM->pixel[y] = (float*)malloc(sizeof(float)*tx); }
    IM->tx = tx;
    IM->ty = ty;
    return IM;
  }

void escreve_imagem_valores(char *nome, imagem_valores_t *I)
  {
    FILE* arq = fopen(nome, "wt");
    assert(arq != NULL); 
    fprintf(arq, "tx = %d\n", I->tx);
    fprintf(arq, "ty = %d\n", I->ty);
    int x, y;
    for(y = 0; y < I->ty; y++)
      { if (y > 0) { fprintf(arq, "\n"); }
        for(x = 0; x < I->tx; x++)
          { float *px = &(I->pixel[y][x]);
            fprintf(arq, "%d %d %f \n", x, y, (*px));
          }
      }
    fclose(arq);
  }

imagem_valores_t *le_imagem_valores(char *nome)
  {
    //fprintf(stderr, "NOME: %s \n",nome);
    FILE* arq = fopen(nome, "rt");
    
    assert(arq != NULL);
    int tx, ty;
    fscanf(arq, "tx = %d\n", &tx);
    fscanf(arq, "ty = %d\n", &ty);
    imagem_valores_t *IM = cria_imagem_valores(tx, ty);
    long int N = tx*ty;
    long int i = 0;
    for (i = 0; i < N; i++)
      { /* L�posi�o e valor de um pixel: */
        int x, y;
        float u;
        int nread = fscanf(arq,"%d %d %f", &x, &y, &u);
        /* Testa EOF normal: */
        if (nread == 0) { break; }
        /* Testa EOF inesperado: */
        assert(nread == 3); 
        
        /* Testa seq�cia dos pixels: */
	//fprintf(stderr, "%d %d %d %d\n",x,y,tx,ty);
        assert(x == i % tx);
        assert(y == i / tx);
        
        /* Guarda pixel na imagem: */
        IM->pixel[y][x] = u;
      }
    fclose(arq);
    return IM;
  }
