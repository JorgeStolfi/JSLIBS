/* Veja imagem_vetores.h */
/* Last edited on 2010-03-17 17:35:45 by stolfi */ 

#include <imagem_vetores.h>

#include <assert.h>
#include <stdio.h>
#include <stdlib.h>
#include <affirm.h>

imagem_vetores_t *cria_imagem_vetores(int tx, int ty)
  {
    imagem_vetores_t *IM = (imagem_vetores_t*)malloc(sizeof(imagem_vetores_t));
    int y;
    
    IM->pixel = (vetor_t**)malloc(sizeof(vetor_t*)*ty);
    for(y = 0; y < ty; y++)
      { IM->pixel[y] = (vetor_t*)malloc(sizeof(vetor_t)*tx); }
    IM->tx = tx;
    IM->ty = ty;
    return IM;
  }

void escreve_imagem_vetores(char *nome, imagem_vetores_t *I)
  {
    FILE* arq = fopen(nome, "wt");
    assert(arq != NULL); 
    fprintf(arq, "tx = %d\n", I->tx);
    fprintf(arq, "ty = %d\n", I->ty);
    int x, y;
    for(y = 0; y < I->ty; y++)
      { if (y > 0) { fprintf(arq, "\n"); }
        for(x = 0; x < I->tx; x++)
          { vetor_t *px = &(I->pixel[y][x]);
            fprintf(arq, "%d %d", x, y);
            fprintf(arq, " %f %f %f \n", px->x, px->y, px->z);
          }
      }
    fclose(arq);
  }

imagem_vetores_t *le_imagem_vetores(char *nome)
  {
    FILE* arq = fopen(nome, "rt");
    assert(arq != NULL);
    int tx, ty;
    demand(fscanf(arq, "tx = %d\n", &tx) == 1,"Arquivo IVET invalido - erro TX");
    demand(fscanf(arq, "ty = %d\n", &ty) == 1,"Arquivo IVET invalido - erro TY");
    imagem_vetores_t *IM = cria_imagem_vetores(tx, ty);
    long int N = tx*ty;
    long int i = 0;
    for (i = 0; i < N; i++)
      {
        /* L�posi�o e normal de um pixel: */
        int x, y;
        float ux,uy,uz;
        int nread = fscanf(arq,"%d %d %f %f %f", &x, &y, &ux, &uy, &uz);
        /* Testa EOF normal: */
	//printf("Lendo pixel %d %d\n",x,y);
        if (nread == 0) { break; }
        /* Testa EOF inesperado: */
        assert(nread == 5); 
        
        /* Testa seq�cia dos pixels: */
        assert(x == i % tx);
        assert(y == i / tx);
        
        /* Guarda pixel na imagem: */
        IM->pixel[y][x].x = ux;
        IM->pixel[y][x].y = uy;
        IM->pixel[y][x].z = uz;
      }
    fclose(arq);
    return IM;
  }

