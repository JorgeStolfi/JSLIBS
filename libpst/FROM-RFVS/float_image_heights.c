#include <float_image_heights.h>
#include <assert.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>



#define TAM_CORTE 10
  /* N� vale a pena recurs� abaixo deste tamanho (acho). */


void elimina_termo_linear(float_image_t* IZ, float_image_t* IDX, float_image_t*IDY, float_image_t*IW){
// 	int tx = IDX->tx; assert(tx == IDY->tx);assert(tx == IW->tx);
//     	int ty = IDX->ty; assert(ty == IDY->ty);assert(ty == IW->ty);
// 
// 	assert(IZ->tx == tx +1);
// 	assert(IZ->ty == ty +1);

	int tx = IDX->sz[1]; assert(tx == IDY->sz[1]);assert(tx == IW->sz[1]);
    	int ty = IDX->sz[2]; assert(ty == IDY->sz[2]);assert(ty == IW->sz[2]);

	assert(IZ->sz[1] == tx +1);
	assert(IZ->sz[2] == ty +1);

	int x,y;
	double sum_w_dx,sum_w_dy,sum_w;
	sum_w_dx = sum_w_dy = sum_w = 0;
	for(x = 0; x < tx; x++){
		for(y = 0; y < ty; y++){
			/*derivadas dadas da normal*/
// 			double ndx = IDX->pixel[y][x];
// 			double ndy = IDY->pixel[y][x];
			double ndx = float_image_get_sample(IDX,0,x,y);
			double ndy = float_image_get_sample(IDY,0,x,y);
			/*derivadas calculadas das alturas*/
// 			double cdx = (IZ->pixel[y+1][x+1] + IZ->pixel[y][x+1] - IZ->pixel[y+1][x] - IZ->pixel[y][x])/2.0;
			double iz11,iz01,iz10,iz00;
			iz11 = float_image_get_sample(IZ,0,x+1,y+1);
			iz01 = float_image_get_sample(IZ,0,x,y+1);
			iz10 = float_image_get_sample(IZ,0,x+1,y);
			iz00 = float_image_get_sample(IZ,0,x,y);
			double cdx = (iz11 + iz10 - iz01 - iz00)/2.0;
			//double cdy = (IZ->pixel[y+1][x+1] - IZ->pixel[y][x+1] + IZ->pixel[y+1][x] - IZ->pixel[y][x])/2.0;
			double cdy = (iz11 - iz10 + iz01 - iz00)/2.0;
			//double w = IW->pixel[y][x];
			double w = float_image_get_sample(IW,0,x,y);
			sum_w_dx+= w*(ndx - cdx);
			sum_w_dy+= w*(ndy - cdy);
			sum_w+= w;
			
		}
	}
	if(sum_w == 0){
		fprintf(stderr,"AVISO: Peso total do mapa de pesos igual a zero !\n");
		return ;
	}

	double Cx,Cy;
	Cx = sum_w_dx/sum_w;
	Cy = sum_w_dy/sum_w;
	for(x = 0; x < tx+1; x++){
		for(y = 0; y < ty+1; y++){
// 			IZ->pixel[y][x]+= Cx*x + Cy*y;
			double val = float_image_get_sample(IZ,0,x,y) + Cx*x + Cy*y;
			float_image_set_sample(IZ,0,x,y,val);
		}
	}
	

}

void ajusta_pesos_e_normais(float_image_t* IW,float_image_t* IN, double expPeso){
	int i,j;
	for(i = 0; i < IW->sz[1];i++){
		for(j = 0; j < IW->sz[2];j++){
			//IW->pixel[j][i] = pow(IW->pixel[j][i],expPeso);
			double val = float_image_get_sample(IW,0,i,j);
			if( val < 0.0)  val = 0;
			val = pow(val,expPeso);
			if( val < 1.0e-10){
				 val = 1.0e-10;
				 float_image_set_sample(IN,0,i,j,0);
				 float_image_set_sample(IN,1,i,j,0);
				 float_image_set_sample(IN,2,i,j,1.0);
			}
			float_image_set_sample(IW,0,i,j,val);
			
		}
	}
}

float_image_t *calcula_alturas_recursiva
  ( float_image_t *IDX, 
    float_image_t *IDY, 
    float_image_t *IW,
    int prof,
    char *prefDebug
  )
  {
    /* Pega tamanho das imagens e verifica compatibilidade: */
//     int tx = IDX->tx; assert(tx == IDY->tx);assert(tx == IW->tx);
//     int ty = IDX->ty; assert(ty == IDY->ty);assert(ty == IW->ty);

     int tx = IDX->sz[1]; assert(tx == IDY->sz[1]);assert(tx == IW->sz[1]);
     int ty = IDX->sz[2]; assert(ty == IDY->sz[2]);assert(ty == IW->sz[2]);
	

    int ind = 2*prof+2; /* Indenta�o dos coment�ios. */

    fprintf(stderr, "%*sIn�io do n�el %d (IDX,IDY,IW = %d�%d)...\n", ind, "", prof, tx, ty);
    char *nomeArq;

    /* Escreve matrizes de derivadas, se foi pedido: */
    if (prefDebug != NULL)
      {
        nomeArq = NULL; char *nomeArq = jsprintf("%s-%02d-dZdX.ival", prefDebug, prof);
        fprintf(stderr, "%*sEscrevendo %s...\n", ind, "", nomeArq);
        //escreve_imagem_valores(nomeArq, IDX);
	FILE* arq_IDX = fopen(nomeArq,"wt");
	assert(arq_IDX != NULL);
	float_image_write(arq_IDX,IDX);
	fclose(arq_IDX);
        free(nomeArq);
        
        nomeArq = NULL; char *nomeArq = jsprintf("%s-%02d-dZdY.ival", prefDebug, prof);
        fprintf(stderr, "%*sEscrevendo %s...\n", ind, "", nomeArq);
        //escreve_imagem_valores(nomeArq, IDY);
	FILE* arq_IDY = fopen(nomeArq,"wt");
	assert(arq_IDY != NULL);
	float_image_write(arq_IDY,IDY);
	fclose(arq_IDY);
        free(nomeArq);
      }
    
    /* Cria a imagem de alturas: */
    float_image_t* IZ;
    
    /* Decide se deve usar o m�odo multi-escala: */
    int pequena = ((tx < TAM_CORTE) && (ty < TAM_CORTE));
    int atomica = ((tx == 1) || (ty == 1));
    if (pequena || atomica)
      { /* Inicializa {IZ} com zeros: */
        //IZ = cria_imagem_valores(tx+1, ty+1);
	IZ = float_image_new(1,tx+1,ty+1);
        int x, y;
        for(y = 0; y < IZ->sz[2]; y++) 
          for(x = 0; x < IZ->sz[1]; x++) 
            {
		// IZ->pixel[y][x] = 0; 
		float_image_set_sample(IZ,0,x,y,0.0);
	    }
      }
    else
      { /* Inicializa {IZ} com por multi-escala. */
      
        /* Reduz imagens de derivadas pela metade: */
        fprintf(stderr, "%*sReduzindo imagens de derivadas...\n", ind, "");
        float_image_t *JDX = reduz_derivadas(IDX);
        float_image_t *JDY = reduz_derivadas(IDY);
	float_image_t *JW = reduz_pesos(IW);
        
        /* Calcula imagem de alturas reduzida: */
        float_image_t *JZ = calcula_alturas_recursiva(JDX, JDY,JW, prof+1, prefDebug);
        
        /* Expande imagem de alturas para tamanho original: */
        fprintf(stderr, "%*sExpandindo alturas para %d %d...\n", ind, "", tx+1, ty+1);
        IZ = expande_alturas(JZ, tx+1, ty+1);
      }

    /* Monta sistema linear: */
    long int N = (tx+1)*(ty+1);
    fprintf(stderr, "%*sCriando o sistema %ld %ld...\n", ind, "", N, N);
    sistema_t* S = cria_sistema(N);
    preenche_sistema(IDX, IDY,IW, S);

    if (prefDebug != NULL)
      {
        nomeArq = NULL; char *nomeArq = jsprintf("%s-%02d.sist", prefDebug, prof);
        fprintf(stderr, "%*sEscrevendo %s...\n", ind, "", nomeArq);
        escreve_sistema(S, nomeArq, IZ->sz[1], IZ->sz[2]);
        free(nomeArq);
      }

    /* Resolve o sistema, obtendo as alturas: */
    fprintf(stderr, "%*sResolvendo Sistema...\n", ind, "");
    long int max_iter = 1000;
    double tol = 0.00001;
    double *Z = (double*)malloc(sizeof(double)*S->N);
    int para = 0; /* 1 significa solu�o paralela, 0 sequencial. */
    int szero = 1; /* 1 significa ajustar soma em zero, 0 deixar livre. */
    extrai_valores_da_imagem(IZ, Z);
    resolve_sistema(S, Z, max_iter, tol, para, szero);
    coloca_valores_na_imagem(Z, IZ);
    
    /* Escreve matriz parcial de alturas, se foi pedido: */
    if (prefDebug != NULL)
      {
        nomeArq = NULL; char *nomeArq = jsprintf("%s-%02d-alturas.ival", prefDebug, prof);
        fprintf(stderr, "%*sEscrevendo %s...\n", ind, "", nomeArq);
        //escreve_imagem_valores(nomeArq, IZ);
	FILE* arq_IZ = fopen(nomeArq,"wt");
	assert(arq_IZ != NULL);
	float_image_write(arq_IZ,IZ);
	fclose(arq_IZ);
        free(nomeArq);
      }
    int kx,ky;
    kx = IZ->sz[1];
    ky = IZ->sz[2];
    fprintf(stderr, "%*sFim do n�el %d (IZ = %d %d)...\n", ind, "", prof, kx, ky);

    return IZ;
  }


void preenche_sistema_interpolacao(float_image_t* ID,float_image_t*IW,sistema_t* S){
/* Obtem o tamanho das imagens de derivadas e verifica consist�cia: */
    int tx = ID->sz[1];assert(tx == IW->sz[1]);
    int ty = ID->sz[2];assert(ty == IW->sz[2]);
    
    /* Verifica tamanho do sistema: */
    assert(S->N == (tx*ty));
    
    auto long int indZ(int y, int x); 
      /* Calcula o �dice do píxel {[y][x]} na imagem {IZ}. */
    
    long int indZ(int y, int x) { return x + y*tx; }
    
    /* Garante que o sistema �represent�el: */
    assert(MAX_COEFFS >= 5);

    /* Enumera v�tices da imagem {ID}*/
    int x, y;
    for(y = 0; y < ty; y++)
      {
        for(x = 0; x < tx; x++)
          {
            /* Calcula �dice {i} da equa�o associada derivada ID(x,y) */
            long int i = indZ(y, x); 
            
            /* Pega a equa�o {i} e limpa a mesma: */
            equacao_t *eqi = &(S->eq[i]);
            int k;
            for (k = 0; k < MAX_COEFFS; k++) 
              { eqi->coef[k] = 0.0; eqi->ind[k] = -1; }
            eqi->indep = 0.0;

	   float peso = float_image_get_sample(IW,0,x,y);
	   float d = float_image_get_sample(ID,0,x,y);
	   if( isnan(d) ){ 
		peso = 0;
		d = 0;
	   }
	   /*Monta equação da média*/
	   int num_vizinhos = 0;	
	   auto void acrescenta_termo(int x1,int y1);
	   auto void acrescenta_termo(int x1,int y1){
		if( (x1 >= 0) && (x1 < tx) && (y1 >=0) && (y1 < ty)){
			num_vizinhos++;
			eqi->ind[num_vizinhos] = indZ(y1,x1);
			eqi->coef[num_vizinhos] = + 1.00 - peso;
		}
	   }	

	   acrescenta_termo(x-1,y);
	   acrescenta_termo(x+1,y);
	   acrescenta_termo(x,y-1);
	   acrescenta_termo(x,y+1);
           eqi->ind[0] = i;
           eqi->coef[0] =  -num_vizinhos;
           eqi->indep = -num_vizinhos*d*peso;

 	}
     }
      
	
}

void interpola_derivadas_recursivo(float_image_t *IDX, 
    float_image_t *IDY, 
    float_image_t *IW,
    int prof,
    char *prefDebug
  ){
	    /* Pega tamanho das imagens e verifica compatibilidade: */
//     int tx = IDX->tx; assert(tx == IDY->tx);assert(tx == IW->tx);
//     int ty = IDX->ty; assert(ty == IDY->ty);assert(ty == IW->ty);

     int tx = IDX->sz[1]; assert(tx == IDY->sz[1]);assert(tx == IW->sz[1]);
     int ty = IDX->sz[2]; assert(ty == IDY->sz[2]);assert(ty == IW->sz[2]);
	

    int ind = 2*prof+2; /* Indenta�o dos coment�ios. */

    fprintf(stderr, "%*sInicio do nivel %d (IDX,IDY,IW = %d %d)...\n", ind, "", prof, tx, ty);
    char *nomeArq;

    /* Escreve matrizes de derivadas, se foi pedido: */
    if (prefDebug != NULL)
      {
        nomeArq = NULL; char *nomeArq = jsprintf("%s-%02d-dZdX.ival", prefDebug, prof);
        fprintf(stderr, "%*sEscrevendo %s...\n", ind, "", nomeArq);
        //escreve_imagem_valores(nomeArq, IDX);
	FILE* arq_IDX = fopen(nomeArq,"wt");
	assert(arq_IDX != NULL);
	float_image_write(arq_IDX,IDX);
	fclose(arq_IDX);
        free(nomeArq);
        
        nomeArq = NULL; char *nomeArq = jsprintf("%s-%02d-dZdY.ival", prefDebug, prof);
        fprintf(stderr, "%*sEscrevendo %s...\n", ind, "", nomeArq);
        //escreve_imagem_valores(nomeArq, IDY);
	FILE* arq_IDY = fopen(nomeArq,"wt");
	assert(arq_IDY != NULL);
	float_image_write(arq_IDY,IDY);
	fclose(arq_IDY);
        free(nomeArq);
      }
    

    /* Decide se deve usar o m�odo multi-escala: */
    int pequena = ((tx < TAM_CORTE) && (ty < TAM_CORTE));
    int atomica = ((tx == 1) || (ty == 1));
    if (pequena || atomica)
      {
    
	
      }
    else
      { /* Inicializa {IZ} com por multi-escala. */
      
        /* Reduz imagens de derivadas pela metade: */
        fprintf(stderr, "%*sReduzindo imagens de derivadas...\n", ind, "");
        float_image_t *JDX = reduz_derivadas(IDX);
        float_image_t *JDY = reduz_derivadas(IDY);
	float_image_t *JW = reduz_pesos(IW);
        
        /* Calcula imagem de alturas reduzida: */
        interpola_derivadas_recursivo(JDX, JDY,JW, prof+1, prefDebug);
        
        /* Expande imagem de alturas para tamanho original: */
        fprintf(stderr, "%*sExpandindo alturas para %d %d...\n", ind, "", tx+1, ty+1);
        float_image_t* KDX = expande_derivadas(JDX,tx, ty);
	float_image_t* KDY = expande_derivadas(JDY,tx, ty);
	float_image_set_channel(IDX,0,KDX,0);
	float_image_set_channel(IDY,0,KDY,0);
	float_image_free(KDX);
	float_image_free(KDY);
      }

    /* Monta sistema linear: */
    long int N = (tx*ty);
    fprintf(stderr, "%*sCriando o sistema %ld %ld...\n", ind, "", N, N);
    sistema_t* S = cria_sistema(N);
    int eixo;
    for(eixo = 0; eixo < 2; eixo++){
	float_image_t* ID = (eixo ==0? IDX: IDY);
	char* nome_eixo = (eixo == 0? "X":"Y");
	preenche_sistema_interpolacao(ID,IW, S);
	
	if (prefDebug != NULL)
	{
		nomeArq = NULL; char *nomeArq = jsprintf("%s-%02d-%s.sist", prefDebug, prof,nome_eixo);
		fprintf(stderr, "%*sEscrevendo %s...\n", ind, "", nomeArq);
		escreve_sistema(S, nomeArq, ID->sz[1], ID->sz[2]);
		free(nomeArq);
	}
	
	/* Resolve o sistema, interpolando as derivadas: */
	fprintf(stderr, "%*sResolvendo Sistema...\n", ind, "");
	long int max_iter = 1000;
	double tol = 0.00001;
	double *Z = (double*)malloc(sizeof(double)*S->N);
	int para = 0; /* 1 significa solu�o paralela, 0 sequencial. */
	int szero = 0; /* 1 significa ajustar soma em zero, 0 deixar livre. */
	extrai_valores_da_imagem(ID,Z);
	resolve_sistema(S, Z, max_iter, tol, para, szero);
	coloca_valores_na_imagem(Z, ID);
    }
    /* Escreve matriz parcial de alturas, se foi pedido: */
    if (prefDebug != NULL)
      {
       nomeArq = NULL; char *nomeArq = jsprintf("%s-%02d-dZdX-int.ival", prefDebug, prof);
        fprintf(stderr, "%*sEscrevendo %s...\n", ind, "", nomeArq);
        //escreve_imagem_valores(nomeArq, IDX);
	FILE* arq_IDX = fopen(nomeArq,"wt");
	assert(arq_IDX != NULL);
	float_image_write(arq_IDX,IDX);
	fclose(arq_IDX);
        free(nomeArq);
        
        nomeArq = NULL; char *nomeArq = jsprintf("%s-%02d-dZdY-int.ival", prefDebug, prof);
        fprintf(stderr, "%*sEscrevendo %s...\n", ind, "", nomeArq);
        //escreve_imagem_valores(nomeArq, IDY);
	FILE* arq_IDY = fopen(nomeArq,"wt");
	assert(arq_IDY != NULL);
	float_image_write(arq_IDY,IDY);
	fclose(arq_IDY);
        free(nomeArq);
      }
    
    fprintf(stderr, "%*sFim do nível %d ...\n", ind, "", prof);

}


void escreve_sistema(sistema_t *S, char *nome, int tx, int ty)
  {
    FILE* arq = fopen(nome, "wt");
    assert(arq != NULL); 

    long int i;
    for(i = 0; i < S->N; i++)
      {
        /* Pega a equa�o {i}: */
        equacao_t *eqi = &(S->eq[i]);
        
        /* Calcula indices do pixel correspondente �inc�nita {Z[i]}: */
        int xi = i%tx; int yi = i/tx;
         
        /* Imprime �dices do pixel principal: */
        fprintf(arq, "eq[%d][%d]:", xi, yi);

        /* Imprime os termos lineares da equa�o: */
        int k;
        for(k = 0; k < MAX_COEFFS; k++)
          { /* Pega mais uma inc�nita {Z[j]} que entra na equa�o {i}: */
            long int j = eqi->ind[k];
            /* Calcula �dices do pixel correspondente: */
            int xj = j%tx; int yj = j/tx;
            double coef_j = eqi->coef[k];
            if (k > 0) { fprintf(arq, " + "); }
            fprintf(arq, "%f*Z[%d][%d]", coef_j, xj, yj);
          }
        /* Imprime o termo independente da equa�o: */
        fprintf(arq, " = %f", S->eq[i].indep);
        fprintf(arq, "\n");
      }
    fclose(arq);
  }

void extrai_valores_da_imagem(float_image_t *IZ, double *Z)
  {
    long int N = IZ->sz[1] * IZ->sz[2];
    long int i;
    for(i = 0; i < N; i++)
      { int x = i%IZ->sz[1];
        int y = i/IZ->sz[1];
        //Z[i] = IZ->pixel[y][x];
	Z[i] = float_image_get_sample(IZ,0,x,y);
      }
  }



void coloca_valores_na_imagem(double *Z, float_image_t *IZ)
  {
    long int N = IZ->sz[1] * IZ->sz[2];
    long int i;
    for(i = 0; i < N; i++)
      { int x = i%IZ->sz[1];
        int y = i/IZ->sz[1];
        //IZ->pixel[y][x] = Z[i];
	float_image_set_sample(IZ,0,x,y,Z[i]);
      }
  }

float_image_t *calcula_derivadas(float_image_t *IN, char eixo, float_image_t* IW){
    assert((eixo == 'x') || (eixo == 'y'));

    int tx = IN->sz[1]; 
    int ty = IN->sz[2];
    //float_image_t *ID = cria_imagem_valores(tx, ty);
    float_image_t *ID = float_image_new(1,tx,ty);
    
    int x,y;
    for(y =0; y < ID->sz[2]; y++)
      { for (x = 0; x < ID->sz[1]; x++)
	   { 
		r3_t v;
		v.c[0] = float_image_get_sample(IN,0,x,y);
		v.c[1] = float_image_get_sample(IN,1,x,y);
		v.c[2] = float_image_get_sample(IN,2,x,y);
		//vetor_t *v = &(IN->pixel[y][x]);
            //float *dpix = &(ID->pixel[y][x]);
	    float dpix = float_image_get_sample(ID,0,x,y);
            /* A normal deve estar apontando para cima: */
	    if( r3_norm(&v) == 0) {
		float_image_set_sample(ID,0,x,y,NAN);
		float_image_set_sample(IW,0,x,y,1.0e-10);
	    }else{
		double ZMIN = 0.10; //Minumum valid value for z coordinate of a normal
		if (v.c[2] <= ZMIN)
		{ fprintf(stderr, "normal invalida %+8.5f %+8.5f %+8.5f\n", v.c[0], v.c[1], v.c[2]); 
			
			v.c[2] = 1.0e-10;
			r3_dir(&v,&v);
			r3_scale(sqrt(1 - ZMIN*ZMIN),&v,&v);
			v.c[2] = ZMIN;
			r3_dir(&v,&v);
			float_image_set_sample(IN,0,x,y,v.c[0]);
			float_image_set_sample(IN,1,x,y,v.c[1]);
			float_image_set_sample(IN,2,x,y,v.c[2]);
			float_image_set_sample(IW,0,x,y,1.0e-10);
	//                 float *wpix = &(IW->pixel[y][x]);
	//                 (*wpix) = 1.0e-10;
		}
		/* Calcula derivada a partir da normal: */
		if(eixo == 'x')
		{
			//(*dpix) = - v->x / v->z;
			dpix = -v.c[0]/v.c[2];
			float_image_set_sample(ID,0,x,y,dpix);
		}
		else if (eixo == 'y')
		{ 
			//(*dpix) = - v->y / v->z;
			dpix = -v.c[1]/v.c[2];
			float_image_set_sample(ID,0,x,y,dpix);
		}
	   }
          }
      }
    return ID;
  }

float_image_t *reduz_derivadas(float_image_t *ID)
  {
    int txI = ID->sz[1]; int txJ = (txI+1)/2;
    int tyI = ID->sz[2]; int tyJ = (tyI+1)/2;
    //float_image_t *JD = cria_imagem_valores(txJ, tyJ);
    float_image_t *JD = float_image_new(1,txJ, tyJ);
    int xJ, yJ;
    for(yJ = 0; yJ < tyJ; yJ++)
      { for(xJ = 0; xJ < txJ; xJ++)
          { /* Pixel da imagem reduzida: */
            //float *pxJ = &(JD->pixel[yJ][xJ]);
	    /* Pixels da imagem original (preenchendo bordas): */
            int xI = 2*xJ, yI = 2*yJ;
            //double d00 = ID->pixel[yI+0][xI+0];
	    double d00 = float_image_get_sample(ID,0,xI+0,yI+0);
            //double d01 = ((xI+1 < txI) ? ID->pixel[yI+0][xI+1] : d00);
	    double d01 = ((xI+1 < txI) ? float_image_get_sample(ID,0,xI+1,yI+0) : NAN);
            //double d10 = ((yI+1 < tyI) ? ID->pixel[yI+1][xI+0] : d00);
	    double d10 = ((yI+1 < tyI) ? float_image_get_sample(ID,0,xI+0,yI+1) : NAN);
            //double d11 = ((xI+1 < txI) && (yI+1 < tyI) ? ID->pixel[yI+1][xI+1] : d00);
	    double d11 = ((xI+1 < txI) && (yI+1 < tyI) ? float_image_get_sample(ID,0,xI+1,yI+1) : NAN);
	    /* Tanto {Z} quanto {X} e {Y} encolheram, logo a derivada n� muda: */
	   /*calcula a média , ignorando NANs*/
	   double soma = 0;
	   int ns = 0;
 	   if(!(isnan(d00)) ){ soma+=d00; ns++;}
	   if(!(isnan(d01)) ){ soma+=d01; ns++;}
	   if(!(isnan(d10)) ){ soma+=d10; ns++;}
	   if(!(isnan(d11)) ){ soma+=d11; ns++;}
	   double avg = (ns == 0? NAN: soma/(double)ns);
           // (*pxJ) = (d00 + d01 + d10 + d11)/4;
	    float_image_set_sample(JD,0,xJ,yJ,avg);
          }
      }
    return JD;
  }

float_image_t *reduz_pesos(float_image_t *ID)
  {
    int txI = ID->sz[1]; int txJ = (txI+1)/2;
    int tyI = ID->sz[2]; int tyJ = (tyI+1)/2;
    //float_image_t *JD = cria_imagem_valores(txJ, tyJ);
    float_image_t *JD = float_image_new(1,txJ, tyJ);
    int xJ, yJ;
    for(yJ = 0; yJ < tyJ; yJ++)
      { for(xJ = 0; xJ < txJ; xJ++)
          { /* Pixel da imagem reduzida: */
            /* Pixels da imagem original (preenchendo bordas): */
            int xI = 2*xJ, yI = 2*yJ;
            double d00 = float_image_get_sample(ID,0,xI+0,yI+0);
            double d01 = ((xI+1 < txI) ? float_image_get_sample(ID,0,xI+1,yI+0) : d00);
            double d10 = ((yI+1 < tyI) ? float_image_get_sample(ID,0,xI+0,yI+1) : d00);
            double d11 = ((xI+1 < txI) && (yI+1 < tyI) ? float_image_get_sample(ID,0,xI+1,yI+1) : d00);
            /* Tanto {Z} quanto {X} e {Y} encolheram, logo a derivada n� muda: */
            double mean =  4.0/(1.0/d00 + 1.0/d01 + 1.0/d10 + 1.0/d11);
	    
	    float_image_set_sample(JD,0,xJ,yJ, mean);
          }
      }
    return JD;
  }

float_image_t *expande_alturas(float_image_t *JZ, int tx, int ty)
  {
    int txJ = JZ->sz[1]; assert(tx/2 + 1 == txJ);
    int tyJ = JZ->sz[2]; assert(ty/2 + 1 == tyJ);
    //float_image_t *IZ = cria_imagem_valores(tx, ty);
    float_image_t *IZ = float_image_new(1,tx, ty);
    int xI, yI;
    for(yI = 0; yI < ty; yI++)
      { for(xI = 0; xI < tx; xI++)
          { /* Pixel da imagem expandida: */
          //  float *pxI = &(IZ->pixel[yI][xI]);
            /* Pixels relevantes da imagem original: */
            int xJ0 = xI/2, xJ1 = (xI+1)/2;
            int yJ0 = yI/2, yJ1 = (yI+1)/2;
            assert(xJ1 < txJ); assert(yJ1 < tyJ);
            //double d00 = JZ->pixel[yJ0][xJ0];
	    double d00 = float_image_get_sample(JZ,0,xJ0,yJ0);
            //double d01 = JZ->pixel[yJ0][xJ1];
	    double d01 = float_image_get_sample(JZ,0,xJ1,yJ0);
            //double d10 = JZ->pixel[yJ1][xJ0];
	    double d10 = float_image_get_sample(JZ,0,xJ0,yJ1);
            //double d11 = JZ->pixel[yJ1][xJ1];
	    double d11 = float_image_get_sample(JZ,0,xJ1,yJ1);
            /* Alturas devem ser ampliadas tamb�: */
            //(*pxI) = (d00 + d01 + d10 + d11)/2;
	    float_image_set_sample(IZ,0,xI,yI,(d00 + d01 + d10 + d11)/2);
          }
      }
    return IZ;
  }

float_image_t *expande_derivadas(float_image_t *JD, int tx, int ty)
  {
    int txJ = JD->sz[1]; assert((tx +1)/2 == txJ);
    int tyJ = JD->sz[2]; assert((ty +1)/2 == tyJ);
    //float_image_t *IZ = cria_imagem_valores(tx, ty);
    float_image_t *ID = float_image_new(1,tx, ty);
    int xI, yI;
    for(yI = 0; yI < ty; yI++)
      { for(xI = 0; xI < tx; xI++)
          { /* Pixel da imagem expandida: */
          //  float *pxI = &(IZ->pixel[yI][xI]);
            /* Pixels relevantes da imagem original: */
            int xJ0 = xI/2;
            int yJ0 = yI/2;
            //double d00 = JZ->pixel[yJ0][xJ0];
	    double d00 = float_image_get_sample(JD,0,xJ0,yJ0);
            //double d01 = JZ->pixel[yJ0][xJ1];
	   
            /* Alturas devem ser ampliadas tamb�: */
            //(*pxI) = (d00 + d01 + d10 + d11)/2;
	    float_image_set_sample(ID,0,xI,yI,d00);
          }
      }
    return ID;
  }

void preenche_sistema(float_image_t *IDX, float_image_t *IDY,float_image_t *IW, sistema_t* S)
  {
    /* 
      Cada equa�o e cada inc�nita do sistema {S} correponde a um
      elemento da imagem de alturas {IZ}, ou seja, a um v�tice
      (canto de pixel) das matrizes {IDX} e {IDY}. 
      
      Por defini�o, o v�tice {IZ[y][x]} tem coordenadas {(x,y)}.
      Portanto o pixel {IDX[y][x]} �o quadrado de lado 1 cujo canto
      inferior esquerdo �o v�tice de �dices {[y][x]}, isto � cuja
      diagonal vai de {(x,y)} a {(x+1,y+1)}.
      
      Em geral, a equa�o de {S} associada ao v�tice {[y][x]} tem a
      forma {Lz[y][x] = Ld[y][x]}, onde {Lz,Ld} s� duas estimativas
      do laplaciano {d2z/dx2 + d2z/dy2} da fun�o altura no v�tice
      {[y][x]}. A estimativa {Lz} �obtida da imagem de alturas {IZ}
      por diferen�s finitas entre a altura {IZ[y][x]} no v�tice e as
      alturas de seus quatro vizinhos mais pr�imos, e portanto �uma
      fun�o linear dessas cinco inc�nitas.
      
      A estimativa {Ld[y][x]} �obtida a partir das derivadas
      fornecidas {IDX} e {IDY}, calculadas nos centros dos quatro
      pixels {[y][x]} dessas imagens incidentes ao v�tice {[y][x]};
      isto � {IDX[y-dy][x-dx]} e {IDY[y-dy][x-dx]} onde {dy} e {dx}
      s� 0 ou 1. Note-se que {Ld} �independente das inc�nitas.
      
      Na verdade, esta descri�o vale apenas para os v�tices {[y][x]}
      que est� cercados por quatro pixels nas imagens {IDX} e {IDY}.
      Para v�tices com �dices {[0][x]}, ao longo da borda inferior
      dessas imagens (menos os dois v�tices extremos {[0][0]} e
      {[0][tx]}, a equa�o �{Dz[0][x] = Dd[0][x]}, onde {Dz} e {Dd}
      s� duas estimativas para a derivada vertical da altura no ponto
      {(x,0.5)}. A estimativa {Dz} �calculada por diferen�s finitas
      entre as alturas nos dois v�tice mais pr�imos, isto �
      {Dz[0][x] = IZ[1][x] - IZ[0][x]}. A estimativa {Dd} �calculada
      a partir de {IDX[y][x-dx]} e {IDY[y][x-dx]} onde {dx} �0 ou 1.
      As equa�es ao longo das outras tr� bordas s� an�ogas.
      
      No v�tice {[0][0]} (canto inferior esquerdo das imagens
      {IDX} e {IDY}), a equa�o �tamb� {Dz[y][x] = Dd[y][x]}, mas as
      estimativas se referem �derivada da altura no ponto {(0.5,0.5)}
      e ao longo do vetor {(1,1)}. Ou seja {Dz[0][0] = IZ[1][1] -
      IZ[0][0]}, e {Dd[0][0] = IDX[0][0] + IDY[0][0]}. As equa�es
      para os outros tr� v�tices ({[0][tx]}, {[ty][0]} e {[ty][tx]})
      s� an�ogas.
      
      Este sistema �resolvido pelo m�odo de Gauss-Jordan. Este
      sistema �indeterminado pois existe uma solu�o n� nula para o
      sistema homog�eo {M Z = 0}, que �todas as alturas iguais a 1.
      (Note que toda linha de {M} tem soma zero, pois �uma estimativa
      da derivada ou do laplaciano da altura.)
    */

    /* Obtem o tamanho das imagens de derivadas e verifica consist�cia: */
    int tx = IDX->sz[1]; assert(tx == IDY->sz[1]);assert(tx == IW->sz[1]);
    int ty = IDX->sz[2]; assert(ty == IDY->sz[2]);assert(ty == IW->sz[2]);
    
    /* Verifica tamanho do sistema: */
    assert(S->N == (tx+1)*(ty+1));
    
    auto long int indZ(int y, int x); 
      /* Calcula o �dice do v�tice {[y][x]} na imagem {IZ}. */
    
    long int indZ(int y, int x) { return x + y*(tx+1); }
    
    /* Garante que o sistema �represent�el: */
    assert(MAX_COEFFS >= 5);

    /* Enumera v�tices das imagens {IDX,IDY}, i.e. elems da imagem de alturas: */
    int x, y;
    for(y = 0; y <= ty; y++)
      {
        for(x = 0; x <= tx; x++)
          {
            /* Calcula �dice {i} da equa�o associada �altura {IZ[y][x]}: */
            long int i = indZ(y, x); 
            
            /* Pega a equa�o {i} e limpa a mesma: */
            equacao_t *eqi = &(S->eq[i]);
            int k;
            for (k = 0; k < MAX_COEFFS; k++) 
              { eqi->coef[k] = 0.0; eqi->ind[k] = -1; }
            eqi->indep = 0.0;

	    /*pega células vizinhas refletindo as bordas*/
	    int xm = ((x > 0) ? x-1: 0);
	    int xp = ((x < tx) ? x: tx-1);
	    int ym = ((y > 0) ? y-1: 0);
	    int yp = ((y < ty) ? y: ty-1);
	    /*Obtem os dados relevantes nas células*/
	    double dxmm = float_image_get_sample(IDX,0,xm,ym);
	    double dxmp = float_image_get_sample(IDX,0,xm,yp);
 	    double dxpm = float_image_get_sample(IDX,0,xp,ym);
            double dxpp = float_image_get_sample(IDX,0,xp,yp);

	    double dymm = float_image_get_sample(IDY,0,xm,ym);
	    double dymp = float_image_get_sample(IDY,0,xm,yp);
 	    double dypm = float_image_get_sample(IDY,0,xp,ym);
            double dypp = float_image_get_sample(IDY,0,xp,yp);

	    double wmm = float_image_get_sample(IW,0,xm,ym);
	    double wmp = float_image_get_sample(IW,0,xm,yp);
 	    double wpm = float_image_get_sample(IW,0,xp,ym);
            double wpp = float_image_get_sample(IW,0,xp,yp);

	    /*Calcula derivadas e pesos nas arestas*/
	    double wdxmo = (dxmm*wmm + dxmp*wmp);
	    double wdxpo = (dxpm*wpm + dxpp*wpp);
            double wdyom = (dymm*wmm + dypm*wpm);
	    double wdyop = (dymp*wmp + dypp*wpp);

	    double wmo = wmm+wmp;
	    double wpo = wpm+wpp;
	    double wom = wmm+wpm;
	    double wop = wmp+wpp;
	    
            /*elimina as arestas que não existem*/
	    if( x ==0 ){wmo = 0; wdxmo = 0;}
	    if( x == tx ){wpo = 0; wdxpo = 0;}
	    if( y ==0 ){wom = 0; wdyom = 0;}
	    if( y == ty ){wop = 0; wdyop = 0;}
	   double woo =  wmo+wpo+wom+wop;
	   assert(woo > 0);
	   /*monta a equação*/
	   k = 0;
           double indep =0;
	   eqi->ind[k] = i;                eqi->coef[k] = -woo; k++;
	   if( wmo != 0) { eqi->ind[k] = indZ(y,x-1); eqi->coef[k] = wmo; indep-=wdxmo;  k++;}
	   if( wpo != 0) { eqi->ind[k] = indZ(y,x+1); eqi->coef[k] = wpo; indep+=wdxpo;  k++;}
	   if( wom != 0) { eqi->ind[k] = indZ(y-1,x); eqi->coef[k] = wom; indep-=wdyom;  k++;}
	   if( wop != 0) { eqi->ind[k] = indZ(y+1,x); eqi->coef[k] = wop; indep+=wdyop;  k++;}
	   eqi->indep = indep;

          
      }
  }
}
