/* Last edited on 2017-06-22 18:14:28 by stolfilocal */

#define _GNU_SOURCE // Para obter "asprintf"
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <assert.h>
//#include <values.h>
#include <float.h>
#include <vetorn.h>
#include <float_image.h>
#include <float_image_io_pnm.h>
#include <hash.h>
#include <normais.h>
#include <sym_eigen.h>
#include <i2.h>
#include <affirm.h>

struct ListaPosicao{
	int posicao[2];
	double dist;
	struct ListaPosicao* proximo;
};

typedef struct ListaPosicao listapos;

struct Bucket{
	double* centro; // assinatura média do bucket
	double raio; // maior distancia, tab chamada de "Rô"[i,j]
	struct ListaPonteiro* itens;

};

struct ListaPonteiro{
       int ponteiro;
       struct ListaPonteiro* proximo;
};



typedef struct ListaPonteiro listapont;

typedef struct Bucket bucket;



struct BucketGrid
{
  int tam_grid;
  int num_luzes;

  double* u; /* Direção principal da nuvem de pontos {[0..num_luzes-1]}. */
  double* v; /* Direção principal da nuvem de pontos {[0..num_luzes-1]}. */
  bucket** buckets; /* Matriz de buckets {[0..tam_grid-1][0..tam_grid-1]}. */
  double bu, bv; /* Posição do baricentro projetado em {u,v} */
  int64_t** m_statistic_find; /* Numero de querys que são hasheadas para o bucket {[iu][iv]} */
  int64_t** m_statistic_eval; /* Numero de distancias euclideanas calculadas no bucket [x][y] */
  int64_t** m_statistic_eval_h; /* Numero de distancias euclideanas calculadas no bucket [x][y] por bucket mapeado*/
  int64_t** m_statistic_scan; /* Numero de vezes que o bucket é examinado nas buscas */
  int64_t** m_statistic_scan_h; /* Numero de buckets scaneados por mapeamento*/
  int64_t** m_statistic_hashed; /* Número de normais mapeadas para um buckets*/
  int64_t max_count_st; /* Máximo de {m_statistic[iu][iv]} sobre todos os buckets. */
  int tamanho_lista_preordenada; /* Índices dos buckets, em ordem de distância crescente. */
  listapos* lista_preordenada; /* Índices dos buckets, em ordem de distância crescente. */
  double R; /* Raio da grade de buckets no plano {u,v}. */
  double* baricentro; /* Centro da grade de buckets. */
  int64_t** circle_map; /**/
};

int procura_no_bucket(double so[], bucket* bitem, Tabela *tab, double *menor_distancia, int *num_euclidP);
void calcula_sistema_de_coordenadas(Tabela *tab, double u[], double v[], double bar[], double *bu, double *bv);
void criaMatrizCovariancia2(Tabela*tab, double M[], double bar[]);

double produtoEscalar(double v[], double u[], int tam);
void flushMatriz(FILE *arq, double M[], int num_luzes, char *nome);
listapont* insereLista(int ponteiro, listapont* antigo);
void flushTabela(Tabela* tab);
void normalizaVetor(double* vet, int tam);

i2_t mapeiaHash
  ( const double ass[],
    double u[],
    double v[],
    const double baricentro[],
    int num_luzes,
    int tam_grid,
    double R
    );

int tam_grid_bucket(int num_linhas);
double cobrinha ( double k);
double calculaDistanciaBucketMinima(int i, int j, int iq, int jq, double tam_bucket);

bucketGrid* CriaBucketGrid(Tabela* tab,int gridsize)
{
  bucketGrid* grid = (bucketGrid*)malloc(sizeof(bucketGrid));
  int num_linhas  = get_num_linhas(tab);
  int num_luzes = get_num_luzes(tab);
  int tam_grid = gridsize;
  if(gridsize == -1) tam_grid  = tam_grid_bucket(num_linhas);

  grid->num_luzes = num_luzes;
  fprintf(stderr,"GRID SIZE %dx%d\n",tam_grid,tam_grid);
  /* Aloca e inicializa as estatísticas. */
  grid->m_statistic_find = (int64_t**)malloc(sizeof(int64_t*)*tam_grid);
  grid->m_statistic_eval = (int64_t**)malloc(sizeof(int64_t*)*tam_grid);
  grid->m_statistic_eval_h = (int64_t**)malloc(sizeof(int64_t*)*tam_grid);
  grid->m_statistic_scan = (int64_t**)malloc(sizeof(int64_t*)*tam_grid);
  grid->circle_map = (int64_t**)malloc(sizeof(int64_t*)*tam_grid);
  grid->m_statistic_scan_h = (int64_t**)malloc(sizeof(int64_t*)*tam_grid);
  grid->m_statistic_hashed = (int64_t**)malloc(sizeof(int64_t*)*tam_grid);
  int i, j;
  for(i = 0; i < tam_grid; i++){
    grid->m_statistic_find[i] = (int64_t*) malloc(sizeof(int64_t)* tam_grid);
    grid->m_statistic_eval[i] = (int64_t*) malloc(sizeof(int64_t)* tam_grid);
    grid->m_statistic_eval_h[i] = (int64_t*) malloc(sizeof(int64_t)* tam_grid);
    grid->m_statistic_scan[i] =  (int64_t*) malloc(sizeof(int64_t)* tam_grid);
    grid->circle_map[i]	= (int64_t*) malloc(sizeof(int64_t)* tam_grid);
    grid->m_statistic_scan_h[i] = (int64_t*) malloc(sizeof(int64_t)* tam_grid);
    grid->m_statistic_hashed[i] = (int64_t*) malloc(sizeof(int64_t)* tam_grid);
    for(j = 0; j < tam_grid;j++){
      grid->m_statistic_find[i][j] = 0;
      grid->m_statistic_eval[i][j] = 0;
      grid->m_statistic_eval_h[i][j] = 0;
      grid->m_statistic_scan[i][j] = 0;
      grid->circle_map[i][j] = 0;
      grid->m_statistic_scan_h[i][j] = 0;
      grid->m_statistic_hashed[i][j] = 0;
    }
  }
  fprintf(stderr,  "Alocadas matrizes de estatistica\n");

  /* Calcula sistema de coordenadas da grade: */
  double *u = (double*)malloc(sizeof(double)*num_luzes);
  double *v = (double*)malloc(sizeof(double)*num_luzes);
  double *bar = (double*)malloc(sizeof(double)*num_luzes);
  calcula_sistema_de_coordenadas(tab, u, v, bar, &grid->bu, &grid->bv);
  /* fprintf(stderr,  "U = ");
     for(t = 0; t < num_luzes;  t++) fprintf(stderr,  "%f ",u[t]);
     fprintf(stderr,  "\n");

     fprintf(stderr,  "V = ");
     for(t = 0; t < num_luzes;  t++) fprintf(stderr,  "%f ",v[t]);
     fprintf(stderr,  "\n\n"); */
  grid->baricentro = bar;
  grid->u = u;
  grid->v = v;
  fprintf(stderr,  "void calcula_sistema_de_coordenadas(Tabela *tab, double u[], double v[], double bar[], double *bu, double *bv)Calculado o sistema de coordenadas da grade.\n");

  /* Aloca e inicializa buckets: */
  grid->tam_grid = tam_grid;
  grid->buckets = (bucket**)notnull(malloc(sizeof(bucket*)*tam_grid),"NOMEM for bucket grid");
  for(i = 0; i < tam_grid; i++){

    grid->buckets[i] = (bucket*)notnull(malloc(sizeof(bucket)*tam_grid),"NOMEM for bucket grid - 2");
    for(j = 0; j < tam_grid; j++){
      grid->buckets[i][j].itens = NULL;
      grid->buckets[i][j].centro = (double*)notnull(malloc(sizeof(double)*num_luzes),"NOMEM for centroid");

      grid->buckets[i][j].raio = 0;
      int k;
      for(k = 0 ;k < num_luzes;k++){
	grid->buckets[i][j].centro[k] = 0;
      }
    }
  }
  fprintf(stderr,  "Buckets criados.\n");

  /* Calcula raio {grid->R} da grade. */
  fprintf(stderr,  "\n");
  //flushTabela(tab);
  grid->R = -1;
  for(i = 0; i < num_linhas; i++){
    const double *go = get_intdir(tab,i);
    double dgo[num_luzes];
    int t;
    for (t = 0; t < num_luzes; t++) dgo[t] = go[t] -  bar[t];
    double su = fabs(produtoEscalar(dgo,u,num_luzes));
    double sv = fabs(produtoEscalar(dgo,v,num_luzes));
    // fprintf(stderr,  "linha %5d SU = %f SV = %f\n",i,su,sv);
    if(grid->R < su) grid->R = su;
    if(grid->R < sv) grid->R = sv;
  }
  fprintf(stderr,  "R: %f\n",grid->R);


  /* Enchendo os buckets com linhas da tabela {tab}, calcula centro de cada bucket. */
  for(i = 0; i < num_linhas; i++){

    const double *go = get_intdir(tab,i);
    i2_t indice = mapeiaHash(go,u,v,bar,num_luzes,tam_grid,grid->R);
    int x = indice.c[0];
    int y = indice.c[1];

    grid->buckets[x][y].itens = insereLista(i,grid->buckets[x][y].itens);
    for(j = 0; j < num_luzes; j++){

      grid->buckets[x][y].centro[j] += (go[j]);
    }
  }


  int delta_t = 3;
  for(i = 0 ; i < tam_grid; i++){
    for(j = 0; j < tam_grid; j++){
	int ix,iy;
	int test = 0;
	for(ix = -delta_t; ix <= delta_t; ix++){
		if(((ix + i)>= 0) && ((ix + i) < tam_grid)){
			for(iy = -delta_t; iy <=delta_t; iy++){
				if(((iy+j)>= 0) && ((iy+j) < tam_grid)){
					int x,y;
					x = i+ix;
					y = j + iy;
					if(grid->buckets[x][y].itens != NULL) test = 1;
				}
			}
		}
	}
	grid->circle_map[i][j] = test;
	//if(test == 0) fprintf(stderr,"INVALID BUCKET AT %d %d\n",i,j);
	//else fprintf(stderr,"VALID BUCKET AT %d %d\n",i,j);
  }
 }
 // plota_bucket_map("map.plt", grid,0);

  for(i = 0 ; i < tam_grid; i++){
    for(j = 0; j < tam_grid; j++){
      /* Conta pontos no bucket: */
      int compr_lista = 0;
      listapont* p;
      for(p = grid->buckets[i][j].itens; p != NULL; p = p->proximo) { compr_lista++; }
      int k;
      if(compr_lista > 0){
	for(k = 0; k < num_luzes;k++){
	  grid->buckets[i][j].centro[k] /= compr_lista;
	}
      }
      /* Calcula raio do bucket: */

      double maior_dist = -1;
      listapont* maior_entrada = NULL;
      for(p = grid->buckets[i][j].itens; p != NULL; p = p->proximo) {
	const double *go = get_intdir(tab,p->ponteiro);

	double atual_dist = dist_euclid(grid->buckets[i][j].centro,go,num_luzes);

	if(maior_dist < atual_dist) {
	  maior_dist = atual_dist;
	  maior_entrada = p;
	}
      }


	/*double maior_dist = -1;
      listapont* maior_entrada = NULL;
      for(p = grid->buckets[i][j].itens; p != NULL; p = p->proximo) {
	const double *go = get_intdir(tab,p->ponteiro);
	double SU,SV;
	int indice;
	double sig[num_luzes];
	for(indice = 0; indice < num_luzes;indice++){ sig[indice] = go[indice] - grid->baricentro[indice];}
	SU = produtoEscalar(sig,grid->u,num_luzes);
	SV = produtoEscalar(sig,grid->v,num_luzes);
	double sigma = ((2.0*grid->R)/((float)grid->tam_grid))*sqrt(2);
	double atual_dist = sqrt((SU*SU) + (SV*SV));

	if(maior_dist < atual_dist) {
	  maior_dist = atual_dist;
	  maior_entrada = p;
	}
      }	*/

      if( maior_dist > 0.95){
	 const double *maior_go = get_intdir(tab,maior_entrada->ponteiro);
	int kk;
        fprintf(stderr,"\nMAIOR_GO: ");
	for(kk = 0; kk < num_luzes;kk++){
	 	fprintf(stderr, " %8.6f",maior_go[kk]);
         }
	fprintf(stderr,"\nCENTRO  : ");
        for(kk = 0; kk < num_luzes;kk++){
	 	fprintf(stderr, " %8.6f",grid->buckets[i][j].centro[kk]);
         }
	fprintf(stderr,"\n");
      }
      grid->buckets[i][j].raio = maior_dist;
    }
  }

  fprintf(stderr,  "Preenchidos os buckets.\n");

  /* Cria lista preordenadas de posições relativas de buckets: */
  criaListaPreordenada(grid);
  fprintf(stderr,  "Criada a lista de posições.\n");
  fprintf(stderr,  "Grade de buckets construída!");

  return grid;
}

void calcula_sistema_de_coordenadas(Tabela *tab, double u[], double v[], double bar[], double *bu, double *bv)
{
  int num_luzes = get_num_luzes(tab);
  double d[num_luzes];
  double e[num_luzes];

  double R[num_luzes*num_luzes];
  double MM[num_luzes*num_luzes];
  double Mbar[num_luzes];
  criaMatrizCovariancia2(tab, MM, Mbar);
  syei_tridiagonalize(num_luzes, MM, d, e, R);
  int p;
  syei_trid_eigen(num_luzes, d, e, R, &p, 1);
  //  fprintf(stderr,  "Autovetores encontrados - %d \n",p);
  if(p< 2){
    fprintf(stderr,  "Erro ao calcular autovalores");
    exit(1);
  }
  int i;
  for(i = 0;i < num_luzes;i++) {
    u[i] = R[(p-1)*num_luzes + i];
    v[i] = R[(p-2)*num_luzes + i];
    bar[i] = Mbar[i];
  }
  (*bu) = produtoEscalar(bar,u,num_luzes);
  (*bv) = produtoEscalar(bar,v,num_luzes);
}

int64_t** acessaMatriz_Statistic_Euclid(bucketGrid* bg){
	return bg->m_statistic_eval;
}

int64_t** acessaMatriz_Statistic_Scan(bucketGrid* bg){
	return bg->m_statistic_scan;
}

double produtoEscalar(double* v, double* u, int tam){
      int i;
      double soma = 0;
      for(i = 0 ; i< tam;i++){
	//   fprintf(stderr,  "U[%d]:%f V[%d]:%f\n",i,u[i],i,v[i]);
            soma = soma + ((v[i])*(u[i]));
      }
      return soma;

}

listapont* insereLista(int ponteiro, listapont* antigo){
           listapont* grid;
           grid = (listapont*) malloc(sizeof(listapont));
           grid->ponteiro = ponteiro;
           grid->proximo = antigo;
           return grid;
}

void flushMatriz(FILE *arq, double M[], int num_luzes, char *nome) {
     int i,j;
     for(i = 0; i < num_luzes;i++){
           for(j = 0; j < num_luzes ; j++){
	     fprintf(arq,"%010f ",M[i*num_luzes + j]);
           }
           fprintf(arq,"\n");
     }
     fclose(arq);
}

void flushTabela(Tabela* tab){
     int luz,i;
     char nome[200];
     FILE* arq;
     int num_luzes = get_num_luzes(tab);
     sprintf(nome,"out/tabela.txt");
     arq = fopen(nome,"wt");
     for(i = 0; i < get_num_linhas(tab); i++){
		fprintf(arq,"LINHA [%d]:\n",i);
		const double *go = get_intdir(tab,i);
           	for(luz = 0; luz < num_luzes; luz++){ fprintf(arq,"%.6lf ",go[luz]); }
		fprintf(arq,"Magnitude: %.6lf\n",get_intmag(tab,i));
		fprintf(arq,"\n");
           }

	fclose(arq);
}


void flushBucketGrid(FILE* arq, bucketGrid* bg)
{
  int i,j;
  fprintf(arq,"TAMANHO: %d \n",bg->tam_grid);
  fprintf(arq,"BARICENTRO \n");
  for(j = 0; j < bg->num_luzes; j++){
    fprintf(arq,"B[%4d]: %f \n",j,bg->baricentro[j]);
  }
  fprintf(arq,"Vetor U\n");
  for(j = 0; j < bg->num_luzes; j++){
    fprintf(arq,"U[%4d]: %f \n",j,bg->u[j]);
  }

  fprintf(arq,"Vetor V\n");
  for(j = 0; j < bg->num_luzes; j++){
    fprintf(arq,"V[%4d]: %f \n",j,bg->v[j]);
  }

  int totdist = 0;
  fprintf(arq,"Buckets \n");
  for(i = 0; i < bg->tam_grid;i++){
    for(j = 0; j < bg->tam_grid; j++){
      bucket b = bg->buckets[i][j];
      if (b.itens != NULL){
	fprintf(arq,"  Bucket[%d][%d]:\n",i,j);
	fprintf(arq,"    COLISOES: %lld\n",bg->m_statistic_find[i][j]);
	fprintf(arq,"    QTD DIST EUCLID: %lld\n",bg->m_statistic_eval[i][j]);
	totdist+= bg->m_statistic_eval[i][j];
	fprintf(arq,"    RAIO: %.6f\n",b.raio);
	fprintf(arq,"    Centro:\n      ");
	int k;
	for(k = 0 ; k < bg->num_luzes; k++) {
	  fprintf(arq," [%d]:%f ", k, b.centro[k]);
	}
	int num_itens = 0;
	listapont* p;
	for (p = b.itens; p != NULL; p = p->proximo) { num_itens++;	}
	fprintf(arq,"\n    NUM ITENS: %5d", num_itens);
	fprintf(arq,"\n    ITENS: ");
	for (p = b.itens; p != NULL; p = p->proximo) {
	  fprintf(arq," %d", p->ponteiro);
	}
	fprintf(arq,"\n");
      }
    }
  }

  fprintf(arq, "  BARICENTRO UV: %f %f \n", bg->bu, bg->bv);
  fprintf(arq, "  RAIO DA GRADE: %f\n", bg->R);
  fprintf(arq, "  TOTAL DE DIST EUCLID: %d\n", totdist);
  fflush(arq);
}

void criaMatrizCovariancia2(Tabela*tab, double M[], double bar[])
{
  int num_linhas  = get_num_linhas(tab);
  int num_luzes = get_num_luzes(tab);
  //A matriz P está implicitamente contida na Tabela
  int i,j;
  /* Calculamos o baricentro {Mbar} de todas as entradas: */
  for (j = 0; j < num_luzes; j++) { bar[j] = 0; }
  for(i = 0; i < num_linhas; i++) {
    const double *go = get_intdir(tab,i);
    for(j = 0; j < num_luzes; j++) {  bar[j] += go[j]; }
  }
  for (j = 0; j < num_luzes; j++) { bar[j] /= num_linhas; fprintf(stderr,  "  B[%4d]; %f ", j, bar[j]); }
  fprintf(stderr,  "\n");
  fprintf(stderr,  "Matriz Covariancia - Baricentro Calculado.\n");

  /* Agora calculamos { M[r,s] = SUM{ tab[k,r]*tab[k,s] : k = 0..num_linhas-1 } } */
  int r, s;
  for(r = 0; r< num_luzes;r++){
    for(s = 0; s < num_luzes;s++){
      M[(r*num_luzes)+ s] = 0;
    }
  }
  double dgo[num_luzes];
  for(i = 0; i< num_linhas; i++) {
    const double *go = get_intdir(tab,i);
    for(j = 0; j < num_luzes;j++) {  dgo[j] = go[j] - bar[j]; }
    for(r = 0; r < num_luzes; r++) {
      for(s = 0; s < num_luzes; s++) {
	M[(r*num_luzes)+ s] += dgo[r]*dgo[s];
      }
    }
  }
  fprintf(stderr,  "Matriz Covariancia calculada.\n");
}

void normalizaVetor(double* vet, int tam){
	int i;
	double soma = 0.0001;
	for(i =0 ; i< tam;i++){
		soma = soma + (vet[i] * vet[i]);
	}
        soma = sqrt(soma);
	for(i =0 ; i< tam;i++){
		vet[i] = vet[i]/soma;
	}
}


i2_t mapeiaHash
  ( const double ass[],
    double u[],
    double v[],
    const double baricentro[],
    int num_luzes,
    int tam_grid,
    double R
    )
{
  int N = tam_grid;
  i2_t indice;
  double s[num_luzes];
  int i;
  for (i = 0; i < num_luzes;i++) s[i] = ass[i] - baricentro[i];
  double su = produtoEscalar(s,u,num_luzes);
  double sv = produtoEscalar(s,v,num_luzes);
  assert(fabs(su) <= 1.0001);
  assert(fabs(sv) <= 1.0001);
  indice.c[0] = (int)floor((su/R + 1)*(N/2.0));
  if(indice.c[0] < 0) indice.c[0] = 0;
  if(indice.c[0] > (N-1)) indice.c[0]= N-1;
  indice.c[1] = floor((sv/R + 1)*(N/2.0));
  if(indice.c[1] < 0) indice.c[1] = 0;
  if(indice.c[1] > (N-1)) indice.c[1]= N-1;
  // fprintf(stderr,  "XR: %f %d YR: %f %d \n",su,indice.c[0],sv,indice.c[1]);
  return indice;
}


void showBucketsPPM(char* prefix,bucketGrid* bg, int num_linhas){
    int i,j;
    char nomearq[500];
    /* Acha maximo tamanho de qualquer bucket list: */
    int max_count = 1;
    int min_count = num_linhas;
    for( i= 0; i< bg->tam_grid;i++){
      for(j = 0; j< bg->tam_grid;j++){
        int count = 0;
        listapont *lista;
        for(lista = bg->buckets[i][j].itens; lista != NULL; lista = lista->proximo){ count++; }
        if(count > max_count) max_count = count;
        if(count < min_count) min_count = count;
      }
    }
    float_image_t  *im = float_image_new(3, bg->tam_grid, bg->tam_grid);
    for(i = 0; i < bg->tam_grid; i++){
      for(j = 0; j < bg->tam_grid; j++){
        int count = 0;
        listapont *lista;
        for(lista = bg->buckets[i][j].itens;lista != NULL; lista = lista->proximo){ count++; }
        double rel_count = ((double)count/(double)max_count);
        int c;
        for (c = 0; c < 3; c++) { float_image_set_sample(im, i,j, c, rel_count); }
      }
    }
    sprintf(nomearq,"%sBucketGridSizes.ppm",prefix);
    bool_t isMask = FALSE; /* Assume pixels have a smooth distribution. */
    float_image_write_pnm_named(nomearq, im, isMask, 1.000, 0.000, FALSE, TRUE, FALSE);
    float_image_free(im);
}

void showBucketsEPS(char* prefix,bucketGrid* bg, int num_linhas){
    char nomearq[500];
    sprintf(nomearq,"%sBucketGridSizes",prefix);
    plota_bucket_sizes(nomearq,bg,1);
    //fprintf(stderr,  "Plotou Bucket Sizes");
    sprintf(nomearq,"%sBucketGridRaios",prefix);
    plota_bucket_raios(nomearq,bg,1);
    //fprintf(stderr,  "Plotou Bucket Raios");
    sprintf(nomearq,"%sBucketGridDesvios",prefix);
    plota_bucket_desvios(nomearq,bg,1);
    //fprintf(stderr,  "Plotou Bucket Desvios");
}

void showBucketsPLT(char* prefix,bucketGrid* bg, int num_linhas){
    char nomearq[500];
    sprintf(nomearq,"%sBucketGridSizes",prefix);
    plota_bucket_sizes(nomearq,bg,0);
    //fprintf(stderr,  "Plotou Bucket Sizes");
    sprintf(nomearq,"%sBucketGridRaios",prefix);
    plota_bucket_raios(nomearq,bg,0);
    //fprintf(stderr,  "Plotou Bucket Raios");
    sprintf(nomearq,"%sBucketGridDesvios",prefix);
    plota_bucket_desvios(nomearq,bg,0);
    //fprintf(stderr,  "Plotou Bucket Desvios");
}

void showBucketsRAW(char* prefix,bucketGrid* bg, int num_linhas){
    char nomearq[500];
    sprintf(nomearq,"%sBucketGridSizes",prefix);
    plota_bucket_sizes(nomearq,bg,2);
    //fprintf(stderr,  "Plotou Bucket Sizes");
    sprintf(nomearq,"%sBucketGridRaios",prefix);
    plota_bucket_raios(nomearq,bg,2);
    //fprintf(stderr,  "Plotou Bucket Raios");
    sprintf(nomearq,"%sBucketGridDesvios",prefix);
    plota_bucket_desvios(nomearq,bg,2);
    //fprintf(stderr,  "Plotou Bucket Desvios");
}


void showBucketsData(char* prefix,bucketGrid* bg,int num_linhas){
	char* filename = NULL;
	char *filename = jsprintf("%sBucketData.txt",prefix);
	FILE* arq;
	arq = fopen(filename,"wt");
	if(arq == NULL){
		fprintf(stderr,"showBucketData %s FAILED !\n",filename);
		return;
	}
	int i, j;
	//first plot N
	fprintf(arq,"#N = %d\n",bg->tam_grid);
	//plot number of lights
	fprintf(arq,"#Lights = %d\n",bg->num_luzes);
	//plot Bu,Bv
	fprintf(arq,"#(bu,bv) = %9.6f %9.6f \n",bg->bu,bg->bv);
	//plot radius of grid
	fprintf(arq,"#R = %f\n",bg->R);
	//plot Bucket Size
	fprintf(arq,"#BucketSize = %f\n",(2.0*bg->R)/((float)bg->tam_grid));
	//plot baricenter
	fprintf(arq,"#B = ");
	for(i = 0; i < bg->num_luzes;i++){
		fprintf(arq,"%9.6f ",bg->baricentro[i]);
	}
	fprintf(arq,"\n");
	//plot U
	fprintf(arq,"#U = ");
	for(i = 0; i < bg->num_luzes;i++){
		fprintf(arq,"%9.6f ",bg->u[i]);
	}
	fprintf(arq,"\n");
	//plot V
	fprintf(arq,"#V = ");
	for(i = 0; i < bg->num_luzes;i++){
		fprintf(arq,"%9.6f ",bg->v[i]);
	}
	fprintf(arq,"\n");

	//now plot RAW data
	fprintf(arq,"# i j vmin umin vmax umax entries raio (centro) find scan  eval scanned_h eval_h hashed\n");
	//hashed - number of bucket access
	//scanned - how many times where bucket was examined
	//compr - how many euclidean distances in a bucket
	for(i = 0; i< bg->tam_grid;i++){
    		for(j = 0; j< bg->tam_grid;j++){
			//coordinates of the BG
			fprintf(arq,"%04d    %04d ",i,j);
			//now umin, vmin, umax,vmax
			double umin,vmin,umax,vmax;
			double su = (bg->R)*(((2*i)/(float)(bg->tam_grid)) -1);
			double sv = (bg->R)*(((2*j)/(float)(bg->tam_grid)) -1);
			umin = su;
			umax = su + (2.0*bg->R)/((float)bg->tam_grid);
			vmin = sv;
			vmax = sv + (2.0*bg->R)/((float)bg->tam_grid);
			int tableEntries = 0;
			listapont *listap;
			for(listap = bg->buckets[i][j].itens;listap != NULL; listap = listap->proximo){
				tableEntries++;
			}
			//plot umin,vmin,umax,vmax, ents
			fprintf(arq,"%9.6f %9.6f %9.6f %9.6f %04d ",umin,vmin,umax,vmax, tableEntries);
			//plot raio
			fprintf(arq,"%9.6f ",bg->buckets[i][j].raio);
			int k;
			//plot (raio - centro)
			for(k = 0; k < bg->num_luzes; k++){
				//fprintf(arq,"%9.6f ",bg->buckets[i][j].raio - bg->buckets[i][j].centro[i]);
				fprintf(arq,"%9.6f ",bg->buckets[i][j].centro[k]);
			}
			//plot find scanned eval mapped_h eval_h hashed
			fprintf(arq,"%06lld %06lld %06lld %06lld %06lld %06lld",
				bg->m_statistic_find[i][j],
				bg->m_statistic_scan[i][j],
				bg->m_statistic_eval[i][j],
				bg->m_statistic_scan_h[i][j],
				bg->m_statistic_eval_h[i][j],
				bg->m_statistic_hashed[i][j]
			);
			fprintf(arq,"\n");
		}
	}
	fclose(arq);
}


void plota_bucket_raios(char* filename, bucketGrid* bg,int isEPS){
   /* Start {gnuplot} and set up a pipe {wp->pp} to it: */
    int hSize = 600;
    int vSize = 400;

    char *execCommand = jsprintf("gnuplot -geometry '%dx%d'  >& /dev/null",hSize, vSize);
    FILE *pp;
    if(isEPS == 1){
    	pp = popen(execCommand, "w");
    }else{
	pp = fopen(filename,"wt");
    }
    if (pp == NULL)
      { fprintf(stderr,  "failed starting gnuplot");
        exit(-1);
      }
    /* Initialize fields of {wp} with sensible defaults: */
    if(isEPS != 2){
    	fprintf(pp, "set terminal postscript eps enhanced  \"Times-Roman\" 12\n");

	    fprintf(pp, "set output \"%s.eps\"\n", filename);

    	fprintf(pp,"set hidden3d \n");
    	fprintf(pp,"set view 45,45\n");
    	fprintf(pp,"set size 2,2 \n");
    	fprintf(pp,"set nokey \n");
    	fprintf(pp, "splot '-' with lines\n");
    }
    /* Output the data to be visualized: */
    int i, j;
    for( i= 0; i< bg->tam_grid;i++){
    		for(j = 0; j< bg->tam_grid;j++){
      			double rplot = (bg->buckets[i][j].itens ? bg->buckets[i][j].raio : 0);

			fprintf(pp, "%d\t%d\t%10.6f\n", i,j,rplot);
    		}
		fputc('\n', pp);
   }
   if(isEPS != 2) fputs("e\n", pp);

   fflush(pp);
   if(!(isEPS == 1)) fclose(pp);

}

void plota_bucket_desvios(char* filename, bucketGrid* bg,int isEPS){
  int hSize = 600;
  int vSize = 400;
  int num_luzes = bg->num_luzes;


  char *execCommand = jsprintf("gnuplot -geometry '%dx%d'  >& /dev/null",hSize, vSize);
  FILE *pp ;
  if(isEPS == 1) pp =  popen(execCommand, "w");
  else pp = fopen(filename,"wt");

  if (pp == NULL) { fprintf(stderr,  "failed starting gnuplot"); exit(-1); }
  /* Initialize fields of {wp} with sensible defaults: */
  if(isEPS != 2){
  	fprintf(pp, "set terminal postscript eps enhanced  \"Times-Roman\" 12\n");

  	fprintf(pp, "set output \"%s.eps\"\n", filename);

  	fprintf(pp,"set hidden3d \n");
  	fprintf(pp,"set view 45,45\n");
  	fprintf(pp,"set size 2,2 \n");
  	fprintf(pp,"set nokey \n");
  	fprintf(pp, "splot '-' with lines\n");
  }
  double desloc[num_luzes]; /* Vetor deslocamento do centro do bucket ao plano {u,v}. */
  /* Output the data to be visualized: */
  int i, j,it;
  for (i = 0; i< bg->tam_grid; i++){
    for (j = 0; j < bg->tam_grid; j++) {
      bucket *bitem = &(bg->buckets[i][j]);
      double dplot = 0;

      if (bitem->itens == NULL) continue;
      for(it = 0; it < num_luzes; it++) {
	desloc[it] = bitem->centro[it] - bg->baricentro[it];
      }
      double su,sv;
      su = produtoEscalar(desloc,bg->u,num_luzes);
      sv = produtoEscalar(desloc,bg->v,num_luzes);
      for(it = 0; it < num_luzes; it++) { desloc[it] -= (bg->u[it]*su) + (bg->v[it]*sv); }
      dplot = sqrt(produtoEscalar(desloc,desloc,num_luzes));



      fprintf(pp, "%d\t%d\t%10.6f\n", i,j,dplot);
    }
    fputc('\n', pp);
  }
  if(isEPS != 2) fputs("e\n", pp);

  fflush(pp);
  if(!(isEPS == 1)) fclose(pp);

}

void plota_bucket_sizes(char* filename, bucketGrid* bg,int isEPS){
    /* Start {gnuplot} and set up a pipe {wp->pp} to it: */
    listapont* lista;
    int hSize = 600;
    int vSize = 400;
  // fprintf(stderr,  "ANTES DE PLOTAR");
    char *execCommand = jsprintf("gnuplot -geometry '%dx%d'  >& /dev/null",hSize, vSize);
    FILE *pp;
    if(isEPS == 1)  pp = popen(execCommand, "w");
    else pp = fopen(filename,"wt");

    if (pp == NULL)
      { fprintf(stderr,  "failed starting gnuplot");
        exit(-1);
      }
    if(isEPS != 2){
    	/* Initialize fields of {wp} with sensible defaults: */
     	fprintf(pp, "set terminal postscript eps enhanced  \"Times-Roman\" 12\n");
   	// fprintf(pp,"set terminal x11 \n");
    	fprintf(pp, "set output \"%s.eps\"\n", filename);
    	//mwplot_set_axis_labels(wp, "X", "Y", "Z");
    	//mwplot_set_ranges(wp, NAN, NAN, NAN, NAN, NAN, NAN);
    	fprintf(pp,"set hidden3d \n");
    	fprintf(pp,"set view 45,45\n");
    	fprintf(pp,"set size 2,2 \n");
    	fprintf(pp,"set nokey\n" );
    	fprintf(pp, "splot '-' with lines\n");
    }
    /* Output the data to be visualized: */
    int i, j,count;
    for( i= 0; i< bg->tam_grid;i++){
    		for(j = 0; j< bg->tam_grid;j++){
      			count = 0;
			for(lista = bg->buckets[i][j].itens;lista != NULL; lista = lista->proximo){
					count++;
      				}

      			fprintf(pp, "%d\t%d\t%d\n", i,j,count);

    		}
		fputc('\n', pp);
   }
   if(isEPS != 2) fputs("e\n", pp);
//   fputs("pause -1\n",pp);
   fflush(pp);
   if(!(isEPS == 1)) fclose(pp);
 // fprintf(stderr,  "DEPOIS DE PLOTAR");
	//scanf("%d",&i);

}

void plota_bucket_map(char* filename, bucketGrid* bg,int isEPS){
    /* Start {gnuplot} and set up a pipe {wp->pp} to it: */

    int hSize = 600;
    int vSize = 400;
  // fprintf(stderr,  "ANTES DE PLOTAR");
    char *execCommand = jsprintf("gnuplot -geometry '%dx%d'  >& /dev/null",hSize, vSize);
    FILE *pp;
    if(isEPS == 1)  pp = popen(execCommand, "w");
    else pp = fopen(filename,"wt");

    if (pp == NULL)
      { fprintf(stderr,  "failed starting gnuplot");
        exit(-1);
      }
    if(isEPS != 2){
    	/* Initialize fields of {wp} with sensible defaults: */
     	fprintf(pp, "set terminal postscript eps enhanced  \"Times-Roman\" 12\n");
   	// fprintf(pp,"set terminal x11 \n");
    	fprintf(pp, "set output \"%s.eps\"\n", filename);
    	//mwplot_set_axis_labels(wp, "X", "Y", "Z");
    	//mwplot_set_ranges(wp, NAN, NAN, NAN, NAN, NAN, NAN);
    	fprintf(pp,"set hidden3d \n");
    	fprintf(pp,"set view 45,45\n");
    	fprintf(pp,"set size 2,2 \n");
    	fprintf(pp,"set nokey\n" );
    	fprintf(pp, "splot '-' with lines\n");
    }
    /* Output the data to be visualized: */
    int i, j;
    for( i= 0; i< bg->tam_grid;i++){
    		for(j = 0; j< bg->tam_grid;j++){
      			fprintf(pp, "%d\t%d\t%lld\n", i,j,bg->circle_map[i][j]);

    		}
		fputc('\n', pp);
   }
   if(isEPS != 2) fputs("e\n", pp);
//   fputs("pause -1\n",pp);
   fflush(pp);
   if(!(isEPS == 1)) fclose(pp);
 // fprintf(stderr,  "DEPOIS DE PLOTAR");
	//scanf("%d",&i);

}



void showBucketGridStatistic(char* prefix,bucketGrid* bg){
	int i,j;
	float_image_t *  im;
	im = float_image_new(3, bg->tam_grid, bg->tam_grid);
	char nomearq[500];
        for( i = 0; i < bg->tam_grid; i++){
		for(j = 0; j < bg->tam_grid; j++){
			double rel_count = ((double)bg->m_statistic_find[i][j]/(double)bg->max_count_st)*255.0;
			if(rel_count < 0) rel_count = 0;
      			float_image_set_sample(im,0,i,j,rel_count);
                        float_image_set_sample(im,1,i,j,rel_count);
                        float_image_set_sample(im,2,i,j,rel_count);
		}
	}
	sprintf(nomearq,"%sBucketStatistic.ppm",prefix);
        //salva_arquivo("BucketGrid.ppm",im);
        bool_t isMask = TRUE; /* Map 1 to 0, 0 to {maxval}. */
        float_image_write_pnm_named(nomearq, im, isMask, 1.000, 0.000, FALSE, TRUE, FALSE);
	//salva_arquivo("BucketStatistic.ppm",im);
}

int tam_grid_bucket(int num_linhas){

	FILE* arq;
	int saida;
	double razao;
	double result;
	arq = fopen("out/BucketSize.txt","rt");
	// R = B*B/T
	if(arq == NULL){
		return sqrt(num_linhas)*2;
	}
	//AGORA LEMOS A RAZÃO
	fscanf(arq,"%lf",&razao);
	if(razao <= 0) {
		saida= 1;
	}
	else{
		//result = sqrt(num_linhas/razao);
		result = sqrt(num_linhas/razao);
		saida = (int)floor(result + 0.5);
	}
	fclose(arq);
	fprintf(stderr,  "Tamanho do Bucket Grid: %d\n",saida);
	return saida;



//return sqrt(num_linhas)*2;
}


double cobrinha ( double k){ //MUDAR ESTE NOME !!!
		double maior;
		maior = k;
		if( k < 0) maior = -k;
		maior--;
		if(maior < 0) maior = 0;
		return maior;
}

double calculaDistanciaBucketMinima(int i, int j, int iq, int jq, double tam_bucket){
	double resultado;


	resultado = tam_bucket*sqrt( (cobrinha(i - iq)*cobrinha(i - iq) ) + (cobrinha(j - jq)*cobrinha(j-jq)) );
	return resultado;
}


void criaListaPreordenada(bucketGrid* bg){
	int R,i,j,cont;
	R = bg->tam_grid;
	listapos* vetor_posicoes;
	vetor_posicoes = (listapos*) (malloc(sizeof(listapos)*(R*R*4)));
	if(vetor_posicoes == NULL){
		fprintf(stderr,  "ERRO - Grade não alocada !!!");
	}
	cont = 0;
        for (i = -R; i < R; i++){
		for(j = -R; j < R; j++){
			vetor_posicoes[cont].posicao[0] = i;
			vetor_posicoes[cont].posicao[1] = j;
			vetor_posicoes[cont].dist = sqrt(cobrinha(i)*cobrinha(i)  + cobrinha(j)*cobrinha(j));
		//	fprintf(stderr,  "Dist %d %d %f\n", i,j,vetor_posicoes[cont].dist);
			cont++;
		}
	}
	fprintf(stderr,  "GRADE : %d Cont:%d\n", R*R*4, cont);
	fprintf(stderr,  "Antes da ordenação !!!!");

	auto int compara_lista (const void * a, const void * b);

	int compara_lista (const void * a, const void * b)
	{
  		listapos posA = *(listapos*)a;
		listapos posB = *(listapos*)b;
		//distA = sqrt(cobrinha(posA.posicao[0])*cobrinha(posA.posicao[0])  + cobrinha(posA.posicao[1])*cobrinha(posA.posicao[1]));
		//distB = sqrt(cobrinha(posB.posicao[0])*cobrinha(posB.posicao[0])  + cobrinha(posB.posicao[1])*cobrinha(posB.posicao[1]));
		//se a > b retorna 1
		// se a = b retorna 0
		// se a < b retorna -1
		if( posA.dist > posB.dist){
			return 1;
		}
		if(posA.dist < posB.dist){
			return -1;
		}
		//desempate pela dist euclideana
		double dA,dB;
		dA = sqrt(posA.posicao[0]*posA.posicao[0] + posA.posicao[1]*posA.posicao[1]);
		dB = sqrt(posB.posicao[0]*posB.posicao[0] + posB.posicao[1]*posB.posicao[1]);
		if(dA > dB ) return 1;
		if(dB > dA ) return -1;

		return  0;

	}
	qsort(vetor_posicoes,(size_t)(R*R*4),sizeof(listapos),(void*)compara_lista);
	fprintf(stderr,  "Lista Preordenada Ordenada !");
	/*for(i = 0; i < 40; i++){
		fprintf(stderr,  " %d %d %f\n", vetor_posicoes[i].posicao[0],  vetor_posicoes[i].posicao[1],vetor_posicoes[i].dist);
	}*/
	bg->lista_preordenada = NULL;
	bg->tamanho_lista_preordenada = R*R*4;
	/*for(i = (R*R*4)-1; i >= 0 ; i--){
		listapos* grid	= (listapos*) malloc(sizeof(listapos));
		grid->posicao[0] = vetor_posicoes[i].posicao[0];
		grid->posicao[1] = vetor_posicoes[i].posicao[1];
		grid->proximo = bg->lista_preordenada;
		bg->lista_preordenada = grid;
	}*/
	bg->lista_preordenada = vetor_posicoes;

	//free(vetor_posicoes);

}

int localiza_normal_hash(bucketGrid* bg, Tabela* tab, double SO[], double* dist, double *albedo,int* n_euclid_evalsP, int* n_scansP)
{
  int num_luzes = get_num_luzes(tab);
  double tam_bucket = 2.0/(double)bg->tam_grid;

  /* Normaliza o vetor de observações {SO} obtendo a assinatura {so}: */
  double so[num_luzes];
  double Smag;
  extrai_assinatura(SO,so,&Smag,num_luzes);

  /* Mapeia {so} para um bucket inicial da grade: */
  i2_t indice = mapeiaHash(so, bg->u, bg->v, bg->baricentro, num_luzes, bg->tam_grid, bg->R);
  if(n_euclid_evalsP != NULL) *n_euclid_evalsP = 0;
  if(n_scansP != NULL) *n_scansP = 0;
  int x = indice.c[0];
  int y = indice.c[1];
  int xmelhor = -1;
  int ymelhor = -1;
  double menor_distancia = INF;
  int i;
  int melhor_linha = -1;
   if(bg->circle_map[x][y] == 0){
	//printf("ERROR\n\n");
	 //return -1;
   }
   bg->m_statistic_hashed[x][y] = bg->m_statistic_hashed[x][y] +1;
  /*fprintf(stderr,  "ASSINATURA so:");
    for (i = 0 ; i < num_luzes; i++){ fprintf(stderr,  " %f", so[i]); }
    fprintf(stderr,  "\n");*/

  /* Varre buckets em ordem de distância crescente: */
  for(i = 0; i < bg->tamanho_lista_preordenada; i++){
    int xtemp = x + bg->lista_preordenada[i].posicao[0];
    int ytemp = y + bg->lista_preordenada[i].posicao[1];

    if((xtemp < 0) || (ytemp < 0) || (xtemp >= bg->tam_grid) || (ytemp >= bg->tam_grid)){
      /* Saiu do bucket grid, então pula este bucket */
      continue;
    }

    /* Olhamos para este bucket: */
    (bg->m_statistic_scan[xtemp][ytemp])++;
    (bg->m_statistic_scan_h[x][y])++;
     if(n_scansP != NULL)  (*n_scansP)++;

    /* Teste para parar de procurar: */
    double bucket_dist = bg->lista_preordenada[i].dist; /* Dist deste bucket ao bucket original. */
    if((melhor_linha != -1) && ((bucket_dist * tam_bucket) > menor_distancia)) {
      /* A partir daqui não vale a pena buscar: */
      break;
    }

    /* Procura no bucket: */
    bucket *bitem = &(bg->buckets[xtemp][ytemp]);
    int num_euclid = 0;
    int melhor_no_bucket = -1;
    if(bitem->itens != NULL){
    	melhor_no_bucket = procura_no_bucket(so, bitem, tab, &menor_distancia, &num_euclid);
    }
    bg->m_statistic_eval[xtemp][ytemp] = bg->m_statistic_eval[xtemp][ytemp] + num_euclid;
    bg->m_statistic_eval_h[x][y] = bg->m_statistic_eval_h[x][y] + num_euclid;
    if(n_euclid_evalsP != NULL) *n_euclid_evalsP += num_euclid;
    if (melhor_no_bucket != -1) {
      melhor_linha = melhor_no_bucket;
      xmelhor = xtemp;
      ymelhor = ytemp;
    }

  }

  if (melhor_linha != -1) {
    /* Achamos mais um neste bucket: */
    assert(xmelhor != -1);
    assert(ymelhor != -1);
    //fprintf(stderr,  "Melhor %d %d\n", xmelhor, ymelhor);
    bg->m_statistic_find[xmelhor][ymelhor] =  bg->m_statistic_find[xmelhor][ymelhor] + 1;
    if(bg->m_statistic_find[xmelhor][ymelhor] > bg->max_count_st){
      bg->max_count_st = bg->m_statistic_find[xmelhor][ymelhor];
    }
  }
  /* Devolve a distância entre as assinaturas: */
  (*dist) = menor_distancia;
  if ((*dist) > 2.0) { (*dist) = 2.0; }
  (*albedo) = calcula_albedo(tab, melhor_linha, SO);
  return melhor_linha;
}

int procura_no_bucket(double so[], bucket* bitem, Tabela *tab, double *menor_distancia, int *num_euclidP)
{
  int num_luzes = get_num_luzes(tab);
  if (bitem->itens == NULL) { /* O bucket está vazio, passe adiante: */ return -1; }
  int melhor_linha = -1; /* Melhor linha encontada neste bucket, ou -1. */

  double raio = bitem->raio;
  double distancia_media = dist_euclid(so, bitem->centro, num_luzes);
  (*num_euclidP)++;
  listapont *listap;
  for (listap = bitem->itens; listap != NULL; listap = listap->proximo) {
    /* Teste para parar e procurar neste bucket: */
    if ((*menor_distancia)  < (distancia_media - raio )) {
      /* A distância de qualquer ponto no bucket é grande demais para valer a pena: */
      return melhor_linha;
    }

    /* Verifica item: */
    int indice_atual = listap->ponteiro;
    const double *go = get_intdir(tab, indice_atual);
    double distancia = dist_euclid(so,go,num_luzes);
    (*num_euclidP)++;
    if (distancia < (*menor_distancia)) {
      (*menor_distancia) = distancia;
      melhor_linha = indice_atual;
    }
  }
  return melhor_linha;
}

int get_tam_grid(bucketGrid* bg){
	return (bg == NULL ? 0 : bg->tam_grid);
}


void LiberaBucketGrid(bucketGrid* bg){
	free(bg->u);
	free(bg->v);
	free(bg->baricentro);

	int i ;
	for(i = 0; i < bg->tam_grid	; i++){
		free(bg->m_statistic_find[i]);
		free(bg->m_statistic_eval[i]);
		free(bg->m_statistic_eval_h[i]);
		free(bg->m_statistic_scan[i]);
		free(bg->m_statistic_scan_h[i]);
		free(bg->m_statistic_hashed[i]);
		free(bg->circle_map[i]);
	}
	free(bg->m_statistic_find);
	free(bg->m_statistic_eval);
	free(bg->m_statistic_eval_h);
	free(bg->m_statistic_scan);
	free(bg->m_statistic_scan_h);
	free(bg->m_statistic_hashed);
	free(bg->circle_map);

	free(bg->lista_preordenada);



	int j;
	for(i = 0; i < bg->tam_grid	; i++){
		for(j = 0; j < bg->tam_grid	; j++){
			free(bg->buckets[i][j].centro);
			listapont* aux = bg->buckets[i][j].itens;
			while(aux != NULL){
				listapont* temp = aux;
				aux = aux->proximo;
				free(temp);
			}
		}
		free(bg->buckets[i]);
	}
	free(bg->buckets);


}