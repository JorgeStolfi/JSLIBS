
#ifndef hash_H
#define hash_H

#include <tabela.h>
#include <stdint.h>
#include <float_image.h>

typedef struct MatrizCovariancia matcov;
typedef struct BucketGrid bucketGrid;

bucketGrid* CriaBucketGrid(Tabela* tab,int gridzise);
void criaListaPreordenada(bucketGrid* bg);
void showBucketsPPM(char* prefix,bucketGrid* bg, int num_linhas);
void showBucketsEPS(char* prefix,bucketGrid* bg, int num_linhas);
void showBucketsPLT(char* prefix,bucketGrid* bg, int num_linhas);
void showBucketsRAW(char* prefix,bucketGrid* bg, int num_linhas);
void showBucketsData(char* prefix,bucketGrid* bg,int num_linhas);
void showBucketGridStatistic(char* prefix,bucketGrid* bg);

int localiza_normal_hash(bucketGrid* bg, Tabela* tab, double SO[], double* dist, double *albedo,int* n_euclid_evalsP, int* n_scansP);
/* Localiza a entrada na tabela {tab} cujo vetor de observações é mais
  similar ao vetor de observações dado {SO[0..n-1]} (a menos de um
  fator de escala); onde {n} é o número de luzes.  Devolve o índice
  dessa entrada como resultado, e coloca em {*dist} a distância enrte
  as assinatura de {SO} e a dessa entrada. Usa hash bucket grid para 
  acelerar a busca. Se {n_scansP} e/ou {n_euclid_evalsP} 
  forem não NULL, devolve número de buckets escaneados em {n_scansP} e 
  número de distâncias euclidianas calculadas em {n_euclid_evalsP} */

void flushBucketGrid(FILE* arq, bucketGrid* bg);
int get_tam_grid(bucketGrid* bg);

void plota_bucket_sizes(char* filename, bucketGrid* bg,int isEPS);

void plota_bucket_desvios(char* filename, bucketGrid* bg,int isEPS);
  /* Plota a distancia entre o centroide do bucket  e o plano  da grade. */

void plota_bucket_map(char* filename, bucketGrid* bg,int isEPS);

void plota_bucket_raios(char* filename, bucketGrid* bg,int isEPS);

int64_t** acessaMatriz_Statistic_Euclid(bucketGrid* bg);
int64_t** acessaMatriz_Statistic_Scan(bucketGrid* bg);

void LiberaBucketGrid(bucketGrid* bg);

#endif
