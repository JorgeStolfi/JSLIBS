#ifndef FAST_HASH_H
#define FAST_HASH_H


#include <r3.h>
#include <tabela.h>

struct Fast_Hash_t{
    int N; /*Number of cells on the grid*/
    int nLights;
   // int dimensions; /*number of simensions of the array and so on*/
    r3_t** hash_table; /*Stores the normals for each cell*/
    double** albedo_table; /*stores the albedo for each cell*/
    double** weight_table; /*stores the albedo for each cell*/
    double* u; /* Direção principal da nuvem de pontos {[0..num_luzes-1]}. */
    double* v; /* Direção principal da nuvem de pontos {[0..num_luzes-1]}. */
    double bu, bv; /* Posição do baricentro projetado em {u,v} */
    double R; /* Raio da grade de buckets no plano {u,v}. */
    double* baricenter; /* Centro da grade de buckets. */
};


typedef struct Fast_Hash_t fast_hash_t ;

fast_hash_t* create_fasthash(int N,Tabela* tab,double sigma,int normal_degree, int albedo_degree);
r3_t fast_hash_compute_normal( fast_hash_t* fh, const double SO[], double sigma,double omg0,double omg1, double *albedo,double* logProb );
void SaveFastHash(FILE* arq, fast_hash_t* fh);
fast_hash_t* LoadFastHash(FILE* arq);
float_image_t* FastHashAlbedoToFNI(fast_hash_t* fh);
float_image_t* FastHashNormalsToFNI(fast_hash_t* fh);
float_image_t* FastHashWeightsToFNI(fast_hash_t* fh);
void PrintFastHash(FILE* arq, fast_hash_t* fh);

#endif