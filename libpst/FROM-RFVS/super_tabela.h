#ifndef super_tabela_H
#define super_tabela_H

#include <tabela.h>
#include <hash.h>
#include <normais.h>
#include <float_image.h>
#include <r2.h>
#include <r3.h>


typedef struct SuperTabela SuperTabela;

/*Output Functions of ST*/

Tabela* stGetTable(SuperTabela* st);
/*Returns the stored table for SuperTable {st}*/
bucketGrid* stGetBucketGrid(SuperTabela* st);
/*Returns the stored Bucket Grid for SuperTable {st}*/
int stGetK(SuperTabela* st);
/*Return number of light sources used in SuperTable {st}*/
int stGetNumLuzes(SuperTabela* st);
/*Return number of light sources of test subject of SuperTable {st} */
const int* stGetIndLuz(SuperTabela* st);
/*Return the vector contaning indices for the ligh sources used in SuperTable {st}
 The size of the vector is aways equal of the number got by stGetK*/


/*Creation and Destruction functions of ST strtucture*/

SuperTabela* stCreate(
    int num_luzes, 
    int canal, 
    int subsetSize,
    int ind_luz[],
    int gridsize,
    Tabela* mainTable
);
/*Creates and initializate every internal structure inside a SuperTable structure and return a pointer to it
  It needs a array {G} of {num_luzes} gauge images of the original test subject to build internal Table and a  mask image {M} 
  to avoid interpolation of unwanted pixels of gauge images. It uses only intensities of the selected channel {canal} and
   only images contained within the array {ind_luz} with size {k}. The function needs also a list of sampling points
  {ponto[0..num_pontos-1]} with their respective surface normals  {normal[0..num_pontos-1]}. 
   The internal Bucket Grid is built with size {gridsize} , when {gridsize} have value -1 the size of the grid
  will be automaticaly set to SQRT(L)*2 where L is the number of lines of the internal Table
  
*/

void stFree(SuperTabela** st);
/* Frees an SuperTable structure ... IMPORTANT, it wont free internal Table and Bucket Grid.*/

int stSearchNormal(
    SuperTabela* st,
    const double SO[], 
    double *SO_mag,
    double *dist, 
    double *albedo,
    int  *matchedIndex,
    int* n_euclid_evalsP,
    int* n_scansP
);
/*Realizes an search using dist_aplha using the inner table and bucket grid inside the SuperTable {st}
   TODO - Not finished yet, returns {dist} = eucliden distance of best matching signatures, {SO_mag} = 
  magnitude of scene observation vector restricted to the selected lights,
   Se {n_scansP} e/ou {n_euclid_evalsP} 
  forem não NULL, devolve número de buckets escaneados em {n_scansP} e 
  número de distâncias euclidianas calculadas em {n_euclid_evalsP} */

void stGenerateSubsets(int** subsets,int subset_size,int num_luzes,int num_subsets);
void stGenerateSubsets_deprecated(int** subsets,int subsetSize,int m,int M);
void LiberaSuperTabela( SuperTabela* st);
#endif
