#define _GNU_SOURCE
#include <stdio.h>
#include <stdlib.h>
#include <hash.h>
#include <tabela.h>
#include <r3.h>
#include <i2.h>
#include <rn.h>
#include <float_image.h>
#include <assert.h>
#include <jsfile.h>

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

i2_t mapeiaHash
  ( const double ass[], 
    double u[], 
    double v[], 
    const double baricentro[],
    int num_luzes,
    int tam_grid,
    double R
    );

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

float_image_t* generate_hash_normal_map(bucketGrid* bg,Tabela* tab){
  int NX,NY;
  NX = NY = bg->tam_grid;
  float_image_t* nm = float_image_new(3,NX,NY);
  int x,y;
  for(x = 0; x < NX; x++){
    for(y = 0; y < NY; y++){
      bucket b = bg->buckets[x][y];
      r3_t avg_normal = (r3_t){{0,0,0}};
      int count = 0;
      listapont* l;
      for(l = b.itens; l != NULL; l = l->proximo){
	r3_t tab_normal = get_normal(tab,l->ponteiro);
	r3_add(&avg_normal,&tab_normal,&avg_normal);
	count++;
      }
      if(count != 0){
	r3_scale(1.0/(float)count,&avg_normal,&avg_normal);
      }
      float_image_set_sample(nm,0,x,y,avg_normal.c[0]);
      float_image_set_sample(nm,1,x,y,avg_normal.c[1]);
      float_image_set_sample(nm,2,x,y,avg_normal.c[2]);
    }
  }
  return nm;
}

float_image_t** generate_hash_sector(bucketGrid* bg,double* SO,double interval){
  assert(interval > 0);
  int NL = bg->num_luzes;
  int NX,NY;
  NX = NY = bg->tam_grid;
  float_image_t** sec_images = (float_image_t**)malloc(sizeof(float_image_t*)*(NL+1));
  int l;
  for(l = 0; l < NL; l++){
    sec_images[l] = float_image_new(1,NX,NY);
    int ix,iy;
    for(ix = 0; ix < NX; ix++){
      for(iy = 0; iy < NY; iy++){
	float_image_set_sample(sec_images[l],0,ix,iy,-1);
      }
    }
    double pos;
    for(pos = 0;pos <= 1; pos+=interval){
      double SO_ex[NL];
      double so[NL];
      rn_copy(NL,SO,SO_ex);
      SO_ex[l] = pos;
      rn_dir (NL,SO_ex,so);
      i2_t indice = mapeiaHash(so, bg->u,bg->v,bg->baricentro,NL,bg->tam_grid,bg->R);
      int x = indice.c[0];
      int y = indice.c[1];
      float_image_set_sample(sec_images[l],0,x,y,pos);
    }
  }
  
  sec_images[NL] = float_image_new(1,NX,NY);
  int ix,iy;
  for(ix = 0; ix < NX; ix++){
    for(iy = 0; iy < NY; iy++){
      double max = -1;
      int i;
      for(i = 0; i < NL;i++){
	double val = float_image_get_sample(sec_images[i],0,ix,iy);
	if(val > max) max = val;
      }
      float_image_set_sample(sec_images[NL],0,ix,iy,max);
    }
  }
  
  i2_t indice = mapeiaHash(SO, bg->u,bg->v,bg->baricentro,NL,bg->tam_grid,bg->R);
  float_image_set_sample(sec_images[NL],0,indice.c[0],indice.c[1],2);
  
  return sec_images;
}

void SaveFNI(char* filename, float_image_t* img){
  FILE* arq = open_write(filename,TRUE);
  float_image_write(arq,img);
  fclose(arq);
}

int main(int argc,char** argv){
  if(argc < 4){
    fprintf(stderr,"program usage:\nbucket_grid_map <in_file> <out_prefix> <res>");
    return 1;
  }
  char* tab_filename = argv[1];
  char* out_prefix = argv[2];
  int res;
  sscanf(argv[3],"%d",&res);
  Tabela* tab = NULL;
  LoadTable(tab_filename,&tab);
  bucketGrid* bg = CriaBucketGrid(tab,res);
  assert(bg != NULL);
  float_image_t* normal_map = generate_hash_normal_map(bg,tab);
  char* normal_map_filename = NULL;
  char *normal_map_filename = jsprintf("%s-N.fni",out_prefix);
  SaveFNI(normal_map_filename,normal_map);
  showBucketsData(out_prefix,bg,get_num_linhas(tab));
  
//   double sig[] = {0.3822187,0.3127932,0.2915799,0.2955589,0.2423424,0.3001726,0.2546457,0.2541331,0.2315283,0.2179556,0.3014908,0.3367895 };
//   //double sig[] = {0.3585156, 0.2812013 , 0.2978339 ,0.3427969, 0.3261642 , 0.3491941 ,0.2360555 , 0.2418313 , 0.2183628 , 0.1916409 , 0.2555029 , 0.3065341 };
//   float_image_t** sec_maps =  generate_hash_sector(bg,sig,0.01);
//   int i;
//   for ( i = 0 ; i  < 13; i++){
//     normal_map_filename = NULL;
//     char *normal_map_filename = jsprintf("%s_sector_%02d.fni",out_prefix,i);
//     SaveFNI(normal_map_filename,sec_maps[i]);
//   }
  return 0;
}

