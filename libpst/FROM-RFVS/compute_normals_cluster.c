
#define _GNU_SOURCE
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <float_image.h>
#include <float_pnm_image_io.h>
#include "normais.h"
#include "imagem_vetores.h"
#include <least_squares_nd.h>
#include <jsfile.h>
#include "tabela.h"
#include "hash.h"
#include "time.h"
#include <r3.h>
#include <r2.h>
#include "r3x3.h"
#include <rmxn.h>
#include <argparser.h>
#include <assert.h>
#include <string.h>
#include <sys/times.h>
#include <unistd.h>
#include <normals_weights.h>

/*
Lê {o->nLights} imagens de uma cena lambertiana cada uma com uma luz pontal no infinito com direção conhecida,
com  a mesma intensidade fixa sem luz ambiente ou secundária, tais que em cada pixel há pelo menos 3 imagens
onde o pixel está iluminado. Determina os mapas de normais e de albedo.
*/

#define PROG_NAME "compute_normals_cluster"


#define PROG_HELP \
  PROG_NAME " \\\n" \
  "	-prefix {FILE_PREFIX} \\\n" \
  "	-channels {RGB} \\\n" \
  "	-nLights {NLIGHTS} \\\n" \
  "	-nClusters {NLIGHTS} \\\n" \
  "	-sceneImages {Image0 Image1...} \\\n" \
  "     -lightDirections {file0 file1..}   \\\n" \
  "     -clusterDirections {file0 file1..}   \\\n" \
  "	-clusterr {r} \\\n" \
  "     [-K {K Glossiness} ]\\\n" \
  "     [-Amax {Amax} ] \\\n" \
  "     [-usingBestThree ] \\\n" \
  "	[- gamma {GAMMA} ] \\\n" \
  "    " argparser_help_info_HELP ""

#define PROG_INFO \
  "NAME\n" \
  "  Etc. etc..\n" \
  "\n" \
  "SYNOPSIS\n" \
  PROG_HELP "\n" \
  "\n" \
  "OPTIONS" \
  "  Etc. etc.."


#define PROGRESS_METER 0


struct transform_data_t{
  r2_t offset;
  r2_t scale;
};

typedef struct  transform_data_t transform_data_t;

struct options_t{
  char* prefix;
  char* channels;
  int nLights;
  int nClusters;

  char** sceneImages;
  char** lightDirections;
  char** clusterDirections;
  double gamma;

  double Amax;
  bool_t usingBestThree;
  double clusterr;
  double K;
  
   
  };

typedef struct options_t options_t;

void besta(void);

void besta(void) { 
  return;
}


  
double* ConstroiLeastSquares(r3_t dir_luz[] , int n);
double* ConstroiLeastSquares(r3_t dir_luz[] , int n){
   double* LS = rmxn_alloc(3,n);
   double* L = rmxn_alloc(n,3);
   double* LT = rmxn_alloc(3,n);
   int i, j;
   for(i = 0; i < n; i++){
     
     for(j = 0; j < 3; j++){
      int ind = j + (i*3);
      int indT = i + (j*n);
      L[ind] = dir_luz[i].c[j];
      LT[indT] = dir_luz[i].c[j];
     }
   }
   
   double* inv = rmxn_alloc(3,3);
   rmxn_mul (3,n,3, LT, L,inv);
   rmxn_inv(3,inv,inv);
   rmxn_mul (3,3,n,inv, LT,LS);
   free(L);
   free(LT);
   free(inv); 
   return LS;
}

bool_t TestCoplanar(r3_t v0, r3_t v1, r3_t v2);
bool_t TestCoplanar(r3_t v0, r3_t v1, r3_t v2){
  r3x3_t A;
	A.c[0][0] = v0.c[0];
	A.c[0][1] = v0.c[1];
	A.c[0][2] = v0.c[2];
	
	A.c[1][0] = v1.c[0];
	A.c[1][1] = v1.c[1];
	A.c[1][2] = v1.c[2];
	
	A.c[2][0] = v2.c[0];
	A.c[2][1] = v2.c[1];
	A.c[2][2] = v2.c[2];
	double det = r3x3_det(&A);
	if(det == 0) return TRUE;
	return  FALSE;
}

void ShowCoplanar(r3_t* lightDirections, int n);
void ShowCoplanar(r3_t* lightDirections, int n){
  int i;
  fprintf(stderr,"Light Directions\n");
  for(i = 0; i < n; i++){
    fprintf(stderr,"[%d]: ",i);
    r3_gen_print (stderr,&(lightDirections[i]), "%3.4lf","( "," , "," )\n");
  }
  int j,k;
  for(i = 0; i < n; i++){
    for(j = i+1 ; j < n; j++){
      for(k = j+1; k < n; k++){
	r3x3_t A;
	A.c[0][0] = lightDirections[i].c[0];
	A.c[0][1] = lightDirections[i].c[1];
	A.c[0][2] = lightDirections[i].c[2];
	
	A.c[1][0] = lightDirections[j].c[0];
	A.c[1][1] = lightDirections[j].c[1];
	A.c[1][2] = lightDirections[j].c[2];
	
	A.c[2][0] = lightDirections[k].c[0];
	A.c[2][1] = lightDirections[k].c[1];
	A.c[2][2] = lightDirections[k].c[2];
	
	double det = r3x3_det(&A);
	fprintf(stderr,"[%d,%d,%d]: - %lf ",i,j,k,det);
	if(det == 0){
	  fprintf(stderr,"COPLANAR");
	}
	fprintf(stderr,"\n");
	
      }
    }
  }
  
  
}



void FindNLargest(double X[], int nX, int I[], int nI);
  /* Reurns in {I[0..nI-1]} the indies of the {nI} largest elements
    of {X[0..nX-1]}.  Assumes that {nI <= nX}. */

void FindNLargest(double X[], int nX, int I[], int nI) {
        int m = 0;    // As {m} maiores são {I[0..m-1]}.
	int i;
	for(i = 0; i < nX; i++) {
          /* Determina o lugar {j+1} de {i} na lista {I[0..m-1]}, empurrando os demais: */
          int j;
          for (j = m-1; (j >= 0) && (X[i] > X[I[j]]); j--){
            /* Empurra {I[j]} para cima, ou joga fora: */
            if (j < nI-1) { I[j+1] = I[j]; }
          }
          /* Agora {I[j+1]} é o lugar certo para {i}: */
          I[j+1] = i;
          /* Ganahmos mais um, a menos que já tenhamos o suficiente: */
          if (m < nI) { m++; }
	}
}

void Find3Best(double X[], int nX, int I[], int nI,r3_t luz_dir[]) ;
void Find3Best(double X[], int nX, int I[], int nI,r3_t luz_dir[]) {
  double SO[nX];
  rn_copy(nX,X,SO);
  int tries = 0;
  do{
    FindNLargest(SO, nX, I, nI);
    r3_t d0 = luz_dir[I[0]];
    r3_t d1 = luz_dir[I[1]];
    r3_t d2 = luz_dir[I[2]];
    if(!TestCoplanar(d0,d1,d2) ){
      break;
    }else{
      SO[I[0]] = 0; //remove one of the coplanar lights
    }
    tries++;
  }while(tries < (nX -3) );
//   if(tries >= (nX  - 3)){
//     fprintf(stderr,"NO MORE TRIES\nNO MORE TRIES\n");
//   }
}

void ConstroiSistema(double SO[], r3_t luz_dir[], int num_luzes, r3x3_t *A, r3_t *b, double *peso, int imax[]);
/* Examina o vetor de observação {SO[0..num_luzes-1]}, supostamente gerado por 
  luzes pontuais no infinito nas direções {luz_dir[0..num_luzes-1]}.
  Escolhe as três luzes que deixam o pixel mais brilhante, e monta o sistema {*A,*b}  com essas luzes.  Devolve 
  em {*peso} um coef de confiabilidade calculado a partir do determinante do sistema.
  Devolve em {imax[0..2]} os índices das três luzes escolhidas. */

void ConstroiSistema(double SO[], r3_t luz_dir[], int num_luzes, r3x3_t *A, r3_t *b, double *peso, int imax[]){
  //FindNLargest(SO, num_luzes, imax, 3);
  Find3Best(SO, num_luzes, imax, 3,luz_dir);
  /* Verifica se as 3 maiores são todas positivas: */
  if (SO[imax[2]] <= 0.0) {
    // Não temos luzes suficientes para o cálculo:
    (*peso) = 0.0;
    r3x3_zero(A);
    r3_zero(b);
    return;
  }
  /* Preenchemos o sistema com dados das 3 melhores luzes: */
  int col, row;
  for (row = 0; row < 3; row++) {
    int i = imax[row];
    for (col = 0; col < 3; col++) {
      A->c[row][col] = luz_dir[i].c[col];
    }
    b->c[row] = SO[i];
  }
  double smalldet = 1.0e-7;
  double det = fabs(r3x3_det(A));
  (*peso) = 1.0 - exp(-det/smalldet);
}


void computeGO(double* GO,r3_t snp,r3_t* luz_dir,int num_luzes);
void computeGO(double* GO,r3_t snp,r3_t* luz_dir,int num_luzes){
  int i;
  for(i = 0; i < num_luzes; i++){
     double val = r3_dot(&(luz_dir[i]),&snp);
     if(val < 0 ) val = 0;
     GO[i] = val;
  }
}

double computeProb(double* SO,double* GO, int num_luzes);
double computeProb(double* SO,double* GO, int num_luzes){
  double so[num_luzes];
  double go[num_luzes];
  double Smag = rn_dir(num_luzes,SO,so);
  rn_dir(num_luzes,GO,go);
  double dist_eu  = dist_euclid(so,go,num_luzes);
  double sigma = 0.2;
  
  double csi = (dist_eu*Smag/sigma);
  double prob = -0.5*(csi*csi) + 2.0*log(Smag);
  
  
  return prob;
}


double user_cpu_time_usec(void);

double user_cpu_time_usec(void){
  struct tms buf;
  (void)times(&buf);
  return(1000000.0 * ((double) buf.tms_utime)/((double)sysconf(_SC_CLK_TCK)));
}


options_t* parse_args(int argc, char** argv);
options_t* parse_args(int argc, char** argv){
  options_t* o = (options_t*)malloc(sizeof(options_t));
  argparser_t *pp = argparser_new(stderr, argc, argv);
  argparser_set_help(pp, PROG_HELP);
  argparser_set_info(pp, PROG_INFO);
  argparser_process_help_info_options(pp);
  
  argparser_get_keyword(pp, "-prefix");
  o->prefix = argparser_get_next(pp);
  
  argparser_get_keyword(pp, "-channels");
  o->channels = argparser_get_next(pp);

  argparser_get_keyword(pp, "-nLights");
  o->nLights = argparser_get_next_int(pp,3,10000);

  argparser_get_keyword(pp, "-nClusters");
  o->nClusters = argparser_get_next_int(pp,3,10000);

  
  argparser_get_keyword(pp, "-sceneImages");
  o->sceneImages = (char**)malloc(sizeof(char*)*(o->nClusters)*(o->nLights));
  int i;
  for(i = 0; i < (o->nLights)*(o->nClusters); i++){
    o->sceneImages[i] = argparser_get_next(pp);
  }
  
  o->lightDirections = NULL;
  argparser_get_keyword(pp, "-lightDirections");
  o->lightDirections = (char**)malloc(sizeof(char*)*(o->nClusters)*(o->nLights));
  for(i = 0; i < (o->nClusters* o->nLights) ; i++){
    o->lightDirections[i] = argparser_get_next(pp);
  }
 
  o->clusterDirections = NULL;
  argparser_get_keyword(pp, "-clusterDirections");
  o->clusterDirections = (char**)malloc(sizeof(char*)*(o->nClusters));
  for(i = 0; i < o->nClusters; i++){
    o->clusterDirections[i] = argparser_get_next(pp);
  }
  
  argparser_get_keyword(pp, "-clusterr");
  o->clusterr = argparser_get_next_int(pp,0,10000);

 
  o->gamma = 1.0;
  if(argparser_keyword_present(pp, "-gamma")){
    o->gamma = argparser_get_next_double(pp,0,INF);
  }
 
 o->K = 5.0;
  if(argparser_keyword_present(pp, "-K")){
    o->K = argparser_get_next_double(pp,0,INF);
  }
 
 o->Amax = 1.0;
  if(argparser_keyword_present(pp, "-Amax")){
    o->Amax = argparser_get_next_double(pp,0,INF);
  }
  //now are the not so essential

   
   o->usingBestThree = argparser_keyword_present(pp, "-usingBestThree");
   //    if(o->UsingAngleWeights){
//      argparser_get_keyword_next(pp, "r");
//      o->clusterr = argparser_get_next_double(pp,-1000.0,10000);
//      argparser_get_keyword_next(pp, "R");
//      o->clusterR = argparser_get_next_double(pp,-1000.0,10000);
//    }
  
  argparser_finish(pp);
  return o;
}

float_image_t* ReadFNI(char* filename);
float_image_t* ReadFNI(char* filename){
  FILE* arq = open_read(filename,TRUE);
  float_image_t* fim = float_image_read(arq);
  fclose(arq);
  return fim;
}

void WriteFNI(char* filename,float_image_t* fim);
void WriteFNI(char* filename,float_image_t* fim){
  FILE* arq = open_write(filename,TRUE);
  float_image_write(arq,fim);
  fclose(arq);
}



int main(int argc,char** argv){
  fprintf(stderr, "ARGCS %d\n",argc );
  
  options_t* o = parse_args(argc,argv);

  
   
  fprintf(stderr,"Reconstructing with:");
  if(o->lightDirections != NULL){
    fprintf(stderr,"Original static light source\n");
    
  }
     
  /* Lê as imagens da cena e define o tamanho {nx,ny}: */
  //Imagem *S[num_luzes];
  float_image_t  *S[o->nClusters][o->nLights];
  
  int nx = -1, ny = -1,nc;
  int i,k;
  
  int count = 0;
  for(k = 0; k < o->nClusters; k++){
    for(i = 0; i < o->nLights; i++){
      fprintf(stderr, "Abrindo arquivo[%d][%d] %s ... \n",k,i,o->sceneImages[count]);
    // Imagem* im;
      //pnm_image_read(nomes_arquivos_img[i],(&im));
      float_image_t *im = float_pnm_image_read(o->sceneImages[count],FALSE, o->gamma, 0.0,  TRUE,TRUE,FALSE);
      if (i == 0) 
	{nx = im->sz[1]; ny = im->sz[2];nc = im->sz[0]; }
      else
	{ if ((nx != im->sz[1]) || (ny != im->sz[2])) 
	  { fprintf(stderr, "Imagem S[%d] com tamanho inconsistente!\n", i); exit(1); }
	}
      S[k][i] = im;
      count++;
    }
  }
  fprintf(stderr, "Imagens da cena lidas.\n");
  /**/
  
  r3_t luz_dir[o->nClusters][o->nLights];
  r3_t avg_light_dir[o->nClusters];
  r3_t cluster_dir[o->nClusters];
  
  fprintf(stderr,"Reading Direction files \n");
  
  count = 0;
  for(k = 0; k < o->nClusters; k++){
    fprintf(stderr,"Cluster %d\n",k);
    r3_zero(&(avg_light_dir[k]));
    for(i = 0; i < o->nLights;i++){
      char* nome_arq_dir;
      FILE* arq_dir;
      nome_arq_dir = o->lightDirections[count];
      arq_dir = open_read(nome_arq_dir,TRUE);
      double dx,dy,dz;
      int test_read;
      test_read = fscanf(arq_dir,"%lf %lf %lf",&dx,&dy,&dz);
      if(test_read != 3){
	fprintf(stderr,"Error reading file - %d numbers found\n",test_read);
	return 1;
      }
      luz_dir[k][i].c[0] = dx;
      luz_dir[k][i].c[1] = dy;
      luz_dir[k][i].c[2] = dz;
      
      r3_mix_in (1.0,&(luz_dir[k][i]),&(avg_light_dir[k]));
      
      fprintf(stderr,"[%d][%d] Direction: %+8.5f %+8.5f %+8.5f\n",k,i,dx,dy,dz);
      count++;
    }
    r3_dir(&(avg_light_dir[k]),&(avg_light_dir[k]));
    fprintf(stderr,"Average direction:  %+8.5f %+8.5f %+8.5f\n",avg_light_dir[k].c[0],avg_light_dir[k].c[1],avg_light_dir[k].c[2]);
  }
  
  fprintf(stderr,"Reading Cluster files \n");
  for(k = 0; k < o->nClusters; k++){
    char* nome_arq_dir;
    FILE* arq_dir;
    nome_arq_dir = o->clusterDirections[k];
    arq_dir = open_read(nome_arq_dir,TRUE);
    double dx,dy,dz;
    int test_read;
    test_read = fscanf(arq_dir,"%lf %lf %lf",&dx,&dy,&dz);
    if(test_read != 3){
      fprintf(stderr,"Error reading file - %d numbers found\n",test_read);
      return 1;
    }
    cluster_dir[k].c[0] = dx;
    cluster_dir[k].c[1] = dy;
    cluster_dir[k].c[2] = dz;

    fprintf(stderr,"[%d] Direction: %+8.5f %+8.5f %+8.5f\n",k,dx,dy,dz);
  }
  
  
  fprintf(stderr, "Arquivos de direção lidos.\n");
  
  
  //ShowCoplanar(luz_dir,o->nLights);
  /* Lê imagem com normais de referência, se houver: */
  
  /* Criando imagens de saída: */
  float_image_t* imagem_normal_canal = float_image_new(3,nx,ny);
  float_image_t* imagem_albedo_canal = float_image_new(1,nx,ny);
  
  float_image_t** imagem_normais_cluster = (float_image_t**)malloc(sizeof(float_image_t*)*o->nClusters);
  float_image_t** imagem_albedo_cluster = (float_image_t**)malloc(sizeof(float_image_t*)*o->nClusters);
  float_image_t** imagem_walb_cluster = (float_image_t**)malloc(sizeof(float_image_t*)*o->nClusters);
  float_image_t** imagem_wnrm_cluster = (float_image_t**)malloc(sizeof(float_image_t*)*o->nClusters);
  float_image_t** imagem_wsmp_cluster = (float_image_t**)malloc(sizeof(float_image_t*)*o->nClusters);
  float_image_t** imagem_whigh_cluster = (float_image_t**)malloc(sizeof(float_image_t*)*o->nClusters);

  for(k = 0; k < o->nClusters; k++){
    imagem_normais_cluster[k] = float_image_new(3,nx,ny);
    imagem_albedo_cluster[k] = float_image_new(1,nx,ny);
    imagem_walb_cluster[k] = float_image_new(1,nx,ny);
    imagem_wnrm_cluster[k] = float_image_new(1,nx,ny);
    imagem_wsmp_cluster[k] = float_image_new(1,nx,ny);
    imagem_whigh_cluster[k] = float_image_new(1,nx,ny);
  }
  
    /*Total */
  
  int x,y;
  int total = nx*ny;
  fprintf(stderr, "Total de iteracoes a executar: %d\n", total);
  fprintf(stderr, "Processando\n");
  
;
  
  /* Loop sobre os canais: */
  int canal;
  for (canal = 0; canal < nc; canal++) {
 
    /* Devemos proessar este canal? */
    int processa_canal = strchr(o->channels,"RGB"[canal]) != NULL;
    if(processa_canal) fprintf(stderr,"Computing channel %d\n",canal);
    if(!processa_canal) continue;
    /* Abre arquivo da tabela de normais do canal: */
        
    fprintf(stderr,"---------------------------------------------------------------------");

    int k;
    for(k = 0; k < o->nClusters; k++){
      double* LS;
      LS = ConstroiLeastSquares(luz_dir[k],o->nLights);
      /*Loop sobre pixels da cena: */
      int x, y;
      fprintf(stderr,"\n");
      time_t last_tempo = time(NULL);
      for(y = 0 ;y < ny; y++){
	for(x = 0; x < nx; x++){
	  /* Decide se pixel deve ser debugado: */
	  //int debug = (((x == 14) && (y == 209)) )||  (((x == 12) && (y == 207)));
	  int debug = FALSE;
	  if (debug){ fprintf(stderr,"%d %d\n",x,y); }
    	        
	  double albedo = 0.0;        /* Albedo calculado do ponto {p} da cena, ou 0.0. */
	  double peso = 0.0;	    /* Confiabilidade na normal calculada */
	  r3_t snp = (r3_t){{0,0,0}}; /* Normal estimada do ponto {p} da cena, ou (0,0,0). */
	  /* Cria arquivo para debug dos cálculo de nrmal do pixel: */
	  FILE *arq_debug_pixel = NULL;
	  if (debug) { 
	    besta();
	    char *nome_arq_debug_pixel = NULL;
	    char *nome_arq_debug_pixel = jsprintf("%s_debug_%d_%04d_%04d.txt", o->prefix, canal, x, y);
	    arq_debug_pixel = fopen(nome_arq_debug_pixel, "wt");
	    free(nome_arq_debug_pixel);
	  } 
	  /* Extrai o vetor de observação {SO} deste pixel: */
	  double SO[o->nLights]; 
	  for (i = 0; i < o->nLights; i++){
	    SO[i] = float_image_get_sample(S[k][i], canal, x, y);
	  }
	  double Smag = rn_norm(o->nLights,SO);
	  /*If it is too dark, skip */
	  if(Smag < 0.01){
	    peso = 0.0;
	    albedo = 0.0;
	  }else{
	    /* Determina a normal {snp} para este pixel: */
	    double n_res[3];
	    if(o->usingBestThree){
	      r3_t b;
	      double tp;
	      int imax[3];
	      r3x3_t A,A_inv;
	      ConstroiSistema(SO, luz_dir[k], o->nLights, &A, &b, &tp, imax);
	      r3x3_inv(&A, &A_inv);
	      r3_t res;
	      r3x3_map_col(&A_inv, &b, &res);
	      n_res[0] = res.c[0];
	      n_res[1] = res.c[1];
	      n_res[2] = res.c[2];
	     }else{
	      rmxn_mul (3, o->nLights, 1, LS, SO,n_res);
	     }
	     snp.c[0] = n_res[0];
	     snp.c[1] = n_res[1];
	     snp.c[2] = n_res[2];
	     // Calcula o albedo e normaliza a normal:
	     albedo = r3_dir(&snp,&snp);
	     // if (albedo != 0.0) { r3_scale(1.0/albedo, &snp, &snp); }
	     // Ajusta o peso proporcionalmente ao quadrado do albedo:
	     if (debug) { 
	      int i;
	      fprintf(arq_debug_pixel,"%+16.14f %+16.14f %+16.14f   ",snp.c[0],snp.c[1],snp.c[2]);
	      for(i = 0; i < o->nLights;i++) fprintf(arq_debug_pixel,"%16.14f ",SO[i]);
	      fprintf(arq_debug_pixel,"\n");
	      fclose(arq_debug_pixel); 
	     }
	  } 
	  /* Grava normal no mapa de normais: */
	  float_image_set_sample(imagem_normais_cluster[k],0,x,y,snp.c[0]);
	  float_image_set_sample(imagem_normais_cluster[k],1,x,y,snp.c[1]);
	  float_image_set_sample(imagem_normais_cluster[k],2,x,y,snp.c[2]);
	  
	  /* Armazena o albedo em {imagem_albedo} (gama visual): */
          float_image_set_sample(imagem_albedo_cluster[k], 0, x, y, albedo);
	  /*compute initial values for weights*/
	  double walb,wnrm,wsmp,whigh;
	  walb = computeWalb(albedo,o->Amax);
	  wnrm = computeWnrm(snp,avg_light_dir[k],o->clusterr);
	  wsmp = computeWsmp(SO,o->nLights);
	  //whigh = computeWhigh(snp,cluster_dir[k],o->K,o->clusterr);
	  whigh = 1.0;
	  
	  float_image_set_sample(imagem_walb_cluster[k],0,x,y,walb);
	  float_image_set_sample(imagem_wnrm_cluster[k],0,x,y,wnrm);
	  float_image_set_sample(imagem_wsmp_cluster[k],0,x,y,wsmp);
	  float_image_set_sample(imagem_whigh_cluster[k],0,x,y,whigh);
	  
        
	{
	  fprintf(stderr,"\033[1A");
	  int contador = x + (y*nx);
	  int total_compute_pixels = (nx*ny);
	  double normals_per_sec = contador/(double)(time(NULL) - last_tempo  + 0.001);
	  double total_secs_remaining = (total_compute_pixels - contador)/normals_per_sec;
	  int total_seconds = (int)floor(total_secs_remaining);
	  int hour = total_seconds/(60*60);
	  int min = (total_seconds - (hour*60*60))/60;
	  int sec = (total_seconds - (hour*60*60) - (min*60));
	  fprintf(stderr,"[%d][%d][%9d] of [%9d] - %6.6f%% - %6.6f n/s   - %02d h %02d m %02d s       \n",
		canal,k,contador,total_compute_pixels, contador*100.0/total_compute_pixels,normals_per_sec,hour,min,sec);
	}
	

	
	
	//finished X
	}
	//finished Y
      }
      
     //finished this cluster 
    }
    /*We will save here the results of the first processing*/
    for(k = 0; k < o->nClusters;k++){
      char* arq_name = NULL;
      char *arq_name = jsprintf("%s_%d_C%d_normals.fni",o->prefix,canal,k);
      WriteFNI(arq_name,imagem_normais_cluster[k]);
      free(arq_name);
      
      arq_name = NULL;
      char *arq_name = jsprintf("%s_%d_C%d_albedo.fni",o->prefix,canal,k);
      WriteFNI(arq_name,imagem_albedo_cluster[k]);
      free(arq_name);
      
      arq_name = NULL;
      char *arq_name = jsprintf("%s_%d_C%d_walbP.fni",o->prefix,canal,k);
      WriteFNI(arq_name,imagem_walb_cluster[k]);
      free(arq_name);
      
      arq_name = NULL;
      char *arq_name = jsprintf("%s_%d_C%d_wnrmP.fni",o->prefix,canal,k);
      WriteFNI(arq_name,imagem_wnrm_cluster[k]);
      free(arq_name);
      
      arq_name = NULL;
      char *arq_name = jsprintf("%s_%d_C%d_wsmpP.fni",o->prefix,canal,k);
      WriteFNI(arq_name,imagem_wsmp_cluster[k]);
      free(arq_name);
      
      arq_name = NULL;
      char *arq_name = jsprintf("%s_%d_C%d_whighP.fni",o->prefix,canal,k);
      WriteFNI(arq_name,imagem_whigh_cluster[k]);
      free(arq_name);
    }
    /*We saved now  post-process the normal*/
    
    fprintf(stderr,"\nFinished Cluster processing, interpolating normals\n");
   int  contador = 0;
    
    for(y = 0 ;y < ny; y++){
	for(x = 0; x < nx; x++){
	  r3_t avg_norm;
	  r3_t prv_norm;
	  double A;
	  double diff = 1.0;
	  double epsilon = 10e-10;
	  int num_iter = 0;
	  int MAX_ITER = 200;
	  double in_SO[o->nClusters][o->nLights];
	  for(k = 0; k < o->nClusters; k++){
	    int i;
	    for(i = 0; i < o->nLights; i++){
	      in_SO[k][i] = float_image_get_sample(S[k][i],canal,x,y);
	    }
	  }
	  
	  
	 while((diff > epsilon) && (num_iter < MAX_ITER) ){
	    /*Compute the normal from the weights*/
	    prv_norm = avg_norm;
	    double sW = 0;
	    r3_zero(&avg_norm);
	    for(k = 0; k < o->nClusters; k++){
	      
	      r3_t norm;
	      norm.c[0] = float_image_get_sample(imagem_normais_cluster[k],0,x,y);
	      norm.c[1] = float_image_get_sample(imagem_normais_cluster[k],1,x,y);
	      norm.c[2] = float_image_get_sample(imagem_normais_cluster[k],2,x,y);
	
	      
	      double walb = float_image_get_sample(imagem_walb_cluster[k],0,x,y);
	      double wnrm = float_image_get_sample(imagem_wnrm_cluster[k],0,x,y);
	      double wsmp = float_image_get_sample(imagem_wsmp_cluster[k],0,x,y);
	      double whigh = float_image_get_sample(imagem_whigh_cluster[k],0,x,y);
	      double albedo = float_image_get_sample(imagem_albedo_cluster[k],0,x,y);
	      
     
	      double weights = walb*wnrm*wsmp*whigh;
	      r3_mix_in (albedo*weights, &(norm), &avg_norm);
	      sW+= weights*albedo;
	    }
	    
	    if(sW != 0){
	      r3_scale(1/sW,&avg_norm,&avg_norm);
	    
	      A = r3_dir(&avg_norm,&avg_norm);
	    }else{
	      A = 0;
	      r3_zero(&avg_norm);
	    }
	    
	    
	    /*update weights*/
	    for(k = 0; k < o->nClusters; k++){
	      //float_image_set_sample(imagem_walb_cluster[k],0,x,y,computeWalb(A,o->Amax));
	      float_image_set_sample(imagem_wnrm_cluster[k],0,x,y,computeWnrm(avg_norm,avg_light_dir[k],o->clusterr));
	    //  float_image_set_sample(imagem_wsmp_cluster[k],0,x,y,computeWsmp(in_SO[k],o->nLights));
	      float_image_set_sample(imagem_whigh_cluster[k],0,x,y,computeWhigh(avg_norm,cluster_dir[k],o->K,o->clusterr));  
	    }
	    /*compute new normal*/
	    
	    if(num_iter ==0){
	      diff = 1.0;
	    }else{
	      diff = r3_dist(&prv_norm,&avg_norm);
	    }
	    num_iter++;
    
	  }
	  
	  //CORRECTS NORMALS WITH ABSURD VALUES
	  if(avg_norm.c[2] < 0) {
	    r3_zero(&avg_norm);
	  }
	  
	  float_image_set_sample(imagem_normal_canal,0,x,y,avg_norm.c[0]);
	  float_image_set_sample(imagem_normal_canal,1,x,y,avg_norm.c[1]);
	  float_image_set_sample(imagem_normal_canal,2,x,y,avg_norm.c[2]);
	  float_image_set_sample(imagem_albedo_canal,0,x,y,A);
	 
	  {
	    fprintf(stderr,"\033[1A");
	    int contador = x + (y*nx);
	    int total_compute_pixels = (nx*ny);
	    fprintf(stderr,"[%d][%d][%9d] of [%9d] - %6.6f%%       \n",
		  canal,k,contador,total_compute_pixels, contador*100.0/total_compute_pixels);
	  }
	contador++;
	  
	}
    }
    
    char* arq_name = NULL;
    char *arq_name = jsprintf("%s_%d_normals.fni",o->prefix,canal);
    WriteFNI(arq_name,imagem_normal_canal);
    free(arq_name);
    
    arq_name = NULL;
    char *arq_name = jsprintf("%s_%d_albedo.fni",o->prefix,canal);
    WriteFNI(arq_name,imagem_albedo_canal);
    free(arq_name);
    
    /*Save all the final weights... uff!*/
    
    for(k = 0; k < o->nClusters;k++){
      arq_name = NULL;
      char *arq_name = jsprintf("%s_%d_C%d_walbF.fni",o->prefix,canal,k);
      WriteFNI(arq_name,imagem_walb_cluster[k]);
      free(arq_name);
      
      arq_name = NULL;
      char *arq_name = jsprintf("%s_%d_C%d_wnrmF.fni",o->prefix,canal,k);
      WriteFNI(arq_name,imagem_wnrm_cluster[k]);
      free(arq_name);
      
      arq_name = NULL;
      char *arq_name = jsprintf("%s_%d_C%d_wsmpF.fni",o->prefix,canal,k);
      WriteFNI(arq_name,imagem_wsmp_cluster[k]);
      free(arq_name);
      
      arq_name = NULL;
      char *arq_name = jsprintf("%s_%d_C%d_whighF.fni",o->prefix,canal,k);
      WriteFNI(arq_name,imagem_whigh_cluster[k]);
      free(arq_name);
    }
    
    
    
    fprintf(stderr,"Channel %d processed\n",canal);
   
   //finished this channel
  }
   
  
  
  
  
  fprintf(stderr, "Concluido!\nO programa rodou com sucesso!\n");

  return 0;

}


