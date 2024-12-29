
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

#define PROG_NAME "compute_normals_simple2"


#define PROG_HELP \
  PROG_NAME " \\\n" \
  "	-prefix {FILE_PREFIX} \\\n" \
  "	-channels {RGB} \\\n" \
  "	-nLights {NLIGHTS} \\\n" \
  "	-sceneImages {Image0 Image1...} \\\n" \
  "     [-lightsDirections {file0 file1..} ]  \\\n" \
  "     [-lightFunctions \\\n" \
  "     	albedo ALBEDOFILE \\\n" \
  "     	azimuth AZIMUTHFILE \\\n" \
  "     	elevation ELEVATIONFILE ]\\\n" \
  "	[-transform {OFFSETX} {OFFSETY} {SCALEX} {SCALEY} ]\\\n" \
  "	[- gamma {GAMMA} ] \\\n" \
  "	[-invertAlbedo] \\\n" \
  "	[-scaleAlbedo {ALBEDO SCALE} ] \\\n" \
  "	[-generateWhite] \\\n" \
  "	[-generateKili] \\\n" \
  "     [-UsingRANSAC \\\n" \
  "     	V {VARIANCE} \\\n" \
  "     	E {E} ]\\\n" \
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
  char** sceneImages;
  char** lightDirections;
  char* albedo_funcfile;
  char* azimuth_funcfile;
  char* elevation_funcfile;
  transform_data_t* transform;
  double gamma;
  
  bool_t generateWhite;
  bool_t generateKili;
  bool_t invertAlbedo;
  double scaleAlbedo;
  
  bool_t usingBestThree;
  bool_t UsingAngleWeights;
  double clusterr;
  double clusterR;
  
  bool_t UsingRANSAC;
  double V,E;
};

typedef struct options_t options_t;

void besta(void);

void besta(void) { 
  return;
}

void processa_arq_entrada
( FILE *arq, 
  int* num_luzes,
  double *gamma,
  double* bias,
  char*** nomes_arquivos_dir,
  char*** nomes_arquivos_img
  );

  
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


void computeGO(double* GO,r3_t snp,r3_t* luz_dir,int num_luzes){
  int i;
  for(i = 0; i < num_luzes; i++){
     double val = r3_dot(&(luz_dir[i]),&snp);
     if(val < 0 ) val = 0;
     GO[i] = val;
  }
}

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


void AnglesToDir(double ele,double azi, double* dir){
     dir[0] = -sin(azi)*cos(ele);
     dir[1] = -cos(azi)*cos(ele);
     dir[2] = sin(ele);
//   dir[0] = -sin(azi)*cos(ele);
//   dir[1] = -cos(azi)*cos(ele);
//   dir[2] = sin(ele);
}

void DirToAngles(double* dir, double* aziP, double* eleP){
  
  double ele,azi;
  double x,y,z;
  x = -dir[0];
  y = dir[1];
  z = dir[2];
  
  ele = asin(z);
  azi = atan(x/y);
  
  if(ele < (M_PI/2.0)) {
    if(y < 0) azi+=M_PI;
  }else{
    if(y > 0) azi+=M_PI;
  }
  
  if(azi < 0) azi+=2*M_PI;
  
  *aziP= azi;
  *eleP = ele;
}

void UpdateLightSourceValues(r3_t* luz_dir,
			     r3_t* original_dir,
			     int num_luzes,
			     poly_function_t** az_func,
			     poly_function_t** el_func,
			     poly_function_t** int_func,
			     double x, double y,
			     bool_t invert
){
  
  int l;
  
  for(l = 0; l < num_luzes; l++){
    double azimuth,elevation,intensity;
    double X[] = {x,y};
    
    if(int_func != NULL){
      intensity = EvaluatePvalue(int_func[l],X);
      if(invert){ intensity = 1.0/intensity;}
    }else{
      intensity = r3_norm(&original_dir[l]);
    }
    r3_t dir;
    if( (az_func != NULL) && (el_func != NULL)){
      //REMEMBER, THE VALUES OF AZIMUTH AND ELEVATION ARE IN RADIANS
      azimuth = EvaluatePvalue(az_func[l],X);
      elevation = EvaluatePvalue(el_func[l],X);
      AnglesToDir(elevation,azimuth,dir.c);
    }else{
      dir = original_dir[l];
    }
    
    r3_dir(&dir,&dir);
    r3_scale(intensity,&dir,&dir);
    luz_dir[l] = dir;
  }
  
  
  
}

poly_function_t*** ReadPolyData(char* filename){
  FILE* arq = open_read(filename,TRUE);
  poly_function_t*** p;
  int num_lights;
  fscanf(arq,"%d",&num_lights);
  p = (poly_function_t***)malloc(sizeof(poly_function_t**)*num_lights);
  int l;
  for(l = 0; l < num_lights;l++){
    int num_channels;
    fscanf(arq,"%d",&num_channels);
    p[l] = (poly_function_t**)malloc(sizeof(poly_function_t*)*num_channels);
    int c;
    for(c = 0; c < num_channels; c++){
	p[l][c] = (poly_function_t*)malloc(sizeof(poly_function_t));
	poly_function_t* polym = p[l][c];
	fscanf(arq,"%d",&(polym->dimensions));
	fscanf(arq,"%d",&(polym->num_coefs));
	polym->coefs = (double**)malloc(sizeof(double*)*(polym->num_coefs));
	polym->weights = (double*)malloc(sizeof(double)*(polym->num_coefs));
	int icf;
	for(icf = 0; icf < (polym->num_coefs);icf++){
	  polym->coefs[icf] = (double*)malloc(sizeof(double)*(polym->dimensions));
	  int d;
	  for(d = 0; d < (polym->dimensions); d++){
	    fscanf(arq,"%lf",&(polym->coefs[icf][d]));
	  }
	  fscanf(arq,"%lf",&(polym->weights[icf]));
	}
    }
  }
  fclose(arq);
  return p;
}



double user_cpu_time_usec(void);

double user_cpu_time_usec(void){
  struct tms buf;
  (void)times(&buf);
  return(1000000.0 * ((double) buf.tms_utime)/((double)sysconf(_SC_CLK_TCK)));
}



void processa_arq_entrada
( FILE *arq, 
  int* num_luzes,
  double *gamma,
  double *bias,
  char*** nomes_arquivos_dir,
  char*** nomes_arquivos_img
  )
{
   
  fscanf(arq,"%d",num_luzes);
    
  int numero_arquivos = *num_luzes;
  fprintf(stderr, "%d arquivos a serem processados\n",numero_arquivos);

  //aloca vetor de arquivos !
  *nomes_arquivos_dir = (char**) malloc(sizeof(char*)*(numero_arquivos));
  *nomes_arquivos_img = (char**) malloc(sizeof(char*)*(numero_arquivos));

  //aloca vetor de nomes de imagens
  int ind;
  for(ind = 0; ind < (numero_arquivos); ind++){
    (*nomes_arquivos_img)[ind] = (char*)malloc(sizeof(char)*400);
  }

  //aloca vetor de nomes de direcoes
  for(ind = 0; ind < (numero_arquivos); ind++){
    (*nomes_arquivos_dir)[ind] = (char*)malloc(sizeof(char)*400);
  }
  
  for(ind = 0; ind < (numero_arquivos);ind++){
  	fscanf(arq,"%s",(*nomes_arquivos_dir)[ind]);
  }

  //le os nomes de arquivo de imagem
  for(ind = 0; ind < (numero_arquivos);ind++){
    fscanf(arq,"%s",(*nomes_arquivos_img)[ind]);
  }

  //Captura a cor dp gabarito, 1,1,1 significa Branco
  //Le gamma das imagens
  fscanf(arq,"%lf",gamma);
  (*bias) = ((*gamma) == 1.000 ? 0.000 : VIEW_BIAS); /* Hack... */
  fprintf(stderr, "Gamma das imagens:%lf \n",*gamma);

  
  fprintf(stderr, "Direcoes:\n");
  for(ind = 0; ind < (numero_arquivos);ind++){
    fprintf(stderr, "%s\n",(*nomes_arquivos_dir)[ind]);
  }
  fprintf(stderr, "Imagens:\n");
  for(ind = 0; ind < (numero_arquivos);ind++){
    fprintf(stderr, "%s\n",(*nomes_arquivos_img)[ind]);
  }

  //fprintf(stderr,  "Gabaritos:%s %s %s \n",nomes_de_arquivo[0],nomes_de_arquivo[1],nomes_de_arquivo[2]);
  //fprintf(stderr,  "Imagens:%s %s %s \n",nomes_arquivos_img[0],nomes_arquivos_img[1],nomes_arquivos_img[2]);
}


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
  fprintf(stderr,"We have %s channels !\n",o->channels);

  argparser_get_keyword(pp, "-nLights");
  o->nLights = argparser_get_next_int(pp,3,10000);


  argparser_get_keyword(pp, "-sceneImages");
  o->sceneImages = (char**)malloc(sizeof(char*)*(o->nLights));
  int i;
  for(i = 0; i < o->nLights; i++){
    o->sceneImages[i] = argparser_get_next(pp);
  }
  
  o->lightDirections = NULL;
  if(argparser_keyword_present(pp, "-lightDirections")){
    o->lightDirections = (char**)malloc(sizeof(char*)*(o->nLights));
    for(i = 0; i < o->nLights; i++){
      o->lightDirections[i] = argparser_get_next(pp);
    }
  }

  o->albedo_funcfile =  NULL;
  o->azimuth_funcfile =  NULL;
  o->elevation_funcfile =  NULL;
  

  if(argparser_keyword_present(pp, "-lightFunctions")){

      argparser_get_keyword_next(pp, "albedo");
      o->albedo_funcfile = argparser_get_next(pp);
      fprintf(stderr,"AF:%s\n",o->albedo_funcfile);
      if(strcmp("NONE",o->albedo_funcfile) == 0){
	o->albedo_funcfile = NULL;
      }
      argparser_get_keyword_next(pp, "azimuth");
      o->azimuth_funcfile = argparser_get_next(pp);
      fprintf(stderr,"ZF:%s\n",o->azimuth_funcfile);
      if(strcmp("NONE",o->azimuth_funcfile) == 0){
	o->azimuth_funcfile = NULL;
      }
      
      argparser_get_keyword_next(pp, "elevation");
      o->elevation_funcfile = argparser_get_next(pp);
      fprintf(stderr,"EF:%s\n",o->elevation_funcfile);
      if(strcmp("NONE",o->elevation_funcfile) == 0){
	o->elevation_funcfile = NULL;
      }
      
  }

  bool_t testFuncs =   (o->albedo_funcfile == NULL) || (o->azimuth_funcfile == NULL) || (o->elevation_funcfile == NULL);

  if((o->lightDirections == NULL) && (testFuncs) ){
    argparser_error(pp, "Insuficient data for reconstruction");
    
  }
  

  o->transform = NULL;
  if(argparser_keyword_present(pp, "-transform")){
    o->transform = (transform_data_t*)malloc(sizeof(transform_data_t));
    o->transform->offset.c[0] = argparser_get_next_double(pp,-INF,INF);
    o->transform->offset.c[1] = argparser_get_next_double(pp,-INF,INF);
    o->transform->scale.c[0] = argparser_get_next_double(pp,0,INF);
    o->transform->scale.c[1] = argparser_get_next_double(pp,0,INF);
  }
 
  o->gamma = 1.0;
  if(argparser_keyword_present(pp, "-gamma")){
    o->gamma = argparser_get_next_double(pp,0,INF);
  }
 
  //now are the not so essential
   o->invertAlbedo = argparser_keyword_present(pp, "-invertAlbedo");
   o->generateWhite = argparser_keyword_present(pp, "-generateWhite");
   o->generateKili = argparser_keyword_present(pp, "-generateKili");
  
   o->scaleAlbedo = 1.0;
   if(argparser_keyword_present(pp, "-scaleAlbedo")){
    o->scaleAlbedo = argparser_get_next_double(pp,0,INF);
   }
   
   o->UsingRANSAC = argparser_keyword_present(pp, "-UsingRANSAC");
   if(o->UsingRANSAC){
     argparser_get_keyword_next(pp, "V");
     o->V = argparser_get_next_double(pp,0.0,10000);
     argparser_get_keyword_next(pp, "E");
     o->E = argparser_get_next_double(pp,0.0,10000);
     fprintf(stderr,"Using Ransac (%lf,%lf)\n",o->V,o->E);
   }
   
   o->usingBestThree = argparser_keyword_present(pp, "-usingBestThree");
   o->UsingAngleWeights = argparser_keyword_present(pp, "-usingAngleWeights");
   if(o->UsingAngleWeights){
     argparser_get_keyword_next(pp, "r");
     o->clusterr = argparser_get_next_double(pp,-1000.0,10000);
     argparser_get_keyword_next(pp, "R");
     o->clusterR = argparser_get_next_double(pp,-1000.0,10000);
   }
  
  argparser_finish(pp);
  return o;
}

void computeNormalsRANSAC(double SO[] ,r3_t* luz_dir,int nLights,double n_res[] ,double* ransac_weight,double V, double E,double *num_valids){
  int nComb =  nLights*(nLights-1)*(nLights-2)/6.0;
  int valids[nComb];
  r3_t normals[nComb];
  double alb[nComb];
  int i,j,k;
  
  r3x3_t A,A_inv;
  r3_t B;
  r3_t snp;
  int count = 0;
  int count_valids = 0;
  
  for(i = 0; i < nLights; i++){
    for(j = (i+1); j < nLights; j++){
      for(k = (j+1); k < nLights;k++){
	r3_t light_vector[3];
	light_vector[0] = luz_dir[i];
	light_vector[1] = luz_dir[j];
	light_vector[2] = luz_dir[k];
	double c_SO[3];
	c_SO[0] =  SO[i];
	c_SO[1] =  SO[j];
	c_SO[2] =  SO[k];
	
	int row, col;
	for(row = 0; row < 3; row++){
	  B.c[row] = c_SO[row];
	  for(col = 0; col < 3; col++){
	    A.c[row][col] = light_vector[row].c[col];
	  }
	}
	r3x3_inv(&A, &A_inv);
        r3x3_map_col(&A_inv, &B, &snp);
	alb[count] = r3_dir(&snp,&snp);
	normals[count] = snp;
	double nz = snp.c[2];
	if( (nz < 0) || (nz*light_vector[0].c[2] < 0) || (nz*light_vector[1].c[2] < 0) || (nz*light_vector[2].c[2] < 0) ){
	  valids[count] = 0;	  
	}else{
	  valids[count] = 1;
	  count_valids++;
	}
	
	count++;
	assert(count <= nComb);
      }
    }
  }
  
  *num_valids = count_valids;
  if(count_valids == 0){
    n_res[0] = n_res[1] = n_res[2] = 0;
    *ransac_weight = 0;
  }else if(count_valids == 1){
    for(i = 0; i < nComb; i++){
      if(valids[i] == 1) break;
    }
    double po = 0.25;
    double Vfinal = po*1 + (1-po)*V;
    
    n_res[0] = alb[i]*normals[i].c[0];
    n_res[1] = alb[i]*normals[i].c[1];
    n_res[2] = alb[i]*normals[i].c[2];
    *ransac_weight = 1/Vfinal;
  }else if(count_valids == 2){
    r3_t normal_valid[2];
    double albs[2];
    int count_norm = 0;
    
    for(i = 0; i < nComb; i++){
      if(valids[i] == 1) {
	normal_valid[count_norm] = normals[i];
	albs[count_norm] = alb[i];
	count_norm++;
      }
    }
    r3_t avg_norm;
    r3_add(&(normal_valid[0]),&(normal_valid[1]),&avg_norm);
    r3_dir(&avg_norm,&avg_norm);
    double albedo = (albs[0]+albs[1])/2.0;
    double Vfinal = V + ( 1 - r3_dot(&(normal_valid[0]),&(normal_valid[1])));
    n_res[0] = albedo*avg_norm.c[0];
    n_res[1] = albedo*avg_norm.c[1];
    n_res[2] = albedo*avg_norm.c[2];
    *ransac_weight = 1/Vfinal;
  }else{
    double Va[nComb];
    double P[nComb];
    double weights[nComb];
    double previousW[nComb];
    double A[nComb];
    double albedo;
    
    r3_t avg_norm;
  
    //init weights
    for(i = 0; i < nComb; i++){
      weights[i] = 0;
      if(valids[i]){
	weights[i] = 1.0;
	P[i] = 1.0/(double)count_valids;
      }
    }

    
    double diff = 10;
    int iter = 0;
    while((diff > 10e-6) && (iter < 200)){
     
      r3_zero(&avg_norm);
      albedo = 0;
      //compute average normal
      for(i = 0; i < nComb; i++){
	if(valids[i]){
	  r3_t scl_nrm;
	  r3_scale(weights[i],&(normals[i]),&scl_nrm);
	  r3_add(&scl_nrm,&avg_norm,&avg_norm);
	  albedo+= alb[i];
	}
      }
      
      r3_dir(&avg_norm,&avg_norm);
      //compute difference angle;
      for(i = 0; i < nComb; i++){
	if(valids[i]){
	  A[i] = r3_dot(&(normals[i]),&avg_norm);
	}
      }
      //recompute weights
      for(i = 0; i < nComb; i++){
	if(valids[i]){
	  Va[i] = (P[i])*1 + (1 - P[i])*V;
	  previousW[i] = weights[i];
	  weights[i] = 1/Va[i];
	  P[i] = 1/(1 + E*exp(-(1-cos(A[i]))/V));
	  if(!((weights[i] >= 0) )){
	    fprintf(stderr,"W: %lf P %lf\n",weights[i],P[i]);
	  }
	  assert((P[i] >= 0) && (P[i] <= 1));
	  assert((weights[i] >= 0) );
	}  
      }
      
      diff = rn_dist(nComb,weights,previousW);
      iter++;
    }
    
    //usar maior dos pesos ou a média ?
    double maxW = 0;
    double sumW = 0;
    for(i = 0; i < nComb; i++){
	if(valids[i]){
	  sumW+= weights[i];
	  if(maxW < weights[i]) maxW = weights[i];
	}
    }
    *ransac_weight = sumW;
    //r3_dir(&avg_norm,&avg_norm);
    r3_scale(albedo,&avg_norm,&avg_norm);
    n_res[0] = avg_norm.c[0];
    n_res[1] = avg_norm.c[1];
    n_res[2] = avg_norm.c[2];
      
  }
}


void computeNormalsRANSAC2(double SO[] ,r3_t* luz_dir,int nLights,double n_res[] ,double* ransac_weight,double V,double *num_valids){
  int nComb =  nLights*(nLights-1)*(nLights-2)/6.0;
  int valids[nComb];
  r3_t normals[nComb];
  double alb[nComb];
  int i,j,k;
  
  r3x3_t A,A_inv;
  r3_t B;
  r3_t snp;
  int count = 0;
  int count_valids = 0;
   
  for(i = 0; i < nLights; i++){
    for(j = (i+1); j < nLights; j++){
      for(k = (j+1); k < nLights;k++){
	r3_t light_vector[3];
	light_vector[0] = luz_dir[i];
	light_vector[1] = luz_dir[j];
	light_vector[2] = luz_dir[k];
	double c_SO[3];
	c_SO[0] =  SO[i];
	c_SO[1] =  SO[j];
	c_SO[2] =  SO[k];
	
	int row, col;
	for(row = 0; row < 3; row++){
	  B.c[row] = c_SO[row];
	  for(col = 0; col < 3; col++){
	    A.c[row][col] = light_vector[row].c[col];
	  }
	}
	r3x3_inv(&A, &A_inv);
        r3x3_map_col(&A_inv, &B, &snp);
	alb[count] = r3_dir(&snp,&snp);
	
	normals[count] = snp;
	double nz = snp.c[2];
	if( (nz < 0) || (nz*light_vector[0].c[2] < 0) || (nz*light_vector[1].c[2] < 0) || (nz*light_vector[2].c[2] < 0) ){
	  valids[count] = 0;	  
	}else if(TestCoplanar(light_vector[0],light_vector[1],light_vector[2])){
	  valids[count] = 0;
	}else{
	  valids[count] = 1;
	  count_valids++;
	}
	
	count++;
	assert(count <= nComb);
      }
    }
  }
  
  int choosen_i = 0;
  int ncons[nComb];
  for(i  = 0; i < nComb; i++){
    ncons[i] = 0;
    for(j  = 0; j < nComb; j++){
      double angle = r3_dot(&(normals[i]),&(normals[j]));
      if((1 - angle) < (5*V)){
	ncons[i]  = ncons[i]+1;
      }
    }
    if(ncons[i] > ncons[choosen_i]){
      choosen_i = i;
    }
  }
  
  double albedo = 0;
  r3_t avg_norm;
  *ransac_weight = (ncons[choosen_i]-1)/V;
  *num_valids = ncons[choosen_i];
  r3_zero(&avg_norm);
  for(j  = 0; j < nComb; j++){
    double angle = r3_dot(&(normals[choosen_i]),&(normals[j]));
    if((1 - angle) < (5*V)){
      r3_add(&(normals[j]),&avg_norm,&avg_norm);
      albedo+=alb[j];
    }
  }
  assert(albedo > 0);
  assert(ncons[choosen_i] > 0);
  
  
  albedo = albedo/(double)ncons[choosen_i];
  r3_scale(1/(double)ncons[choosen_i],&avg_norm,&avg_norm);
  r3_dir(&avg_norm,&avg_norm);
  assert(r3_norm(&avg_norm) > 0);
  n_res[0] = albedo*avg_norm.c[0];
  n_res[1] = albedo*avg_norm.c[1];
  n_res[2] = albedo*avg_norm.c[2];
  
}



int main(int argc,char** argv){
  fprintf(stderr, "ARGCS %d\n",argc );
  
  options_t* o = parse_args(argc,argv);

  
  
  poly_function_t*** int_functions = NULL;
  poly_function_t*** az_functions = NULL;
  poly_function_t*** el_functions  = NULL;
  bool_t UsingPolyDataFile = FALSE;
  
  fprintf(stderr,"Reconstructing with:");
  if(o->lightDirections != NULL){
    fprintf(stderr,"Original static light source\n");
    
  }
  if(o->albedo_funcfile != NULL){
    UsingPolyDataFile = TRUE;
    fprintf(stderr,"Variable albedo\n");
    int_functions = ReadPolyData(o->albedo_funcfile); //it needs to be in order [channels][lights]
  }
  
  if(o->azimuth_funcfile != NULL){
    fprintf(stderr,"Variable azimuth\n");
    UsingPolyDataFile = TRUE;
    az_functions = ReadPolyData(o->azimuth_funcfile);
  }
  
  if(o->elevation_funcfile != NULL){
    fprintf(stderr,"Variable elevation\n");
     el_functions = ReadPolyData(o->elevation_funcfile);
  }
    


  /* Lê as imagens dos gabaritos: */
     
  /* Lê as imagens da cena e define o tamanho {nx,ny}: */
  //Imagem *S[num_luzes];
  float_image_t  *S[o->nLights];
  
  int nx = -1, ny = -1,nc;
  int i;
  
  for(i = 0; i < o->nLights; i++){
    fprintf(stderr, "Abrindo arquivo[%d] %s ... \n",i,o->sceneImages[i]);
   // Imagem* im;
    //pnm_image_read(nomes_arquivos_img[i],(&im));
    float_image_t *im = float_pnm_image_read(o->sceneImages[i],FALSE, o->gamma, 0.0,  TRUE,TRUE,FALSE);
    if (i == 0) 
      {nx = im->sz[1]; ny = im->sz[2];nc = im->sz[0]; }
    else
      { if ((nx != im->sz[1]) || (ny != im->sz[2])) 
	{ fprintf(stderr, "Imagem S[%d] com tamanho inconsistente!\n", i); exit(1); }
      }
     S[i] = im;
  }
  fprintf(stderr, "Imagens da cena lidas.\n");
  /**/
  
  r3_t luz_dir[o->nLights];
  r3_t original_dir[o->nLights];
  
  if(o->lightDirections != NULL){
    for(i = 0; i < o->nLights;i++){
      char* nome_arq_dir;
      FILE* arq_dir;
      nome_arq_dir = o->lightDirections[i];
      fprintf(stderr,"Abrindo arquivo de direcao %s ...\n",nome_arq_dir);
      arq_dir = fopen(nome_arq_dir,"rt");
      if(arq_dir == NULL){
	fprintf(stderr,"Nao conseguiu abrir arquivo !\n");
	return 1;
      }
      double dx,dy,dz;
      int test_read;
      test_read = fscanf(arq_dir,"%lf %lf %lf",&dx,&dy,&dz);
	  
      if(test_read != 3){
	fprintf(stderr,"Error reading file - %d numbers found\n",test_read);
	return 1;
      }
      luz_dir[i].c[0] = dx;
      luz_dir[i].c[1] = -dy;
      luz_dir[i].c[2] = dz;
      fprintf(stderr,"Direction: %+8.5f %+8.5f %+8.5f\n",dx,dy,dz);
    //  free(nome_arq_dir);
      original_dir[i] = luz_dir[i];
    }
    fprintf(stderr, "Arquivos de direção lidos.\n");
  }
  
  //ShowCoplanar(luz_dir,o->nLights);
  /* Lê imagem com normais de referência, se houver: */
  
  /* Criando imagens de saída: */
  float_image_t  *fi_normais_canal[nc];
  for(i = 0; i < nc; i++){
    fi_normais_canal[i] = float_image_new(3,nx,ny);
  }
  float_image_t  *imagem_albedo = float_image_new(nc, nx, ny);
  float_image_t  *imagem_peso = float_image_new(3, nx, ny);
  float_image_t  *imagem_peso_canal = float_image_new(1, nx, ny);
  
  fprintf(stderr, "Imagens de saída criadas.\n");

  int x,y;
  /* Abre arquivos de saída: */
    
  int total = nx*ny;
  fprintf(stderr, "Total de iteracoes a executar: %d\n", total);
  fprintf(stderr, "Processando\n");
  
  //tabela com normais calculadas e pesos, indexadas por canal
  r3_t **normais_calculadas = malloc(nc*sizeof(r3_t *));
  double ** pesos = (double**) malloc(sizeof(double*)* nc);
  
  /* Loop sobre os canais: */
  int canal;
  for (canal = 0; canal < nc; canal++) {
 
    /* Devemos proessar este canal? */
    int processa_canal = strchr(o->channels,"RGB"[canal]) != NULL;
    if(processa_canal) fprintf(stderr,"Computing channel %d\n",canal);
    if(!processa_canal) continue;
    /* Abre arquivo da tabela de normais do canal: */
    char *nome_arq_normais_canal = NULL;
    char *nome_arq_normais_canal = jsprintf("%s_%d_normals.fni", o->prefix, canal);
    FILE* arq_normais_canal = fopen(nome_arq_normais_canal,"wt");
  //  fprintf(arq_normais_canal,"tx = %d\n",nx);
   // fprintf(arq_normais_canal,"ty = %d\n",ny);

    /* Abre arquivo da pesos do canal: */
    char *nome_arq_pesos_canal = NULL;
    char *nome_arq_pesos_canal = jsprintf("%s_%d_weights.fni", o->prefix, canal);
    FILE* arq_pesos_canal = fopen(nome_arq_pesos_canal,"wt");
   // fprintf(arq_pesos_canal,"tx = %d\n",nx);
   // fprintf(arq_pesos_canal,"ty = %d\n",ny);
    char* nome_arq_valid_lights = NULL;
    char *nome_arq_valid_lights = jsprintf("%s_%d_valids.fni", o->prefix, canal);
    FILE* arq_valid_lights = fopen(nome_arq_valid_lights,"wt");
  

    fprintf(stderr,"---------------------------------------------------------------------");

    /* Aloca tabelas de pesos e normais calculadas para este canal: */
    normais_calculadas[canal] = malloc(nx*ny*sizeof(r3_t));
    float_image_t* fi_valid_lights = float_image_new(1,nx,ny);
    pesos[canal] = (double*)malloc(sizeof(double)*nx*ny);
    /*build LS  system*/
    double* LS;
    LS = ConstroiLeastSquares(luz_dir,o->nLights);
        
    /* Loop sobre pixels da cena: */
    
    int x, y;
    fprintf(stderr,"\n");
    double Lmax = -1;
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

        if (processa_canal) {
	  
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
    	    //SO[i] = get_intensity(x, y, S[i], canal, gamma);
	    SO[i] = float_image_get_sample(S[i], canal, x, y);
    	  }
	  double Smag = rn_norm(o->nLights,SO);
	  /*If it is too dark, skip */
	  if(Smag < 0.01){
	    peso = 0.0;
	    albedo = 0.0;
	  }else{
	    
	    /* Determina a normal {snp} para este pixel: */
	    if( UsingPolyDataFile){
	      poly_function_t** el_func = NULL;
	      if(o->elevation_funcfile != NULL) el_func = el_functions[canal];
	      poly_function_t** az_func = NULL;
	      if(o->azimuth_funcfile != NULL) az_func = az_functions[canal];
	      poly_function_t** int_func = NULL; 
	      if(o->albedo_funcfile != NULL) int_func = int_functions[canal];
	      double ix = x;
	      double iy = y;
	      if(o->transform != NULL){
		  ix = ix/o->transform->scale.c[0];
		  iy = iy/o->transform->scale.c[1];
		  ix = ix + o->transform->offset.c[0]; 
		  iy = iy + o->transform->offset.c[1];
		  
	      }
	      UpdateLightSourceValues(luz_dir,original_dir,o->nLights,az_func,el_func,int_func,ix,iy,o->invertAlbedo);
	      
	      if(!o->UsingRANSAC)  LS = ConstroiLeastSquares(luz_dir,o->nLights);
	      
	    }
	    double n_res[3];
	    double ransac_weight;
	    double num_valids;
	    if(o->UsingRANSAC){
	      computeNormalsRANSAC2(SO,luz_dir,o->nLights,n_res,&ransac_weight,o->V,&num_valids);
	      float_image_set_sample(fi_valid_lights,0,x,y,num_valids);
	    }if(o->usingBestThree){
	    
	      r3_t b;
	      double tp;
	      int imax[3];
	      r3x3_t A,A_inv;
	      ConstroiSistema(SO, luz_dir, o->nLights, &A, &b, &tp, imax);
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
	  //   printf("%f %f %f\n,",n_res[0],n_res[1],n_res[2]);
    /* Computes the matrix product {M = A * B}, where {A} has size
      {m x p} and {B} has size {p x n}. The matrix {M} must be disjoint
      from {A} and {B}. */
	    
	    
	      // Resolve o sistema
	      // Calcula o albedo e normaliza a normal:
	      albedo = r3_dir(&snp,&snp)*(o->scaleAlbedo);
	    // if (albedo != 0.0) { r3_scale(1.0/albedo, &snp, &snp); }
	      // Ajusta o peso proporcionalmente ao quadrado do albedo:
	      if(!o->UsingRANSAC){
		double GO[o->nLights];
		computeGO(GO,snp,luz_dir,o->nLights);
		
		peso = computeProb(SO,GO,o->nLights);
		if((peso != +INF) && (peso != -INF) && (!isnan(peso) ) ){
		  if(peso > Lmax) Lmax = peso;	
		}
	      }else{
		peso = ransac_weight;
	      }

	    if (debug) { 
	      int i;
	      fprintf(arq_debug_pixel,"%+16.14f %+16.14f %+16.14f   ",snp.c[0],snp.c[1],snp.c[2]);
	      for(i = 0; i < o->nLights;i++) fprintf(arq_debug_pixel,"%16.14f ",SO[i]);
	      fprintf(arq_debug_pixel,"\n");
	    
	      fclose(arq_debug_pixel); 
	    }
	  } 
	  /* Registra as luzes escolhidas na imagem de seleção de luzes: */
         
	}
    	
	/* Grava normal no mapa de normais: */
        //fprintf(arq_normais_canal,"%d %d %f %f %f\n", x, y, snp.c[0], -snp.c[1], snp.c[2]);
	float_image_set_sample(fi_normais_canal[canal],0,x,y,snp.c[0]);
	float_image_set_sample(fi_normais_canal[canal],1,x,y,-snp.c[1]);
	float_image_set_sample(fi_normais_canal[canal],2,x,y,snp.c[2]);
    	
        /* Grava peso no mapa de pesos: */
       // fprintf(arq_pesos_canal,"%d %d %f\n", x, y, peso);
        //float_image_set_sample(imagem_peso, 0, x, y, peso);
    	
        /* Salva normal e peso para imagens médias de canais: */
	/* !!! Deveria levar em conta o albedo da cena no canal !!! */
        int ip = x + y*nx;
        normais_calculadas[canal][ip] = snp; 
        pesos[canal][ip] = peso;
	
	

        /* Armazena o albedo em {imagem_albedo} (gama visual): */
        double albedo_dimming = 1.0; /* Fator de redução para caso do albedo ser maior que 1.0. */
        float_image_set_sample(imagem_albedo, canal, x, y, albedo_dimming*albedo);
        //int albedo_int = sample_from_intensity(albedo_dimming*albedo, get_maxval(imagem_albedo), VIEW_GAMMA);
    	//set_sample(imagem_albedo,x,y,canal,albedo_int);

        /* Armazena peso na {imagem_peso} (escala linear): */
        //int peso_int = sample_from_intensity(peso, get_maxval(imagem_peso), 1.0);
        //set_sample(imagem_peso, x,y,canal, (short int) peso_int);
	float_image_set_sample(imagem_peso, canal, x, y, peso);
	float_image_set_sample(imagem_peso_canal, 0, x, y, peso);
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
	  fprintf(stderr,"[%d][%9d] of [%9d] - %6.6f%% - %6.6f n/s   - %02d h %02d m %02d s       \n",
		canal,contador,total_compute_pixels, contador*100.0/total_compute_pixels,normals_per_sec,hour,min,sec);
	}
      }
    }
    
    //calcula mapa de pesos
    //First misc weights
    
    //Now official weights
    
    if (processa_canal) {
      for(y = 0 ;y < ny; y++){
	for(x = 0; x < nx; x++){
	   float weight;
	  if( !o->UsingRANSAC) {
	      if(!o->UsingAngleWeights){
		float TINY_LOGPROB = log(1.0e-10);
		float TINY_DOT = 0.01;
		
		float prob = float_image_get_sample(imagem_peso_canal,0,x,y);
	      
		if(prob == +INF) weight = 1.0;
		else if( prob == -INF) weight = 0.0;
		else if( isnan(prob) ) weight = 0.0;
		else if( (prob - Lmax) < TINY_LOGPROB ) weight = 0.0;
		else{
		  weight = exp(prob - Lmax);
		}
		
		r3_t norm = (r3_t){{float_image_get_sample(fi_normais_canal[canal],0,x,y), float_image_get_sample(fi_normais_canal[canal],1,x,y), float_image_get_sample(fi_normais_canal[canal],2,x,y)}};
		r3_t view_dir  = (r3_t){{0,0,1}};
		double dot = r3_dot(&norm,&view_dir);
	      if(dot < TINY_DOT){
		weight = 0.0 ;
	      }
	     }else{
	       //weights for the 12-light PS
	       r3_t avg_luz_dir;
	       r3_zero(&avg_luz_dir);
	       int ii;
	       for(ii = 0; ii < o->nLights; ii++){
		 r3_add(&luz_dir[ii],&avg_luz_dir,&avg_luz_dir);
	       }
	       
	       r3_dir(&avg_luz_dir,&avg_luz_dir);
	       r3_t norm = (r3_t){{float_image_get_sample(fi_normais_canal[canal],0,x,y), float_image_get_sample(fi_normais_canal[canal],1,x,y), float_image_get_sample(fi_normais_canal[canal],2,x,y)}};
	       weight = computeWnrm(norm,avg_luz_dir,o->clusterr);
	     }
	  }else{
	    weight = float_image_get_sample(imagem_peso_canal,0,x,y);
	    assert(weight >= 0);
	   
	  }
	  float_image_set_sample(imagem_peso_canal,0,x,y,weight);
	  float_image_set_sample(imagem_peso,canal,x,y,weight);
	}
      }
      
    }

    float_image_write(arq_normais_canal,fi_normais_canal[canal]);
    float_image_write(arq_valid_lights,fi_valid_lights);
    fclose(arq_valid_lights);
    /* Fecha arquivo de normais do canal: */
    fprintf(arq_normais_canal,"\n"); fclose(arq_normais_canal);

    /* Fecha arquivo de pesos do canal: */
    float_image_write(arq_pesos_canal,imagem_peso_canal);
    fprintf(arq_pesos_canal,"\n"); fclose(arq_pesos_canal);
  
    
   
  }
  
  /*Add Ki/Li computing */
  if(o->generateKili){
    for(canal = 0; canal < nc; canal++){
      
      int processa_canal = strchr(o->channels,"RGB"[canal]) != NULL;
      if(!processa_canal) continue;
      for(i = 0; i < o->nLights;i++){
	char *nome_arq_Li = NULL;
	//char *nome_arq_Li = jsprintf("%s_%d_G%d_L%02d.pgm",o->prefix,canal,ind_gab,i);
	char *nome_arq_Li = jsprintf("%s_%d_L%02d.pgm",o->prefix,canal,i);
	char *nome_arq_Ki = NULL;
	//char *nome_arq_Ki = jsprintf("%s_%d_G%d_K%02d.pgm",o->prefix,canal,ind_gab,i);
	char *nome_arq_Ki = jsprintf("%s_%d_K%02d.fni",o->prefix,canal,i);
	fprintf(stderr,"Gerando imagens %s %s \n",nome_arq_Li,nome_arq_Ki);
	float_image_t* fi_Li = float_image_new(1,nx,ny);
	float_image_t* fi_Ki = float_image_new(1,nx,ny);
	for(x = 0; x < nx; x++){
	  for(y = 0; y < ny; y++){
	    
	    if( UsingPolyDataFile){
	      poly_function_t** el_func = NULL;
	      if(o->elevation_funcfile != NULL) el_func = el_functions[canal];
	      poly_function_t** az_func = NULL;
	      if(o->azimuth_funcfile != NULL) az_func = az_functions[canal];
	      poly_function_t** int_func = NULL; 
	      if(o->albedo_funcfile != NULL) int_func = int_functions[canal];
	      double ix, iy;
	      ix = x;
	      iy = y;
	      if(o->transform != NULL){
		ix = ix/o->transform->scale.c[0];
		iy = iy/o->transform->scale.c[1];
		ix = ix + o->transform->offset.c[0]; 
		iy = iy + o->transform->offset.c[1];
		
	      }
	      UpdateLightSourceValues(luz_dir,original_dir,o->nLights,az_func,el_func,int_func,ix,iy,o->invertAlbedo);
	    }
	    double go[o->nLights];      /* Assinatura normalizada de {q}, ou (0..). */
	    double Gmag;
	    r3_t norm ;
	    norm.c[0] = float_image_get_sample(fi_normais_canal[canal],0,x,y);
	    norm.c[1] = float_image_get_sample(fi_normais_canal[canal],1,x,y);
	    norm.c[2] = float_image_get_sample(fi_normais_canal[canal],2,x,y);
	    if(r3_norm(&norm) < 0.01){
	      float_image_set_sample(fi_Li,0,x,y,0);
	      float_image_set_sample(fi_Ki,0,x,y,0);
	      continue;
	    }
	    double GO[o->nLights];
	    
	    computeGO(GO,norm,luz_dir,o->nLights);
	    
	    Gmag = rn_dir(o->nLights,GO,go);
	    if(Gmag == 0){
	      rn_zero(o->nLights,go);
	    }
	    double SO[o->nLights];
	    double so[o->nLights];
	    double Smag;
	    int i_light;
	    for(i_light = 0; i_light < o->nLights; i_light++){
	      SO[i_light] = float_image_get_sample(S[i_light], canal, x, y);	
	    }
	    extrai_assinatura(SO, so,&Smag,o->nLights);
	    double albedo  = EstAlbedo_00(so,Smag,go,Gmag,o->nLights,0.2,0.2,0.2);
	    
	    double Li = go[i]*Gmag;
	    double Ki = SO[i]/(albedo*Li);
// 	    if((isnan(albedo)) || (isnan(Li)) || (isnan(Ki)) ){
// 	      Li = Ki = 0;
// 	    }
	    float_image_set_sample(fi_Li,0,x,y,Li);
	    float_image_set_sample(fi_Ki,0,x,y,Ki);
	    }
	}
	
    
	float_pnm_image_write(nome_arq_Li, fi_Li,FALSE, VIEW_GAMMA, VIEW_BIAS,TRUE,TRUE,FALSE);
	float_image_free(fi_Li);
	free(nome_arq_Li);
	FILE* arq_Ki = fopen(nome_arq_Ki,"wt");
	float_image_write(arq_Ki, fi_Ki);
	fclose(arq_Ki);
	float_image_free(fi_Ki);
	free(nome_arq_Ki);
      }
      //FILE* arq_Li = fopen(nome_arq_Li,"wt");
    

    }
  }
  
  /*This is a test thing, where i compute the value of a white field for the lighting source */
  if(o->generateWhite){
    for (canal = 0; canal < nc; canal++){
      float_image_t* white_field_image[o->nLights];
      float_image_t* light_field_image[o->nLights];
      for(i = 0 ; i  < o->nLights; i++){
	white_field_image[i]  = float_image_new(1,nx,ny);
	light_field_image[i]  = float_image_new(3,nx,ny);
      }
      for(x = 0; x < nx; x++){
	for(y = 0; y < ny; y++){
	  if( UsingPolyDataFile){
	    poly_function_t** el_func = NULL;
	    if(o->elevation_funcfile != NULL) el_func = el_functions[canal];
	    poly_function_t** az_func = NULL;
	    if(o->azimuth_funcfile != NULL) az_func = az_functions[canal];
	    poly_function_t** int_func = NULL; 
	    if(o->albedo_funcfile != NULL) int_func = int_functions[canal];
	    UpdateLightSourceValues(luz_dir,original_dir,o->nLights,az_func,el_func,int_func,x,y,o->invertAlbedo);
	  }
	  double GO[o->nLights];
	  r3_t plain_dir = (r3_t){{0,0,1.0}};
	  computeGO(GO,plain_dir,luz_dir,o->nLights);
	  for(i = 0;i < o->nLights; i++){
	    float_image_set_sample(white_field_image[i],0,x,y,GO[i]);
	    float_image_set_sample(light_field_image[i],0,x,y,luz_dir[i].c[0]);
	    float_image_set_sample(light_field_image[i],1,x,y,luz_dir[i].c[1]);
	    float_image_set_sample(light_field_image[i],2,x,y,luz_dir[i].c[2]);
	  }
	}
      }
      for(i = 0;i < o->nLights; i++){
	char* white_field_filename;
	char *white_field_filename = jsprintf("%s_%d_L%d_white.pgm",o->prefix,canal,i);
	float_pnm_image_write(white_field_filename, white_field_image[i],FALSE, VIEW_GAMMA, VIEW_BIAS,TRUE,TRUE,FALSE);
	char* light_field_filename;
	char *light_field_filename = jsprintf("%s_%d_L%d_light.fni",o->prefix,canal,i);
	FILE* light_field_file = open_write(light_field_filename,TRUE);
	float_image_write(light_field_file,light_field_image[i]);
	fclose(light_field_file);
	float_image_free(white_field_image[i]);
	float_image_free(light_field_image[i]);
      }
    }
  }
  
  
  
  /* Grava e libera a imagem de albedo: */
  char *nome_imagem_albedo = NULL;
  char *nome_imagem_albedo = jsprintf("%s_albedo.ppm", o->prefix);
  fprintf(stderr, "gravando %s ...\n", nome_imagem_albedo);
  float_pnm_image_write(nome_imagem_albedo, imagem_albedo,FALSE, VIEW_GAMMA, VIEW_BIAS,TRUE,TRUE,FALSE);
  
  char *nome_imagem_albedo = jsprintf("%s_albedo.fni", o->prefix);
  //fprintf(stderr, "gravando %s ...\n", nome_imagem_albedo);
  FILE* arq_imagem_albedo = open_write(nome_imagem_albedo,TRUE);
  float_image_write(arq_imagem_albedo, imagem_albedo);
  fclose(arq_imagem_albedo);
  

  /* Grava e libera a imagem de peso por canal: */
  char *nome_imagem_peso = NULL;
  char *nome_imagem_peso = jsprintf("%s_weights.ppm", o->prefix);
  fprintf(stderr, "gravando %s ...\n", nome_imagem_peso);
  float_pnm_image_write(nome_imagem_peso, imagem_peso,FALSE, VIEW_GAMMA, VIEW_BIAS,TRUE,TRUE,FALSE);

  /* Abre arquivo da tabela de médias de normais: */
//   char *nome_arq_normais_media = NULL;
//   char *nome_arq_normais_media = jsprintf("%s%d_normals.fni", nome_arq_prefixo_tabelas, 3);
//   FILE* arq_normais_media = fopen(nome_arq_normais_media,"wt");
  //fprintf(arq_normais_media,"tx = %d\n",nx);
  //fprintf(arq_normais_media,"ty = %d\n",ny);
  
  /* Abre arquivo da tabela de médias de pesos: */
//   char *nome_arq_pesos_media = NULL;
//   char *nome_arq_pesos_media = jsprintf("%s%d_weights.fni", nome_arq_prefixo_tabelas, 3);
//   FILE* arq_pesos_media = fopen(nome_arq_pesos_media,"wt");
//   //fprintf(arq_pesos_media,"tx = %d\n",nx);
//   //fprintf(arq_pesos_media,"ty = %d\n",ny);
//   
//   /* Cria imagem de médias de pesos: */
//   float_image_t *imagem_peso_media = float_image_new(1, nx, ny);
//   float_image_t* fi_normais_media = float_image_new(3, nx, ny);
//   /* Calcula imagens e grava tabelas com médias de normais e pesos dos três canais: */
//   for (y = 0; y < ny; y++) {
//     for (x = 0; x < nx; x++) {
// 
//       int ip = x + nx*y;
//      
//       /* Combina as normais calculadas para os três canais: */
//       /* Deveria levar em conta os respectivos dist_busca. */
//       r3_t normal_media = calcula_media_das_normais
//         ( pesos[0][ip], normais_calculadas[0][ip],
//           pesos[1][ip], normais_calculadas[1][ip],
//           pesos[2][ip], normais_calculadas[2][ip]
//         );
//       //fprintf(arq_normais_media,"%d %d %f %f %f\n", x, y, normal_media.c[0], -normal_media.c[1], normal_media.c[2]);
// 	float_image_set_sample(fi_normais_media,0,x,y,normal_media.c[0]);
// 	float_image_set_sample(fi_normais_media,1,x,y,normal_media.c[1]);
// 	float_image_set_sample(fi_normais_media,2,x,y,normal_media.c[2]);
//       
//       /* Combina os pesos de confiabilidade para os três canais: */
//       /* Deveria levar em conta os respectivos dist_busca. */
//       double peso_media = calcula_media_dos_pesos
//         ( pesos[0][ip], normais_calculadas[0][ip],
//           pesos[1][ip], normais_calculadas[1][ip],
//           pesos[2][ip], normais_calculadas[2][ip],
//           normal_media
//         );
//       
//       //fprintf(arq_pesos_media,"%d %d %f\n", x, y, peso_media);
//       float_image_set_sample(imagem_peso_media, 0, x, y, peso_media);
//     }
//   }
//   
//   /* Grava a imagem de médias de pesos: */
//   char *nome_imagem_peso_media = NULL;
//   char *nome_imagem_peso_media = jsprintf("%s-peso_media.ppm", nome_arq_prefixo_imagens);
//   fprintf(stderr, "gravando %s ...\n", nome_imagem_peso_media);
//   float_pnm_image_write(nome_imagem_peso_media, imagem_peso_media, 1.000, 0.000,TRUE,TRUE,FALSE);
// 
//   float_image_write(arq_pesos_media,imagem_peso_media);
//   float_image_write(arq_normais_media,fi_normais_media);
//   /* Fecha os arquivos de médias: */
//   fclose(arq_normais_media);
//   fclose(arq_pesos_media);
  
  fprintf(stderr, "Concluido!\nO programa rodou com sucesso!\n");

  return 0;

}


