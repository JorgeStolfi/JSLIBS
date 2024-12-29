#define PROG_NAME "compute_normals_from_table"

#define _GNU_SOURCE
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <i2.h>
#include <rn.h>
#include <frgb.h>
#include <values.h>
#include <float_image.h>
#include <float_pnm_image_io.h>
#include <argparser.h>
#include "normais.h"
#include "imagem_vetores.h"
#include "tabela.h"
#include "super_tabela.h"
#include "hash.h"
#include "kdtree.h"
#include <time.h>
#include <stdint.h>
#include <jsfile.h>
#include <string.h>
#include <estlogprob.h>
#include <assert.h>
#include <limits.h>
#include <sys/times.h>
#include <unistd.h>
#include <sys/time.h>
#include <sys/resource.h>
#include <fast_hash.h>
#include <analytic_inversion.h>
#include <analytic_alpha.h>

#define PROG_HELP \
  PROG_NAME " \\\n" \
  "    -nLights {NUM} \\\n" \
  "    -prefix {FILE_PREFIX} \\\n" \
  "    [ -refNormalMap {FILE_NAME} ] \\\n" \
  "    [ -channels {R | G | B| RG | RB | GB | RGB } \\\n" \
  "    -gamma {NUM} [ -gray ] \\\n" \
  "    -sceneImages {FILENAME}...  \\\n" \
  "    [ -rectangle {XMIN] {XMAX] {YMIN] {YMAX} ] \\\n" \
  "    [ -sigma {NUM} \\\n" \
  "    [ -logProbFunction {NUM} ] \\\n" \
  "    [ -albedoFunction {NUM} ] \\\n" \
  "    { -brute | -hash | -bruteAndHash } \\\n" \
  "    { -nsig | -alpha | -nsigOrAlpha } \\\n" \
  "    -tableFiles {FILENAME}... \\\n" \
  "    [ -mask {FILENAME} ] \\\n" \
  "    [ -debug {NP} {H[1]} {V[1]} ... {H[NP]} {V[NP]} ] \\\n" \
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

#define PROGRESS_METER 1

typedef enum { NSIG_APENAS = 0, NSIG_OR_ALPHA = 1, ALPHA_APENAS = 2  } alphaOption_t;
typedef enum { BRUTA_APENAS = 1, HASH_APENAS = 2, BRUTA_E_HASH = 3 , SUPER_TABLE = 4, KDTREE = 5, RANSAC_DEBURRO = 6, FAST_HASH = 7, FAST_RANSAC = 8 , ANALY_INVERSION = 9, ANALY_ALPHA = 10} opcao_hash_t;

typedef struct options_t {
	/* General parameters: */
  	int nLights;
  	char *prefix;
  	char *refNormalMap;          /* Reference normal map, or NULL if not specified. */
  	int channel;              /* Channel number to be processed 0-2 (RGB)*/
  	int hMin, hMax, vMin, vMax;  /* Rectangle of interest, or {[-BIG .. +BIG]^2} (PNM indices). */
	bool_t computeKiLi;
	bool_t computeDist;
   	/* Input image parameters: */
  	double gamma;
	char* mask;
  	double bias;
  	bool_t gray;
  	char** scene_image_name;  /* Scene image filenames are {scene_image_name[0..nLights-1]} */
  	/* Observation matching parameters: */
  	double sigma;
  	double omg0;
  	double omg1;
  	int logProbFunction;
  	int albedoFunction;
  	int hashOption;
  	int alphaOption;
  	/* Bucket Grid only*/
  	int gridSize;
  	int showBucketsEPS;
  	int showBucketsPPM;
  	int showBucketsPLT;
  	int showBucketData;
  	int performanceTestOnly;
  	/*Super Table Options*/
  	int numSubsets;
  	int subsetSize;
  	int** subsets;
	/*Analytic Alpha options*/
	int nSteps;
	double tolAlbedo;
	double tolNormal;
  	/* Signature table parameters: */
  	/* Debugging */
  	int nDebug;
  	i2_t *debug;  /* PNM-style indices of pixels to debug. */
  	/*Gauge Parameters*/
	int nGauges;
  	frgb_t* gaugeAlbedo;
	char** gaugeTag;
	/*for RANSAC*/
	double RansacProbOutlier;
	int RansacFastNum;
	int** RansacLightSeq;
	  
} options_t;

options_t *parse_options(int argc, char **argv);

double user_cpu_time_usec(void);

double user_cpu_time_usec(void){
	struct tms buf;
	(void)times(&buf);
	return(1000000.0 * ((double) buf.tms_utime)/((double)sysconf(_SC_CLK_TCK)));
}

void findMinMaxMean(float_image_t* img, double* pmin, double* pmax, double* pmean,double* prms);
void findMinMaxMean(float_image_t* img, double* pmin, double* pmax, double* pmean,double* prms){
	int NC,NX,NY;
	NC = img->sz[0];
	NX = img->sz[1];
	NY = img->sz[2];
	float total = NC*NX*NY;
	int c,x,y;
	
	double min = +INF;
	double max = -INF;
	double mean = 0;
	double rms = 0;
	for(c = 0; c < NC; c++){
		for(x = 0; x < NX; x++){
			for(y = 0; y < NY; y++){
				float val = float_image_get_sample(img,c,x,y);
				if(min > val) min = val;
				if(max < val) max = val;
				mean+= val;
				rms+=val*val;
			}
		}
	}
	mean = mean/(double)total;
	rms =  sqrt(rms/(double)total);
	*pmean = mean;
	*pmax = max;
	*pmin = min;
	*prms = rms;

}

options_t *parse_options(int argc, char **argv)
{
	argparser_t *pp = argparser_new(stderr, argc, argv);
	argparser_set_help(pp, PROG_HELP);
  	argparser_set_info(pp, PROG_INFO);
  	argparser_process_help_info_options(pp);
	
  	options_t *o = (options_t *)malloc(sizeof(options_t));
	
  	if (argparser_keyword_present(pp, "-logProbFunction")) {
    		o->logProbFunction = argparser_get_next_int(pp, 0, 999);
  	} else {
    		o->logProbFunction = 0;
  	}
  	fprintf(stderr, "  -logProbFunction %d \\\n", o->logProbFunction);

  	if (argparser_keyword_present(pp, "-albedoFunction")) {
    		o->albedoFunction = argparser_get_next_int(pp, 0, 999);
  	} else {
    		o->albedoFunction = 0;
  	}
  	fprintf(stderr, "  -albedoFunction %d \\\n", o->albedoFunction);

  	if (argparser_keyword_present(pp, "-sigma")) {
    		o->sigma = argparser_get_next_double(pp, 0.0001, 10.0000);
  	} else {
    		o->sigma = 0.05;
  	}
  	fprintf(stderr, "  -sigma %lf \\\n", o->sigma);

  	if (argparser_keyword_present(pp, "-omega")) {
    		o->omg0 = argparser_get_next_double(pp, 0.0001, 10.0000);
    		o->omg1 = argparser_get_next_double(pp, 0.0001, 10.0000);
  	} else {
    		o->omg0 = 0.05;
    		o->omg1 = 0.005;
  	}
  	fprintf(stderr, "  -omega %lf %lf \\\n", o->omg0, o->omg1);

  	if (argparser_keyword_present(pp, "-channel")){
		char* channels = argparser_get_next(pp);
		char* allowed_channels = "RGB";
		int ind_chan;
		o->channel = -1;
		for(ind_chan = 0; ind_chan < strlen(allowed_channels); ind_chan++){
			//OBS: Only teh FIRST channel will be selected. Others are ignored	
			if( strchr(channels,allowed_channels[ind_chan]) != NULL ){ 
				if(o->channel == -1){
					o->channel = ind_chan;
				}
			}
		}
		if(o->channel == -1){
			argparser_error(pp, "Must specify an valid channel {R,G,B} !");
		}
	}else{
		//Default - Channel 0
		o->channel = 0;
	}
    		
  	
  
  	fprintf(stderr, "  -channel %d \\\n", o->channel);

  	if (argparser_keyword_present(pp, "-rectangle")) {
    		o->hMin = argparser_get_next_int(pp, INT_MIN, INT_MAX);
    		o->hMax = argparser_get_next_int(pp, o->hMin, INT_MAX);
    		o->vMin = argparser_get_next_int(pp, INT_MIN, INT_MAX);
    		o->vMax = argparser_get_next_int(pp, o->vMin, INT_MAX);
    		fprintf(stderr, "  -rectangle %d %d  %d %d \\\n", o->hMin, o->hMax, o->vMin, o->vMax);
  	} else {
    		o->hMin = o->vMin = INT_MIN;
    		o->hMax = o->vMax = INT_MAX;
  	}

	if (argparser_keyword_present(pp, "-nGauges")) {
		o->nGauges = argparser_get_next_int(pp, 1, 100);
	} else {
		o->nGauges = 1;
	}	
	

	if (argparser_keyword_present(pp, "-gaugeTag")) {
		//end of party , time to read gauge tags all again !
		o->gaugeTag = (char**)malloc(sizeof(char*)*(o->nGauges));
		int ii;
		for(ii = 0; ii < o->nGauges; ii++){
		  o->gaugeTag[ii] = argparser_get_next(pp);
		}
	}else{
		o->gaugeTag = (char**)malloc(sizeof(char*));
		o->gaugeTag[0] = "0";
	}
	

  	argparser_get_keyword(pp, "-nLights");
  	o->nLights = argparser_get_next_int(pp, 3, 1000);
  	fprintf(stderr, "  -nLights %d \\\n", o->nLights);


	argparser_get_keyword(pp, "-prefix");
  	o->prefix = argparser_get_next(pp);
  	fprintf(stderr, "  -prefix %s \\\n", o->prefix);
  
  	argparser_get_keyword(pp, "-sceneImages");
  	fprintf(stderr, "Imagens:\n");
  	o->scene_image_name = (char**) malloc(sizeof(char*)*(o->nLights));
  	int ind;
  	for (ind = 0; ind < (o->nLights); ind++){
    		o->scene_image_name[ind] = argparser_get_next(pp);
    		fprintf(stderr, "    S[%02d] = %s\n", ind, o->scene_image_name[ind]);
  	}

	o->computeKiLi = FALSE;
  	if (argparser_keyword_present(pp, "-computeKiLi")){
		o->computeKiLi = TRUE;
	}
  
	o->computeDist = FALSE;
  	if (argparser_keyword_present(pp, "-computeDist")){
		o->computeDist = TRUE;
	}

  	if (argparser_keyword_present(pp, "-gridSize")) {
    		o->gridSize = argparser_get_next_int(pp, 1, 1000000);
  	} else {
    		o->gridSize = -1;
  	}
  	fprintf(stderr,"   -gridSize %d \\\n",o->gridSize);

  
  	if(argparser_keyword_present(pp,"-showBucketsPPM")){
		o->showBucketsPPM = 1;
		fprintf(stderr,"    -showBucketsPPM \\\n");
  	} 
  	else{
		o->showBucketsPPM = 0;
  	}
  	if(argparser_keyword_present(pp,"-showBucketsEPS")){
		o->showBucketsEPS = 1;
		fprintf(stderr,"    -showBucketsEPS \\\n");
  	} 
  	else{
		o->showBucketsEPS = 0;
  	}
  
  	if(argparser_keyword_present(pp,"-showBucketsPLT")){
		o->showBucketsPLT = 1;
		fprintf(stderr,"    -showBucketsPLT \\\n");
  	} 
  	else{
		o->showBucketsPLT = 0;
  	}
  
  	o->showBucketData = 1;
  
	o->mask = NULL;
	if (argparser_keyword_present(pp,"-mask")){
		o->mask = argparser_get_next(pp);
	}

  	if(argparser_keyword_present(pp,"-performanceTestOnly") ){
		o->performanceTestOnly = 1;
		fprintf(stderr,"    -performanceTestOnly\\\n");
   	}else{
		o->performanceTestOnly = 0;
   	}

  	if(argparser_keyword_present(pp,"-numSubsets") ){
		o->numSubsets = argparser_get_next_int(pp, 1, 1000000);
		fprintf(stderr,"    -numSubsets %d\\\n",o->numSubsets);
		argparser_get_keyword(pp,"-subsetSize");
		o->subsetSize = argparser_get_next_int(pp, 3, 100000);
		fprintf(stderr,"    -subsetSize %d\\\n",o->subsetSize);
		int random_subsets = !argparser_keyword_present(pp,"-subsets");
		if(random_subsets){
			fprintf(stderr,"USING Random Subsets\n");
		}else{
			fprintf(stderr,"USING User-Defined Subsets\n");
		}
		o->subsets = (int**)malloc(sizeof(int*)*(o->numSubsets));
		srand ( time(NULL) );
		for(ind = 0; ind < o->numSubsets;ind++){
			o->subsets[ind] = (int*)malloc(sizeof(int)*(o->subsetSize));
		}
		if(random_subsets){
			stGenerateSubsets(o->subsets,o->subsetSize,o->nLights,o->numSubsets);
		}else{
			for(ind = 0; ind < o->numSubsets; ind ++){
				int ind2;
				for(ind2 = 0; ind2 < o->subsetSize; ind2++){
					o->subsets[ind][ind2] = argparser_get_next_int(pp,0, 100000);
				}
			}
		}
		int n_uso[o->nLights];
		for(ind = 0; ind < o->nLights; ind++){
			n_uso[ind] = 0;
		}
		for(ind = 0; ind < o->numSubsets; ind ++){
			fprintf(stderr,"Subset[%d]: ",ind);
			int ind2;
			for(ind2 = 0; ind2 < o->subsetSize; ind2++){
				fprintf(stderr,"%03d ",o->subsets[ind][ind2]);
				n_uso[o->subsets[ind][ind2]]++;
			}
			fprintf(stderr,"\n");
		}
		fprintf(stderr,"Light Usage:\n");
		for(ind = 0; ind < o->nLights; ind++){
			fprintf(stderr,"Light[%03d]: %5d\n",ind,n_uso[ind]);
		}
	}else{
		o->numSubsets = 0;
		o->subsetSize = 0;
		o->subsets = NULL;
  	}
	

	o->gaugeAlbedo = (frgb_t*) malloc(sizeof(frgb_t)*o->nGauges);

	int ind_gab;
	int test_albedo = argparser_keyword_present(pp, "-gaugeAlbedo");
	for(ind_gab = 0; ind_gab < o->nGauges; ind_gab++){
		if (test_albedo) {
    			o->gaugeAlbedo[ind_gab].c[0] = argparser_get_next_double(pp, 0.001, 1000.000);
    			o->gaugeAlbedo[ind_gab].c[1] = argparser_get_next_double(pp, 0.001, 1000.000);
    			o->gaugeAlbedo[ind_gab].c[2] = argparser_get_next_double(pp, 0.001, 1000.000);
  		} else {
    			o->gaugeAlbedo[ind_gab].c[0] = o->gaugeAlbedo[ind_gab].c[1] = o->gaugeAlbedo[ind_gab].c[2] = 1.000;
  		}
  		fprintf(stderr, "  -gaugeAlbedo %5.3lf %5.3lf %5.3lf \\\n", o->gaugeAlbedo[ind_gab].c[0],o->gaugeAlbedo[ind_gab].c[1],o->gaugeAlbedo[ind_gab].c[2]);
	}
  	
	argparser_get_keyword(pp, "-gamma");
  	o->gamma = argparser_get_next_double(pp, 0.100, 9.000);
  	o->bias = (o->gamma == 1.000 ? 0.000 : VIEW_BIAS); /* Hack... */
  	fprintf(stderr, "Gamma das imagens:%lf bias: %lf\n", o->gamma, o->bias);

  	o->gray = argparser_keyword_present(pp, "-gray");

  	if (argparser_keyword_present(pp, "-refNormalMap")) {
    		o->refNormalMap = argparser_get_next(pp);
    		fprintf(stderr, "  -refNormalMap %s \\\n", o->refNormalMap);
  	} else {
    		o->refNormalMap = NULL;
  	}

  	/* Analisa opção de hash: */
  	o->hashOption = BRUTA_APENAS;
  	if (argparser_keyword_present(pp, "-brute")){
    		o->hashOption = BRUTA_APENAS;
    		fprintf(stderr, "Força bruta apenas.\n");
  	} else if (argparser_keyword_present(pp, "-hash")){
    		o->hashOption = HASH_APENAS;
    		fprintf(stderr, "Hash (bucket grid) apenas.\n");
  	} else if (argparser_keyword_present(pp, "-bruteAndHash")){
    		o->hashOption = BRUTA_E_HASH;
    		fprintf(stderr, "Força bruta e hash, comparando os resultados.\n");
  	} else if (argparser_keyword_present(pp, "-superTable")){
    		o->hashOption = SUPER_TABLE;
    		fprintf(stderr, "Utilizando SuperTable.\n");
	} else if (argparser_keyword_present(pp, "-kdTree")){
    		o->hashOption = KDTREE;
    		fprintf(stderr, "Utilizando KDTree.\n");
  	} else if (argparser_keyword_present(pp, "-ransacBurro")){
    		o->hashOption = RANSAC_DEBURRO;
    		fprintf(stderr, "Utilizando Ransac de Burro.\n");
		//o->V = argparser_get_next_double(pp,0,10000);
		o->RansacProbOutlier = argparser_get_next_double(pp,0,10000);
	} else if (argparser_keyword_present(pp, "-fastHash")){
		o->hashOption = FAST_HASH;
    		fprintf(stderr, "Utilizando Fast Hash.\n");
	} else if (argparser_keyword_present(pp, "-fastRansac")){
    		o->hashOption = FAST_RANSAC;
    		fprintf(stderr, "Utilizando Fast Ransac.\n");
		//o->V = argparser_get_next_double(pp,0,10000);
		o->RansacProbOutlier = argparser_get_next_double(pp,0,10000);
		int limit_sets =  o->nLights*(o->nLights-1)*(o->nLights-2)/6.0;
		//o->RansacFastNum = argparser_get_next_int(pp,limit_sets,limit_sets);
		
		
		if(argparser_keyword_present(pp, "-ransacSeq")){
		  o->RansacFastNum = argparser_get_next_int(pp,1,limit_sets);
		  o->RansacLightSeq = (int**)malloc(sizeof(int*)*(o->RansacFastNum));
		  int i;
		  fprintf(stderr,"RANSAC USER DEFINED SEQUENCE\n");
		  for(i = 0; i < o->RansacFastNum; i++){
		      o->RansacLightSeq[i] = (int*)malloc(sizeof(int)*3);
		      o->RansacLightSeq[i][0] = argparser_get_next_int(pp,0,o->nLights -1);
		      o->RansacLightSeq[i][1] = argparser_get_next_int(pp,0,o->nLights -1);
		      o->RansacLightSeq[i][2] = argparser_get_next_int(pp,0,o->nLights -1);
		      
		  }
		}else{
		  fprintf(stderr,"RANSAC AUTOSEQUENCE\n");
		  o->RansacFastNum = argparser_get_next_int(pp,1,limit_sets);
		  o->RansacLightSeq = (int**)malloc(sizeof(int*)*(o->RansacFastNum));
		  int i,j,k;
		  int count = 0;
		  for(i = 0; i < o->nLights; i++){
		    for(j = i+1; j < o->nLights; j++){
		      for(k = j+1; k < o->nLights; k++){
			o->RansacLightSeq[count] = (int*)malloc(sizeof(int)*3);
			o->RansacLightSeq[count][0] = i;
			o->RansacLightSeq[count][1] = j;
			o->RansacLightSeq[count][2] = k;
			count++;
			if(count == o->RansacFastNum) break;
		      }
		      if(count == o->RansacFastNum) break;
		    }
		    if(count == o->RansacFastNum) break;
		  }
		}
		
		fprintf(stderr,"RANSAC SEQUENCE\n");
		int i;
		for(i = 0; i < o->RansacFastNum; i++){
		  fprintf(stderr,"[%03d]: %02d %02d %02d\n",i,o->RansacLightSeq[i][0],o->RansacLightSeq[i][1],o->RansacLightSeq[i][2]);
		}
		fprintf(stderr,"END SEQUENCE\n");
		
		
	}else if( argparser_keyword_present(pp, "-analyticInversion")){
	  o->hashOption = ANALY_INVERSION;
	  fprintf(stderr, "Utilizando Analytic Inversion.\n");
	}else if( argparser_keyword_present(pp, "-analyticAlpha")){
	  o->hashOption = ANALY_ALPHA;
	  fprintf(stderr, "Utilizando Analytic Alpha Inversion.\n");
	  o->tolAlbedo = 0.001;
	  if(argparser_keyword_present(pp, "-tolAlbedo")){
	    o->tolAlbedo = argparser_get_next_double(pp,0.0,1.0);
	  }
	  o->tolNormal = M_PI/180.0;
	  if(argparser_keyword_present(pp, "-tolNormal")){
	    o->tolNormal = argparser_get_next_double(pp,0.0,M_PI);
	  }
	  o->nSteps = 10;
	  if(argparser_keyword_present(pp, "-nSteps")){
	    o->nSteps = argparser_get_next_int(pp,0,10000);
	  }
	}else {
    		argparser_error(pp, "Must specify a hash/brute option");
  	}

  	/* Analisa opção de estimador de probabilidade usada na busca: */
  	o->alphaOption = NSIG_APENAS;
  	if (argparser_keyword_present(pp, "-nsig")) {
    		o->alphaOption = NSIG_APENAS;
    		fprintf(stderr, "Usando apenas comparação de assinaturas normalizadas.\n");
  	} else if (argparser_keyword_present(pp, "-alpha")) {
    		o->alphaOption = ALPHA_APENAS;
    		fprintf(stderr, "Usando apenas análise de Bayes.\n");
  	} else if (argparser_keyword_present(pp, "-nsigOrAlpha")) {
    		o->alphaOption = NSIG_OR_ALPHA;
    		fprintf(stderr, "Usando comp. de assinaturas, se erro grande usa análise de Bayes.\n");
  	} else {
    		argparser_error(pp, "Must specify an alpha/nsig option");
  	}

	/* Hash só pode ser usado com comparação de assinaturas normalizadas: */
  	if ((o->hashOption != BRUTA_APENAS) && (o->alphaOption == ALPHA_APENAS)) {
    		argparser_error(pp, "Hash é inútil quando se usa análise de Bayes apenas");
  	}

  	/* Opções de debug: */
  	if (argparser_keyword_present(pp, "-debug")) {
    		o->nDebug = argparser_get_next_int(pp, 0, INT_MAX);
    		o->debug = (i2_t*) malloc(o->nDebug * sizeof(i2_t));
    		int id;
    		for (id = 0; id < o->nDebug; id++) {
      			o->debug[id].c[0] = argparser_get_next_int(pp, 0, INT_MAX);
      			o->debug[id].c[1] = argparser_get_next_int(pp, 0, INT_MAX);
      			fprintf(stderr, "  debugging pixel [%d %d] (PNM indices)\n", o->debug[id].c[0], o->debug[id].c[1]);
    		}
  	} else {
    		o->nDebug = 0;
    		o->debug = NULL;
  	}
  	argparser_finish(pp);
	return o;
}

int main(int argc, char** argv){
	options_t *o = parse_options(argc, argv);

	/* Lê as imagens da cena e define o tamanho {nx,ny}: */
  	float_image_t  *S[o->nLights];
  	int nx = -1, ny = -1;
  	int i;
  	for(i = 0; i < o->nLights; i++){
    		fprintf(stderr, "Abrindo arquivo[%d] %s ... \n",i,o->scene_image_name[i]);
    		float_image_t *im = float_pnm_image_read(o->scene_image_name[i],FALSE, o->gamma, o->bias, TRUE,TRUE,FALSE);
    		//assert(im->sz[0]== 3);
    		if (i == 0){
			nx = im->sz[1]; ny = im->sz[2];
		}else{
			if ((nx != im->sz[1]) || (ny != im->sz[2])){
				fprintf(stderr, "Imagem S[%d] com tamanho inconsistente!\n", i); exit(1);
			}
      		}
		if(o->channel >= im->sz[0]){
			fprintf(stderr, "Imagem S[%d] não apresenta canal %d!\n", i,o->channel); exit(1);
		}
     		S[i] = im;
     	
  	}
  	fprintf(stderr, "Imagens da cena Lidas.\n");

  	/* Ajusta o retângulo a calcular: */
  	if (o->hMin < 0)      { o->hMin = 0; }
  	if (o->hMax > nx - 1) { o->hMax = nx - 1; }
  	assert(o->hMin <= o->hMax);
  	if (o->vMin < 0)      { o->vMin = 0; }
  	if (o->vMax > ny - 1) { o->vMax = ny - 1; }
  	assert(o->vMin <= o->vMax);

  	/* Compute the toal number of pixel to be computed: */
  	int nx_comp = o->hMax - o->hMin + 1;
  	int ny_comp = o->vMax - o->vMin + 1;
  	int total_compute_pixels = nx_comp*ny_comp;
  	fprintf(stderr, "Total de iteracoes a executar: %d\n", total_compute_pixels);
  	fprintf(stderr, "Processando\n");
	
	/* Lê imagem com normais de referência, se houver: */
  	//imagem_vetores_t *nref = NULL;
	float_image_t* nref = NULL;
  	if (o->refNormalMap != NULL) {
    		//nref = le_imagem_vetores(o->refNormalMap);
		FILE* arq_nref = fopen(o->refNormalMap,"rt");
		assert(arq_nref != NULL);
		nref = float_image_read(arq_nref);
		fclose(arq_nref);
  	}

	/*Lê imagem de mascara se houver*/
	float_image_t* imagem_mask = NULL;
	if((o->mask != NULL) && (strcmp(o->mask,"NONE") != 0) ){
		imagem_mask = float_pnm_image_read(o->mask,TRUE,1.0, 0.0, TRUE,TRUE,FALSE);
	}
	
	
	
	/* Criando imagens de saída: */
  	float_image_t  *imagem_gab_select = float_image_new(1, nx, ny);
  	float_image_t  *imagem_albedo = float_image_new(1, nx, ny);
	float_image_t  *imagem_euclid_evals = NULL;
	float_image_t  *imagem_scans = NULL;
	float_image_t* fi_count_valids  = NULL;
	

	if((o->hashOption == HASH_APENAS) || (o->hashOption == BRUTA_E_HASH) || (o->hashOption == KDTREE) ){
		imagem_euclid_evals = float_image_new(1, nx, ny);
		imagem_scans = float_image_new(1, nx, ny);
	}
	
	if(o->hashOption == FAST_RANSAC){
	  fi_count_valids = float_image_new(1, nx, ny);
	}
        fprintf(stderr, "Imagens de saída criadas.\n");

	
	
	
	
	
  	/* Tabelas de normais e pesos, indexadas por {[canal][x + y*nx]}: */
  	r3_t *normais_calculadas;
	
  	double *logPrSG_busca; 
	
 	estima_log_prob_S_G_t *estLogPrSG = escolhe_estLogPrSG(o->logProbFunction);
  	estima_albedo_t *estAlbedo = escolhe_estAlbedo(o->albedoFunction);

  	/* Marca tempo de processamento real: */
  	time_t tempo_start = time(NULL);

  	struct tm *timeinfo = localtime(&tempo_start);

  	/* Loop sobre os canais: */
  	
	int canal = o->channel;
  	int x,y;	
	int ind_gab;
	int64_t stats_scan;
	int64_t stats_scan_gab[o->nGauges];
	int64_t stats_euclid;
	int64_t stats_euclid_gab[o->nGauges];
	int linha;
	int linha_gab[o->nGauges];
	int linha_bruta;
  	int linha_bruta_gab[o->nGauges];
	int linha_hash;
  	int linha_hash_gab[o->nGauges];
	/* Abre arquivo da tabela de normais do canal: */
	
	char *nome_arq_normais_canal = NULL;
	char *nome_arq_normais_canal = jsprintf("%s_%d_normals.fni", o->prefix, canal);
	FILE* arq_normais_canal = fopen(nome_arq_normais_canal,"wt");
	char* nome_arq_gab_select_canal = NULL;
	char *nome_arq_gab_select_canal = jsprintf("%s_%d_gabselect.txt", o->prefix, canal);
	FILE* arq_gab_select_channel = fopen(nome_arq_gab_select_canal,"wt");
	float_image_t* fi_normais_canal = float_image_new(3,nx,ny);
	//Abre imagens com normais para cada gabarito
	float_image_t* fi_normais_canal_gab[o->nGauges];
	int* indice_normais_gab[o->nGauges];
	for (ind_gab = 0; ind_gab < o->nGauges; ind_gab++){
		fi_normais_canal_gab[ind_gab] = float_image_new(3,nx,ny);
		indice_normais_gab[ind_gab] = malloc(nx*ny*sizeof(int));
	}

	/* Abre arquivos de saída: */
  	char *nome_im_gab_select = NULL;
  	char *nome_im_gab_select = jsprintf("%s_%d_gab_select.pgm",o->prefix,canal);

	/* Abre arquivo da tabela de erros de normais por canal: */
    	FILE* arq_err_normais_canal = NULL;
	float_image_t* fi_err_normais_canal = NULL;
    	if(nref != NULL){
		char *nome_arq_err_normais_canal = NULL;
		char *nome_arq_err_normais_canal = jsprintf("%s_%d_norm_errors.fni", o->prefix, canal);
		arq_err_normais_canal = fopen(nome_arq_err_normais_canal,"wt");
		fi_err_normais_canal = float_image_new(1,nx,ny);
	}
	fprintf(stderr,"---------------------------------------------------------------------");
	fprintf(stderr, "INICIO DO PROCESSAMENTO DO CANAL %d: %s", canal, asctime(timeinfo));
	/* Aloca tabelas de pesos e normais calculadas para este canal: */
    	normais_calculadas = malloc(nx*ny*sizeof(r3_t));
	logPrSG_busca = malloc(nx*ny*sizeof(double));
	/* Inicializa as estatísticas de busca deste canal: */
	int stats_brut_gab[o->nGauges];
	int stats_hash_gab[o->nGauges];
	int stats_brut = 0;
	int stats_hash = 0;
	stats_scan = 0;
	stats_euclid = 0;
	for(ind_gab = 0; ind_gab < o->nGauges; ind_gab++){
    		stats_brut_gab[ind_gab] = 0;
    		stats_hash_gab[ind_gab] = 0;
		stats_scan_gab[ind_gab] = 0;
		stats_euclid_gab[ind_gab] = 0;
	}
	/* Carrega a tabelas {tab_vector} para este canal: */
    	Tabela** tab = NULL;
	int num_linhas[o->nGauges];
	tab = (Tabela**)malloc(sizeof(Tabela*)*o->nGauges);
	for(ind_gab = 0; ind_gab < o->nGauges;  ind_gab++){
		char *nome_arq_table = NULL;
		tab[ind_gab] = NULL;
		//char *nome_arq_table = jsprintf("%s_%d_G%d_TableData.txt", o->prefix, canal,ind_gab);
		char *nome_arq_table = jsprintf("%s_%d_G%s_TableData.txt", o->prefix, canal,o->gaugeTag[ind_gab]);
		fprintf(stderr,"Loading signature table %s G%d...",nome_arq_table,ind_gab);
		LoadTable(nome_arq_table,&(tab[ind_gab]));
		if(tab != NULL){
			num_linhas[ind_gab] = get_num_linhas(tab[ind_gab]);
			fprintf(stderr,"OK !\n");
		}else{
			fprintf(stderr,"FAILED !\n");
		}
		assert(tab[ind_gab] != NULL);
	}
	fprintf(stderr, "Tabela carregadas.\n");
	


	/* Gera a tabela de hash, se solicitado: */
	bucketGrid* bg[o->nGauges];
	int tam_grid[o->nGauges];
	for(ind_gab = 0; ind_gab < o->nGauges; ind_gab++){
		tam_grid[ind_gab] = 0;
		bg[ind_gab] = NULL;
		if ((o->hashOption != SUPER_TABLE) && (o->hashOption != BRUTA_APENAS) && (o->hashOption != KDTREE)) {
			bg[ind_gab] = CriaBucketGrid(tab[ind_gab],o->gridSize);
			fprintf(stderr, "Bucket grid gerada.\n");
			tam_grid[ind_gab] = get_tam_grid(bg[ind_gab]);
		}
	}
	
	/*Gera Super Tabelas se Solicitado*/
   	SuperTabela* ST_set[o->nGauges][o->numSubsets];
	float_image_t* imagem_subset = NULL;
   	if (o->hashOption == SUPER_TABLE) {
		imagem_subset = float_image_new(1,nx,ny);
		for(ind_gab = 0; ind_gab < o->nGauges; ind_gab++){
			int ind;
			for(ind = 0; ind < o->numSubsets; ind++){
				fprintf(stderr,"CREATING SUPERTABLE %d G%d\n",ind,ind_gab);
				ST_set[ind_gab][ind] = stCreate(
					o->nLights, 
    					canal, 
    					o->subsetSize,
    					o->subsets[ind],
    					o->gridSize,
					tab[ind_gab]
				);
			}

   		}
	}
	
	fast_hash_t* FH_set[o->nGauges];
	if (o->hashOption == FAST_HASH) {
	    for(ind_gab = 0; ind_gab < o->nGauges; ind_gab++){
	      fprintf(stderr,"LOADING FASTHASH G%d\n",ind_gab);
	      char* fh_name = NULL;
	      char *fh_name = jsprintf("%s_%d_G%s_FastHash.txt", o->prefix, canal,o->gaugeTag[ind_gab]);
	      FILE* fh_file = open_read(fh_name,TRUE);
	      FH_set[ind_gab] = LoadFastHash(fh_file);
	      free(fh_name);
	      PrintFastHash(stderr, FH_set[ind_gab]);
	      fclose(fh_file);
	    }
	}
	
	fast_hash_t* FH_ransac_set[o->nGauges][o->RansacFastNum];
	if (o->hashOption == FAST_RANSAC) {
	  for(ind_gab = 0; ind_gab < o->nGauges; ind_gab++){
	      int ind_set;
	      for(ind_set = 0; ind_set < o->RansacFastNum; ind_set++){
		fprintf(stderr,"LOADING FASTHASH SUBSET G%d\n",ind_gab);
		char* fh_name = NULL;
		char *fh_name = jsprintf("%s_%d_S%d_G%s_FastHash.txt", o->prefix, canal,ind_set,o->gaugeTag[ind_gab]);
		FILE* fh_file = open_read(fh_name,TRUE);
		FH_ransac_set[ind_gab][ind_set] = LoadFastHash(fh_file);
		free(fh_name);
		PrintFastHash(stderr, FH_ransac_set[ind_gab][ind_set]);
		fclose(fh_file);
	      }
	  }
	}
	
	kdtree_t* kdTree[o->nGauges];
	if(o->hashOption == KDTREE){
		fprintf(stderr,"Generating KdTrees\n");
		kdTree[ind_gab] = NULL;
		for(ind_gab = 0; ind_gab < o->nGauges; ind_gab++){
			fprintf(stderr,"Generating KDTree %d...",ind_gab);
			kdTree[ind_gab] = buildKdTree(tab[ind_gab]);
			fprintf(stderr,"OK\n");
		}
		fprintf(stderr,"Tree generation finished\n");
	}
	
	analytic_inversion_t* ai = NULL;
	if (o->hashOption == ANALY_INVERSION || o->hashOption == ANALY_ALPHA) {
	      fprintf(stderr,"LOADING ANALYTIC MODEL\n");
	      char* an_name = NULL;
	      char *an_name = jsprintf("%s_%d_G%s_Imodel.txt", o->prefix, canal,o->gaugeTag[0]);
	      FILE* an_file = open_read(an_name,TRUE);
	      ai = analytic_inversion_read(an_file);
	      assert(ai->m == o->nLights);
	      free(an_name);
	      fclose(an_file);
	 }
	
	
	/* Estrutura para debug de pixels: */
	normal_debug_opts_t *dbopt = notnull(malloc(sizeof(normal_debug_opts_t)), "no mem");
	/* Loop sobre pixels da cena: */
	int contador = 0;
	time_t last_tempo = time(NULL);
	fprintf(stderr,"\n");
	/* Marca o CPU time para o canal: */
	double usec_busca_total = user_cpu_time_usec() ;	
	if(o->performanceTestOnly){
		for(y = 0 ;y < ny; y++){
			for(x = 0; x < nx; x++){
				int h = x;   /* Indice horizontal estilo PPM. */
        			int v = y;   /* Indice vertical estilo PPM. */	
        			bool_t process_pixel = ((h >= o->hMin) && (h <= o->hMax) && (v >= o->vMin) && (v <= o->vMax));
				if(imagem_mask != NULL){
					float val = float_image_get_sample(imagem_mask,0,x,y);
					if(val < 1.0) process_pixel = FALSE;
				}
				if((!process_pixel )) continue;
				double SO[o->nLights];
				int preto = 1;
				double albedo = 0.0;        /* Albedo calculado do ponto {p} da cena, ou 0.0. */
		        	double logPrSG_gab[ind_gab];
				for (i = 0; i < o->nLights; i++){
					SO[i] = float_image_get_sample(S[i], canal, x, y);
					if (SO[i] > 0.0) { preto = 0; }
				}
				for(ind_gab = 0; ind_gab < o->nGauges ; ind_gab++){
					if ((o->hashOption == BRUTA_APENAS) || (o->hashOption == BRUTA_E_HASH)){
						if (o->alphaOption == NSIG_OR_ALPHA) {
							linha_bruta_gab[ind_gab] = localiza_alternativa (
								tab[ind_gab], SO, 
								estLogPrSG, estAlbedo, o->sigma, o->omg0, o->omg1, 
								&(logPrSG_gab[ind_gab]), &albedo, 
								NULL
							);
						} else {
							linha_bruta_gab[ind_gab] = localiza_simples (
								tab[ind_gab], SO, 
								estLogPrSG, estAlbedo, o->sigma, o->omg0, o->omg1,
								&(logPrSG_gab[ind_gab]), &albedo, 
								NULL
							);
						}
						(stats_brut_gab[ind_gab])++;
					}
					if ((o->hashOption == HASH_APENAS) || (o->hashOption == BRUTA_E_HASH)){
						if (preto) {
							linha_hash_gab[ind_gab] = 0;
							(stats_hash_gab[ind_gab])++;
							logPrSG_gab[ind_gab] = -INF;
							albedo = 0;
						}
						else {
							linha_hash_gab[ind_gab] = localiza_normal_hash(bg[ind_gab], tab[ind_gab], SO, &(logPrSG_gab[ind_gab]), &albedo,NULL,NULL);
							(stats_hash_gab[ind_gab])++;
							stats_hash++;
						}
					}
					if( o->hashOption == KDTREE){
						linha_hash_gab[ind_gab] = localiza_normal_kdtree(kdTree[ind_gab], tab[ind_gab], SO, &(logPrSG_gab[ind_gab]), &albedo,NULL);
					}
					if(o->hashOption == ANALY_INVERSION){
					  double albe;
					  double logpro ;
					  r3_t normal_analy;
					  double G_vec[o->nLights];
					  analytic_inversion_compute_normal(SO,ai, &normal_analy, &albe, G_vec, &logpro);
					  logPrSG_gab[ind_gab] = logpro;
					  linha_hash_gab[ind_gab] = 0;
					  (stats_brut_gab[ind_gab])++;
					}
				}
			}
		}
		usec_busca_total = user_cpu_time_usec() - usec_busca_total ;
	}
	
	r3_t rnp;
    	for(y = 0 ;(y < ny) && (!o->performanceTestOnly); y++){
      		for(x = 0; x < nx; x++){
        		/* Decide se pixel deve ser processado: */
        		int h = x;            /* Indice horizontal estilo PPM. */
        		int v = y;   /* Indice vertical estilo PPM. */
        		bool_t process_pixel = ((h >= o->hMin) && (h <= o->hMax) && (v >= o->vMin) && (v <= o->vMax));
			if(imagem_mask != NULL){
					float val = float_image_get_sample(imagem_mask,0,x,y);
					if(val < 1.0) process_pixel = FALSE;
			}
			if((!process_pixel )) continue;
			/* Decide se pixel deve ser debugado: */
        		int id;
			bool_t debug = FALSE;
			for (id = 0; id < o->nDebug; id++) {
        			if ((h == o->debug[id].c[0]) && (v == o->debug[id].c[1])) { debug = TRUE; }
			}
			if (debug){ fprintf(stderr,"debugging pixel [%d %d] = PNM [%d %d]\n",x,y,h,v); }
			/* Extrai o vetor de observação {SO} deste pixel: */
			double SO[o->nLights];
			int preto = 1;
			for (i = 0; i < o->nLights; i++){
				SO[i] = float_image_get_sample(S[i], canal, x, y);
				if (SO[i] > 0.0) { preto = 0; }
			}
			double so[o->nLights];
			double Smag;
			extrai_assinatura(SO, so, &Smag, o->nLights);
			rnp = (r3_t){{ 0,0,0 }};
			if(nref != NULL){
				rnp.c[0] = float_image_get_sample(nref, 0, x,y);
				rnp.c[1] = float_image_get_sample(nref, 1, x,y);
				rnp.c[2] = float_image_get_sample(nref, 2, x,y);
			}
			double albedo = 0.0;        /* Albedo calculado do ponto {p} da cena, ou 0.0. */
			double albedo_gab[o->nGauges]; 
       			r3_t snp = (r3_t){{0,0,0}}; /* Normal estimada do ponto {p} da cena, ou (0,0,0). */
			r3_t snp_gab[o->nGauges];
			double logPrSG = -INF;
			double logPrSG_gab[ind_gab];
			
			for(ind_gab = 0; ind_gab < o->nGauges; ind_gab++){
				albedo_gab[ind_gab] = 0.0;
		        	snp_gab[ind_gab] = (r3_t){{0,0,0}};
       				logPrSG_gab[ind_gab] = -INF;      /* Maximo {log(Pr(S|G)} para {G} na tabela, ou {-oo} */
       				if ( process_pixel) {
					/* Extrai a normal de referência {rnp} para este pixel: */
					/* Cria arquivo para debug dos cálculo de nrmal do pixel: */
					char* debug_pixel_prefix = NULL;
					FILE *arq_debug_pixel = NULL;
					if (debug) {
						char *nome_arq_debug_pixel = NULL;
						char *debug_pixel_prefix = jsprintf("%s_%d_G%s_debug_%04d_%04d",o->prefix, canal,o->gaugeTag[ind_gab], h, v);
						//char *nome_arq_debug_pixel = jsprintf("%s_debug_%d_%04d_%04d_G%d.txt", o->prefix, canal, h, v,ind_gab);
						//char *nome_arq_debug_pixel = jsprintf("%s_debug_%d_%04d_%04d_G%s.txt", o->prefix, canal, h, v,o->gaugeTag);
						char *nome_arq_debug_pixel = jsprintf("%s_%d_G%s_debug_%04d_%04d.txt", o->prefix, canal,o->gaugeTag[ind_gab], h, v);
						arq_debug_pixel = fopen(nome_arq_debug_pixel, "wt");
						(*dbopt) = (normal_debug_opts_t) {
							.c = canal,
							.hp = h,
							.vp = v,
							.arq_debug = arq_debug_pixel,
							.prefixo = o->prefix,
							.mapa_gabarito = TRUE,
							.normal_ref = rnp,
							.gera_plots_q = TRUE
						};
          				}
					/*Ticks antes da busca*/
					/* Localiza {SO} na tabela, determinando a linha {linha} da mesma (ponto {q} do gabarito): */
					if ((o->hashOption == BRUTA_APENAS) || (o->hashOption == BRUTA_E_HASH)){
						if (o->alphaOption == NSIG_OR_ALPHA) {
							linha_bruta_gab[ind_gab] = localiza_alternativa (
							tab[ind_gab], SO, 
							estLogPrSG, estAlbedo, o->sigma, o->omg0, o->omg1, 
							&(logPrSG_gab[ind_gab]), &(albedo_gab[ind_gab]), 
							(debug ? dbopt : NULL)
							);
						} else {
							linha_bruta_gab[ind_gab] = localiza_simples (
							tab[ind_gab], SO, 
							estLogPrSG, estAlbedo, o->sigma, o->omg0, o->omg1,
							&(logPrSG_gab[ind_gab]), &(albedo_gab[ind_gab]), 
							(debug ? dbopt : NULL)
							);
						}
						(stats_brut_gab[ind_gab])++;
					}
					if ((o->hashOption == HASH_APENAS) || (o->hashOption == BRUTA_E_HASH)){
						if (preto) {
							linha_hash_gab[ind_gab] = 0;
							logPrSG_gab[ind_gab] = -INF;
							albedo_gab[ind_gab] = 0;
						}
						else {
							int n_euclid_evals,n_scans;
							linha_hash_gab[ind_gab] = localiza_normal_hash(bg[ind_gab], tab[ind_gab], SO, &(logPrSG_gab[ind_gab]), &(albedo_gab[ind_gab]),&n_euclid_evals,&n_scans);
							if(imagem_euclid_evals != NULL) float_image_set_sample(imagem_euclid_evals,0,x,y,n_euclid_evals);
							if(imagem_scans != NULL) float_image_set_sample(imagem_scans,0,x,y,n_scans);
							(stats_hash_gab[ind_gab])++;
						}
					}
					if(o->hashOption == KDTREE){
						double n_euclid_evals;
						linha_hash_gab[ind_gab] = localiza_normal_kdtree(kdTree[ind_gab], tab[ind_gab], SO, &(logPrSG_gab[ind_gab]), &(albedo_gab[ind_gab]),&n_euclid_evals);
						if(imagem_euclid_evals != NULL) float_image_set_sample(imagem_euclid_evals,0,x,y,n_euclid_evals);
						if(imagem_scans != NULL) float_image_set_sample(imagem_scans,0,x,y,1);
						(stats_hash_gab[ind_gab])++;
					}
					if(o->hashOption == RANSAC_DEBURRO){
					  double albe;
					  double logpro;
					  snp_gab[ind_gab] = compute_normal_by_RANSAC(SO,tab[ind_gab],o->sigma,o->RansacProbOutlier,o->gaugeAlbedo->c[canal],&albe,&logpro);
					  logPrSG_gab[ind_gab] = logpro;
					  albedo_gab[ind_gab] = albe;
					  linha_hash_gab[ind_gab] = 0;
					}
					if(o->hashOption == FAST_HASH){
					  double albe;
					  double logpro ;
					  snp_gab[ind_gab] = fast_hash_compute_normal(FH_set[ind_gab],SO, o->sigma,o->omg0,o->omg1,&albe,&logpro);
					  logPrSG_gab[ind_gab] = logpro;
					  albedo_gab[ind_gab] = albe;
					  linha_hash_gab[ind_gab] = 0;
					}
					if(o->hashOption == ANALY_INVERSION){
					  double albe;
					  double logpro ;
					  r3_t normal_analy;
					  double G_vec[o->nLights];
					  analytic_inversion_compute_normal(SO,ai, &normal_analy, &albe, G_vec, &logpro);
					  snp_gab[ind_gab] = normal_analy;
					  logPrSG_gab[ind_gab] = logpro;
					  albedo_gab[ind_gab] = albe;
					  linha_hash_gab[ind_gab] = 0;
					  (stats_brut_gab[ind_gab])++;
					}
					if(o->hashOption == ANALY_ALPHA){
					  double albe;
					  double logpro ;
					  r3_t normal_analy;
					  double G_vec[o->nLights];
					 
					  analytic_alpha_compute_normal(SO,o->nLights, ai,TRUE,&normal_analy,&albe,&logpro,G_vec,o->logProbFunction,o->sigma,o->omg0,o->omg1,o->nSteps,o->tolAlbedo,o->tolNormal,get_view_dir(tab[0]), debug_pixel_prefix);
					  snp_gab[ind_gab] = normal_analy;
					  logPrSG_gab[ind_gab] = logpro;
					  albedo_gab[ind_gab] = albe;
					  linha_hash_gab[ind_gab] = 0;
					}
					
					if(o->hashOption == FAST_RANSAC){
					  double albe;
					  double logpro;
					  int valids = 0;
					  snp_gab[ind_gab] = compute_normal_by_fast_RANSAC(SO,o->nLights,FH_ransac_set[ind_gab],o->RansacLightSeq, o->RansacFastNum,o->sigma,o->RansacProbOutlier,o->gaugeAlbedo->c[canal],&albe,&logpro,&valids);
					  logPrSG_gab[ind_gab] = logpro;
					  albedo_gab[ind_gab] = albe;
					  linha_hash_gab[ind_gab] = 0;
					  float_image_set_sample(fi_count_valids,0,x,y,valids);
					}
					if(o->hashOption == SUPER_TABLE){
						int ind_st;
						double best_prob= -MAXDOUBLE;
						double best_albedo;
						Tabela* best_tab;
						int best_line = -1;
						double test_distance;
						double test_albedo;
						int test_line;
						double SO_test_mag;
						int tableIndex;
						int tot_euclid_evals = 0;
						int tot_scans = 0;
						int best_st = 0;
						for(ind_st = 0; ind_st < o->numSubsets; ind_st++){
							int n_euclid_evals,n_scans;
							test_line = stSearchNormal(
    							ST_set[ind_gab][ind_st],
    							SO,
							&SO_test_mag, 
    							&test_distance,
    							&test_albedo,
							&tableIndex,
							&n_euclid_evals,
							&n_scans
							);
							tot_euclid_evals+=n_euclid_evals;
							tot_scans+=n_scans;
							//printf("TEST %f\n",test_distance);
							if(tableIndex != -1){
								const double* go = get_intdir(tab[ind_gab],tableIndex);
								double Gmag = get_intmag(tab[ind_gab],tableIndex);
								//double test_prob = EstLogPrSG_09(so, Smag,go,Gmag,o->nLights, o->sigma,o->omg0,o->omg1);
								double test_prob = estLogPrSG(so, Smag,go,Gmag,o->nLights, o->sigma,o->omg0,o->omg1);
								if(best_prob <  test_prob){
									best_prob = test_prob;
									best_albedo = test_albedo;
									best_line = tableIndex;
									best_tab = stGetTable(ST_set[ind_gab][ind_st]);
									best_st = ind_st;
								}
								if(debug){
									r3_t norm,norm2;
									norm = get_normal(tab[ind_gab], tableIndex);
									norm2 = get_normal(stGetTable(ST_set[ind_gab][ind_st]), test_line);
									fprintf(dbopt->arq_debug,"[%03d]: %04d (%6.3f ,%6.3f ,%6.3f ) ",ind_st,tableIndex,norm.c[0],norm.c[1],norm.c[2]);
									fprintf(dbopt->arq_debug," (%6.3f ,%6.3f ,%6.3f ) ",norm2.c[0],norm2.c[1],norm2.c[2]);
									int inder;
									const int* ind_luz = stGetIndLuz(ST_set[ind_gab][ind_st]);
									fprintf(dbopt->arq_debug,"SO { ");
									for(inder = 0; inder < o->subsetSize; inder++){
										fprintf(dbopt->arq_debug,"%f ",SO[ind_luz[inder]]);
									}
									fprintf(dbopt->arq_debug," }  ");
									fprintf(dbopt->arq_debug,"GO { ");
									for(inder = 0; inder < o->subsetSize; inder++){
										fprintf(dbopt->arq_debug,"%f ",go[ind_luz[inder]]*Gmag);
									}
									fprintf(dbopt->arq_debug," } Dist: %6.4f\n\n",test_prob);
								}
							}
						}
						float_image_set_sample(imagem_subset,0,x,y,best_st);
						if(imagem_euclid_evals != NULL) float_image_set_sample(imagem_euclid_evals,0,x,y,tot_euclid_evals);
						if(imagem_scans != NULL) float_image_set_sample(imagem_scans,0,x,y,tot_scans);
						albedo_gab[ind_gab] = best_albedo;
						logPrSG_gab[ind_gab] = best_prob;
						linha_hash_gab[ind_gab] = best_line;
  					}
					if ((o->hashOption == BRUTA_E_HASH) && (linha_hash != linha_bruta)){
    	    					fprintf(stderr, "BUSCA HASH DIFERENTE DE BUSCA BRUTA\n");
    	    					fprintf(stderr, "HASH = "); print_linha(tab[ind_gab],linha_hash_gab[ind_gab]);
    	    					fprintf(stderr, "BRUT = "); print_linha(tab[ind_gab],linha_bruta_gab[ind_gab]);
    	  				}
    	  				linha_gab[ind_gab] = ((o->hashOption == BRUTA_APENAS) || (o->hashOption == BRUTA_E_HASH) ? linha_bruta_gab[ind_gab] : linha_hash_gab[ind_gab]);
					indice_normais_gab[ind_gab][x + y*nx] = linha_gab[ind_gab];
					/* Obtém o vetor de observações {GO} do ponto {q} do gabarito: */
					double GO[o->nLights]; /* Vetor de observação do ponto {q} no gabarito, ou (0..). */
					const double *go;      /* Assinatura normalizada de {q}, ou (0..). */
					double Gmag;           /* Magnitude de {GO}, ou 0. */
					if(linha >= 0){
						go = get_intdir(tab[ind_gab], linha_gab[ind_gab]);  /* Assinatura do ponto {q} no gabarito. */
						Gmag = get_intmag(tab[ind_gab], linha_gab[ind_gab]);  /* Magnitude do vetor de observações do ponto {q}. */
						for (i = 0; i < o->nLights; i++){ GO[i] = go[i] * Gmag; }
   	  				}else{
            					for (i = 0; i < o->nLights; i++){ GO[i] = 0.0; }
						go = GO; /* Dirty trick. */
						Gmag = 0.0;
					}
					if((o->hashOption != RANSAC_DEBURRO) && (o->hashOption != ANALY_ALPHA) && (o->hashOption != ANALY_INVERSION) && (o->hashOption != FAST_HASH) && (o->hashOption != FAST_RANSAC) ){
					  /* Considera que a normal do ponto é a normal do gabarito: */
					  snp_gab[ind_gab] = get_normal(tab[ind_gab],linha_gab[ind_gab]); /* Normal do ponto {q} do gabarito. */
					}
       	  				/* Dados para estudar os alphas: */
					/* Estima log verossimilhança {log(Pr(S|G)} para a linha escolhida: */
					if (arq_debug_pixel != NULL){ fclose(arq_debug_pixel);}
				
				  if(debug_pixel_prefix != NULL) free(debug_pixel_prefix);
       				}/*fim do processamento do pixel*/
			}
			/* aqui entra a decisão dos dados que serão passados para o arq se saída principal*/
			int selected_gab = 0;
			if( (o->nGauges > 1) ){
				//desempatar com menor distância de log_estprob
				double best_logPrSG;
				int is_selected = 0;
				for(ind_gab = 0; ind_gab < o->nGauges; ind_gab++){
					double test_prob;
					if((o->hashOption == HASH_APENAS) || (o->hashOption == BRUTA_E_HASH) ){
						// we need compute dist_alpha
						const double* go = get_intdir(tab[ind_gab],linha_gab[ind_gab]);
						double Gmag = get_intmag(tab[ind_gab],linha_gab[ind_gab]);
						test_prob = estLogPrSG(so, Smag,go,Gmag,o->nLights, o->sigma,o->omg0,o->omg1);
					}else{
						test_prob = logPrSG_gab[ind_gab]; 
					}
					if(is_selected){
						if(best_logPrSG < test_prob){
							best_logPrSG = test_prob;
							selected_gab = ind_gab;
						}
					}else{
						is_selected = 1;
						best_logPrSG = test_prob;
						selected_gab = ind_gab;
					}
				}
			}else{
				if((o->hashOption == HASH_APENAS) || (o->hashOption == BRUTA_E_HASH) ){
					const double* go = get_intdir(tab[0],linha_hash_gab[0]);
					double Gmag = get_intmag(tab[0],linha_hash_gab[0]);
					logPrSG_gab[0] = estLogPrSG(so, Smag,go,Gmag,o->nLights, o->sigma,o->omg0,o->omg1);	
				}	
			}
			
			/*gather gauge data inside the correspondent structures*/
			double gsamp = ((double)selected_gab)/((double)o->nGauges);
			float_image_set_sample(imagem_gab_select,0,x,y,gsamp );
			fprintf(arq_gab_select_channel,"%d %d %d\n",x,y,selected_gab);
			for(ind_gab = 0; ind_gab < o->nGauges; ind_gab++){
				float_image_set_sample(fi_normais_canal_gab[ind_gab],0,x,y,snp_gab[ind_gab].c[0]);
				float_image_set_sample(fi_normais_canal_gab[ind_gab],1,x,y,snp_gab[ind_gab].c[1]);
				float_image_set_sample(fi_normais_canal_gab[ind_gab],2,x,y,snp_gab[ind_gab].c[2]);	
			}
			/*pic the final data*/
			logPrSG = logPrSG_gab[selected_gab];
			snp = snp_gab[selected_gab];
			albedo = albedo_gab[selected_gab];
			linha_bruta = linha_bruta_gab[selected_gab];
			linha_hash = linha_hash_gab[selected_gab];
			linha = linha_gab[selected_gab];
			/* Armazena o albedo em {imagem_albedo}: */
			double albedo_dimming = 0.90; /* Fator de redução para caso do albedo ser maior que 1.0. */
			float_image_set_sample(imagem_albedo, 0, x, y, albedo_dimming*albedo);
			/* Grava normal no mapa de normais: */
			float_image_set_sample(fi_normais_canal,0,x,y,snp.c[0]);
			float_image_set_sample(fi_normais_canal,1,x,y,snp.c[1]);
			float_image_set_sample(fi_normais_canal,2,x,y,snp.c[2]);
			//fprintf(arq_normais_canal,"%d %d %f %f %f\n", x, y, snp.c[0], snp.c[1], snp.c[2]);
			if(nref != NULL){
				r3_t enp ;
				r3_sub(&snp,&rnp,&enp);
				double err_val = r3_norm(&enp);
				//fprintf(arq_err_normais_canal,"%d %d %f %f %f\n", x, y, enp.c[0], enp.c[1], enp.c[2]);
				float_image_set_sample(fi_err_normais_canal,0,x,y,err_val);
	
			}
			/* Salva normal e log verossimilhança para imagens médias de canais: */
			int ip = x + y*nx;
       			normais_calculadas[ip] = snp;
        		logPrSG_busca[ip] = logPrSG;
			/*perfumaria - A.K.A - contador de progresso*/
			if (process_pixel) {
       				/* Imprime relatório de progresso: */
  				contador++;
				if (PROGRESS_METER) {
					if ((contador % 10) == 0) {
						fprintf(stderr,"\033[1A");
						double normals_per_sec = contador/(double)(time(NULL) - last_tempo  + 0.001);
						double total_secs_remaining = (total_compute_pixels - contador)/normals_per_sec;
						int total_seconds = (int)floor(total_secs_remaining);
						int hour = total_seconds/(60*60);
						int min = (total_seconds - (hour*60*60))/60;
						int sec = (total_seconds - (hour*60*60) - (min*60));
						fprintf(stderr,"[%d][%9d] of [%9d] - %6.6f%% - %6.6f n/s   - %02d h %02d m %02d s       \n",
						canal,contador,total_compute_pixels, contador*100.0/total_compute_pixels,normals_per_sec,hour,min,sec
						);
					}
				} else {
					int step = nx*ny/10;
					if ((contador % step == 0) || (contador == nx*ny)) {
						fprintf(stderr, "canal %d, %6d pixels - %5.5lf%%\n", canal, contador, 100*((double)contador)/(nx*ny));
					}
				}
			}
		}/*fim eixo X*/
    	}/*fim eixo Y*/



	if(! o->performanceTestOnly ) usec_busca_total = user_cpu_time_usec() - usec_busca_total ;
    	/* Fecha arquivo de normais do canal: */
	//fprintf(arq_normais_canal,"\n");
	float_image_write(arq_normais_canal,fi_normais_canal);
	for(ind_gab = 0; ind_gab < o->nGauges; ind_gab++){
		char *nome_arq_normais_canal_temp= NULL;
    		//char *nome_arq_normais_canal_temp = jsprintf("%s_%d_G%d_normals.fni", o->prefix, canal,ind_gab);
		char *nome_arq_normais_canal_temp = jsprintf("%s_%d_G%s_normals.fni", o->prefix, canal,o->gaugeTag[ind_gab]);
		FILE* arq_normais_temp = fopen(nome_arq_normais_canal_temp,"wt");	
		float_image_write(arq_normais_temp,fi_normais_canal_gab[ind_gab]);
		fclose(arq_normais_temp);
	}
	fclose(arq_normais_canal);
	if(nref != NULL){
		double min,max,mean,rms;
		findMinMaxMean(fi_err_normais_canal,&min, &max,&mean,&rms);
		
		double radmin = 2.0*asin(min/2.0), degmin = (180.0/M_PI)*radmin;
		double radmax = 2.0*asin(max/2.0), degmax = (180.0/M_PI)*radmax;
		double radmean = 2.0*asin(mean/2.0), degmean = (180.0/M_PI)*radmean;
		double radrms = 2.0*asin(rms/2.0), degrms = (180.0/M_PI)*radrms;
		float_image_write(arq_err_normais_canal,fi_err_normais_canal);
		float_image_free(fi_err_normais_canal);
		fclose(arq_err_normais_canal);
		char* nome_arq_norm_err_temp;
		char *nome_arq_norm_err_temp = jsprintf("%s_%d_G%s_norm_error_stats.txt", o->prefix, canal,o->gaugeTag[ind_gab]);
		FILE* arq_norm_err_temp = fopen(nome_arq_norm_err_temp,"wt");
		fprintf(arq_norm_err_temp,"#Unit MIN MAX MEAN RMS\n");
		fprintf(arq_norm_err_temp,"Abs %lf %lf %lf %lf\n",min,max,mean,rms);
		fprintf(arq_norm_err_temp,"Rad %lf %lf %lf %lf\n",radmin,radmax,radmean,radrms);
		fprintf(arq_norm_err_temp,"Deg %lf %lf %lf %lf\n",degmin,degmax,degmean,degrms);
		fclose(arq_norm_err_temp);
		free(nome_arq_norm_err_temp);
		
	}
	fclose(arq_gab_select_channel);
	/* Escreve estatísticas da busca para este canal: */
	stats_brut = 0;
	stats_hash = 0;
	for(ind_gab = 0; ind_gab < o->nGauges; ind_gab++){
		stats_brut+= stats_brut_gab[ind_gab];
		stats_hash+= stats_hash_gab[ind_gab];
	}
	stats_brut= stats_brut/(float)o->nGauges;
	stats_hash= stats_hash/(float)o->nGauges;
	double fpercent = (stats_brut*100.0/total_compute_pixels)/3.0;
    	fprintf(stderr, "Normais calculadas por forca bruta - %d ( %3.2lf%% ) \n",stats_brut,fpercent);
    	fpercent = (stats_hash*100.0/total_compute_pixels)/3.0;
    	fprintf(stderr, "Normais calculadas por hash - %d ( %3.2lf%% ) \n",stats_hash,fpercent);
	stats_euclid = 0;
    	stats_scan = 0;
    	int j;
	/* Anexa estatísticas de buckets no arquivo: */
	for(ind_gab = 0; ind_gab < o->nGauges; ind_gab++){
    		char *nome_arq_estat = NULL;
    		//char *nome_arq_estat = jsprintf("%s_G%d_times_%02d.txt", o->prefix, ind_gab,canal);
		//char *nome_arq_estat = jsprintf("%s_G%s_times_%02d.txt", o->prefix, o->gaugeTag,canal);
		if(bg[ind_gab] != NULL){
			char *nome_arq_estat = jsprintf("%s_%d_G%s_times.txt", o->prefix,canal, o->gaugeTag[ind_gab]);
    			FILE* arq_estat = fopen(nome_arq_estat,"wt");
      			for(i = 0; i < get_tam_grid(bg[ind_gab]); i++){
    				for(j = 0; j < tam_grid[ind_gab]; j++){
    					stats_euclid_gab[ind_gab] += acessaMatriz_Statistic_Euclid(bg[ind_gab])[i][j];
    					stats_scan_gab[ind_gab] += acessaMatriz_Statistic_Scan(bg[ind_gab])[i][j];
    				}
      			}
      			double buckets =  tam_grid[ind_gab]*tam_grid[ind_gab];
      			double entries_per_bucket = ((double)num_linhas[ind_gab])/buckets;
      			int queries = stats_brut_gab[ind_gab] + stats_hash_gab[ind_gab];
      			double usec_per_query = usec_busca_total/(queries);
      			double euclids_per_query = ((double) stats_euclid_gab[ind_gab])/(queries);
      			double bucket_scans_per_query = ((double) stats_scan_gab[ind_gab])/(queries);
      			fprintf(arq_estat,"#Queries N T entries_per_bucket t_tot t_avg d_med b_med \n");
      			fprintf(arq_estat," %d %d %d", queries, tam_grid[ind_gab], num_linhas[ind_gab]);
      			fprintf(arq_estat," %lf %lf %lf", entries_per_bucket, usec_busca_total, usec_per_query);
      			fprintf(arq_estat,"  %lf %lf\n", euclids_per_query, bucket_scans_per_query);
			fclose(arq_estat);
    			free(nome_arq_estat);
		}
    		if(o->hashOption == BRUTA_APENAS){
			char *nome_arq_estat = jsprintf("%s_%d_G%s_times.txt", o->prefix,canal, o->gaugeTag[ind_gab]);
    			FILE* arq_estat = fopen(nome_arq_estat,"wt");
			fprintf(arq_estat,"#Queries N T entries_per_bucket t_tot t_avg d_med b_med \n");
			int queries = stats_brut_gab[ind_gab] + stats_hash_gab[ind_gab];
			double usec_per_query = usec_busca_total/(queries);
			double euclids_per_query = num_linhas[ind_gab];
      			double bucket_scans_per_query = 1.0;
      			fprintf(arq_estat," %d %d %d", queries, 0, num_linhas[ind_gab]);
      			fprintf(arq_estat," %lf %lf %lf", euclids_per_query, usec_busca_total, usec_per_query);
      			fprintf(arq_estat,"  %lf %lf\n", euclids_per_query, bucket_scans_per_query);
			fclose(arq_estat);
    			free(nome_arq_estat);
		}
		/* Escreve estatísticas da tabela de hash deste canal: */
    		if((o->hashOption != BRUTA_APENAS) && (o->hashOption != SUPER_TABLE) && (o->hashOption != KDTREE) ){
      			char *nome_arq_grid_show = NULL;
      			//char *nome_arq_grid_show = jsprintf("%s_G%d_statistics_%d_", o->prefix,ind_gab, canal);
			char *nome_arq_grid_show = jsprintf("%s_%d_G%s_", o->prefix,canal,o->gaugeTag[ind_gab]);
		       	if(o->showBucketsEPS)   showBucketsEPS(nome_arq_grid_show,bg[ind_gab],get_num_linhas(tab[ind_gab]));
      			if(o->showBucketsPPM)   showBucketsPPM(nome_arq_grid_show,bg[ind_gab],get_num_linhas(tab[ind_gab]));
      			if(o->showBucketsPLT)   showBucketsPLT(nome_arq_grid_show,bg[ind_gab],get_num_linhas(tab[ind_gab]));
      			if(o->showBucketData)   showBucketsData(nome_arq_grid_show,bg[ind_gab],get_num_linhas(tab[ind_gab]));
			//showBucketGridStatistic(nome_arq_grid_show,bg[ind_gab]);
			char *nome_arq_grid_stats = NULL;
      			//char *nome_arq_grid_stats = jsprintf("%s_G%d_bucket_grid_stats_%d.txt", o->prefix,ind_gab, canal);
			char *nome_arq_grid_stats = jsprintf("%s_%d_G%s_bucket_grid_stats.txt", o->prefix,canal,o->gaugeTag[ind_gab]);
      			FILE* arq_grid_stats = fopen(nome_arq_grid_stats, "wt");
      			flushBucketGrid(arq_grid_stats, bg[ind_gab]);
      			fclose(arq_grid_stats); free(nome_arq_grid_stats);
			LiberaBucketGrid(bg[ind_gab]);
    		}else if((o->hashOption == SUPER_TABLE) ){
			char *nome_arq_grid_show = NULL;
			int ind;
			for(ind =0 ; ind < o->numSubsets; ind++){
				//char *nome_arq_grid_show = jsprintf("%s_G%d_statistics_S%d_%d_", o->prefix,ind_gab,ind, canal);
				char *nome_arq_grid_show = jsprintf("%s_%d_G%s_S%d", o->prefix, canal,o->gaugeTag[ind_gab],ind);
				bucketGrid* buck = stGetBucketGrid(ST_set[ind_gab][ind]);
				if(o->showBucketsEPS)   showBucketsEPS(nome_arq_grid_show,buck,get_num_linhas(tab[ind_gab]));
      				if(o->showBucketsPPM)   showBucketsPPM(nome_arq_grid_show,buck,get_num_linhas(tab[ind_gab]));
      				if(o->showBucketsPLT)   showBucketsPLT(nome_arq_grid_show,buck,get_num_linhas(tab[ind_gab]));
				if(o->showBucketData)   showBucketsData(nome_arq_grid_show,buck,get_num_linhas(tab[ind_gab]));
				fprintf(stderr,"liberado ST %d G%d ",ind, ind_gab);
				LiberaSuperTabela(ST_set[ind_gab][ind]);
			}
			char *nome_arq_subset = NULL;
			char *nome_arq_subset = jsprintf("%s_%d_G%s_subset.fni",o->prefix,canal,o->gaugeTag[ind_gab]);
			FILE* arq_subset = open_write(nome_arq_subset,TRUE);
			float_image_write(arq_subset,imagem_subset);
			fclose(arq_subset);
		}else{
			fprintf(stderr, "TOT DIST EUCLIDIANAS %d",total_compute_pixels*num_linhas[ind_gab]);
		}


		/*Escreve estatísticas da KdTree*/
		if((o->hashOption == KDTREE)){
			char *nome_arq_tree_stats = NULL;
			kdTree[ind_gab]->stats.usec_time = usec_busca_total;
      			//char *nome_arq_tree_stats = jsprintf("%s_G%d_tree_stats_%d.txt", o->prefix,ind_gab, canal);
			char *nome_arq_tree_stats = jsprintf("%s_%d_G%s_tree_stats.txt", o->prefix,canal,o->gaugeTag[ind_gab]);
			FILE* arq_tree_stats = fopen(nome_arq_tree_stats, "wt");
			printfTreeStats(arq_tree_stats,kdTree[ind_gab]);
			fclose(arq_tree_stats);
		}
		/*Calcula distancia do go e so*/
		if(o->computeDist){
			char* nome_arq_dist;
			char *nome_arq_dist = jsprintf("%s_%d_G%s_SGdist.fni",o->prefix,canal,o->gaugeTag[ind_gab]);
			float_image_t* fi_dist = float_image_new(1,nx,ny);
			for(x = 0; x < nx; x++){
				for(y = 0; y < ny; y++){
					int ip = x + y*nx;
					int linha = indice_normais_gab[ind_gab][ip];
					const double* go_or = get_intdir(tab[ind_gab],linha);
					double go[o->nLights];
					double SO[o->nLights];
					double so[o->nLights];
					int i_light;
					for(i_light = 0; i_light < o->nLights; i_light++){
							SO[i_light] = float_image_get_sample(S[i_light], canal, x, y);	
							go[i_light] = go_or[i_light];
					}
					double Smag;
					extrai_assinatura(SO, so,&Smag,o->nLights);	
					double dist = rn_dist (o->nLights, so, go);
					float_image_set_sample(fi_dist,0,x,y,dist);
				
				}
			}
			FILE* arq_dist = open_write(nome_arq_dist,FALSE);
			float_image_write(arq_dist,fi_dist);
			fclose(arq_dist);
			float_image_free(fi_dist);
			free(nome_arq_dist);
		}
		/* Calcula imagems LI e KI*/
		if( o->computeKiLi ){
			for(i = 0; i < o->nLights;i++){
				char *nome_arq_Li = NULL;
				//char *nome_arq_Li = jsprintf("%s_%d_G%d_L%02d.pgm",o->prefix,canal,ind_gab,i);
				char *nome_arq_Li = jsprintf("%s_%d_G%s_L%02d.pgm",o->prefix,canal,o->gaugeTag[ind_gab],i);
				char *nome_arq_Ki = NULL;
				//char *nome_arq_Ki = jsprintf("%s_%d_G%d_K%02d.pgm",o->prefix,canal,ind_gab,i);
				char *nome_arq_Ki = jsprintf("%s_%d_G%s_K%02d.fni",o->prefix,canal,o->gaugeTag[ind_gab],i);
				fprintf(stderr,"Gerando imagens %s %s \n",nome_arq_Li,nome_arq_Ki);
				float_image_t* fi_Li = float_image_new(1,nx,ny);
				float_image_t* fi_Ki = float_image_new(1,nx,ny);
				for(x = 0; x < nx; x++){
					for(y = 0; y < ny; y++){
						double *go;      /* Assinatura normalizada de {q}, ou (0..). */
						double Gmag;
						int ip = x + y*nx;
						if(o->hashOption != ANALY_INVERSION && o->hashOption != ANALY_ALPHA){
						  int linha = indice_normais_gab[ind_gab][ip];
						  go = get_intdir(tab[ind_gab],linha);
						  Gmag = get_intmag(tab[ind_gab],linha);
						}else{
						  assert(ai != NULL);
						  double G[o->nLights];
						  double dx = float_image_get_sample(fi_normais_canal_gab[ind_gab],0,x,y);
						  double dy = float_image_get_sample(fi_normais_canal_gab[ind_gab],1,x,y);
						  double dz = float_image_get_sample(fi_normais_canal_gab[ind_gab],2,x,y);
						  r3_t normal = (r3_t){{dx,dy,dz}};
						  go = (double*)malloc(sizeof(double)*(o->nLights));
						  analytic_inversion_compute_gauge_OV(ai, normal, G);
						  Gmag = rn_dir(o->nLights,G,go);
						  
						}
						double SO[o->nLights];
						double so[o->nLights];
						double Smag;
						int i_light;
						for(i_light = 0; i_light < o->nLights; i_light++){
							SO[i_light] = float_image_get_sample(S[i_light], canal, x, y);	
						}
						extrai_assinatura(SO, so,&Smag,o->nLights);
						double albedo  = estAlbedo(so,Smag,go,Gmag,o->nLights,o->sigma,o->omg0,o->omg1);
						
						double Li = go[i]*Gmag;
// 						double Ki = ((albedo*Li <= 0) ? 0 : SO[i]/(albedo*Li));
						double Ki = SO[i] - (albedo*Li);
						float_image_set_sample(fi_Li,0,x,y,Li);
						float_image_set_sample(fi_Ki,0,x,y,Ki);
						if(o->hashOption == ANALY_INVERSION || o->hashOption == ANALY_ALPHA){
						  free(go);
						}
					}
				}
				
				//FILE* arq_Li = fopen(nome_arq_Li,"wt");
				float_pnm_image_write(nome_arq_Li, fi_Li,FALSE, VIEW_GAMMA, VIEW_BIAS,TRUE,TRUE,FALSE);
				//fclose(arq_Li);
				float_image_free(fi_Li);
				free(nome_arq_Li);
				FILE* arq_Ki = fopen(nome_arq_Ki,"wt");
				float_image_write(arq_Ki, fi_Ki);
				fclose(arq_Ki);
				float_image_free(fi_Ki);
				free(nome_arq_Ki);

			}
		}
		/*free things that wont be usefull anymore*/
		float_image_free(fi_normais_canal_gab[ind_gab]);	
		free(indice_normais_gab[ind_gab]);
		
	}
   
	/*Grava imagem do logPrSG */
	char *nome_arq_logPrSg = NULL;
	//char *nome_arq_logPrSg = jsprintf("%s_%d_logPrSG.fni",o->prefix,canal);
	char *nome_arq_logPrSg = jsprintf("%s_%d_G%s_logPrSG.fni",o->prefix,canal,o->gaugeTag[0]);
	float_image_t* fi_logPrSG = float_image_new(1,nx,ny);
	float Lmax = -1;
	for(x = 0; x < nx; x++){
		for(y = 0; y < ny; y++){
			int ip = x + y*nx;
       			float prob = logPrSG_busca[ip];
			float_image_set_sample(fi_logPrSG,0,x,y,prob);
			if((prob != +INF) && (prob != -INF) && (!isnan(prob) ) ){
				if(prob > Lmax) Lmax = prob;	
			}
		}
	}
	FILE* arq_logPrSG = fopen(nome_arq_logPrSg,"wt");
	float_image_write(arq_logPrSG,fi_logPrSG);
	fclose(arq_logPrSG);
	float_image_free(fi_logPrSG);

	char *nome_arq_pesos = NULL;
	//char *nome_arq_pesos = jsprintf("%s_%d_weights.fni",o->prefix,canal);
	char *nome_arq_pesos = jsprintf("%s_%d_G%s_weights.fni",o->prefix,canal,o->gaugeTag[0]);		
	float_image_t* fi_pesos = float_image_new(1,nx,ny);
	float TINY_LOGPROB = log(1.0e-10);
	float TINY_DOT = 0.01;


	for(x = 0; x < nx; x++){
		for(y = 0; y < ny; y++){
			int ip = x + y*nx;
			float prob = logPrSG_busca[ip];
			float weight;
			if(prob == +INF) weight = 1.0;
			else if( prob == -INF) weight = 0.0;
			else if( isnan(prob) ) weight = 0.0;
			else if( (prob - Lmax) < TINY_LOGPROB ) weight = 0.0;
			else weight = exp(prob - Lmax);
			r3_t norm = (r3_t){{float_image_get_sample(fi_normais_canal,0,x,y), float_image_get_sample(fi_normais_canal,1,x,y), float_image_get_sample(fi_normais_canal,2,x,y)}};
			r3_t view_dir  = get_view_dir(tab[0]);
			double dot = r3_dot(&norm,&view_dir);
			if(dot < TINY_DOT){
				 weight = 0.0 ;
// 				fprintf(stderr,"Pixel [%04d,%04d] Dot = %f View_Dir = ",x,y,dot);
// 				r3_print(stderr,&view_dir);
// 				fprintf(stderr," Norm = ");
// 				r3_print(stderr,&norm);
// 				fprintf(stderr,"\n");
			}
			float_image_set_sample(fi_pesos,0,x,y,weight);
		}
	}
	FILE* arq_pesos = fopen(nome_arq_pesos,"wt");
	float_image_write(arq_pesos,fi_pesos);
	fclose(arq_pesos);
	float_image_free(fi_pesos);
	float_image_free(imagem_mask);

	
	
	if(imagem_euclid_evals != NULL){
		char *nome_arq_euclid_evals = NULL;
		char *nome_arq_euclid_evals = jsprintf("%s_%d_G%s_euclid_evals.fni",o->prefix,canal,o->gaugeTag[0]);	
		FILE* arq_euclid_evals = open_write(nome_arq_euclid_evals,TRUE);
		float_image_write(arq_euclid_evals, imagem_euclid_evals);
		fclose(arq_euclid_evals);
		free(nome_arq_euclid_evals);
		float_image_free(imagem_euclid_evals);
	}
	
	if(fi_count_valids != NULL){
	  char *nome_arq_count_valids = NULL;
	  char *nome_arq_count_valids = jsprintf("%s_%d_G%s_count_valids.fni",o->prefix,canal,o->gaugeTag[0]);
	  FILE* arq_count_valids = open_write(nome_arq_count_valids,TRUE);
	  float_image_write(arq_count_valids, fi_count_valids);
	  fclose(arq_count_valids);
	  free(nome_arq_count_valids);
	  float_image_free(fi_count_valids);
	}

	if(imagem_scans != NULL){
		char *nome_arq_scans = NULL;
		char *nome_arq_scans = jsprintf("%s_%d_G%s_scans.fni",o->prefix,canal,o->gaugeTag[0]);	
		FILE* arq_scans = open_write(nome_arq_scans,TRUE);
		float_image_write(arq_scans, imagem_scans);
		fclose(arq_scans);
		free(nome_arq_scans);
		float_image_free(imagem_scans);
	}

	for (ind_gab = 0; ind_gab < o->nGauges; ind_gab++){
		fprintf(stderr,"\nLiberando Tabela G%d\n",ind_gab);
		LiberaTabela(tab[ind_gab]);
	}
	/* Fecha os arquivos de médias: */
	/* Grava e libera as imagens resultado: */
	int ind;
	if(!o->performanceTestOnly){
		float_pnm_image_write(nome_im_gab_select, imagem_gab_select,FALSE, 1.000, 0.000, TRUE,TRUE,TRUE);
		for (ind = 0; ind < o->nLights; ind++){ float_image_free(S[ind]); }
		char *nome_im_albedo = NULL;
		//char *nome_im_albedo = jsprintf("%s_%d_albedo.pgm", o->prefix,canal);
		char *nome_im_albedo = jsprintf("%s_%d_G%s_albedo.pgm", o->prefix,canal,o->gaugeTag[0]);
		float_pnm_image_write(nome_im_albedo, imagem_albedo,FALSE, VIEW_GAMMA, VIEW_BIAS,TRUE,TRUE,FALSE);
		char *nome_im_albedo = jsprintf("%s_%d_G%s_albedo.fni", o->prefix,canal,o->gaugeTag[0]);
		FILE* arq_albedo = open_write(nome_im_albedo,TRUE);
		float_image_write(arq_albedo,imagem_albedo);
	}
	fprintf(stderr, "Concluido!\nO programa rodou com sucesso!\n");

	time_t tempo_fim = time(NULL);
	timeinfo = localtime(&tempo_fim);
	fprintf(stderr, "FIM DO PROCESSAMENTO: %s",asctime(timeinfo));
	
	double usec_total = difftime(tempo_fim,tempo_start);
	fprintf(stderr, "Executado em %6.3f segs.\n",usec_total/1000000);
	
	return 0;
}
