#define PROG_NAME "compute_normals"

#define _GNU_SOURCE
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <i2.h>
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
#include <time.h>
#include <jsfile.h>
#include <string.h>
#include <assert.h>
#include <limits.h>
#include <sys/times.h>
#include <unistd.h>
#include <sys/time.h>
#include <sys/resource.h>
#include <lighting_models.h>



#define PROG_HELP \
  PROG_NAME " \\\n" \
  "    -nLights {NUM} \\\n" \
  "    -prefix {FILE_PREFIX} \\\n" \
  "    [ -channels {R | G | B| RG | RB | GB | RGB } \\\n" \
  "    -gamma {NUM} [ -gray ] \\\n" \
  "    -nGauges {NUM} \\\n" \
  "    [ -gaugePos {NUM} {NUM} ... ] \\\n" \
  "    -gaugeCenter {NUM} {NUM}... \\\n" \
  "    -gaugeRadius {NUM}... \\\n" \
  "    [ -gaugeStretch {STRX} {STRY}... ] \\\n" \
  "    [ -gaugeViewDir {VDIRX} {VDIRY} {VDIRZ}... ] \\\n" \
  "    [ -gaugeAlbedo {NUM} {NUM} {NUM}... ] \\\n" \
  "    -gaugeImages {FILENAME}...  \\\n" \
  "    [ -tableSize {NUM} ] \\\n" \
  "    [ -yflip  ] \\\n" \
  "    [ -debug {NP} {H[1]} {V[1]} ... {H[NP]} {V[NP]} ] \\\n" \
  "    [-saveSamplePoints] \\\n" \
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

typedef struct gaugeData_t{
	r2_t gaugeCenter;
	double gaugeRadius;
	r2_t gaugeStretch; /*extra stretching of gauge image*/
	r3_t gaugeViewDir;
	frgb_t gaugeAlbedo;
	
} gauge_data_sec_t;

typedef struct options_t {
	/* General parameters: */
	int nLights;
	char *prefix;
	char *channels;              /* Channels to process (a substring of "RGB"). */
	/* Input image parameters: */
	double gamma;
	double bias;
	bool_t gray;
	/* Observation matching parameters: */
	/* Bucket Grid only*/
	/* Gauge parameters: */
	gauge_data_sec_t* gauge_atributes;
	char** gauge_image_name; /* Gauge image filenames are {gauge_image_name[0..nLights-1]} */
	char** gauge_direction_name; /*files containing directions of light source - if not NULL will be used to generate virtual gauge table*/
	r3_t* gaugeDirections; /*contains light source directions of virtual gauge*/
	
	model_type_t lightingModel;
	char** gauge_lparms_name; /*files containing parameters of lihting models*/
	void*** lightingParams;
	double thetaMax;
	/* Light table parameters: */
	int tableSize;
	/* Debugging */
	int nGauges;  /* will be  used now. */
	r2_t* gaugePos;
	bool_t yflip;
	bool_t rawPixels; //does not interpolate pixels, uses pixels centers as sample points
	char** gaugeTags;
	
	bool_t saveSamplePoints;
	
} options_t;

options_t *parse_options(int argc, char **argv);

options_t *parse_options(int argc, char **argv){
	argparser_t *pp = argparser_new(stderr, argc, argv);
	argparser_set_help(pp, PROG_HELP);
	argparser_set_info(pp, PROG_INFO);
	argparser_process_help_info_options(pp);
	
	options_t *o = (options_t *)malloc(sizeof(options_t));
	
	if (argparser_keyword_present(pp, "-channels")) {
		o->channels = argparser_get_next(pp);
	} else {
		o->channels = "RGB";
	}
	
	fprintf(stderr, "  -channels %s \\\n", o->channels);
	
	
	if (argparser_keyword_present(pp, "-tableSize")) {
		o->tableSize = argparser_get_next_int(pp, 2, 1000);
	} else {
		o->tableSize = 30;
	}
	fprintf(stderr, "  -tableSize %d \\\n", o->tableSize);
	
	if (argparser_keyword_present(pp, "-nGauges")) {
		o->nGauges = argparser_get_next_int(pp, 1, 100);
	} else {
		o->nGauges = 1;
	}	
	
	int ind_gab;
	o->gaugeTags = (char**)malloc(sizeof(char*)*o->nGauges);
	if (argparser_keyword_present(pp, "-gaugeTags")) {
	  for(ind_gab = 0; ind_gab < o->nGauges; ind_gab++){
	    o->gaugeTags[ind_gab] = argparser_get_next(pp);
	  }
	}else{
	  for(ind_gab = 0; ind_gab < o->nGauges; ind_gab++){
	    o->gaugeTags[ind_gab] = NULL;
	    char *o->gaugeTags[ind_gab] = jsprintf("%d",ind_gab);
	  }
	}
	o->gauge_atributes = (gauge_data_sec_t*)malloc(sizeof(gauge_data_sec_t)*o->nGauges);
	
	
	
	
	argparser_get_keyword(pp, "-gaugeCenter");
	for(ind_gab = 0; ind_gab < o->nGauges; ind_gab++){
		o->gauge_atributes[ind_gab].gaugeCenter.c[0] = argparser_get_next_double(pp, 0.0, 100000.0);
		o->gauge_atributes[ind_gab].gaugeCenter.c[1] = argparser_get_next_double(pp, 0.0, 100000.0);
		fprintf(stderr, "  -gaugeCenter %lf %lf \\\n", o->gauge_atributes[ind_gab].gaugeCenter.c[0], o->gauge_atributes[ind_gab].gaugeCenter.c[1]);
	}
	
	
	
	argparser_get_keyword(pp, "-gaugeRadius");
	for(ind_gab = 0; ind_gab < o->nGauges; ind_gab++){
		o->gauge_atributes[ind_gab].gaugeRadius = argparser_get_next_double(pp, 1.0, 100000.0);
		fprintf(stderr, "  -gaugeRadius %lf \\\n", o->gauge_atributes[ind_gab].gaugeRadius);
	}

        /*Read gauge stretch
	OBS: It is important that the stretch vector {gauge_stretch} points towards   the optical axis.
	*/
	
	int test_stretch = argparser_keyword_present(pp, "-gaugeStretch");
	for(ind_gab = 0; ind_gab < o->nGauges; ind_gab++){
		if(test_stretch){
			o->gauge_atributes[ind_gab].gaugeStretch.c[0] = argparser_get_next_double(pp, -100000.0, 100000.0);
			o->gauge_atributes[ind_gab].gaugeStretch.c[1] = argparser_get_next_double(pp, -100000.0, 100000.0);
		
		}else{
			o->gauge_atributes[ind_gab].gaugeStretch = (r2_t){{0,0}};
		}
		fprintf(stderr, "  -gaugeStretch %lf %lf\\\n", o->gauge_atributes[ind_gab].gaugeStretch.c[0],o->gauge_atributes[ind_gab].gaugeStretch.c[1]);
	}
	
	/*Read gauge view dir*/
	int test_viewDir = argparser_keyword_present(pp, "-gaugeViewDir");
	for(ind_gab = 0; ind_gab < o->nGauges; ind_gab++){
		if(test_viewDir){
			r3_t* viewDir = &(o->gauge_atributes[ind_gab].gaugeViewDir);
			viewDir->c[0] = argparser_get_next_double(pp, -1.0, 1.0);
			viewDir->c[1] = argparser_get_next_double(pp, -1.0, 1.0);
			viewDir->c[2] = argparser_get_next_double(pp, -1.0, 1.0);
			double m = r3_dir(viewDir,viewDir);
			if(fabs(m -1) >= 1.0e-6){ 
				fprintf(stderr,"!! WARNING - ViewDir (Gauge %d ) not normalized = ",ind_gab);
				r3_print(stderr,viewDir);
				fprintf(stderr,"\n");
			}
		}else{
			o->gauge_atributes[ind_gab].gaugeViewDir = (r3_t){{0,0,1}};
		}
	}
	
	o->rawPixels = argparser_keyword_present(pp, "-rawPixels");
	o->yflip = argparser_keyword_present(pp, "-yflip");
	
	int test_albedo = argparser_keyword_present(pp, "-gaugeAlbedo");
	for( ind_gab = 0; ind_gab < o->nGauges; ind_gab++){
		if (test_albedo) {
			o->gauge_atributes[ind_gab].gaugeAlbedo.c[0] = argparser_get_next_double(pp, 0.001, 1.000);
			o->gauge_atributes[ind_gab].gaugeAlbedo.c[1] = argparser_get_next_double(pp, 0.001, 1.000);
			o->gauge_atributes[ind_gab].gaugeAlbedo.c[2] = argparser_get_next_double(pp, 0.001, 1.000);
		} else {
			o->gauge_atributes[ind_gab].gaugeAlbedo.c[0] = 1.000; 
			o->gauge_atributes[ind_gab].gaugeAlbedo.c[1] = 1.000; 
			o->gauge_atributes[ind_gab].gaugeAlbedo.c[2] = 1.000;
		}
		fprintf(stderr, "  -gaugeAlbedo %5.3lf %5.3lf %5.3lf \\\n",
			o->gauge_atributes[ind_gab].gaugeAlbedo.c[0],
			o->gauge_atributes[ind_gab].gaugeAlbedo.c[1],
			o->gauge_atributes[ind_gab].gaugeAlbedo.c[2]
		);
	}


	
	argparser_get_keyword(pp, "-nLights");
	o->nLights = argparser_get_next_int(pp, 1, 1000);
	fprintf(stderr, "  -nLights %d \\\n", o->nLights);
	
	argparser_get_keyword(pp, "-prefix");
	o->prefix = argparser_get_next(pp);
	fprintf(stderr, "  -prefix %s \\\n", o->prefix);
	
	

	fprintf(stderr,"    -nGauges %d \\\n",o->nGauges);
	
	argparser_get_keyword(pp, "-gaugeImages");
	o->gauge_image_name = (char**) malloc(sizeof(char*)*(o->nLights*o->nGauges));
	fprintf(stderr, "Gabaritos:\n");
	int ind;
	for (ind = 0; ind < (o->nLights*o->nGauges); ind++){
		o->gauge_image_name[ind] = argparser_get_next(pp);
		fprintf(stderr, "    G[%02d] = %s\n",  ind, o->gauge_image_name[ind]);
	}
	
	if(argparser_keyword_present(pp, "-gaugeDirections") ){
		o->gauge_direction_name = (char**) malloc(sizeof(char*)*(o->nLights*o->nGauges));
		fprintf(stderr, "Gabaritos - DIRECOES:\n");
		for (ind = 0; ind < (o->nLights*o->nGauges); ind++){
			o->gauge_direction_name[ind] = argparser_get_next(pp);
			fprintf(stderr, "    Dir[%02d] = %s\n",  ind, o->gauge_image_name[ind]);
		}
		o->gaugeDirections = (r3_t*)malloc(sizeof(r3_t)*o->nLights);
	}
	else{
		o->gauge_direction_name = NULL;
		o->gaugeDirections = NULL;
	} 
	
	o->thetaMax = M_PI/2.0;
	o->lightingModel = -1;
	o->lightingParams = NULL;
	if(argparser_keyword_present(pp, "-useLightingModel") ){
	  char* lighting_model_name = argparser_get_next(pp);
	  if(strcmp(lighting_model_name,"compact") == 0){
	    o->lightingModel = COMPACT_MODEL;
	  }else if(strcmp(lighting_model_name,"harmonic") == 0){
	    o->lightingModel = HARMONIC_MODEL;
	  }else if(strcmp(lighting_model_name,"stereopoly") == 0){
	    o->lightingModel = STEREOPOLY_MODEL;  
	  }else if(strcmp(lighting_model_name,"radial") == 0){
	    o->lightingModel = RADIALBASIS_MODEL;
	  }else {
	    argparser_error(pp,"Must specify a valid lighting model type !");
	  }
	  o->gauge_lparms_name = (char**) malloc(sizeof(char*)*(o->nLights*o->nGauges));
	  fprintf(stderr, "Gabaritos - PARAMS:\n");
	  for (ind = 0; ind < (o->nLights*o->nGauges); ind++){
	    o->gauge_lparms_name[ind] = argparser_get_next(pp);
	    fprintf(stderr, "    Dir[%02d] = %s\n",  ind, o->gauge_lparms_name[ind]);
	  }
	}
	
	o->gaugePos = (r2_t*)malloc(sizeof(r2_t)*(o->nGauges));
	if (argparser_keyword_present(pp, "-gaugePos")) {
		for (ind =0; ind < o->nGauges; ind++){
		o->gaugePos[ind].c[0] = argparser_get_next_double(pp, -100000.0, +200000.0);
		o->gaugePos[ind].c[1] = argparser_get_next_double(pp, -100000.0, +200000.0);
		fprintf(stderr, "    [%02d] %lf  %lf \\\n", ind, o->gaugePos[ind].c[0], o->gaugePos[ind].c[1]);
		}
	} else {
		if (o->nGauges > 1) { 
			argparser_error(pp, "Must spcify \"-gaugePos\" when there are multiple gauges");
		} else {
			for (ind =0; ind < o->nGauges; ind++){
				o->gaugePos[ind] = (r2_t){{ 0,0 }};
			}
		}
	}
	
	o->saveSamplePoints = argparser_keyword_present(pp,"-saveSamplePoints");
	if(o->saveSamplePoints){
	  fprintf(stderr,"SAVING SAMPLE POINTS\n");
	}
	
	argparser_get_keyword(pp, "-gamma");
	o->gamma = argparser_get_next_double(pp, 0.100, 9.000);
	o->bias = (o->gamma == 1.000 ? 0.000 : VIEW_BIAS); /* Hack... */
	fprintf(stderr, "Gamma das imagens:%lf bias: %lf\n", o->gamma, o->bias);
	
	o->gray = argparser_keyword_present(pp, "-gray");
	
	argparser_finish(pp);
	
	return o;
}

void salva_samplepoints(char* nome_samplepoints_file, r2_t* pontos, long int num_pontos);

void salva_samplepoints(char* nome_samplepoints_file, r2_t* pontos, long int num_pontos){
  FILE* arq = open_write(nome_samplepoints_file,TRUE);
  fprintf(arq,"%ld\n",num_pontos);
  long int i;
  for(i = 0; i < num_pontos; i++){
    fprintf(arq,"%9.6lf %9.6lf\n",pontos[i].c[0],pontos[i].c[1]);
  }
  fclose(arq);
}


int main(int argc,char** argv){
	options_t *o = parse_options(argc, argv);
		
	/* Lê as imagens dos gabaritos: */
	
	float_image_t *** G_list = (float_image_t ***)malloc(sizeof(float_image_t **)*o->nGauges);
	int ind;
	for (ind = 0; ind < o->nGauges; ind++) {
		char** pnames;
		pnames= o->gauge_image_name + (ind*o->nLights);
		G_list[ind]= float_pnm_image_list_read( o->nLights,pnames,FALSE, o->gamma, o->bias, TRUE,TRUE,FALSE);
// 		assert(G_list[ind][0]->sz[0] == 3);
	}
	int NC = G_list[0][0]->sz[0];
	
	
	int ind_gabarito;
	if(o->yflip){
	  //correct the gauge coordinates to current orientation system
	  for( ind_gabarito = 0; ind_gabarito < o->nGauges; ind_gabarito++){
		  r2_t gaugeCenter = o->gauge_atributes[ind_gabarito].gaugeCenter;
		  double gaugeRadius = o->gauge_atributes[ind_gabarito].gaugeRadius;
		  r2_t gaugeStretch = o->gauge_atributes[ind_gabarito].gaugeStretch;
		  r3_t gaugeViewDir = o->gauge_atributes[ind_gabarito].gaugeViewDir;
		  int NX,NY;
		  NX = G_list[ind_gabarito][0]->sz[1];
		  NY = G_list[ind_gabarito][0]->sz[2];
		  
		  gaugeCenter.c[1] = NY - gaugeCenter.c[1] -1;
		  gaugeStretch.c[1] = -gaugeStretch.c[1];
		  gaugeViewDir.c[1] = -gaugeViewDir.c[1];
		  
		  o->gauge_atributes[ind_gabarito].gaugeCenter = gaugeCenter;
		  o->gauge_atributes[ind_gabarito].gaugeRadius = gaugeRadius;
		  o->gauge_atributes[ind_gabarito].gaugeStretch = gaugeStretch;
		  o->gauge_atributes[ind_gabarito].gaugeViewDir = gaugeViewDir;
	  }
	}
	
	fprintf(stderr, "Imagens de gabaritos lidas.\n");
	if(o->gaugeDirections != NULL){
		int i;
		for(i = 0; i < o->nLights;i++){
			FILE* arq_dir;
			fprintf(stderr,"Abrindo arquivo de direcao %s ...\n",o->gauge_direction_name[i]);
			arq_dir = fopen(o->gauge_direction_name[i],"rt");
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
			o->gaugeDirections[i].c[0] = dx;
			o->gaugeDirections[i].c[1] = dy; 
			o->gaugeDirections[i].c[2] = dz;
			fprintf(stderr,"Direction: %+8.5f %+8.5f %+8.5f\n",dx,dy,dz);
			fclose(arq_dir);
		}
	}

	
	approx_model_t* lm = NULL;
	int current_model = o->lightingModel; 
	if( current_model > -1){
	 fprintf(stderr,"Using lighting model\n");
	  lm = create_approx_lighting_model(o->lightingModel);
	  o->lightingParams = (void***)malloc(sizeof(void**)*3);
	  int c;
	  
	  for(c = 0; c < NC; c++){
	    int process_channel = strchr(o->channels, "RGB"[c]) != NULL;
	    if(process_channel){
	      o->lightingParams[c] = (void**)malloc(sizeof(void*)*(o->nLights));
	    }
	  }
	  
	  int  ll;
	  for(ll = 0; ll < o->nLights;ll++){
	    FILE* arq_param = open_read(o->gauge_lparms_name[ll],TRUE);
	    for(c = 0; c < NC ; c++){
	      int process_channel = strchr(o->channels, "RGB"[c]) != NULL;
	      if(process_channel){
		o->lightingParams[c][ll] = lm->read_parameters(arq_param);
	      }
	    }
	    fclose(arq_param);
	  }
	}

	int canal;
	//int ind_gabarito;
	for( ind_gabarito = 0; ind_gabarito < o->nGauges; ind_gabarito++){
  		for (canal = 0; canal < NC; canal++) {
			r2_t *pontos;
			r3_t *normais;
			int num_pontos;
			fprintf(stderr, "Gerando os pontos de amostragem...\n");
			r3_t view_dir = o->gauge_atributes[ind_gabarito].gaugeViewDir;
			/* Gera máscara do gabarito: */
			float_image_t  *M = gera_mascara_do_gabarito(
						G_list[ind_gabarito][0]->sz[1],
						G_list[ind_gabarito][0]->sz[2],
						o->gauge_atributes[ind_gabarito].gaugeCenter,
						o->gauge_atributes[ind_gabarito].gaugeRadius,
						o->gauge_atributes[ind_gabarito].gaugeStretch);
			if( o->rawPixels){
				gera_pontos_com_mascara(
					M,
					o->gauge_atributes[ind_gabarito].gaugeCenter,
					o->gauge_atributes[ind_gabarito].gaugeRadius,
					o->gauge_atributes[ind_gabarito].gaugeStretch,
					&pontos, &normais, &num_pontos,view_dir);
			}else{
				gera_pontos_no_gabarito_eliptico(
						o->tableSize,M_PI/2.0,
						o->gauge_atributes[ind_gabarito].gaugeCenter,
						o->gauge_atributes[ind_gabarito].gaugeRadius,
						o->gauge_atributes[ind_gabarito].gaugeStretch,
				 &pontos, &normais, &num_pontos,view_dir);
			}
	
			if(o->saveSamplePoints){
			  char* nome_samplepoints_file = NULL;
			  char *nome_samplepoints_file = jsprintf("%s_%d_G%s_SamplePoints.txt", o->prefix, canal,o->gaugeTags[ind_gabarito]);
			  salva_samplepoints(nome_samplepoints_file,pontos,num_pontos);
			}
			
			char *nome_arq_M = NULL;
			char *nome_arq_M = jsprintf("%s_G%d_gauge_mask.ppm", o->prefix,ind_gabarito);
			float_pnm_image_write(nome_arq_M, M,FALSE, 1.000, 0.000,TRUE,TRUE,FALSE);
			free(nome_arq_M);
    			/* Devemos proessar este canal? */
    			int process_channel = strchr(o->channels, "RGB"[canal]) != NULL;
			if(process_channel){
				char *nome_arq_table = NULL;
				char *nome_arq_table = jsprintf("%s_%d_G%s_TableData.txt", o->prefix, canal,o->gaugeTags[ind_gabarito]);
				Tabela* tab;
				fprintf(stderr,"---------------------------------------------------------------------");
    				fprintf(stderr, "INICIO DO PROCESSAMENTO DO CANAL %d: \n", canal);
				
				bool_t useRealGauges = !((o->lightingParams != NULL) || (o->gaugeDirections != NULL));
				if(o->gaugeDirections != NULL){
				  fprintf(stderr,"Using DIRECTIONS\n");
				}
				if(o->lightingParams != NULL){
				  fprintf(stderr,"Using MODEL\n");
				}
				
				if(useRealGauges){
					fprintf(stderr,"Using REAL DATA\n");
					tab = cria_tabela
					( G_list[ind_gabarito],
					M,
					o->gauge_atributes[ind_gabarito].gaugeAlbedo.c[canal],
					o->nLights,
					canal,
					pontos,
					normais,
					num_pontos,
					view_dir,
					!o->rawPixels
					);
					
				}else if(o->lightingParams != NULL){
				  
				  fprintf(stderr,"Lighting Model selected - ");
				  //Just for indicate what goes on.
				  if(o->lightingModel == COMPACT_MODEL){
				    fprintf(stderr,"COMPACT\n");
				  }else if(o->lightingModel == HARMONIC_MODEL){
				    fprintf(stderr,"HARMONIC\n");
				  }
				  else if(o->lightingModel == STEREOPOLY_MODEL){
				    fprintf(stderr,"STEREOPOLY\n");
				  }
				  else if( o->lightingModel == RADIALBASIS_MODEL){
				    fprintf(stderr,"RADIAL BASIS\n");
				  }else{
				    fprintf(stderr,"No Model selected !\n");
				    assert(FALSE);
				  }
				  
				  tab = cria_tabela_from_model(o->nLights,o->tableSize,view_dir,o->thetaMax,lm,o->lightingParams[canal]);
				  
				}else if(o->gaugeDirections != NULL){
					fprintf(stderr,"Virtual Table selected\n");
					tab = cria_tabela_virtual
						(o->gaugeDirections,
						o->gauge_atributes[ind_gabarito].gaugeRadius,
						o->gauge_atributes[ind_gabarito].gaugeCenter,
						o->gauge_atributes[ind_gabarito].gaugeAlbedo.c[canal],
						o->nLights,
						canal,
						pontos, 
						normais,
						num_pontos,
						view_dir
					);
				}
				if (tab == NULL) { fprintf(stderr, "Erro na geracao de tabela!\n");  return 1; }
				SaveTable(nome_arq_table,tab,TRUE);
    				fprintf(stderr, "Tabela gerada.\n");
			}
		}
	}
	return 0;
}

