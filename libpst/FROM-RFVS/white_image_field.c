#define PROG_NAME "white_image_field"

#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <rn.h>
#include <r2.h>
#include <rmxn.h>
#include <math.h>
#include <float_image.h>
#include <argparser.h>
#include <float_pnm_image_io.h>
#include <least_squares_nd.h>
#include <jsfile.h>

#define PROG_HELP \
  PROG_NAME " \\\n" \
  "	-prefix {FILE_PREFIX} \\\n" \
  "	-inputImage {INPUT IMAGE} \\\n" \
  "	-maskImage {INPUT IMAGE} \\\n" \
  "	[-nImages {NIMAGES} \\\n" \
  "    -adjustImages {FILENAME}...]  \\\n" \
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

struct options_t{
  char** inputImage;
  char* maskImage;
  char*gaugeData;
  char* polyData;
  int degree;
  char* prefix;
  int* nImages;
  int nLights;
  char*** adjustImages;
  bool_t preserveSignatures;
  bool_t saveChannelAndLight;
  bool_t generateWhite;
};

typedef struct options_t options_t;

struct gaugeData_t{
	r2_t gaugeCenter;
	double gaugeRadius;
};

typedef struct gaugeData_t gaugeData_t;

void ReadGaugeData(char* filename, gaugeData_t** gp,int* nGauges){
  int n,i;
  FILE* arq = open_read(filename,TRUE);
  fscanf(arq,"%d",&n);
  *nGauges = n;
  gaugeData_t* g = (gaugeData_t*)malloc(sizeof(gaugeData_t)*n);
  for(i = 0; i < n; i++){
    fscanf(arq,"%lf %lf %lf",&(g[i].gaugeRadius),&(g[i].gaugeCenter.c[0]),&(g[i].gaugeCenter.c[1]));
    fprintf(stderr,"Read %lf %lf %lf\n",(g[i].gaugeRadius),(g[i].gaugeCenter.c[0]),(g[i].gaugeCenter.c[1]));
  }
  *gp = g;
  fclose(arq);
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

void WritePolyData(FILE* arq,poly_function_t* polym){
 	
  fprintf(arq,"%d ",(polym->dimensions));
  fprintf(arq,"%d\n",(polym->num_coefs));
  int icf;
  for(icf = 0; icf < (polym->num_coefs);icf++){
    int d;
    for(d = 0; d < (polym->dimensions); d++){
      fprintf(arq,"%3.4lf ",(polym->coefs[icf][d]));
    }
    fprintf(arq,"%9.6e\n",(polym->weights[icf]));
  }
}



void GenerateSystemFromInputData(float_image_t* img, float_image_t* msk,gaugeData_t* gauges, int nGauges,int* N,double*** Xp,  double** Fp,int channel){
  int NX,NY;
  
  NX = img->sz[1];
  NY = img->sz[2];

  int count = 0;
  int x,y;
  int n;
  double ** X = (double**) malloc(sizeof(double*)* nGauges);
  double *F = (double*) malloc(sizeof(double)* nGauges);


  for(n = 0;n < nGauges; n++){	
    double radius = gauges[n].gaugeRadius; 
    double cx  =  gauges[n].gaugeCenter.c[0];
    double cy  =  gauges[n].gaugeCenter.c[1];
    int iniX = cx - radius;
    int fimX = cx + radius;
    int iniY = cy - radius;
    int fimY = cy + radius;
  //  fprintf(stderr,"Inspecting inside [%d,%d] to [%d,%d]\n",iniX,iniY,fimX,fimY);
    double valMax = -1;
    int maxX,maxY;
    for(x = iniX; x < fimX; x++){
      for( y = iniY; y < fimY; y++){	
	if( ((x - cx)*(x - cx) + (y - cy)*(y - cy)) < (radius*radius) ){
	  if(float_image_get_sample(msk,0,x,y) > 0.5 ){
	    double val = float_image_get_sample(img,channel,x,y);
	    if(val > valMax){
	      valMax = val;
	      maxX = x;
	      maxY = y;
	    }
	  }
	}
      }
    }
    fprintf(stderr,"Max Point[%d]: at (%d %d) = %lf\n",n,maxX,maxY,valMax);
    X[n] = (double*)malloc(sizeof(double)*2);
    X[n][0] = maxX;
    X[n][1] = maxY;
    F[n] = valMax;

  }



  *N = nGauges;
  

  *Xp = X;
  *Fp = F;
    
  
  

  return;
  
}


void GenerateSystemFromRAWData(float_image_t* img, float_image_t* msk,int* N,double*** Xp,  double** Fp,int channel){
  int NX,NY;
  
  NX = img->sz[1];
  NY = img->sz[2];

  
  
  int count = 0;
  int x,y;
  int n;

  for(x = 0 ; x < NX; x++){
    for(y = 0; y < NY; y++){
      double val = float_image_get_sample(msk,0,x,y);
      if(val > 0.5){
	count++;
      }
    }
  }
  n = count;

  double ** X = (double**) malloc(sizeof(double*)* n);
  double *F = (double*) malloc(sizeof(double)* n);

  count = 0;
  for(x = 0 ; x < NX; x++){
    for(y = 0; y < NY; y++){
      double val = float_image_get_sample(msk,0,x,y);
      if(val > 0.5){
	X[count] = (double*)malloc(sizeof(double)*2);
	X[count][0] = x;
	X[count][1] = y;
	F[count] = float_image_get_sample(img,channel,x,y);
	count++;
      }
    }
  }
    



  *N = n;
  

  *Xp = X;
  *Fp = F;
    
  
  

  return;
  
}


float_image_t* AdjustImage(float_image_t* img, poly_function_t** pv,double factor){
  int NC,NX,NY;
  
  NC = img->sz[0];
  NX = img->sz[1];
  NY = img->sz[2];

  float_image_t* adj_img = float_image_new(NC,NX,NY);
  int bug = 0;
  int c,x,y;
  for(c = 0 ; c < NC; c++){
    poly_function_t* p = pv[c];
    for(x = 0; x < NX; x++){
      for(y = 0; y < NY ; y++){
	double X[] = {x,y};
	double white = EvaluatePvalue(p,X);
	white = white*factor;
	double val = float_image_get_sample(img,c,x,y)/(white*1.1);
	
	if( float_image_get_sample(img,c,x,y) > white) {
	   // val = 1.0;
	    bug = 1;
	}
	float_image_set_sample(adj_img,c,x,y,val);
      }
    }
  }

  if(bug==1) fprintf(stderr,"OVEREXPOSED IMAGE !\n");
  return adj_img;

}


double CorrectImageSet(float_image_t** imgs, int n){
  double max_value;
  max_value = 1.0;
  int x,y,c,l;
  int NC,NX,NY;
  NC = imgs[0]->sz[0];
  NX = imgs[0]->sz[1];
  NY = imgs[0]->sz[2];
  for(l = 0; l < n; l++){
    for( c = 0; c < NC; c++){
      for( x = 0; x < NX; x++){
	for(y = 0; y < NY; y++){
	  double value = float_image_get_sample(imgs[l],c,x,y);
	  if(value > max_value) max_value = value;
	}
      }
    }
  }

  fprintf(stderr,"Maximum Value found was %3.4lf - scaling by %3.4lf \n",max_value,1/max_value);
  if(max_value < 1.0e-10) max_value = 1.0;
  /*Now we will rescale everything by the maxvalue, preserving teh signatures*/
  
  for(l = 0; l < n; l++){
    for( c = 0; c < NC; c++){
      for( x = 0; x < NX; x++){
	for(y = 0; y < NY; y++){
	  double value = float_image_get_sample(imgs[l],c,x,y);
	  value = value/max_value;
	  float_image_set_sample(imgs[l],c,x,y,value);
	}
      }
    }
  }

  return 1.0/max_value;
}



float_image_t* GenerateWhiteImage(int NC,int NX,int NY, poly_function_t** pv){
  float_image_t* adj_img = float_image_new(NC,NX,NY);

  int c,x,y;
  for(c = 0 ; c < NC; c++){
    poly_function_t* p = pv[c];
    for(x = 0; x < NX; x++){
      for(y = 0; y < NY ; y++){
	double X[] = {x,y};
	double val = EvaluatePvalue(p,X);
	float_image_set_sample(adj_img,c,x,y,val);
      }
    }
  }

  return adj_img;

}



options_t* parse_args(int argc, char** argv){
  options_t* o = (options_t*)malloc(sizeof(options_t));
  argparser_t *pp = argparser_new(stderr, argc, argv);
  argparser_set_help(pp, PROG_HELP);
  argparser_set_info(pp, PROG_INFO);
  argparser_process_help_info_options(pp);

  argparser_get_keyword(pp, "-prefix");
  o->prefix = argparser_get_next(pp);
  
  argparser_get_keyword(pp, "-maskImage");
  o->maskImage = argparser_get_next(pp);
  
  o->gaugeData =  NULL;
  if(argparser_keyword_present(pp, "-gaugeData")){
    //argparser_get_keyword(pp, "-gaugeData");
    o->gaugeData = argparser_get_next(pp);
  }
  
  o->polyData = NULL;
  if(argparser_keyword_present(pp, "-polyData")){
    o->polyData = argparser_get_next(pp);
  }
  o->saveChannelAndLight = argparser_keyword_present(pp,"-saveChannelAndLight");

  o->preserveSignatures = argparser_keyword_present(pp, "-preserveSignatures");
  o->generateWhite = argparser_keyword_present(pp, "-generateWhite");
  
  argparser_get_keyword(pp, "-nLights");
  o->nLights = argparser_get_next_int(pp, 1, 100000);
  o->nImages = (int*)malloc(sizeof(int)*(o->nLights));
  o->adjustImages = (char***) malloc(sizeof(char**)*(o->nLights));
  o->inputImage = (char**)malloc(sizeof(char*)*(o->nLights));
  int l;
  for(l = 0; l < o->nLights; l++){
    argparser_get_keyword_next(pp, "inputImage");
    o->inputImage[l] = argparser_get_next(pp);
    argparser_get_keyword_next(pp, "nImages");
    o->nImages[l] = argparser_get_next_int(pp, 0, 100000);
    o->adjustImages[l] = (char**)malloc(sizeof(char*)*(o->nImages[l]));
    if(o->nImages[l] > 0){
      argparser_get_keyword_next(pp, "adjustImages");
      int i;
      for(i = 0 ; i < o->nImages[l];i++){
	o->adjustImages[l][i] = argparser_get_next(pp);
      }
    }
  }


  o->degree = 2;
  if(argparser_keyword_present(pp, "-degree")){
    o->degree = argparser_get_next_int(pp, 1, 100000);
  }

  argparser_finish(pp);
  return o;
}

void SaveFNIImage(char* name, float_image_t* img){
  FILE* arq = open_write(name,TRUE);
  float_image_write(arq,img);
  fclose(arq);
}



int main(int argc, char** argv){
  bool_t FLIP = TRUE;
  options_t* o = parse_args(argc,argv);
  
  float_image_t* msk_image = float_pnm_image_read(o->maskImage,TRUE, 1.0, 0.0, FLIP,TRUE,FALSE); 
  
  int nGauges;
  
  gaugeData_t* gv;
  if(o->gaugeData != NULL){
    ReadGaugeData(o->gaugeData,&gv,&nGauges);
  }
   poly_function_t*** ParamPoly =  NULL;
  if(o->polyData != NULL){
    ParamPoly = ReadPolyData(o->polyData);
  }

  int c;
  float_image_t*** adj_img = (float_image_t***)malloc(sizeof(float_image_t**)*(o->nLights));

  FILE* polyfile;
  char* polyfilename;
  char *polyfilename = jsprintf("%s_PolyData.txt",o->prefix);
  polyfile = open_write(polyfilename,TRUE);
  //fprintf(polyfile,"%d\n",o->nLights);

  poly_function_t*** Computed_Poly;
  Computed_Poly = (poly_function_t***)malloc(sizeof(poly_function_t**)*(o->nLights));
  int l;
  int NumChannels;
  for(l = 0 ; l < o->nLights; l++){
    float_image_t* in_image = float_pnm_image_read(o->inputImage[l],FALSE, 1.0, 0.0, FLIP,TRUE,FALSE);
    int NC,NX,NY;
    NC = in_image->sz[0];
    NumChannels = NC;
    NX = in_image->sz[1];
    NY = in_image->sz[2];
    poly_function_t** pv =  (poly_function_t**)malloc(sizeof(poly_function_t*)*NC);
    Computed_Poly[l] = pv;
    adj_img[l] = (float_image_t**)malloc(sizeof(float_image_t*)*(o->nImages[l]));
    if(o->polyData == NULL){
      for(c = 0; c < NC; c++){
	double** X;
	double* F;
	int N;
      
	if( o->gaugeData != NULL){
	  GenerateSystemFromInputData(in_image, msk_image,gv,nGauges,&N,&X,&F,c);
	}else{
	  GenerateSystemFromRAWData(in_image,msk_image,&N,&X,&F,c);
	}
	pv[c] = LS_PolyFitting(2,N,o->degree,X,F);
	PrintfPolyFunction(stderr,pv[c]);
	free(F);
	int n;
	for(n = 0; n < N; n++){
	  free(X[n]);
	}
	free(X);
      }
    }else{
      for(c = 0; c < NC; c++){
	pv[c] = ParamPoly[l][c];
	PrintfPolyFunction(stderr,pv[c]);
      }
    }

    /*fprintf(polyfile,"%d\n",NC);
    for(c = 0; c < NC; c++){
      WritePolyData( polyfile,pv[c]);
    }*/
     

    if(o->generateWhite){
      float_image_t* white_image =  GenerateWhiteImage(NC,NX,NY,pv);
      char* white_image_filename = NULL;
      char *white_image_filename = jsprintf("%s_L%d_white.ppm", o->prefix,l);
      float_pnm_image_write(white_image_filename, white_image,FALSE, 1.0, 0.0,FLIP,TRUE,FALSE); 
      char *white_image_filename = jsprintf("%s_L%d_white.fni", o->prefix,l);
      SaveFNIImage(white_image_filename,white_image);
    }

    int i;
    for(i = 0; i < o->nImages[l]; i++){
      float_image_t* img = float_pnm_image_read(o->adjustImages[l][i],FALSE, 1.0, 0.0, FLIP,TRUE,FALSE);
      adj_img[l][i] = AdjustImage(img,pv,1.0);
      char* adj_image_filename = NULL;
      char *adj_image_filename = jsprintf("%s_L%d_%d_corrected.ppm", o->prefix,l,i);
      float_pnm_image_write(adj_image_filename, adj_img[l][i],FALSE, 1.0, 0.0,FLIP,TRUE,FALSE); 
      char *adj_image_filename = jsprintf("%s_L%d_%d_corrected.fni", o->prefix,l,i);
      SaveFNIImage(adj_image_filename,adj_img[l][i]);
      float_image_free(img);
      //float_image_free(adj_img);
    }
   
    
  }
  

  if(o->preserveSignatures){
    int num_sets = o->nImages[0];
    for(l = 0; l < o->nLights; l++){
      if(num_sets <  o->nImages[l]) num_sets =  o->nImages[l];
    }

    int n;
    for(n = 0; n < num_sets; n++){
      float_image_t* img_set[o->nLights];
      int i;
      for(i = 0; i < o->nLights; i++){
	img_set[i] = adj_img[i][n];
      }
      CorrectImageSet(img_set,o->nLights);
      for(i = 0; i < o->nLights; i++){
	char* adj_image_filename = NULL;
	char *adj_image_filename = jsprintf("%s_L%d_%d.ppm", o->prefix,i,n);
	float_pnm_image_write(adj_image_filename, adj_img[i][n],FALSE, 1.0, 0.0,FLIP,TRUE,FALSE); 
      }
    }
  }
  
  if(o->saveChannelAndLight){
    fprintf(polyfile,"%d\n",NumChannels);
    for(c = 0; c< NumChannels; c++){
      fprintf(polyfile,"%d\n",o->nLights);
      for(l = 0; l < o->nLights; l++){
	 WritePolyData( polyfile,Computed_Poly[l][c]);
      }
    }
  }else{
    fprintf(polyfile,"%d\n",o->nLights);
    for(l = 0; l < o->nLights; l++){
      fprintf(polyfile,"%d\n",NumChannels);
      for(c = 0; c< NumChannels; c++){
	 WritePolyData( polyfile,Computed_Poly[l][c]);
      }
    }
  }

  fclose(polyfile);
  return 0;
}

