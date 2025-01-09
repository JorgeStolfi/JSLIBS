#define _GNU_SOURCE
#define PROG_NAME "compute_rms"

#define _GNU_SOURCE 

#include <stdio.h>
#include <stdlib.h>
#include <float_image.h>
#include <argparser.h>
#include <pst_img_graph.h>
#include <jsfile.h>
#include <string.h>
#include <float.h>
#include <sys_stats.h>
#include <assert.h>
#include <rn.h>



#define PROG_HELP \
  PROG_NAME " \\\n" \
  "	-ref {REFERENCE_MAP} \\\n" \
  "	-cmp {COMPARISON_MAP} \\\n" \
  "	-out {OUTPUT_MAP} \\\n" \
  "	[-weight {WEIGHT_MAP}] \\\n" \
  "	[-channels {Chanel0,Channel1,Channel2...}] \\\n" \
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
  char* refMap;
  char* cmpMap;
  char* outMap;
  char* wgtMap;
  int numChannels;
  int* channels;
};

typedef struct options_t options_t;

options_t* parse_args(int argc, char** argv);
options_t* parse_args(int argc, char** argv){
  options_t* o = (options_t*)malloc(sizeof(options_t));
  argparser_t *pp = argparser_new(stderr, argc, argv);
  argparser_set_help(pp, PROG_HELP);
  argparser_set_info(pp, PROG_INFO);
  argparser_process_help_info_options(pp);
  
  argparser_get_keyword(pp, "-ref");
  o->refMap = argparser_get_next(pp);
  
  argparser_get_keyword(pp, "-cmp");
  o->cmpMap = argparser_get_next(pp);
  
  argparser_get_keyword(pp, "-out");
  o->outMap = argparser_get_next(pp);
  
  o->wgtMap = NULL;
  if(argparser_keyword_present(pp, "-weight")){
    o->wgtMap = argparser_get_next(pp);
  }
  
  o->numChannels = 0;
  o->channels = NULL;
  if(argparser_keyword_present(pp, "-channels")){
    char* origin_ch = argparser_get_next(pp);
    char* string_ch = origin_ch;
    char aux_str[strlen(origin_ch) + 1];
    /*first, count how much numbers we have*/
    int count=0;
    for(string_ch = origin_ch; string_ch != NULL; string_ch = string_ch){
      count++;
      string_ch = strchr(string_ch,',');
      if(string_ch != NULL){
	string_ch = string_ch +1;
      }
    }
    
    o->numChannels = count;
    o->channels = (int*)malloc(sizeof(int)*count);
    
    string_ch = origin_ch;
    count = 0;
    while(strlen(string_ch) != 0){
      int lenght_num = strcspn(string_ch,",");
      strncpy(aux_str,string_ch,lenght_num);
      aux_str[lenght_num] = '\0';
      if( sscanf(aux_str,"%d",&(o->channels[count])) != 1){
	char* str_error = NULL;
	char *str_error = jsprintf("Incorrect parameter \"%s\" for -channels.",aux_str);
	demand(FALSE,str_error);
      }
      string_ch= string_ch + lenght_num + 1;
      count++;
    }
    
    demand(count == o->numChannels, "Missing channel parameter");
   
  }
  
  argparser_finish(pp);
  
  return o;

}

float_image_t * nmap_map_compare
  ( float_image_t *AZ,
    float_image_t *BZ,
    float_image_t *W,
    double *sEZP
  );

float_image_t * nmap_map_compare
  ( float_image_t *AZ,
    float_image_t *BZ,
    float_image_t *W,
    double *sEZP
  )
  { 
    int NC = AZ->sz[0]; assert(BZ->sz[0] == NC);
    int NX = AZ->sz[1]; assert(BZ->sz[1] == NX);
    int NY = AZ->sz[2]; assert(BZ->sz[2] == NY);
    
    if (W != NULL)
      { assert(W->sz[0] == 1);
        assert(W->sz[1] == NX);
        assert(W->sz[2] == NY);
      }
      
    /* Compute the mean values of {AZ,BZ,EZ}: */
    double sum_W = 0;
    int x, y;
    for(y = 0; y < NY; y++)
      { for(x = 0; x < NX; x++)
          { /* Get relevant samples from original image: */
            double vW = 1.0;
	    if(W != NULL) vW = float_image_get_sample(W, 0, x, y);
            sum_W += vW;
          }
      }
   
    /* Unweighted mean of weights: */
    double avgW = sum_W/(NX*NY);
      
    /* Fill {EZ} and RMS values: */
    float_image_t *EZ = float_image_new(1, NX, NY);
    double sum_EZ2W = 0.0;
    for(y = 0; y < NY; y++)
      { for(x = 0; x < NX; x++)
          { /* Get relevant samples from original image: */
            double vA[NC];
            double vB[NC];
	    int c;
	    for(c = 0; c < NC; c++){
	      vA[c] = float_image_get_sample(AZ,c,x,y);
	      vB[c] = float_image_get_sample(BZ,c,x,y);
	    }
            double vW = 1.0;
	    if(W != NULL) vW = float_image_get_sample(W, 0, x, y);
            double vE = rn_dist(NC,vA,vB); 
                      /* Store difference in error image: */
            float_image_set_sample(EZ, 0, x, y, (float)(sqrt(vW/avgW)*vE));
            /* Accumuate squares: */
         
            sum_EZ2W += vW*vE*vE;
          }
      }
    /* Compute the RMS values and errors: */
    double sEZ = (sum_W == 0 ? 0.0 : sqrt(sum_EZ2W/sum_W));
  
    if (sEZP != NULL) { (*sEZP) = sEZ; }

    return EZ;
  }


float_image_t* height_map_shrink_one_pixel(float_image_t* Z);

float_image_t* height_map_shrink_one_pixel(float_image_t* Z){
	int NC = Z->sz[0];
	int NX = Z->sz[1];
	int NY = Z->sz[2];
	float_image_t* IJ = float_image_new(NC,NX-1,NY-1);
	int c;
	for( c = 0; c < NC; c++)
	{
		int x;
		for(x = 0; x < NX -1; x++)
		{
			int y;
      			for(y = 0; y < NY -1; y++)
			{
				double P00 = float_image_get_sample(Z,c,x,y);
				double P10 = float_image_get_sample(Z,c,x+1,y);
				double P01 = float_image_get_sample(Z,c,x,y+1);
				double P11 = float_image_get_sample(Z,c,x+1,y+1);

				double val = (P00 + P01 + P10 + P11)/4.0;
				float_image_set_sample(IJ,c,x,y,val);
      			}
		}
	}
	return IJ;
}

void add_border_weights(float_image_t* IW);

void add_border_weights(float_image_t* IW){
  int NX = IW->sz[1];
  int NY = IW->sz[2];
  int x,y;
  for(y = 0; y < NY; y++){
    for( x = 0; x < NX; x++){
      if( (x== 0) || (y == 0) ){
	float_image_set_sample(IW,0,x,y,0);
      }
    }
  }
}

float_image_t* readFNI(char* filename);


float_image_t* readFNI(char* filename){
	FILE* arq = open_read(filename,TRUE);
	float_image_t* img = float_image_read(arq);
	fclose(arq);
	return img;
}

void writeFNI(char* filename, float_image_t* img);

void writeFNI(char* filename, float_image_t* img){
	FILE* arq = open_write(filename,TRUE);
	float_image_write(arq,img);
	fclose(arq);
}

int main(int argc, char** argv){

	options_t* o = parse_args(argc,argv);
	char* hmREF_filename = o->refMap;
	char* hmCPR_filename = o->cmpMap;
	char* output_filename = o->outMap;
	char* wm_filename =o->wgtMap;

        /* First read the reference file */
	float_image_t* imgREF = readFNI(hmREF_filename);
	/* First read the compared file */
	float_image_t* imgCPR = readFNI(hmCPR_filename);
	/* Enforce same number of channels */
	if(imgREF->sz[0] != imgCPR->sz[0])
	{
		fprintf(stderr,"ERROR: Comparing maps with different number of channels\n");
		assert(FALSE);
	}


	float_image_t* imgW = NULL;
	if(wm_filename != NULL){
		imgW = readFNI(wm_filename);
	}

	
	int NC_REF = imgREF->sz[0];
	int NC_CPR = imgCPR->sz[0];
	
	assert(NC_REF == NC_CPR);
	
	int NX_REF = imgREF->sz[1];
	int NY_REF = imgREF->sz[2];

	int NX_CPR = imgCPR->sz[1];
	int NY_CPR = imgCPR->sz[2];

	int NX_WGT;
	int NY_WGT;

   	fprintf(stderr,"Height maps sizes: REF[%d x %d] and CPR[%d x %d]\n",NX_REF,NY_REF,NX_CPR,NY_CPR);
	if(imgW != NULL){ 
		NX_WGT = imgW->sz[1];
		NY_WGT = imgW->sz[2];
		fprintf(stderr,"weight map W[%d x %d]\n",NX_WGT,NY_WGT);
	}

	
	if ( (NX_REF != NX_CPR) || (NY_REF !=NY_CPR) )
	{
		/*The reference file may be one pixel bigger than the compared one so we reduce it, using the weights*/
		if ( (NX_REF == (NX_CPR+1))  && (NY_REF == (NY_CPR+1)) )
		{	
			fprintf(stderr,"Shrunk reference map: 1 pixel larger\n");
			float_image_t* imgREF_red = height_map_shrink_one_pixel(imgREF);
			float_image_free(imgREF);
			imgREF = imgREF_red;
			NX_REF--;
			NY_REF--;
			
		}
		else if ( (NX_REF == (NX_CPR-1))  && (NY_REF == (NY_CPR-1)) )
		{	/*It can happens that we want compare the inverse ! */
			fprintf(stderr,"Shrunk compared map: 1 pixel larger\n");
			float_image_t* imgCPR_red = height_map_shrink_one_pixel(imgCPR);
			float_image_free(imgCPR);
			imgCPR = imgCPR_red;
			NX_CPR--;
			NY_CPR--;
			
		}
		else
		{
			fprintf(stderr,"ERROR: Compared maps with incompatible sizes\n");
			assert(FALSE);
		}
	}
	
	if( imgW != NULL ){
		if( (NX_CPR != NX_WGT)  || (NY_CPR != NY_WGT) ){
			fprintf(stderr,"ERROR: Weight map with incompatible size\n");
			assert(FALSE);
		}
	}
	
	/* We will add a border zero, to avoid confusions.. it will create a weight map if needed */
	if(imgW == NULL){
	  imgW = float_image_new(1,NX_CPR,NY_CPR);
	  float_image_fill_channel(imgW,0,1.0);
	}
	add_border_weights(imgW);


	
	int NC = imgREF->sz[0];
	if(o->channels == NULL){
	  NC = 1;
	}else{
	  NC = o->numChannels;
	}
	
	fprintf(stderr,"NREF %d NC %d\n",NC_REF,NC);
	
	
	float_image_t* img_diff = float_image_new(NC,NX_CPR,NY_CPR);
	if( (NC_REF ==1 ) || (NC > 1)){
	  fprintf(stderr,"Comparing single-channel map\n");
	  int c;
	  for( c = 0; c < NC ; c++)
	  {
		  int ic = (o->channels == NULL ? 0 : o->channels[c]);
		  assert((ic >= 0) && (ic < NC_REF));
		  double sAZP,sBZP; /* The standard deviations of the values of {imgREF,imgCPR}, each from its own mean value. */
		  double sEZP; /* The root-mean-square value of the error */
		  double sreP; /* The relative error {sEZ/sMZ} where {sMZ = hypot(imgREF,imgCPR)/sqrt(2)}. */
		  
		  float_image_t* imgREF_channel = float_image_new(1,NX_REF,NY_REF);
		  float_image_t* imgCPR_channel = float_image_new(1,NX_CPR,NY_CPR);

		  float_image_set_channel(imgREF_channel, 0, imgREF, ic);
		  float_image_set_channel(imgCPR_channel, 0, imgCPR, ic);

		  float_image_t* img_diff_channel = pst_height_map_compare(
							  imgREF_channel,
							  imgCPR_channel,
							  imgW,
							  TRUE,
							  &sAZP,
							  &sBZP,
							  &sEZP,
							  &sreP
						  );
		  fprintf(stdout,"%9.6lf %9.6lf %9.6lf %9.6lf\n",sEZP,sreP,sAZP,sBZP);
		  
		  float_image_set_channel(img_diff, ic, img_diff_channel, 0);
				  
		  float_image_free(imgREF_channel);
		  float_image_free(imgCPR_channel);
		  float_image_free(img_diff_channel);
	  }
	}else{
	   fprintf(stderr,"Comparing multi-channel map\n");
	   double sEZP; /* The root-mean-square value of the error */
	   img_diff = nmap_map_compare(imgREF,imgCPR,imgW,&sEZP);
	   fprintf(stdout,"%9.6lf %9.6lf %9.6lf %9.6lf\n",sEZP,0.0,0.0,0.0);
	}
	
	
	writeFNI(output_filename,img_diff);

	return 0;
}
