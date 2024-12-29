#define PROG_NAME "fni_to_pov"

#define PROG_HELP \
  PROG_NAME " \\\n" \
  "    [-in {FNI INPUT IMAGE} ] \\\n" \
  "    [-out {INC OUTPUT IMAGE} ] \\\n" \
  "    [ -texture {PPM TEXTURE IMAGE} ] \\\n" \
  "    [ -mask {FNI or PPM IMAGE} ] \\\n" \
  "    [ -disableTexture ] \\\n" \
  "    [ -viewAzim ${VIEW AZIMUTH} ] \\\n" \
  "    [ -viewElev ${VIEW ELEV} ] \\\n" \
  "    [ -viewDist ${VIEW DIST} ] \\\n" \
  "    [ -stepZ {Z STEPS} ] \\\n" \
  "    [ -scaleZ {B STEPS} ] \\\n" \
  "    [ -channel {CHANNEL} ] \\\n" \
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


#define _GNU_SOURCE
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <limits.h>
#include <float.h>
#include <argparser.h>
#include <unistd.h>
#include <string.h>
#include <assert.h>
#include <jsfile.h>
#include <float_image.h>
#include <float_pnm_image_io.h>
#include <imagem_valores.h>
#include <imagem_vetores.h>


struct options_t {
  char*  in_file;
  char* out_file;
  char* tex_file;
  char* mask_file;
  int disableTexture;
  double scale;
  double viewAzim;
  double viewElev;
  double viewDist;
  double stepZ;
  double focusZ;
  double scaleZ;
  int channel;
};

typedef struct options_t options_t;

options_t* parse_options(int argc,char** argv);

options_t* parse_options(int argc,char** argv){
  argparser_t *pp = argparser_new(stderr, argc, argv);
  argparser_set_help(pp, PROG_HELP);
  argparser_set_info(pp, PROG_INFO);
  argparser_process_help_info_options(pp);
  
  options_t* o = (options_t*)malloc(sizeof(options_t));
  
  o->in_file = NULL;
  if(argparser_keyword_present(pp, "-in")){
    o->in_file = argparser_get_next(pp);
  }
  
  o->out_file = NULL;
  if(argparser_keyword_present(pp, "-out")){
    o->out_file = argparser_get_next(pp);
  }
  
  o->mask_file = NULL;
  if(argparser_keyword_present(pp, "-mask")){
    o->mask_file = argparser_get_next(pp);
  }
  
  o->tex_file = NULL;
  if(argparser_keyword_present(pp, "-texture")){
    o->tex_file = argparser_get_next(pp);
    if(strcmp(o->tex_file,"NONE") == 0){
      o->tex_file = NULL;
    }
  }
  
  
  
  o->disableTexture = argparser_keyword_present(pp, "-disableTexture");
   
  o->channel = 0;
  if(argparser_keyword_present(pp, "-channel")){
    o->channel = argparser_get_next_int(pp,0,10000);
  }
  
  o->scale = 1.0;
  if(argparser_keyword_present(pp, "-scale")){
    o->scale = argparser_get_next_double(pp,0,10000);
  }
  
  o->focusZ = 0.0;
  if(argparser_keyword_present(pp, "-focusZ")){
    o->focusZ = argparser_get_next_double(pp,0,10000);
  }
  
  o->stepZ = 0;
  o->scaleZ = 1;
  if(argparser_keyword_present(pp, "-stepZ")){
    o->stepZ = argparser_get_next_double(pp,-100000,100000);
  }
  
  if(argparser_keyword_present(pp, "-scaleZ")){
    o->scaleZ = argparser_get_next_double(pp,-100000,100000);
  }
  
  o->viewAzim = 120;
  o->viewElev =  30;
  o->viewDist = -1;
  
  if(argparser_keyword_present(pp, "-viewAzim")){
    o->viewAzim = argparser_get_next_double(pp,-100000,100000);
  }
  if(argparser_keyword_present(pp, "-viewElev")){
    o->viewElev = argparser_get_next_double(pp,-100000,100000);
  }
  if(argparser_keyword_present(pp, "-viewDist")){
    o->viewDist = argparser_get_next_double(pp,-100000,100000);
  }
  argparser_finish(pp);
  return o;
}





//void generatePOVMesh(FILE* arq,imagem_valores_t* height_map,  imagem_vetores_t* normal_map,float_image_t *img, int nx,int ny);
float_image_t* float_image_read_file(char* filename);

float_image_t* float_image_read_file(char* filename){
	FILE* arq_IN = stdin;
	
	if(filename != NULL){
	 arq_IN = open_read(filename,TRUE);
	}else{
	  fprintf(stderr,"Reading from standard input\n");
	}
    	
    	float_image_t *IN = float_image_read(arq_IN);
	if(filename != NULL){
	  fclose(arq_IN);
	}
	return IN;
}

float_image_t* scale_mask(float_image_t* mask_img);
float_image_t* scale_mask(float_image_t* mask_img){
  
  int WNX = mask_img->sz[1];
  int WNY = mask_img->sz[2];
  float_image_t* cm_img = float_image_new(1,WNX+1,WNY+1);
  int x,y;
  for(x = 0; x < WNX+1; x++){
    for(y = 0; y < WNY+1; y++){
      double w00,w11,w01,w10;
      w11 = w00 = w10 = w01  = 0;
      
      if( (x > 0) && (y > 0) ){ 
	w00 = float_image_get_sample(mask_img,0,x-1,y-1);
      }
      
      if( (x < WNX) && (y > 0) ){ 
	w10 = float_image_get_sample(mask_img,0,x,y-1);
      }
      
      if( (x > 0) && (y < WNY) ){ 
	w01 = float_image_get_sample(mask_img,0,x-1,y);
      }
      
      if( (x < WNX) && (y < WNY) ){ 
	w11 = float_image_get_sample(mask_img,0,x,y);
      }
      
     /* fprintf(stderr,"[%d,%d] - %lf %lf %lf %lf \n",x,y,w00,w10,w01,w11);
      if((x == 121) && (y == 129)){
	fscanf(stdin,"%lf",&w00);
      }*/
       double M = (1.0/w00) + (1.0/w10) + (1.0/w11) + (1.0/w01);
       double val = 4.0/M;
       float_image_set_sample(cm_img,0,x,y,val);
     }
  }
 
   FILE* arq2 = fopen("teste_mask_orig.fni","wt");
  float_image_write(arq2,mask_img);
  fclose(arq2);
 
  FILE* arq = fopen("teste_mask.fni","wt");
  float_image_write(arq,cm_img);
  fclose(arq);
  return cm_img;
}

void generatePOVMesh(
  FILE* arq,
  float_image_t* fni_map,
  float_image_t *texture,
  float_image_t* mask_img,
  int c,
  options_t* o
  );
  

void generatePOVMesh(
  FILE* arq,
  float_image_t* fni_map,
  float_image_t *texture,
  float_image_t* mask_img,
  int c,
  options_t* o)
  {

	/*
	VERY IMPORTANT !
	Heigh maps have (nx +1)x(ny+1) values
	normal and color maps have (nx)x(xy) values
	*/
	int NC,NX,NY;
	NC = fni_map->sz[0];
	NX = fni_map->sz[1];
	NY = fni_map->sz[2];
	
	
	float_image_t* cm_img = NULL;
	if( mask_img != NULL) {
	  int WNX = mask_img->sz[1];
	  int WNY = mask_img->sz[2];
	  if( (WNX == NX) && (WNY == NY ) ){
	    cm_img = mask_img;
	  }else{
	    cm_img = scale_mask(mask_img);
	  }
	}
	
	double MASK_LIMIT =  0.05;
	
	int x, y;
	double maxHeight = DBL_MIN;
	double minHeight = DBL_MAX;
	fprintf(arq,"#declare scene_mesh = \n");
	fprintf(arq,"union{\n");
	for(x = 0; x < NX-1; x++){
		for(y = 0; y < NY-1; y++){
					
			float R,G,B;
			if(texture != NULL){
				R =  float_image_get_sample(texture, 0, x,y);
				if(texture->sz[0] == 3){
					G =  float_image_get_sample(texture, 1, x,y);
					B =  float_image_get_sample(texture, 2, x,y);
				}else{
					G = B = R;
				}
			}
			else{
				R = 0.900;
				G = 0.800;
				B = 0.500;
			}
			
                        double scale = o->scale;
			double w00,w11,w10,w01;
			w00 = w11 = w10 = w01 = 1;
			if( cm_img != NULL){
			  w00 = float_image_get_sample(cm_img, c, x,y);
			  w01 =  float_image_get_sample(cm_img, c, x,y+1);
			  w10 =  float_image_get_sample(cm_img, c, x+1,y);
			  w11 =  float_image_get_sample(cm_img, c, x+1,y+1);
			}
			
			if( fmax(fmax(w00,w10),fmax(w11,w01)) <= MASK_LIMIT ) {
			  continue;
			}
			
			double h00 = float_image_get_sample(fni_map, c, x,y)*scale;
			double h01 =  float_image_get_sample(fni_map, c, x,y+1)*scale;
			double h10 =  float_image_get_sample(fni_map, c, x+1,y)*scale;
			double h11 =  float_image_get_sample(fni_map, c, x+1,y+1)*scale;
			double hhh = (h00*w00 + h10*w10 + h01*w01 + h11*w11)/(w00+w11+w10+w01);
			
			//a pixel is made of two triangles
			fprintf(arq,"  union{\n");
			if( (w00*w01) > MASK_LIMIT ){
			  fprintf(arq,"     triangle{");
			  fprintf(arq," <%.1f,%.1f,%.5f>",x+0.0,y+0.0,h00);
			  fprintf(arq," <%.1f,%.1f,%.5f>",x+0.0,y+1.0,h01);
			  fprintf(arq," <%.1f,%.1f,%.5f>",x+0.5,y+0.5,hhh);
			  fprintf(arq,"}\n");
			}
			
			if( (w01*w11) > MASK_LIMIT ){
			  fprintf(arq,"     triangle{");
			  fprintf(arq," <%.1f,%.1f,%.5f>",x+0.0,y+1.0,h01);
			  fprintf(arq," <%.1f,%.1f,%.5f>",x+1.0,y+1.0,h11);
			  fprintf(arq," <%.1f,%.1f,%.5f>",x+0.5,y+0.5,hhh);
			  fprintf(arq,"}\n");
			}
			
			if( (w11*w10) > MASK_LIMIT){
			  fprintf(arq,"     triangle{");
			  fprintf(arq," <%.1f,%.1f,%.5f>",x+1.0,y+1.0,h11);
			  fprintf(arq," <%.1f,%.1f,%.5f>",x+1.0,y+0.0,h10);
			  fprintf(arq," <%.1f,%.1f,%.5f>",x+0.5,y+0.5,hhh);
			  fprintf(arq,"}\n");
			}

			if( (w10*w00) > MASK_LIMIT ){
			  fprintf(arq,"     triangle{");
			  fprintf(arq," <%.1f,%.1f,%.5f>",x+1.0,y+0.0,h10);
			  fprintf(arq," <%.1f,%.1f,%.5f>",x+0.0,y+0.0,h00);
			  fprintf(arq," <%.1f,%.1f,%.5f>",x+0.5,y+0.5,hhh);
			  fprintf(arq,"}\n");
			}

			if( !o->disableTexture ){
			  fprintf(arq," texture{ pigment {color rgb <%.3f,%.3f,%.3f> }",R,G,B);
			  fprintf(arq," finish  {diffuse 0.8 ambient 0.2} } ");
			  
			}
			fprintf(arq,"}\n");
			//dumb manner to find maxHeight
			double hmax = fmax(fmax(h00,h01), fmax(h10,h11));
			double hmin = fmin(fmin(h00,h01), fmin(h10,h11));
			if(maxHeight < hmax) maxHeight = hmax;
			if(minHeight > hmin) minHeight = hmin;
		}
	}
	fprintf(arq,"}\n");
	fprintf(arq,"#declare maxX = %d; \n",NX);
	fprintf(arq,"#declare maxY = %d; \n",NY);
	if (o->focusZ == 0.0) o->focusZ = (minHeight + maxHeight)/2.0;
	fprintf(arq,"#declare focusZ = %lf; \n",o->focusZ);
	fprintf(arq,"#declare maxZ = %lf; \n",maxHeight);
	fprintf(arq,"#declare minZ = %lf; \n",minHeight);
	fprintf(arq,"#declare viewAzim = %lf; \n",o->viewAzim);
	fprintf(arq,"#declare viewElev = %lf; \n",o->viewElev);
	fprintf(arq,"#declare viewDist = %lf; \n",o->viewDist);
	if(o->scaleZ == 0 ) o->scaleZ = 1.0;
	if(o->stepZ == 0) o->stepZ = o->scaleZ/20;
	fprintf(arq,"#declare stepZ = %lf; \n",o->stepZ);
	fprintf(arq,"#declare scaleZ = %lf; \n",o->scaleZ);
	

}


float_image_t* read_texture_file(char* filename);

float_image_t* read_texture_file(char* filename){
  FILE* arq = open_read(filename,TRUE);
  char magic_c;
  fscanf(arq,"%c",&magic_c);
  fclose(arq);
  if ( magic_c== 'P'){
    return float_pnm_image_read(filename,FALSE,1,0,TRUE,TRUE,FALSE);
  }else{
    return float_image_read_file(filename);
  }
}

int main(int argc, char** argv){
  options_t* o = parse_options(argc,argv);
  
  float_image_t* fni_img = float_image_read_file(o->in_file);
  float_image_t* mask_img = NULL;
  float_image_t* text_img = NULL;
  if(o->mask_file != NULL) mask_img = float_image_read_file(o->mask_file);
  if(o->tex_file != NULL) text_img = read_texture_file(o->tex_file);

  FILE* arq_out = stdout;
  if(o->out_file == NULL){
    fprintf(stderr,"Writing into standard output\n");
  }else{
    arq_out = open_write(o->out_file,TRUE);
  }
  
  if(o->viewDist == -1){
    o->viewDist = 5.0*hypot(fni_img->sz[1],fni_img->sz[2]);
  }
  fprintf(stderr,"Plotting channel %d\n",o->channel);
  generatePOVMesh(arq_out,fni_img,text_img,mask_img,o->channel,o);
  
  if(o->out_file == NULL){
	fclose(arq_out);
  }
	
  return 0;
  
}
