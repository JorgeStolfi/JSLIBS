#include <stdio.h>
#include <stdlib.h>
#include <jsfile.h>
#include <r3.h>
#include <float_image.h>
#include <string.h>
#include <tabela.h>
#include <float_pnm_image_io.h>


float_image_t* GenerateTableImage(Tabela* tab,int resolution,int magnitude_mode){
  /* Generates a FNI image {resolutionxresolution} where each channel contains:	
      0 - Population
      1 - Average Magnitude
      2 - N - Average Signature for each channel
  */
  int nLights = get_num_luzes(tab);
  int nLines = get_num_linhas(tab);
  float_image_t* img = float_image_new(nLights+2,resolution,resolution);
  
  int l;
  for(l = 0; l < nLines; l++){
    r3_t normal = get_normal(tab,l);
    double Gmag = get_intmag(tab,l);
    const double* go = get_intdir(tab,l);

    int x, y;
    x = ((resolution -1)/2.0)*(normal.c[0] + 1); /* from [-1,1] domain to [0,res-1]*/
    y = ((resolution -1)/2.0)*(normal.c[1] + 1);

    float population = float_image_get_sample(img,0,x,y);
    population = population+ 1;
    float_image_set_sample(img,0,x,y,population);
    
    float magnitude = float_image_get_sample(img,1,x,y);
    magnitude = magnitude+Gmag;
    float_image_set_sample(img,1,x,y,magnitude);

    int i;
    for(i = 0; i < nLights; i++){
      float sig = float_image_get_sample(img,i+2,x,y);
      
      if(magnitude_mode == 1){	
	sig = sig+(go[i]*Gmag);
      }else{
	sig = sig+go[i];
      }
      float_image_set_sample(img,i+2,x,y,sig);
    }

  }
  /* Average what should be averaged */
  int x,y;
  for(x = 0; x < resolution; x++){
    for(y = 0; y < resolution; y++){
      float population = float_image_get_sample(img,0,x,y);
      if(population != 0.0 ){
	int i;
	for(i = 1; i < nLights+ 2; i++){
	  float val = float_image_get_sample(img,i,x,y);
	  val = val/population;
	  float_image_set_sample(img,i,x,y,val);
	}
      }
    }
  }

  return img;
}


int main(int argc, char** argv){
  char* table_filename;
  char* table_img_filename;
  int resolution;
  int mag_or_sig;
  int generatePGM = 0;
  
  if(argc < 4){
    fprintf(stderr,"Program usage:\ntable_to_fni <inTable> <outPrefix> <resolution> <M|MP>\n");
    return 1;
  }

  mag_or_sig = 0; /*It means signature mode*/
  if(argc >= 5){
    if(strcmp("M",argv[4]) == 0){
      mag_or_sig = 1;
      fprintf(stderr,"Magnitude mode\n");
    }
    if(strcmp("MP",argv[4]) == 0){
      mag_or_sig =1;
      generatePGM = 1;
      fprintf(stderr,"Magnitude and PGM mode\n");
    }
  }
  table_filename = argv[1];
  char* prefix  = argv[2];
  char *table_img_filename = jsprintf("%s.fni",prefix);
  sscanf(argv[3],"%d",&resolution);

  Tabela* tab;
  LoadTable(table_filename,&tab);
  
  fprintf(stderr,"Generating FNI image from Table data with Resolution %d...",resolution);
  float_image_t* img = GenerateTableImage(tab,resolution,mag_or_sig);
  FILE* arq = open_write(table_img_filename,FALSE);
  float_image_write(arq,img);
  fclose(arq);

  if(generatePGM == 1){
     int nLights = get_num_luzes(tab);
      int i;
      int NX,NY;
      NX = resolution;
      NY = resolution;
      for(i = 0;i <  nLights; i++){
	char* img_name;

	char *img_name = jsprintf("%s_%d.pgm",prefix,i);
	float_image_t* m_img = float_image_new(1,NX,NY);
	float_image_set_channel(m_img, 0, img,2+ i);
	float_pnm_image_write(img_name,m_img,FALSE,1.0,0.0,TRUE,TRUE,FALSE);
	float_image_free(m_img);
      }
  }
  fprintf(stderr,"OK\n");
  return 0;
}

