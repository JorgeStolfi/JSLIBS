#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <time.h>
#include <assert.h>
#include <float_image.h>
#include <float_pnm_image_io.h>

#define IS_MAX 0
#define IS_MIN 1
#define IS_AVG 2
#define IS_VAR 3

float_image_t* MaxImages(float_image_t** image_list, int N,int mode){
	float_image_t* novo;
	int nc,nx,ny;
	nc = image_list[0]->sz[0];
	nx = image_list[0]->sz[1];
	ny = image_list[0]->sz[2];
	
	int i;
	for(i = 0; i < N; i++){
		if( nc > image_list[i]->sz[0]){ nc = image_list[i]->sz[0];}
		if( nx > image_list[i]->sz[1]){ nx = image_list[i]->sz[1];}
		if( ny > image_list[i]->sz[2]){ ny = image_list[i]->sz[2];}
	}

	novo = float_image_new(nc,nx,ny);
	
	int c,x,y;
	for(c = 0; c < nc; c++){
		for(x = 0; x < nx; x++){
			for( y = 0; y < ny; y++){
				float value = float_image_get_sample(image_list[0],c,x,y);
				double val = 0;
				double vec[N];
				for(i = 0; i < N; i++){
					
					float v2 = float_image_get_sample(image_list[i],c,x,y);
					vec[i] = v2;
					if(mode == IS_MAX){
						if (value < v2){
							value = v2;
						}
					}else if( mode == IS_MIN){
						if (value > v2){
							value = v2;
						}
					}else if((mode == IS_AVG) || (mode == IS_VAR)){
						val+= v2;
					}
				}
				if(mode == IS_AVG){
				  value = val/(double)N;
				}
				if( mode == IS_VAR){
				  double avg = val/(double)N;
				  double var = 0;
				  for(i = 0; i < N; i++){
				    var+= ((vec[i] - avg)*(vec[i] - avg))/(double)N;
				  }
				  value = var;
				}
				float_image_set_sample(novo,c,x,y,value);

			}	
		}
	}
	
	return novo;	

}


int main(int argc, char** argv){
	if(argc < 4){
		fprintf(stderr,"Chamada do programa:\nprograma {MAX|MIN|AVG|VAR} <outfile> <file1,file2,file3...>\n");
		return 1;
	}
	char* text_mode = argv[1];
	int mode = IS_MAX;
	if(strcmp("MAX",text_mode) == 0){
		mode = IS_MAX;
	}else if( strcmp("MIN",text_mode) == 0){
		mode = IS_MIN;
	}else if( strcmp("AVG",text_mode) == 0){
		mode = IS_AVG;
	}else if( strcmp("VAR",text_mode) == 0){
		mode = IS_VAR;
	}else{
	  fprintf(stderr,"No mode selected !\n");
	  return 1;
	}
	char* outfilename = argv[2];
	int N = argc - 3;
	char** filename_list = argv + 3;
    	float_image_t** image_list;
	image_list = float_pnm_image_list_read(N,filename_list,FALSE,1.00,0.0,FALSE,TRUE,FALSE);
    
	float_image_t* outimage = MaxImages(image_list,N,mode);
	
	float_pnm_image_write(outfilename,outimage,FALSE,1.0,0.0,FALSE,TRUE,FALSE);
	FILE* arq;
	arq = fopen("treco.txt","wt");
		float_image_write(arq,outimage);
	fclose(arq);

	return 0;
}
