#include <stdio.h>
#include <jsfile.h>
#include <float_image.h>

void ground_float_image(float_image_t* IZ);

void ground_float_image(float_image_t* IZ){
  long int NC = IZ->sz[0];
  long int NX = IZ->sz[1];
  long int NY = IZ->sz[2];

  float least_values[NC] ; 
  
  int c,x,y;
  for(c = 0; c < NC; c++){
    least_values[c] = float_image_get_sample(IZ,c,0,0);
    for( x = 0; x < NX; x++){
      for(y = 0; y < NY; y++){
	float v = float_image_get_sample(IZ,c,x,y);
	if( v < least_values[c]){ least_values[c] = v; }
      }
    }
  }

  for(c = 0; c < NC; c++){
    fprintf(stderr,"Channel %d : LEAST VALUE: %9.6f\n",c,least_values[c]);
    for( x = 0; x < NX; x++){
      for(y = 0; y < NY; y++){
	float v = float_image_get_sample(IZ,c,x,y);
	float_image_set_sample(IZ,c,x,y,v - least_values[c] + 1);
      }
    }
  }

  return;
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
    if(argc != 3) {
      fprintf(stderr,"Program usage\nground_float_image <infile> <outfile>\n");
      fprintf(stderr,"Normalize image, such that its least value is zero\n");
      return 1;
    }
    
    char* in_filename = argv[1];
    char* out_filename = argv[2];
    float_image_t* IZ = readFNI(in_filename);
    ground_float_image(IZ);
    writeFNI(out_filename,IZ);
    return 0;
}