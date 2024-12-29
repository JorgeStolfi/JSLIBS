#include <stdio.h>
#include <float_image.h>
#include <jsfile.h>
#include <assert.h>

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

float_image_t* discard_3rd_channel(float_image_t* img);

float_image_t* discard_3rd_channel(float_image_t* img){
  int NC = img->sz[0];
  int NX = img->sz[1];
  int NY = img->sz[2];

  assert(NC == 3);
  float_image_t* im = float_image_new(2,NX,NY);
  float_image_set_channel(im, 0,img,0);
  float_image_set_channel(im, 1,img,1);

  return im;
}

int main(int argc, char** argv){
  if(argc != 3){
    fprintf(stderr,"Program usage:\ndiscard_channel <in-image> <out-image>\n");
    return 1;
  }
  
  char* in_filename = argv[1];
  char* out_filename = argv[2];
  
  float_image_t* in_im = readFNI(in_filename);
  float_image_t* out_im = discard_3rd_channel(in_im);
  writeFNI(out_filename,out_im);
  
  return 0;
  
}