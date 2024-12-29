#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <time.h>
#include <assert.h>
#include <float_image.h>
#include <float_pnm_image_io.h>
#include <float_image_paint.h> 
#include <float_image_gradient.h>

float_image_t* computeGradient(float_image_t* img, double noise);

float_image_t* computeGradient(float_image_t* img, double noise){
	float_image_t* grad_img = float_image_gradient_sqr_relative(img,noise, TRUE);
	return grad_img;
}


int main(int argc,char** argv){
	if(argc < 4){
		fprintf(stderr,"Program Usage:\npnmgrad <infile> <outfile> <noise>\n");
		return 1;
	}
	char* inFile = argv[1];
	char* outFile = argv[2];
	double noise = 0.25;
	sscanf(argv[3],"%lf",&noise);

	float_image_t * img = float_pnm_image_read(
						inFile,    /* PPM/PGM/PBM file name (with extension). */
						FALSE,
						1.0,   /* Gamma to use in decoding (1 = linear decoding). */
						0.0,    /* Offset to use in decoding. */
						FALSE,     /* If TRUE, reverses the indexing of rows. */
						TRUE,    /* If TRUE, prints "reading {fname}..." to {stderr}. */
						FALSE  /* If TRUE, prints conversion diagnostics to {stderr}. */
						);
	float_image_t* grad_img = computeGradient(img,noise);
	float_pnm_image_write(
				outFile,        /* PPM/PGM/PBM file name (with extension). */            
    				grad_img, /* Image to write. */
				FALSE,
    				1.0,       /* Gamma to use in encoding (1 = linear encoding). */    
    				0.0,        /* Offset to use in encoding. */                         
				FALSE,     /* If TRUE, reverses the indexing of rows. */
				TRUE,    /* If TRUE, prints "reading {fname}..." to {stderr}. */
				FALSE  /* If TRUE, prints conversion diagnostics to {stderr}. */
				);
	return 0;
}
