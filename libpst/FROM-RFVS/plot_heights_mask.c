#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <limits.h>
#include <float.h>
#include <argparser.h>
#include <unistd.h>
#include <string.h>
#include <assert.h>
#include <float_image.h>
#include <float_pnm_image_io.h>
#include <imagem_valores.h>
#include <imagem_vetores.h>




//void generatePOVMesh(FILE* arq,imagem_valores_t* height_map,  imagem_vetores_t* normal_map,float_image_t *img, int nx,int ny);
float_image_t* float_image_read_file(char* filename);

float_image_t* float_image_read_file(char* filename){
	FILE* arq_IN = fopen(filename,"rt");
    	assert(arq_IN != NULL);
    	float_image_t *IN = float_image_read(arq_IN);
    	fclose(arq_IN);
	return IN;
}

void generatePOVMesh(FILE* arq,float_image_t* height_map,  float_image_t* normal_map,float_image_t *img,float_image_t* mask_img, int nx,int ny,double scale);

void generatePOVMesh(FILE* arq,float_image_t* height_map,  float_image_t* normal_map,float_image_t *img,float_image_t* mask_img, int nx,int ny, double scale){
	/*
	VERY IMPORTANT !
	Heigh maps have (nx +1)x(ny+1) values
	normal and color maps have (nx)x(xy) values

	TODO - Use normal map to produce smooth triangles !
	*/
	int x, y;
	double maxHeight = 0;
	fprintf(arq,"#declare scene_mesh = \n");
	fprintf(arq,"union{\n");
	for(x = 0; x < nx; x++){
		for(y = 0; y < ny; y++){
			if(mask_img != NULL){
			  double mval = float_image_get_sample(mask_img, 0, x,y);
			  if(mval < 0.5) continue ;
			}
			  
			float R,G,B;
			if(img != NULL){
				R =  float_image_get_sample(img, 0, x,y);
				if(img->sz[0] == 3){
					G =  float_image_get_sample(img, 1, x,y);
					B =  float_image_get_sample(img, 2, x,y);
				}else{
					G = B = R;
				}
			}
			else{
				R = 0.900;
				G = 0.800;
				B = 0.500;
			}
			
                        // h{X}{Y} = height at corner {(x+X,y+Y)}
// 			double h00 = height_map->pixel[y][x];
// 			double h01 = height_map->pixel[y+1][x];
// 			double h10 = height_map->pixel[y][x+1];
// 			double h11 = height_map->pixel[y+1][x+1];
			double h00 = float_image_get_sample(height_map, 0, x,y)*scale;
			double h01 =  float_image_get_sample(height_map, 0, x,y+1)*scale;
			double h10 =  float_image_get_sample(height_map, 0, x+1,y)*scale;
			double h11 =  float_image_get_sample(height_map, 0, x+1,y+1)*scale;
			double hhh = (h00 + h10 + h01 + h11)/4.0;
			
			//a pixel is made of two triangles
			fprintf(arq,"  union{\n");
			fprintf(arq,"     triangle{");
			fprintf(arq," <%f,%f,%f>",x+0.0,y+0.0,h00);
			fprintf(arq," <%f,%f,%f>",x+0.0,y+1.0,h01);
			fprintf(arq," <%f,%f,%f>",x+0.5,y+0.5,hhh);
			fprintf(arq,"}\n");
			
			fprintf(arq,"     triangle{");
			fprintf(arq," <%f,%f,%f>",x+0.0,y+1.0,h01);
			fprintf(arq," <%f,%f,%f>",x+1.0,y+1.0,h11);
			fprintf(arq," <%f,%f,%f>",x+0.5,y+0.5,hhh);
			fprintf(arq,"}\n");
			
			fprintf(arq,"     triangle{");
			fprintf(arq," <%f,%f,%f>",x+1.0,y+1.0,h11);
			fprintf(arq," <%f,%f,%f>",x+1.0,y+0.0,h10);
			fprintf(arq," <%f,%f,%f>",x+0.5,y+0.5,hhh);
			fprintf(arq,"}\n");

			fprintf(arq,"     triangle{");
			fprintf(arq," <%f,%f,%f>",x+1.0,y+0.0,h10);
			fprintf(arq," <%f,%f,%f>",x+0.0,y+0.0,h00);
			fprintf(arq," <%f,%f,%f>",x+0.5,y+0.5,hhh);
			fprintf(arq,"}\n");

			
			fprintf(arq," texture{ pigment {color rgb <%f,%f,%f> }",R,G,B);
			fprintf(arq," finish  {diffuse 0.8 ambient 0.2} } ");
			fprintf(arq,"}\n");
			
			//dumb manner to find maxHeight
			double hmax = fmax(fmax(h00,h01), fmax(h10,h11));
			if(maxHeight < hmax) maxHeight = hmax;
		}
	}
	fprintf(arq,"}\n");
	fprintf(arq,"#declare maxX = %d; \n",nx+1);
	fprintf(arq,"#declare maxY = %d; \n",ny+1);
	fprintf(arq,"#declare maxZ = %lf; \n",maxHeight);

}

int main(int argc, char** argv){
	
	if(argc < 5){
		printf("Program usage:\nplot_heights <height_map> <normal_map> <image_file> <mask_file> <output_file> [scale]\n");
		return 1;
	}
	float_image_t* height_map;
	float_image_t* normal_map = NULL;
	float_image_t* img= NULL;
	float_image_t* mask_img= NULL;
	fprintf(stderr,"Opening Height Map\n");
	//height_map = le_imagem_valores(argv[1]);
	height_map = float_image_read_file(argv[1]);
	
	//if(strcmp(argv[2],"NONE") != 0){ fprintf(stderr,"Opening Normal Map\n"); normal_map = le_imagem_vetores(argv[2]);}
	if(strcmp(argv[2],"NONE") != 0){ fprintf(stderr,"Opening Normal Map\n"); normal_map = float_image_read_file(argv[2]);}
	if(strcmp(argv[3],"NONE") != 0){ fprintf(stderr,"Opening Texture\n"); img = float_pnm_image_read(argv[3],FALSE,1,0,TRUE,TRUE,TRUE);}
	if(strcmp(argv[4],"NONE") != 0){ fprintf(stderr,"Opening Mask\n"); mask_img = float_pnm_image_read(argv[4],TRUE,1,0,TRUE,TRUE,TRUE);} 
	fprintf(stderr,"Opening Output File\n");
	FILE* output_file = fopen(argv[5],"wt");
	double scale = 1.0;
	if(argc ==6){
		scale = atof(argv[6]);
	}
	assert(output_file != NULL);
	generatePOVMesh(output_file,height_map,normal_map,img,mask_img, height_map->sz[1] -1,height_map->sz[2]-1,scale);
	fclose(output_file);
	
	return 0;
}
