#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <limits.h>
#include <float.h>
#include <argparser.h>
#include <unistd.h>
#include <string.h>
#include "i2.h"
#include <r2.h>
#include <r3.h>
#include <tabela.h>
#include <jsfile.h>
#include <quad.h>
#include <delaunay.h>
#include <delaunay_plot_POV.h>



void generatePOVSceneTri(char* filename,Tabela* tab,int numLight,double scale,double offset);
double*  buildHeightVector(Tabela* tab, int numLight,double scale,double offset);
void plot_delaunay_POV (quad_arc_t e, delaunay_site_t *st, int nsites, double height[] ,char *filename,double scale,double offset);



double*  buildHeightVector(Tabela* tab, int numLight,double scale,double offset){
	int i;
	int nsites = get_num_linhas(tab);
	double* h = (double*)malloc(sizeof(double)*nsites);
	for(i = 0;i < nsites; i++){
		const double* go = get_intdir(tab,i);
		double gmag = get_intmag(tab,i);
		h[i] = (go[numLight]*gmag*scale) + offset;
	}
	return h;
}

void plot_delaunay_POV (quad_arc_t e, delaunay_site_t *st, int nsites, double height[] ,char *filename,double scale,double offset)
  { 
      
    /* Create Postscript document or EPS figure stream. */
    FILE* wr = open_write(filename,TRUE);
    float minX,minY,maxX,maxY;
    minX = minY = -1;
    maxX = maxY = +1;
    float maxZ = offset;
    float minZ = offset;
    int i;
    for(i = 0; i < nsites; i++){
	maxZ = fmax(maxZ,height[i]);
	minZ = fmin(minZ,height[i]);
    }
    fprintf(wr,"#declare minX = %f;\n",minX);
    fprintf(wr,"#declare maxX = %f;\n",maxX);
    fprintf(wr,"#declare minY = %f;\n",minY);
    fprintf(wr,"#declare maxY = %f;\n",maxY);
    fprintf(wr,"#declare maxZ = %f;\n",maxZ);
    fprintf(wr,"#declare minZ = %f;\n",minZ);
    fprintf(wr,"#declare zeroLevel = %f;\n",offset);
    fprintf(wr,"#declare graph = \n");
    fprintf(wr,"  mesh{\n");
    delaunay_plot_POV_triangles(wr, e,height,"graph_tx_out");
    fprintf(wr,"  }\n");
    fprintf(wr,"#declare skirt = \n");
    fprintf(wr,"  mesh{\n");
    delaunay_plot_POV_skirt(wr,e,height,"skirt_tx_out");
    fprintf(wr,"  }\n");
    fclose(wr);
  }




void generatePOVSceneTri(char* filename,Tabela* tab,int numLight,double scale,double offset){
	//the main idea is use an workspace with 2 decimal digits of precision, so we multiply all value by escala
	FILE* arq;
	arq = open_write(filename,TRUE);
	
	r3_t view_dir = get_view_dir(tab);
        r3x3_t roda_normal = compute_normal_correction_matrix(view_dir);

	int nsites = get_num_linhas(tab);
	delaunay_site_t *st = notnull(malloc(nsites*sizeof(delaunay_site_t)), "no mem");
	int i;
	for (i = 0; i < nsites; i++) {
	    r3_t normal_torta = get_normal(tab,i);
	    r3_t normal_bruta;
            r3x3_map_col(&roda_normal,&(normal_torta),&(normal_bruta));
	    st[i].index = i;
      	    st[i].p.c[0] = normal_bruta.c[0];
            st[i].p.c[1] = normal_bruta.c[1];
        }
	quad_arc_t e = delaunay_build (st, nsites);
	double* height = buildHeightVector(tab,numLight,scale,offset);
	plot_delaunay_POV(e, st, nsites, height, filename,scale,offset);
	free(height);
	fclose(arq);
	free(st);
}

int main(int argc, char** argv){
	if(argc < 5){
		fprintf(stderr,"Program Usage: program <input_file> <first_light> <last_light> <output_prefix> [scale offset]\n");
		return 0;
	} 
	char* tab_filename = argv[1];
	int first_light = atoi(argv[2]);
	int last_light = atoi(argv[3]);
	char* out_prefix = argv[4];
	double scale = 1.0;
	double offset = 0.0;
	if(argc == 6){
		fprintf(stderr,"Program Usage: program <input_file> <first_light> <last_light> <output_prefix> [scale offset]\n");
		return 0;
	}
	if(argc >= 7){
		
		scale=atof(argv[5]);
		offset=atof(argv[6]);
	}
	Tabela* tab = NULL;
	
	LoadTable(tab_filename,&tab);
	if(tab == NULL){
		fprintf(stderr,"ERROR - FILE \"%s\" CANT BE OPENED !\n",argv[2]);
		return 1;
	}
	int l;
	for(l = first_light; l <= last_light;l++){
		char* outputFileName = NULL;
		char *outputFileName = jsprintf("%s_%d.inc",out_prefix,l);
		generatePOVSceneTri(outputFileName,tab,l,scale,offset);
	}
	return 0;
	
}
