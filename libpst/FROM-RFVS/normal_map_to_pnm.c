#include <stdio.h>
#include <math.h>
#include <assert.h>
#include <float_image.h>
#include <float_pnm_image_io.h>

int main(int argc, char** argv){
	float_image_t * im;
	FILE* arq;
// 	long int x,y;
// 	long int i,j;
	double dx,dy,dz;
	if(argc != 3){
		fprintf(stderr, "Chamada incorreta do Programa !\n");
		fprintf(stderr, "Chamada: programa <tabela de entrada> <imagem-saida>\n");
		return 0;
	}
	fprintf(stderr, "Passou\n");
	/*fprintf(stderr, "Digite a dimensâo X e Y da imagem:");
	scanf(" %ld",&x);
	scanf(" %ld",&y);
	sscanf(argv[3],"%d",&x);
	sscanf(argv[4],"%d",&y);*/


	fprintf(stderr, "Abrindo Imagem\n");
	arq = fopen(argv[1],"rt");
	float_image_t* normal_map;
	normal_map = float_image_read(arq);
	int nx,ny;
	nx = normal_map->sz[1];
	ny = normal_map->sz[2];
	/*fscanf(arq, "tx = %ld\n", &x);
	fprintf(stderr, "X = %ld\n",x);
	fscanf(arq, "ty = %ld\n", &y);
	fprintf(stderr, "Y = %ld\n",y);
	*/im = float_image_new(3, nx, ny);
	fprintf(stderr, "Gerando imagem\n");
	//while (1) {
	int ix,iy;
	
	for(ix = 0 ; ix < nx; ix++){
		for(iy = 0; iy < ny; iy++){
			/*int nread = fscanf(arq,"%ld %ld %lf %lf %lf",&i,&j,&dx,&dy,&dz);
			//  printf("%d %d %f %f %f\n",i,j,dx,dy,dz);
			int debug = ((i == 14) && (( y - 1 - j) == 209)) ||  (((i == 12) && ((y - 1 - j) == 207)));
			if(debug) {
				fprintf(stderr,"%4ld %4ld %+19.16f %+19.16f %+19.16f\n",i,j,dx,dy,dz);
			}
			if(nread != 5 ) printf("[%ld][%ld] %lf %lf %lf nread = %d\n",i,j,dx,dy,dz,nread);
			if (nread == 0) { break; }
			assert(nread == 5);
			*/
			dx = float_image_get_sample(normal_map,0,ix,iy);
			dy = float_image_get_sample(normal_map,1,ix,iy);
			dz = float_image_get_sample(normal_map,2,ix,iy);
			float_image_set_sample(im, 0,ix,iy, (dx + 1)/2);
			float_image_set_sample(im, 1,ix,iy, (dy + 1)/2);
			float_image_set_sample(im, 2,ix,iy, dz);
		//	if((i == (x -1)) && (j == (y-1))){ break; }
		}
	}
	fprintf(stderr, "Salvando Imagem\n");
	float_pnm_image_write(argv[2],im, FALSE,VIEW_GAMMA, VIEW_BIAS,TRUE,TRUE,FALSE);
	fprintf(stderr, "Concluido !\n");
	return 0;
}
