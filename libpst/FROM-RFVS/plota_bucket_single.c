#include <stdio.h>
#include "hash.h"
#include "tabela.h"
#include <jsfile.h>
#include <rn.h>
#include <assert.h>

void generateDetals(FILE* arq, bucketGrid* bg);
void generateDetals(FILE* arq, bucketGrid* bg){
	
	// The whole scene ranges from -1 to +1 in X and Y,
	// from about 0 to {max_value} in Z.  The grid proper ranges from {-grid_size*cell_size} to {+grid_size*cell_size}
        // in X and Y; each cell is {cell_size} wide. Other scaling is left to the POV program.
        // The scene is divided into quarters, {grid_sticks_{I}{J}},
        // {grid_balls_{I}{J}}, and {grid_plate_{I}{J}}, for {i,J} in {0,1}.
 	double cell_size = 2.0*(bg->R/bg->tam_grid);
	int grid_size = bg->tam_grid;
        fprintf(arq,"// Scene description\n");
 	fprintf(arq,"#declare grid_size = %d ;\n",grid_size);
	int is,js,ir,jr;
// 	double stick_radius = 0.3*cell_size;
	int Hsize=grid_size/2;
	double hN = ((double)grid_size/2.0);
	for(is = 0; is < 2; is++){
		for(js = 0; js < 2; js++){
			int Imin = is*Hsize;
			int Jmin = js*Hsize;
			int Isize = (is == 0 ? Hsize:grid_size - Hsize);
			int Jsize = (js == 0 ? Hsize:grid_size - Hsize);
			double Zcellmin = 0.0005;                       
			double Zcellmax = 0.002;                       
			fprintf(arq,"#declare grid_cells_%d%d = \n",is,js);
			if (grid_size < 80){
				fprintf(arq,"union{\n");
				for(ir = 0; ir < Isize;ir++){
					for(jr = 0; jr < Jsize;jr++){
						int i = Imin + ir;
						int j = Jmin + jr;
						double Xlo = (i-hN+0.1)*cell_size;
						double Ylo = (j-hN+0.1)*cell_size;
						double Xhi = (i-hN+0.9)*cell_size;
						double Yhi = (j-hN+0.9)*cell_size;
						fprintf(arq,"  box{\n");
						fprintf(arq,"    <%lf,%lf,%lf>,",Xlo,Ylo,Zcellmin);
						fprintf(arq,"    <%lf,%lf,%lf>\n",Xhi,Yhi,Zcellmax);
						fprintf(arq,"  }\n");
					}
				}
				fprintf(arq,"}\n");
			}else{
				double Xgridmin = (Imin-hN)*cell_size;
				double Xgridmax = (Imin+Isize-hN)*cell_size;
				double Ygridmin = (Jmin-hN)*cell_size;
				double Ygridmax = (Jmin+Isize-hN)*cell_size; 
				fprintf(arq,"  box{\n");
				fprintf(arq,"    <%lf,%lf,%lf>,\n",Xgridmin,Ygridmin,Zcellmin);
				fprintf(arq,"    <%lf,%lf,%lf>\n",Xgridmax,Ygridmax,Zcellmax);
				fprintf(arq,"  }\n");
			}
			//draw a thick plate
			double Xplatemin = (is == 0 ? -1 : (Hsize-hN)*cell_size);
			double Xplatemax = (is == 0 ? (Hsize-hN)*cell_size : +1);
			double Yplatemin = (js == 0 ? -1 : (Hsize-hN)*cell_size);
			double Yplatemax = (js == 0 ? (Hsize-hN)*cell_size : +1); 
			double Zplatemin = 0;
			double Zplatemax = 0.001;                       
			fprintf(arq,"#declare grid_plate_%d%d = \n",is,js);
			fprintf(arq,"  box{\n");
			fprintf(arq,"    <%lf,%lf,%lf>,\n",Xplatemin,Yplatemin,Zplatemin);
			fprintf(arq,"    <%lf,%lf,%lf>\n",Xplatemax,Yplatemax,Zplatemax);
			fprintf(arq,"  }\n");
			fprintf(arq,"\n");
		}
	}
}

void generatePOVDots(FILE* arq, bucketGrid* bg, Tabela* tab,char* tag);
void generatePOVDots(FILE* arq, bucketGrid* bg, Tabela* tab,char* tag){
	
	// The whole scene ranges from -1 to +1 in X and Y,
	// from about 0 to {max_value} in Z.  The grid proper ranges from {-grid_size*cell_size} to {+grid_size*cell_size}
        // in X and Y; each cell is {cell_size} wide. Other scaling is left to the POV program.
        // The scene is divided into quarters, {grid_sticks_{I}{J}},
        // {grid_balls_{I}{J}}, and {grid_plate_{I}{J}}, for {i,J} in {0,1}.
	double cell_size = 2.0*bg->R/(double)bg->tam_grid;
	double* u = bg->u;
	double* v = bg->v;
	double bu = bg->bu;
	double bv = bg->bv;
	double* bar = bg->baricentro;
	int n_lights = bg->num_luzes;
	int is,js;
// 	double stick_radius = 0.3*cell_size;
	for(is = 0; is < 2; is++){
		for(js = 0; js < 2; js++){
			fprintf(arq,"#declare grid_sticks_%d%d = union{}\n",is,js);
			fprintf(arq,"#declare %s_%d%d = \n",tag,is,js);
			fprintf(arq,"union{\n");
			int linha;
			int NL = get_num_linhas(tab);
			for(linha = 0; linha < NL; linha++){
					  const double* gsig = get_intdir(tab,linha);
					  double sig[n_lights];
					  int ll;
					  for (ll = 0; ll < n_lights; ll++){
					    sig[ll] = gsig[ll];
					  }
					  double sigu = (rn_dot(n_lights,sig,u) - bu);
					  double sigv = (rn_dot(n_lights,sig,v) - bv);
					  int l;
					  double sum = 0;
					  for(l = 0; l < n_lights; l++){
					    double pi = bar[l] + sigu*u[l] + sigv*v[l];
					    double gi = sig[l];
					    sum+=(pi-gi)*(pi-gi);
					  }
					  double height  = sqrt(sum);
					  double radius = 0.5*cell_size;
					  fprintf(arq,"  sphere{\n");
					  fprintf(arq,"    <%lf,%lf,%lf>,",sigu,sigv,height);
					  fprintf(arq,"    %lf\n",radius);
					  fprintf(arq,"  }\n");
					  
				}
				
			
			
			fprintf(arq,"}\n");
			
		}
	}
}
// void generatePOVDots(FILE* arq, bucketGrid* bg, Tabela* tab,char* tag){
// 	
// 	// The whole scene ranges from -1 to +1 in X and Y,
// 	// from about 0 to {max_value} in Z.  The grid proper ranges from {-grid_size*cell_size} to {+grid_size*cell_size}
//         // in X and Y; each cell is {cell_size} wide. Other scaling is left to the POV program.
//         // The scene is divided into quarters, {grid_sticks_{I}{J}},
//         // {grid_balls_{I}{J}}, and {grid_plate_{I}{J}}, for {i,J} in {0,1}.
// 	double cell_size = 2.0*bg->R/(double)bg->tam_grid;
// 	double* u = bg->u;
// 	double* v = bg->v;
// 	double bu = bg->bu;
// 	double bv = bg->bv;
// 	double* bar = bg->baricentro;
// 	int n_lights = bg->num_luzes;
// 	int grid_size = bg->tam_grid;
// 	int is,js,ir,jr;
// // 	double stick_radius = 0.3*cell_size;
// 	int Hsize=grid_size/2;
// 	for(is = 0; is < 2; is++){
// 		for(js = 0; js < 2; js++){
// 			int Imin = is*Hsize;
// 			int Jmin = js*Hsize;
// 			int Isize = (is == 0 ? Hsize:grid_size - Hsize);
// 			int Jsize = (js == 0 ? Hsize:grid_size - Hsize);
// 			fprintf(arq,"#declare grid_sticks_%d%d = union{}\n",is,js);
// 			fprintf(arq,"#declare %s_%d%d = \n",tag,is,js);
// 			fprintf(arq,"union{\n");
// 			for(ir = 0; ir < Isize;ir++){
// 				for(jr = 0; jr < Jsize;jr++){
// 					int i = Imin + ir;
// 					int j = Jmin + jr;
// 					
// 					listapont* b = bg->buckets[i][j].itens;
// 					while(b != NULL){
// 					  int linha = b->ponteiro;
// 					  const double* gsig = get_intdir(tab,linha);
// 					  double sig[n_lights];
// 					  int ll;
// 					  for (ll = 0; ll < n_lights; ll++){
// 					    sig[ll] = gsig[ll];
// 					  }
// 					  double sigu = (rn_dot(n_lights,sig,u) - bu);
// 					  double sigv = (rn_dot(n_lights,sig,v) - bv);
// 					  int l;
// 					  double sum = 0;
// 					  for(l = 0; l < n_lights; l++){
// 					    double pi = bar[l] + sigu*u[l] + sigv*v[l];
// 					    double gi = sig[l];
// 					    sum+=(pi-gi)*(pi-gi);
// 					  }
// 					  double height  = sqrt(sum);
// 					  double radius = 0.5*cell_size;
// 					  fprintf(arq,"  sphere{\n");
// 					  fprintf(arq,"    <%lf,%lf,%lf>,",sigu,sigv,height);
// 					  fprintf(arq,"    %lf\n",radius);
// 					  fprintf(arq,"  }\n");
// 					  
// 					  b = b->proximo;
// 					}
// 				}
// 			}
// 			
// 			fprintf(arq,"}\n");
// 			
// 		}
// 	}
// }


int main(int argc, char** argv){
  char* tabname;
  char* tabname2 =  NULL;
  char* povname;
  int gridsize;
  
  if(argc <  4 ){
    fprintf(stderr,"program usage <tabfile> <povfile> <gridsize> <table2>\n");
    return 1;
  }
  
  tabname = argv[1];
  povname = argv[2];
  sscanf(argv[3],"%d",&gridsize);
  
  if(argc == 5){
    tabname2 = argv[4];
  }
  
  Tabela* tab;
  LoadTable(tabname,&tab);
  assert(tab != NULL);
  bucketGrid* bg = CriaBucketGrid(tab,gridsize);
  FILE* arq_pov = open_write(povname, TRUE);
  generateDetals(arq_pov, bg);
  generatePOVDots(arq_pov,bg,tab,"grid_balls");
  if( tabname2 != NULL){
    Tabela* tab2;
    LoadTable(tabname2,&tab2);
    generatePOVDots(arq_pov,bg,tab2,"grid_balls2");
  }

  fclose(arq_pov);
  return 0;
}