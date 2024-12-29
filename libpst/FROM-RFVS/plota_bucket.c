#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <limits.h>
#include <float.h>
#include <argparser.h>
#include <unistd.h>
#include <string.h>



#define PROG_NAME "plota_bucket"

#define PROG_HELP \
  PROG_NAME " \\\n" \
  "    -filename {INPUT_FILE}  \\\n" \
  "    -prefix {FILE_PREFIX} \\\n" \
  "    -imageDimensions {WIDTH} {HEIGTH} {DEPTH}\\\n" \
  "    -plotType {GNUPLOT|POVRAY} \\\n" \
  "    -dataType {ENTRIES|RADIUS|FIND|SCANNED|EVAL|DISTRAD|SCANNED_H|} \\\n" \
  "    -HWRatio {HWRATIO} \\\n" \
  "    -cameraInclination {INCLINATION} \\\n" \
  "    -cameraAzimuth {AZIMUTH}  \\\n" \
  "    -cameraDist {DIST}  \\\n" \
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


typedef enum {GNUPLOT = 0, POVRAY = 1} opcao_plot_type_t;
typedef enum { ENTRIES = 0, RADIUS = 1, FIND = 2 , SCANNED =3 , EVAL =4,
	       MI_DIST = 5, DISTRAD = 6,SCANNED_H =7, EVAL_H = 8,HASHED = 9,
	       EVAL_S_PER_SCAN_S = 10, EVAL_H_PER_HASH =11, SCAN_H_PER_HASH = 12, FIND_PER_SCAN_S = 13} opcao_data_t;
char* plot_type_desc[] = {"GNUPLOT","POVRAY"};
char* data_type_desc[] = {"ENTRIES","RADIUS","FIND","SCANNED","EVAL","MI_DIST","DISTRAD","SCANNED_H","EVAL_H","HASHED","EVAL_S_PER_SCAN_S",
			  "EVAL_H_PER_HASH","SCAN_H_PER_HASH", "FIND_PER_SCAN_S" };


struct bucket_grid_data_t{
	int N;
	int nLights;
	double bu,bv;
	double R;
	double bucketSize;
	double* baricenter;
	double* u;
	double* v;
	int** entries; //Number of table entries of the BG
	double** radius; // Radius of an bucket of the BG
	double*** centro; //vector with the distance of the center of bucket minus mi
	double** umin;
	double** umax;
	double** vmin;
	double** vmax;
	int** find; //how much a normals was mapped to an bucket
	int** scanned_h; //how much a normals was mapped to an bucket
	int** scanned; //how much times a bucket was scanned
	int** eval; //number of euclidean distances executed in a bucket
	int** eval_h; //number of euclidean distances executed in a bucket
	int** hashed;
};

typedef struct bucket_grid_data_t bucket_grid_data;

struct plot_config_t{
	//data related to image output
	int h,w;
	int depth;
	char format[4];
	opcao_plot_type_t plot_type;
	//data related to the data to be ploted
	opcao_data_t data_type;
	//data relative to scene perspective
	double HWRatio; //ratio Height/Width of the grid to be projected
	double camera_inclination;
	double camera_azimuth;
	double camera_dist;
	
};

typedef struct plot_config_t plot_config;

typedef double grid_func_t( int x, int y);

void advanceUntilEqual(FILE* arq,char ch);
bucket_grid_data* loadBucketGridData(char* filename);
void generatePOVHistogram(FILE* arq,grid_func_t* func ,int grid_size,double max_value);
void generatePOVPopsicles(FILE* arq,grid_func_t* h_func,grid_func_t* r_func,int grid_size,double cell_size,double max_value);
void generatePOVCameraConfig(FILE* arq,plot_config* pg,int maxValue,int grid_size);
void generatePOVFiles(char* prefix ,grid_func_t* h_func,grid_func_t* r_func,int grid_size,double cell_size,plot_config* pg);
void plotScene(char* prefix, plot_config* pg,bucket_grid_data* bg_data);
void getConfig(char* prefix,char* filename, plot_config* pg, int argc, char** argv);


void advanceUntilEqual(FILE* arq,char ch){
	char c;
	c = fgetc(arq);
	while( c != EOF){
		if(c == ch) break;
		c = fgetc(arq);
	}
}


bucket_grid_data* loadBucketGridData(char* filename){
	FILE* arq = fopen(filename,"rt");
	
	if(arq == NULL){
		return NULL;
	}
	bucket_grid_data* bg = (bucket_grid_data*)malloc(sizeof(bucket_grid_data));
	advanceUntilEqual(arq,'=');
	fscanf(arq,"%d",&(bg->N));
	advanceUntilEqual(arq,'=');
	fscanf(arq,"%d",&(bg->nLights));
	advanceUntilEqual(arq,'=');
	fscanf(arq,"%lf %lf",&(bg->bu),&(bg->bv));
	advanceUntilEqual(arq,'=');
	fscanf(arq,"%lf",&(bg->R));
	advanceUntilEqual(arq,'=');
	fscanf(arq,"%lf",&(bg->bucketSize));
	advanceUntilEqual(arq,'=');
	bg->baricenter = (double*)malloc(sizeof(double)*(bg->nLights));
	int i,j;
	for(i = 0; i < bg->nLights; i++){
		fscanf(arq,"%lf",&(bg->baricenter[i]));
	}
	advanceUntilEqual(arq,'=');
	bg->u = (double*)malloc(sizeof(double)*(bg->nLights));
	for(i = 0; i < bg->nLights; i++){
		fscanf(arq,"%lf",&(bg->u[i]));
	}
	advanceUntilEqual(arq,'=');
	bg->v = (double*)malloc(sizeof(double)*(bg->nLights));
	for(i = 0; i < bg->nLights; i++){
		fscanf(arq,"%lf",&(bg->v[i]));
	}
	advanceUntilEqual(arq,'#');
	advanceUntilEqual(arq,'\n');
	
	bg->entries = (int**)malloc(sizeof(int*)*(bg->N));
	bg->scanned_h = (int**)malloc(sizeof(int*)*(bg->N));
	bg->find = (int**)malloc(sizeof(int*)*(bg->N));
	bg->hashed = (int**)malloc(sizeof(int*)*(bg->N));
	bg->scanned = (int**)malloc(sizeof(int*)*(bg->N));
	bg->eval = (int**)malloc(sizeof(int*)*(bg->N));
	bg->eval_h = (int**)malloc(sizeof(int*)*(bg->N));
	bg->radius = (double**)malloc(sizeof(double*)*(bg->N));
	bg->centro = (double***)malloc(sizeof(double**)*(bg->N));
	bg->umin = (double**)malloc(sizeof(double*)*(bg->N));
	bg->umax = (double**)malloc(sizeof(double*)*(bg->N));
	bg->vmin = (double**)malloc(sizeof(double*)*(bg->N));
	bg->vmax = (double**)malloc(sizeof(double*)*(bg->N));

	for( i = 0; i < bg->N;i++){
		bg->entries[i] = (int*)malloc(sizeof(int)*(bg->N));
		bg->find[i] = (int*)malloc(sizeof(int)*(bg->N));
		bg->scanned_h[i] = (int*)malloc(sizeof(int)*(bg->N));
		bg->hashed[i] = (int*)malloc(sizeof(int)*(bg->N));
		bg->scanned[i] = (int*)malloc(sizeof(int)*(bg->N));
		bg->eval[i] = (int*)malloc(sizeof(int)*(bg->N));
		bg->eval_h[i] = (int*)malloc(sizeof(int)*(bg->N));
		bg->radius[i] =  (double*)malloc(sizeof(double)*(bg->N));
		bg->centro[i] =  (double**)malloc(sizeof(double*)*(bg->N));
		bg->umin[i] =  (double*)malloc(sizeof(double)*(bg->N));
		bg->umax[i] =  (double*)malloc(sizeof(double)*(bg->N));
		bg->vmin[i] =  (double*)malloc(sizeof(double)*(bg->N));
		bg->vmax[i] =  (double*)malloc(sizeof(double)*(bg->N));
		for(j = 0; j < bg->N; j++){
			bg->centro[i][j] =  (double*)malloc(sizeof(double)*(bg->nLights));
		}
	}
	int count = 0;
	int total = (bg->N)*(bg->N);
	while( (!feof(arq))   && (count < total)){
		fscanf(arq,"%d %d", &i,&j);
		//printf("I&J: %d %d\n",i,j);
	
		double umin,umax,vmin,vmax;
		fscanf(arq,"%lf %lf %lf %lf",&umin,&vmin,&umax,&vmax);
		
		bg->umin[i][j] = umin;
		bg->umax[i][j] = umax;
		bg->vmin[i][j] = vmin;
		bg->vmax[i][j] = vmax;
		fscanf(arq,"%d",&(bg->entries[i][j]));
		//printf("ENTRIES: %d ",bg->entries[i][j]);
		
		fscanf(arq,"%lf",&(bg->radius[i][j]));
		//printf("RADIUS: %lf ",bg->radius[i][j]);

		int k;
		//printf("MI_DIST: ");
		for(k = 0; k < bg->nLights;k++){
			fscanf(arq,"%lf",&(bg->centro[i][j][k]));
		//	printf("%lf ",bg->centro[i][j][k]);
		}
		//printf("\n");
		
		fscanf(arq,"%d",&(bg->find[i][j]));
		//printf("HASHED: %d ",bg->hashed[i][j]);
		fscanf(arq,"%d",&(bg->scanned[i][j]));
		//printf("SCANNED: %d ",bg->scanned[i][j]);
		fscanf(arq,"%d",&(bg->eval[i][j]));
		//printf("COMPR: %d\n ",bg->compr[i][j]);
		fscanf(arq,"%d",&(bg->scanned_h[i][j]));
		fscanf(arq,"%d",&(bg->eval_h[i][j]));
		fscanf(arq,"%d",&(bg->hashed[i][j]));
		//char c = getchar();
		count++;
	}
	fclose(arq);
	
	return bg;
}


void generatePOVHistogram(FILE* arq,grid_func_t* func,int grid_size,double max_value){
	// The whole scene ranges from {-grid_size} to {+grid_size} in X and Y,
	// from about 0 to {max_value} in Z.  Each cell is 1 unit wide in X and Y.
        // Scaling is left to the POV program.
        // The scene is divided into quarters, {grid_sticks_{I}{J}} and
        // {grid_plate_{I}{J}}, for {i,J} in {0,1}.
	fprintf(arq,"// Scene description\n");
 	fprintf(arq,"#declare grid_size = %d ;\n",grid_size);
	fprintf(arq,"#declare max_value = %f ;\n",max_value);
	int is,js,ir,jr;
	int Hsize = grid_size/2;
	for(is = 0; is < 2; is++){
		for(js = 0; js < 2; js++){
			int Imin = is*Hsize;
			int Jmin = js*Hsize;
			int Isize = (is == 0 ? Hsize:grid_size - Hsize);
			int Jsize = (js == 0 ? Hsize:grid_size - Hsize);
			fprintf(arq,"#declare grid_sticks_%d%d = \n",is,js);
			fprintf(arq,"union{\n");
			for(ir = 0; ir < Isize;ir++){
				for(jr = 0; jr < Jsize;jr++){
					int i = Imin + ir;
					int j = Jmin + jr;
					double value = func(i,j);
					if(isfinite(value) && (value != 0)){
						fprintf(arq,"  box{\n");
						fprintf(arq,"    <%lf,%lf,%lf>,",(double)i,(double)j,0.0);
						fprintf(arq,"    <%lf,%lf,%lf>\n",(double)i+1,(double)j+1,func(i,j));
						fprintf(arq,"  }\n");
					}
				}
			}
			//draw a thick plate
			fprintf(arq,"}\n");

			fprintf(arq,"#declare grid_plate_%d%d = \n",is,js);
			fprintf(arq,"  box{ <%d,%d,0> <%d,%d,%lf> }\n",Imin,Jmin,Imin + Isize,Jmin + Jsize,0.001*max_value);
			fprintf(arq,"\n");
	
		}
		
	}
}

void generatePOVPopsicles(FILE* arq,grid_func_t* h_func,grid_func_t* r_func,int grid_size,double cell_size,double max_value){
	
	// The whole scene ranges from -1 to +1 in X and Y,
	// from about 0 to {max_value} in Z.  The grid proper ranges from {-grid_size*cell_size} to {+grid_size*cell_size}
        // in X and Y; each cell is {cell_size} wide. Other scaling is left to the POV program.
        // The scene is divided into quarters, {grid_sticks_{I}{J}},
        // {grid_balls_{I}{J}}, and {grid_plate_{I}{J}}, for {i,J} in {0,1}.
	
        fprintf(arq,"// Scene description\n");
 	fprintf(arq,"#declare grid_size = %d ;\n",grid_size);
	fprintf(arq,"#declare max_value = %f ;\n",max_value);
	int is,js,ir,jr;
	double stick_radius = 0.3*cell_size;
	int Hsize=grid_size/2;
	double hN = ((double)grid_size/2.0);
	for(is = 0; is < 2; is++){
		for(js = 0; js < 2; js++){
			int Imin = is*Hsize;
			int Jmin = js*Hsize;
			int Isize = (is == 0 ? Hsize:grid_size - Hsize);
			int Jsize = (js == 0 ? Hsize:grid_size - Hsize);
                        fprintf(arq,"#declare grid_sticks_%d%d = \n",is,js);
			fprintf(arq,"union{\n");
			for(ir = 0; ir < Isize;ir++){
				for(jr = 0; jr < Jsize;jr++){
					int i = Imin + ir;
					int j = Jmin + jr;
					double Xct = (i-hN+0.5)*cell_size;
					double Yct = (j-hN+0.5)*cell_size;
					double height = h_func(i,j);
					double radius = r_func(i,j);
					if(isfinite(height) && (isfinite(radius)) && (radius >= 0)){
						if( (height > radius)){
							fprintf(arq,"  cylinder{\n");
							fprintf(arq,"    <%lf,%lf,%lf>,",Xct,Yct,0.0);
							fprintf(arq,"    <%lf,%lf,%lf>,",Xct,Yct,height);
							fprintf(arq,"  %lf\n",stick_radius);
							fprintf(arq,"  }\n");
						}
						
					}
				}
			}
			fprintf(arq,"}\n");
			fprintf(arq,"#declare grid_balls_%d%d = \n",is,js);
			fprintf(arq,"union{\n");
			for(ir = 0; ir < Isize;ir++){
				for(jr = 0; jr < Jsize;jr++){
					int i = Imin + ir;
					int j = Jmin + jr;
					double Xct = (i-hN+0.5)*cell_size;
					double Yct = (j-hN+0.5)*cell_size;
					double height = h_func(i,j);
					double radius = r_func(i,j);
					if(isfinite(height) && (isfinite(radius)) && (radius >= 0)){
						fprintf(arq,"  sphere{\n");
						fprintf(arq,"    <%lf,%lf,%lf>,",Xct,Yct,height);
						fprintf(arq,"    %lf\n",fmax(radius,stick_radius));
						fprintf(arq,"  }\n");
					}
				}
			}
			
			fprintf(arq,"}\n");
			//draw the cells
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

void generatePOVCameraConfig(FILE* arq,plot_config* pg,int maxValue,int grid_size){
	double azim,incl;
	azim = (pg->camera_azimuth)*(M_PI)/180.0;
	incl = (pg->camera_inclination)*(M_PI)/180.0;
	//double camera_dist = 4.5*maxValue;
	/*fprintf(arq,"// Camera configs \n");
	fprintf(arq,"#declare image_rel_width = %d\n",pg->w);
	fprintf(arq,"#declare image_rel_height = %d\n",pg->h);
	fprintf(arq,"#declare camera_focus = < 0,0,0>\n");
	fprintf(arq,"#declare camera_disp = %f * <%f, %f, %f>;\n",camera_dist,camera_dist*cos(azim),camera_dist*cos(incl),camera_dist*sin(incl));
	fprintf(arq,"\n");*/
	double h = pg->camera_dist*grid_size;
	

	fprintf(arq,"//Camera Params\n");
	fprintf(arq,"camera {\n");
  	fprintf(arq,"\tlocation <%f, %f, %f>\n",cos(azim)*h ,sin(azim)*h  ,sin(incl)*h );
  	fprintf(arq,"\tup <0, 0, 1>\n");
  	fprintf(arq,"\tsky <0, 0, 1>\n");
  	fprintf(arq,"\tlook_at <0, 0, 0>\n");
	fprintf(arq,"}\n");
}

void generatePOVFiles(char* prefix , grid_func_t* h_func ,grid_func_t* r_func ,int grid_size,double cell_size,plot_config* pg){
	
	int i,j;
	double max_value = DBL_MIN;
	int maxi,maxj;
	
	for(i = 0; i < grid_size;i++){
		for(j = 0; j < grid_size ;j++){
			double value = h_func(i,j) + r_func(i,j);
			if(max_value < value ){
				max_value = value;
				maxi = i;
				maxj = j;
			}
		}
	}
	fprintf(stderr,"Max Value[%d][%d] %f\n",maxi,maxj,max_value);
	/*FILE* arq_camera;
	char* cameraFileName = NULL;
	char *cameraFileName = jsprintf("%s/camera_params.inc",prefix);
	arq_camera = fopen(cameraFileName,"wt");
	if(arq_camera == NULL){
		fprintf(stderr,"ERROR - FILE \"%s\" CANT BE OPENED !\n",cameraFileName);
		return ;
	}
	//generatePOVCameraConfig(arq_camera,pg,max_value,grid_size);
	fclose(arq_camera);*/
	
	
	
	FILE* arq_scene;
	char* sceneFileName = NULL;
	char *sceneFileName = jsprintf("%s/boxes_%s.inc",prefix,data_type_desc[pg->data_type]);
	arq_scene = fopen(sceneFileName,"wt");
	if(arq_scene == NULL){
		fprintf(stderr,"ERROR - FILE \"%s\" CANT BE OPENED !\n",sceneFileName);
		return ;
	}
	if( pg->data_type != DISTRAD ) generatePOVHistogram(arq_scene,h_func,grid_size,max_value);
	else generatePOVPopsicles(arq_scene,h_func,r_func,grid_size,cell_size,max_value);
	fclose(arq_scene);
	

}

void plotScene(char* prefix, plot_config* pg,bucket_grid_data* bg_data){
	
	grid_func_t* h_func;
	grid_func_t* r_func;
	auto double func_Entries(int x, int y);
	auto double func_Radius(int x, int y);
	auto double func_Find(int x, int y);
	auto double func_Scanned_h(int x, int y);
	auto double func_Hashed(int x, int y);
	//auto double func_BH(int x, int y);
	//auto double func_EH(int x, int y);
	auto double func_Scanned(int x, int y);
	auto double func_Eval(int x, int y);
	auto double func_Eval_h(int x, int y);
	auto double func_MiDist(int x, int y);
	auto double func_Zero(int x, int y);
	auto double func_Scan_h_per_Hash(int x, int y);
	auto double func_Eval_s_per_Scan_s(int x, int y);
	auto double func_Eval_h_per_Hash(int x, int y);
	auto double func_Find_per_Scan_s(int x, int y);

	double func_Zero(int x, int y){
		return 0;
	}

	double func_Entries(int x, int y){
		return bg_data->entries[x][y];
	}
	
	double func_Radius(int x, int y){
		return bg_data->radius[x][y];
	}

	double func_Find(int x, int y){
		return bg_data->find[x][y];
	}

	double func_Scanned_h(int x, int y){
		return bg_data->scanned_h[x][y];
	}
	
	double func_Hashed(int x, int y){
		return bg_data->hashed[x][y];
	}	


	double func_Eval_s_per_Scan_s(int x, int y){
		if( bg_data->scanned[x][y] != 0){
			return bg_data->eval[x][y]/(double)bg_data->scanned[x][y];
		}else{
			return 0;
		}
	}

	double func_Eval_h_per_Hash(int x, int y){
		if( bg_data->hashed[x][y] != 0){
			return bg_data->eval_h[x][y]/(double)bg_data->hashed[x][y];
		}else{
			return 0;
		}
	}

	double func_Scan_h_per_Hash(int x, int y){
		if( bg_data->hashed[x][y] != 0){
			return bg_data->scanned_h[x][y]/(double)bg_data->hashed[x][y];
		}else{
			return 0;
		}
	}	
	
	double func_Find_per_Scan_s(int x, int y){
		if( bg_data->scanned[x][y] != 0){
			return 10.0*(bg_data->find[x][y]/(double)bg_data->scanned[x][y]);
		}else{
			return 0;
		}
	}

	double func_Scanned(int x, int y){
		return bg_data->scanned[x][y];
	}

	double func_Eval(int x, int y){
		return bg_data->eval[x][y];
	}
	double func_Eval_h(int x, int y){
		return bg_data->eval_h[x][y];
	}
	
	double func_MiDist(int x, int y){
		double mi[bg_data->nLights];
		double vec[bg_data->nLights];
		double value;
		int i;
		if(bg_data->entries[x][y] == 0){
			return NAN;
		}
		for(i = 0; i < bg_data->nLights; i++){
			mi[i]= bg_data->centro[x][y][i];
		}
		for(i = 0; i < bg_data->nLights; i++){
			//gets center bucket
			vec[i] = bg_data->u[i]*((bg_data->umin[x][y] + bg_data->umax[x][y])/2.0);
			vec[i]+= bg_data->v[i]*((bg_data->vmin[x][y] + bg_data->vmax[x][y])/2.0);
			vec[i]+= bg_data->baricenter[i];
			
		}
		for(i = 0; i < bg_data->nLights; i++){
			vec[i] = mi[i] - vec[i];
		}
		value = 0;
		for(i = 0; i < bg_data->nLights; i++){
			value += vec[i]*vec[i];
 		}
		return sqrt(value);
	}

	r_func = func_Zero;
	switch(pg->data_type){
		case ENTRIES:
			h_func = &func_Entries;
			break;
		case RADIUS:
			h_func = &func_Radius;
			break;
		case FIND:
			h_func = &func_Find;
			break;
		case SCANNED:
			h_func = &func_Scanned;
			break;
		case EVAL:
			h_func = &func_Eval;
			break;
		case EVAL_H:
			h_func = &func_Eval_h;
			break;
		case SCANNED_H:
			h_func = &func_Scanned_h;
			break;
		case HASHED:
			h_func = &func_Hashed;
			break;
		case MI_DIST:
			h_func = &func_MiDist;
			break;
		case EVAL_S_PER_SCAN_S:
			h_func = &func_Eval_s_per_Scan_s;
			break;
		case EVAL_H_PER_HASH:
			h_func = &func_Eval_h_per_Hash;
			break;
		case SCAN_H_PER_HASH:
			h_func = &func_Scan_h_per_Hash;
			break;
		case FIND_PER_SCAN_S:
			h_func = &func_Find_per_Scan_s;
			break;
		case DISTRAD:
			h_func = &func_MiDist;
			r_func = &func_Radius;
			break;
		default:
			fprintf(stderr,"ERROR - NO DATA TO BE PLOTED !\n");
	}
        double cell_size = 2.0*bg_data->R/(double)bg_data->N;
	switch(pg->plot_type){
		case GNUPLOT:
			break;
		case POVRAY:
			generatePOVFiles(prefix,h_func,r_func,bg_data->N,cell_size,pg);
			break;
		default:
			fprintf(stderr,"ERROR - NO PLOT TYPE DEFINED !\n");
	}
}

void getConfig(char* prefix,char* filename, plot_config* pg, int argc, char** argv){
	
	argparser_t *pp = argparser_new(stderr, argc, argv);
	argparser_set_help(pp, PROG_HELP);
	argparser_set_info(pp, PROG_INFO);
	argparser_process_help_info_options(pp);
	
	argparser_get_keyword(pp, "-filename");
	char *filename_str = argparser_get_next(pp);
	strcpy(filename,filename_str);
	argparser_get_keyword(pp, "-prefix");
	char *prefix_name = argparser_get_next(pp);
	strcpy(prefix,prefix_name);

	
	fprintf(stderr,"Program configuration \n");	

	if (argparser_keyword_present(pp, "-imageDimensions")) {
		pg->w = argparser_get_next_int(pp, 0, INT_MAX);
		pg->h = argparser_get_next_int(pp, 0, INT_MAX);
		pg->depth = argparser_get_next_int(pp, 0, INT_MAX);
	}
	else{
		pg->w = 640;
		pg->h = 480;
		pg->depth = 8; 
	}
	fprintf(stderr,"Image Output -> %dx%dx%d\n",pg->w,pg->h,pg->depth);

	if (argparser_keyword_present(pp, "-plotType")) {
		char* plot_type = argparser_get_next(pp);
		if(strcmp(plot_type,"GNUPLOT") == 0) pg->plot_type = GNUPLOT;
		else if(strcmp(plot_type,"POVRAY") == 0) pg->plot_type = POVRAY;
		else{
			argparser_error(pp, "Must specify an plot type");
		}
	}
	else{
		pg->plot_type = POVRAY;
	}
	fprintf(stderr,"Plot Type - %s\n",plot_type_desc[pg->plot_type]);
	
	
	if (argparser_keyword_present(pp, "-dataType")) {
		char* data_type = argparser_get_next(pp);
		if(strcmp(data_type,"ENTRIES") == 0) pg->data_type = ENTRIES;
		else if(strcmp(data_type,"RADIUS") == 0) pg->data_type = RADIUS;
		else if(strcmp(data_type,"FIND") == 0) pg->data_type = FIND;
		else if(strcmp(data_type,"SCANNED_H") == 0) pg->data_type = SCANNED_H;
		else if(strcmp(data_type,"SCANNED") == 0) pg->data_type = SCANNED;
		else if(strcmp(data_type,"EVAL") == 0) pg->data_type = EVAL;
		else if(strcmp(data_type,"HASHED") == 0) pg->data_type = HASHED;
		else if(strcmp(data_type,"EVAL_H") == 0) pg->data_type = EVAL_H;
		else if(strcmp(data_type,"MI_DIST") == 0) pg->data_type = MI_DIST;
		else if(strcmp(data_type,"DISTRAD") == 0) pg->data_type = DISTRAD;
		else if(strcmp(data_type,"EVAL_S_PER_SCAN_S") == 0) pg->data_type = EVAL_S_PER_SCAN_S;
		else if(strcmp(data_type,"EVAL_H_PER_HASH") == 0) pg->data_type = EVAL_H_PER_HASH;
		else if(strcmp(data_type,"SCAN_H_PER_HASH") == 0) pg->data_type = SCAN_H_PER_HASH;
		else if(strcmp(data_type,"FIND_PER_SCAN_S") == 0) pg->data_type = FIND_PER_SCAN_S;
		else{
			argparser_error(pp, "Must specify an data type to plot");
		}
	}
	else{
		pg->data_type = ENTRIES;
	}
        fprintf(stderr,"Data Type - %s \n",data_type_desc[pg->data_type]);
	
	if (argparser_keyword_present(pp, "-HWRatio")) {
		pg->HWRatio = argparser_get_next_double(pp, 0.001, 1.000);	
	}
	else{
		pg->HWRatio = 1.0/3.0;
	}
	fprintf(stderr,"HW Ratio - %f\n",pg->HWRatio);
	
	if (argparser_keyword_present(pp, "-cameraInclination")) {
		pg->camera_inclination = argparser_get_next_double(pp, 0.0, 90.0);
	}
	else{
		pg->camera_inclination = 45.0;
	}
	fprintf(stderr,"Camera Inclination - %f \n",pg->camera_inclination);

	if (argparser_keyword_present(pp, "-cameraAzimuth")) {
		pg->camera_azimuth = argparser_get_next_double(pp, 0.0, 360.0);
	}
	else{
		pg->camera_azimuth = 0.0;
	}
	fprintf(stderr,"Camera Azimuth - %f \n",pg->camera_azimuth);
	
	if (argparser_keyword_present(pp, "-cameraDist")) {
		pg->camera_dist = argparser_get_next_double(pp, 0.001, 100000.0);
	}
	else{
		pg->camera_dist = 1.25;
	}
	fprintf(stderr,"Camera Dist - %f \n",pg->camera_dist);
	
}

int main(int argc, char** argv){
	
	
	plot_config pg;
	char prefix[500];
	char filename[500];
	
	getConfig(prefix,filename,&pg,argc,argv);
	bucket_grid_data* bg = NULL;
	bg = loadBucketGridData(filename);
	if(bg == NULL){
		fprintf(stderr,"ERROR - FAILED TO GATHER BUCKET GRID DATA\n");
		return 1;
	}	
	plotScene(prefix, &pg,bg);
	
	return 0;
}


