#include <stdio.h>
#include <tabela.h>
#include <float_image.h>
#include <argparser.h>
#include <float_pnm_image_io.h>
#include <assert.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#define IS_PPM 0
#define IS_TABLE 1


#define PROG_NAME "pnm_to_sig"

#define PROG_HELP \
  PROG_NAME " \\\n" \
  "    -nImages {NUM} \\\n" \
  "    -prefix {FILE_PREFIX} \\\n" \
  "    -Images {FILENAME}...  \\\n" \
  "    -nTables {NUM} \\\n" \
  "    -Tables {FILENAME}... \\\n" \
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

typedef struct options_t {
	char* prefix;
	char** Images;
	double gamma;
	int nImages;
	char** Tables;
	int channel;
	int nTables;
	int log;
	
} options_t;

typedef struct statistics_t{
	double total_time;
	double average_time;
	int querys_found;
	int total_querys;
	int euclidian_dists;
} statistics_t;



void PrintSignatureImage(FILE* arq_sig,FILE* arq_mag,FILE* arq_index,float_image_t** im,int num_images,int channel);

void PrintSignatureImage(FILE* arq_sig,FILE* arq_mag,FILE* arq_index,float_image_t** im,int num_images,int channel){
	int nc,nx,ny;
	nc = im[0]->sz[0];
	nx = im[0]->sz[1];
	ny = im[0]->sz[2];
	assert(channel < nc);
	int x,y;
	fprintf(arq_index,"%d\n",nx*ny);
	for(x = 0; x < nx; x++){
		for(y = 0; y < ny; y++){
			double SO[num_images];
			double sMag = 0;
			int i;
			for(i = 0; i < num_images; i++){
				SO[i] = float_image_get_sample(im[i],channel,x,y);
				sMag+= SO[i]*SO[i];
			}
			sMag = sqrt(sMag);
			for(i = 0; i < num_images; i++){
				fprintf(arq_sig,"%4.7lf ",SO[i]/sMag);
			}
			fprintf(arq_sig,"\n");
			fprintf(arq_mag,"%4.7lf\n",sMag);
			fprintf(arq_index,"%d %d\n",x,y);
		}
	}
}


void PrintSignatureTable(FILE* arq_sig, FILE* arq_mag, Tabela* tab);

void PrintSignatureTable(FILE* arq_sig, FILE* arq_mag, Tabela* tab){
	int num_linhas = get_num_linhas(tab);
	int i;
	for(i = 0; i < num_linhas; i++){
		const double* dir = get_intdir(tab,i);
		int j;
		for(j = 0; j < get_num_luzes(tab); j++){
			fprintf(arq_sig,"%4.7lf ",dir[j]);
		}
		fprintf(arq_sig,"\n");
		double gMag = get_intmag(tab,i);
		fprintf(arq_mag,"%4.7lf\n",gMag);
	}

}


options_t *parse_options(int argc, char **argv);
options_t *parse_options(int argc, char **argv){
	argparser_t *pp = argparser_new(stderr, argc, argv);
	argparser_set_help(pp, PROG_HELP);
  	argparser_set_info(pp, PROG_INFO);
  	argparser_process_help_info_options(pp);
	
	options_t* o =  (options_t *)malloc(sizeof(options_t));

	argparser_get_keyword(pp, "-prefix");
	o->prefix = argparser_get_next(pp);
	if (argparser_keyword_present(pp, "-nImages")){
		o->nImages = argparser_get_next_int(pp, 1, 1000);
		argparser_get_keyword(pp, "-gamma");
  		o->gamma = argparser_get_next_double(pp, 0.100, 9.000);		
		o->Images = (char**) malloc(sizeof(char*)*o->nImages);
		argparser_get_keyword(pp, "-Images");
		int i;
		for(i = 0; i < o->nImages; i++){
			o->Images[i] = argparser_get_next(pp);
		}
	}else{
		o->nImages = 0;
		o->Images = NULL;
	}
	if (argparser_keyword_present(pp, "-nTables")){
		o->nTables = argparser_get_next_int(pp, 1, 1000);
		o->Tables = (char**) malloc(sizeof(char*)*o->nTables);
		argparser_get_keyword(pp, "-Tables");
		int i;
		for(i = 0; i < o->nTables; i++){
			o->Tables[i] = argparser_get_next(pp);
		}
	}else{
		o->nTables = 0;
		o->Tables = NULL;
	}
	if (argparser_keyword_present(pp, "-channel")){
		o->channel = argparser_get_next_int(pp, 0, 2);
	}else{
		o->channel = 0;
	}
	o->log = 0;
	if (argparser_keyword_present(pp, "-log")){
		o->log = 1;
		if( (o->nImages == 0) || (o->nTables == 0) ){
			argparser_error(pp, "Must specify image and table processing for log!");
		}
	}
	if( (o->nImages == 0) && (o->nTables == 0) ){
		
		argparser_error(pp, "Must specify image or table processing !");
	}

	return o;	
}

float_image_t* GenerateNormalMapfromLSH(FILE* arq,FILE* index_file,int nx,int ny,Tabela* tab,float_image_t* query_image);
void GenerateNextSignatureRemains(float_image_t* query_image,float_image_t** image_list,int nImages, char* prefix, int c);
void ComputeStatsFromQuery(statistics_t*  stats,float_image_t* query_image);


int main(int argc, char** argv){
	options_t* o;
	o = parse_options(argc,argv);
	float_image_t** image_list = NULL;
	if( o->nImages > 0){
		double bias = 0.0;
		if(o->gamma == 0.0) bias = VIEW_BIAS;
		image_list =  float_pnm_image_list_read(o->nImages,o->Images,FALSE,o->gamma,bias,FALSE,TRUE,FALSE);
		int c;
		c = o->channel;
		

		if(o->log == 0){
			char out_sig_name[500];
			sprintf(out_sig_name,"%s_%d_imgsig.txt",o->prefix,c);
			char out_mag_name[500];
			sprintf(out_mag_name,"%s_%d_imgmag.txt",o->prefix,c);
			char out_index_name[500];
			sprintf(out_index_name,"%s_%d_imgidx.txt",o->prefix,c);
			FILE* out_sig_file = fopen(out_sig_name,"wt");
			if(out_sig_file == NULL){
				fprintf(stderr,"Arquivo de saida nao pode ser aberto %s\n",out_sig_name);
				return 1;
			}
			
			FILE* out_mag_file = fopen(out_mag_name,"wt");
			if(out_mag_file == NULL){
				fprintf(stderr,"Arquivo de saida nao pode ser aberto %s\n",out_mag_name);
				return 1;
			}
			FILE* out_index_file = fopen(out_index_name,"wt");
			if(out_index_file == NULL){
				fprintf(stderr,"Arquivo de saida nao pode ser aberto %s\n",out_index_name);
				return 1;
			}
			puts("RESETED");
			PrintSignatureImage(out_sig_file,out_mag_file,out_index_file,image_list,o->nImages,c);
			int nx = image_list[0]->sz[1];
			int ny = image_list[0]->sz[2];
			float_image_t* query_image;
			FILE* query_file;
			char out_query_name[500];
			sprintf(out_query_name,"%s_%d_query.txt",o->prefix,c);
			query_file = fopen(out_query_name,"w");
			query_image = float_image_new(5,nx,ny);
			float_image_fill_channel(query_image,0,0.0);
			float_image_fill_channel(query_image,1,0.0);
			float_image_fill_channel(query_image,2,0.0);
			float_image_fill_channel(query_image,3,0.0);
			float_image_fill_channel(query_image,4,0.0);
			float_image_write(query_file,query_image);
			fclose(query_file);
			

			fclose(out_sig_file);
			fclose(out_mag_file);
		}
		
		
	}
	
	if(o->nTables > 0){
		int c;
		c = o->channel;
		char out_sig_name[500];
		sprintf(out_sig_name,"%s_%d_tabsig.txt",o->prefix,c);
		char out_mag_name[500];
		sprintf(out_mag_name,"%s_%d_tabmag.txt",o->prefix,c);
		FILE* out_sig_file = fopen(out_sig_name,"wt");
		if(out_sig_file == NULL){
			fprintf(stderr,"Arquivo de saida nao pode ser aberto %s\n",out_sig_name);
			return 1;
		}
		
		FILE* out_mag_file = fopen(out_mag_name,"wt");
		if(out_sig_file == NULL){
			fprintf(stderr,"Arquivo de saida nao pode ser aberto %s\n",out_mag_name);
			return 1;
		}
		Tabela* tab;
		LoadTable(o->Tables[c],&tab);
		PrintSignatureTable(out_sig_file,out_mag_file,tab);
		fclose(out_sig_file);
		fclose(out_mag_file);
			
		if(o->log != 0){
			statistics_t stats;
			FILE* log_file;
			char log_file_name[500];
			sprintf(log_file_name,"%s_search_%d.txt",o->prefix,c);
			log_file = fopen(log_file_name,"rt");
			if(log_file == NULL){
				fprintf(stderr,"Arquivo de entrada nao pode ser aberto %s\n",log_file_name);
				return 1;
			}
			int nx = image_list[0]->sz[1];
			int ny = image_list[0]->sz[2];
			//Open the residual image
			float_image_t* query_image;
			char query_file_name[500];
			sprintf(query_file_name,"%s_%d_query.txt",o->prefix,c);
				
			FILE* query_file = fopen(query_file_name,"r");
			query_image = float_image_read(query_file);
			fclose(query_file);
			assert( (nx == query_image->sz[1]) && (ny == query_image->sz[2]) && ( query_image->sz[0] == 5) );
			//
			char out_index_name[500];
			sprintf(out_index_name,"%s_%d_imgidx.txt",o->prefix,c);
			FILE* index_file = fopen(out_index_name,"rt");
			float_image_t* normal_map = GenerateNormalMapfromLSH(log_file,index_file,nx,ny,tab,query_image);
			fclose(index_file);
			fclose(log_file);
			//Write the residual file bunch of sigs
			GenerateNextSignatureRemains(query_image,image_list,o->nImages,o->prefix,c);
			puts("AQUI");
			ComputeStatsFromQuery(&stats,query_image);
			
			char query_image_name[500];
			sprintf(query_image_name,"%s_querymap_%d.fni",o->prefix,c);
			FILE* query_image_file = fopen(query_image_name,"wt");
		//		float_pnm_image_write(normal_file_name, normal_map,o->gamma, bias,FALSE,TRUE,FALSE);
			float_image_write(query_image_file,query_image);
			fclose(query_image_file);			

			
			char normal_file_name[500];
			sprintf(normal_file_name,"%s_normalmap_%d.fni",o->prefix,c);
			double bias = 0.0;
			if(o->gamma == 0.0) bias = VIEW_BIAS;
			FILE* normal_map_file = fopen(normal_file_name,"wt");
		//		float_pnm_image_write(normal_file_name, normal_map,o->gamma, bias,FALSE,TRUE,FALSE);
			float_image_write(normal_map_file,normal_map);
			fclose(normal_map_file);
			FILE* stats_file;
			char stats_file_name[500];
			sprintf(stats_file_name,"%s_stats_%d.txt",o->prefix,c);
			stats_file = fopen(stats_file_name,"wt");
			if(stats_file == NULL){
				fprintf(stderr,"Arquivo de entrada nao pode ser aberto %s\n",log_file_name);
				return 1;
			}
			fprintf(stats_file,"Total Time: %lf\n",stats.total_time);
			fprintf(stats_file,"Average time: %lf \n", stats.total_time/(float)stats.querys_found);
			fprintf(stats_file,"Querys found %d of %d ( %lf%%)\n",stats.querys_found,stats.total_querys,100.0*stats.querys_found/(float)stats.total_querys);
			fprintf(stats_file,"Total Euclid %d Average %9.6lf \n",stats.euclidian_dists, stats.euclidian_dists/(float)stats.total_querys);
			fclose(stats_file);
			
		}
	}
	
		
		
	

	return 0;
}

float_image_t* GenerateNormalMapfromLSH(FILE* arq,FILE* index_file,int nx,int ny,Tabela* tab,float_image_t* query_image){
	float_image_t* novo = float_image_new(3,nx,ny);
	int x,y;
	int total_lines ;
	fscanf(index_file,"%d",&total_lines);
	int i;

	//first fill query_image 	
	for(i = 0; i < total_lines; i++){
		r3_t normal;
		r3_zero(&normal);
		int test = fscanf(index_file,"%d %d",&x,&y);
		if(test != 2){
			fprintf(stderr,"ERROR AT line channel %d\n",i);
			assert(test ==2);
		}
		char temp[500];
		double query_time;
		int lines_found;
		int dist_euclids;
		fscanf(arq,"%s %lf %s %d %s %d",temp,&query_time,temp,&lines_found,temp,&dist_euclids);
		//printf("Query  %lf Lines Found %d Dists %d",query_time,lines_found,dist_euclids);
		float_image_set_sample(query_image,4,x,y,dist_euclids); //euclidean distances, written whenere we find something or not !
		if(lines_found !=0){
			int line;
			double dist;
			fscanf(arq,"%s %d %s %lf",temp,&line, temp, &dist);
		//	printf("Line %d Distance %lf \n ",line, dist);
		//	normal = get_normal(tab,line);

			
			int already_found = float_image_get_sample(query_image,1,x,y);
			if( already_found == 0){
				float_image_set_sample(query_image,0,x,y,query_time*1000.0); //channel 0 is query time
				float_image_set_sample(query_image,1,x,y,lines_found); //channel 1 is lines_found
				float_image_set_sample(query_image,2,x,y,line); //channel 2 line
				float_image_set_sample(query_image,3,x,y,dist); //channel 2 distance
				
	
			}
			
		}else{
		//	printf("\n");
			int already_found = float_image_get_sample(query_image,1,x,y);
			if( already_found != 0){
				//int line = float_image_get_sample(query_image,2,x,y);
				//normal = get_normal(tab,line);
			}
		}
		//printf("%lf %lf %lf \n",normal.c[0],normal.c[1],normal.c[2]);
// 		float_image_set_sample(novo,0,x,y,normal.c[0]);
// 		float_image_set_sample(novo,1,x,y,normal.c[1]);
// 		float_image_set_sample(novo,2,x,y,normal.c[2]);
	}

	int NX,NY;
	NX = query_image->sz[1];
	NY = query_image->sz[2];

	for(x = 0; x < NX; x++){
		for(y = 0; y < NY; y++){
			int line = float_image_get_sample(query_image,2,x,y);
			r3_t normal = (r3_t){{ 0, 0 ,0}};
			int already_found = float_image_get_sample(query_image,1,x,y);
			if(already_found != 0) 	normal = get_normal(tab,line);

			float_image_set_sample(novo,0,x,y,normal.c[0]);
 			float_image_set_sample(novo,1,x,y,normal.c[1]);
	 		float_image_set_sample(novo,2,x,y,normal.c[2]);
		}
	}

	
	return novo;
}

void GenerateNextSignatureRemains(float_image_t* query_image,float_image_t** image_list,int nImages, char* prefix, int c){
	int x,y;
	int nx = image_list[0]->sz[1];
	int ny = image_list[0]->sz[2];
	
	FILE* arq_index;
	FILE* arq_mag;
	FILE* arq_sig;
	
	char out_sig_name[500];
	sprintf(out_sig_name,"%s_%d_imgsig.txt",prefix,c);
	char out_mag_name[500];
	sprintf(out_mag_name,"%s_%d_imgmag.txt",prefix,c);
	
	
	arq_mag = fopen(out_mag_name,"wt");
	arq_sig = fopen(out_sig_name,"wt");

	int total = 0;
	for(x = 0; x < nx; x++){
		for(y = 0; y < ny; y++){
			int already_found = float_image_get_sample(query_image,1,x,y);
			if( already_found == 0){
				double SO[nImages];
				double sMag = 0;
				int i;
				for(i = 0; i < nImages; i++){
					SO[i] = float_image_get_sample(image_list[i],c,x,y);
					sMag+= SO[i]*SO[i];
				}
				sMag = sqrt(sMag);
				for(i = 0; i < nImages; i++){
					fprintf(arq_sig,"%4.7lf ",SO[i]/sMag);
				}
				fprintf(arq_sig,"\n");
				fprintf(arq_mag,"%4.7lf\n",sMag);
				total++;
			}
		}
	}
	fclose(arq_sig);
	fclose(arq_mag);
	
	char out_index_name[500];
	sprintf(out_index_name,"%s_%d_imgidx.txt",prefix,c);
	arq_index = fopen(out_index_name,"wt");
	fprintf(arq_index,"%d\n",total);
	for(x = 0; x < nx; x++){
		for(y = 0; y < ny; y++){
			int already_found = float_image_get_sample(query_image,1,x,y);
			if( already_found == 0){
				fprintf(arq_index,"%d %d\n",x,y);
			}
		}
	}
	fclose(arq_index);
	
	char out_query_name[500];
	sprintf(out_query_name,"%s_%d_query.txt",prefix,c);
	FILE* arq_query = fopen(out_query_name,"wt");
	float_image_write(arq_query,query_image);
	fclose(arq_query);
}

void ComputeStatsFromQuery(statistics_t*  stats,float_image_t* query_image){
	int x,y;
	int nx = query_image->sz[1];
	int ny = query_image->sz[2];
	stats->querys_found = 0;
	stats->euclidian_dists = 0;
	stats->total_time = 0;
	for(x = 0; x < nx; x++){
		for(y = 0; y < ny; y++){
			stats->euclidian_dists+= float_image_get_sample(query_image,4,x,y);
			int already_found = float_image_get_sample(query_image,1,x,y);
			if( already_found != 0){
				/*	double total_time;
					int querys_found;
					int total_querys;
					int euclidian_dists;
				*/
				stats->querys_found++;
				stats->total_time+= float_image_get_sample(query_image,0,x,y);
				
			}
		}
	}
	stats->total_querys = nx*ny;
	
}
