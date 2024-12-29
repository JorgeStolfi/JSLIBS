#define _GNU_SOURCE
#include <stdio.h>
#include <stdlib.h>
#include <float_image.h>
#include <fast_hash.h>
#include <jsfile.h>
#include <assert.h>
// #include <mcheck.h>

int main(int argc, char** argv){
//   mtrace();
  
  if(argc < 6){
    fprintf(stderr,"program usage:\ngenerate_fasthash <prefix> <res> <sigma> <Normal deg> <Albedo deg>\n");
    return 1;
  }
  
  char* prefix = argv[1];
  int resolution;
  sscanf(argv[2],"%d",&resolution);
  double sigma;
  sscanf(argv[3],"%lf",&sigma);
  int normal_deg;
  sscanf(argv[4],"%d",&normal_deg);
  int albedo_deg;
  sscanf(argv[5],"%d",&albedo_deg);
  
  char* inTabname = NULL;
  char *inTabname = jsprintf("%s_TableData.txt",prefix);
  char* outHashname = NULL;
  char *outHashname = jsprintf("%s_FastHash.txt",prefix);
  char* outHashNormalname = NULL;
  char *outHashNormalname = jsprintf("%s_FastHash-N.fni",prefix);
  char* outHashAlbedoname = NULL;
  char *outHashAlbedoname = jsprintf("%s_FastHash-A.fni",prefix);
  char* outWeightname = NULL;
  char *outWeightname = jsprintf("%s_FastHash-W.fni",prefix);
  
  Tabela* tab;
  LoadTable(inTabname,&tab);
  if(tab == NULL){
     fprintf(stderr,"Cannot load %s \n",inTabname);
     assert(FALSE);
  }
  
  fprintf(stderr,"Generating Fast Hash...");
  fast_hash_t* fh = create_fasthash(resolution,tab,sigma,normal_deg,albedo_deg);
  fprintf(stderr,"OK\n");
  
  FILE* arq = open_write(outHashname,TRUE);
  SaveFastHash(arq,fh);
  fclose(arq);
  
  float_image_t* img_normal = FastHashNormalsToFNI(fh);
  arq = open_write(outHashNormalname,TRUE);
  float_image_write(arq,img_normal);
  fclose(arq);
  
  float_image_t* img_albedo = FastHashAlbedoToFNI(fh);
  arq = open_write(outHashAlbedoname,TRUE);
  float_image_write(arq,img_albedo);
  fclose(arq);
  
  float_image_t* img_weight = FastHashWeightsToFNI(fh);
  arq = open_write(outWeightname,TRUE);
  float_image_write(arq,img_weight);
  fclose(arq);
  
  return 0;
  
}