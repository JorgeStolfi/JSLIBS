#define _GNU_SOURCE
#include <stdio.h>
#include <analytic_inversion.h>
#include <jsfile.h>


int main(int argc, char** argv){
  if(argc < 4){
    fprintf(stderr,"Program usage:\analytic_inversion_test <in_file> <out_prefix> <mapSize>\n");
    return 1;
  }
  char* in_name=argv[1];
  char* ou_prefix=argv[2];
  int mapSize;
  sscanf(argv[3],"%d",&mapSize);
  
  FILE* in_arq = open_read(in_name,TRUE);
  analytic_inversion_t * ai = analytic_inversion_read(in_arq);
  fclose(in_arq);
  r3_t view_dir = (r3_t){{0,0,1}};
  analytic_inversion_generate_SG_plot(ou_prefix, ai, mapSize ,view_dir);
  
 char* inv_filename_map =  NULL;
 char *inv_filename_map = jsprintf("%s_Imapa.fni",ou_prefix);
 FILE* inv_arq_map = open_write(inv_filename_map,TRUE);
 float_image_t* im_map = analytic_inversion_plot_map(ai,mapSize);
 float_image_write(inv_arq_map,im_map);
 fclose(inv_arq_map);
  
  return 0;
  
}