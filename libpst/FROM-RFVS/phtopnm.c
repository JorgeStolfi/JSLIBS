#include <bool.h>
#include <frgb_ops.h>
#include <affirm.h>
#include <argparser.h>
#include <jspnm_image.h>
#include <jspnm.h>
#include <jsfile.h>

#include <jspnm_Pike_F100_image.h>
// pnm_image_t* FNItoPNM(float_image_t* fni_img){
//   int NC,NX,NY;
//   NC = fni_img->sz[0];
//   NX = fni_img->sz[1];
//   NY = fni_img->sz[2];
//   pnm_image_t* pnm_img = pnm_image_new(NX, NY, NC);
//   for(c = 0; c < NC; c++){
//     for(x =  0; x < NX; x++){
//       for(y = 0; y < NY; y++){
// 	double val = float_image_get_sample(fni_img,c,x,y);
// 	pnm_image_set_sample(pnm_img,c, x, y, val);
//       }
//     }	
//   }
//   return pnm_img;
// }

pnm_image_t* ReadPNM(char* filename);
pnm_image_t* ReadPNM(char* filename){
//   FILE* arq = stdin;
//   if(filename != NULL){
//     arq = open_read(filename,TRUE);
//   }else{
//     fprintf(stderr,"Reading from standard input\n");
//   }
  pnm_image_t * img  = pnm_image_read(filename,TRUE);
//   if(filename != NULL){
//     fclose(arq);
//   }
  return img;
}

void WritePNM(char* filename, pnm_image_t* img);
void WritePNM(char* filename, pnm_image_t* img){
  FILE* arq = stdout;
  if(filename != NULL){
    arq = open_write(filename,TRUE);
  }else{
    fprintf(stderr,"Writing into standard output\n");
  }
  bool_t forceplain = FALSE;
  pnm_image_fwrite(arq, img, forceplain);
  if(filename != NULL){
    fclose(arq);
  }
}

int main(int argc, char** argv){
  if(argc < 3){
    fprintf(stderr,"Program usage:\nphtopnm <infile> <outfile> [GRGB | RGGB]\n");
    return 1;
  }
  char* in_file = argv[1];
  char* out_file = argv[2];
  char* chr_debayer = "GRGB";
  if( argc >= 4){
    chr_debayer = argv[3];
  }
  
  pnm_image_t* in_img = ReadPNM(in_file);
  pnm_image_t* out_img = NULL;
  
  if( strcmp(chr_debayer,"RGGB") == 0){
    out_img = pnm_Canon_EOS50D_raw16_debayer(in_img, TRUE,FALSE);
  }else{
    out_img = pnm_Pike_F100_raw16_debayer(in_img, TRUE,FALSE);
  }
  WritePNM(out_file,out_img);
  
  return 0;
}