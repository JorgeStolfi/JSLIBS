#include <stdio.h>
#include <stdlib.h>
#include <float_pnm_image_io.h>
#include <float_image.h>
#include <math.h>
#include <jsfile.h>
#include <pam.h>
#include <assert.h>


float_image_t* Read_PAM(char* filename){
  FILE* arq = fopen(filename,"r");
  assert(filename != NULL);
  struct pam inpam;
//   inpam.comment_p = NULL;
  tuple** tuple_pam = pnm_readpam(arq,&inpam,sizeof(struct pam));
  
  int NX,NY,NC;
  NC = 1;
  NX = inpam.width;
  NY = inpam.height;
  unsigned int maxval = inpam.maxval;
  fprintf(stderr,"MAXVAL = %d\n",maxval);
  
  float_image_t* novo = float_image_new(NC,NX,NY);
  int x,y;
  for(x = 0; x < NX; x++){
    for(y = 0; y < NY; y++){
      tuple pixel_val = tuple_pam[y][x];
      unsigned int ival = pixel_val[0];
      double val = ival/(double)maxval;
      int iy = NY - y - 1 ;
      float_image_set_sample(novo,0,x,iy,val);
    }
  }
  
  fclose(arq);
  return novo;
}

float_image_t* JoinNormals(float_image_t** images){
  int NC,NX,NY;
  NC = images[0]->sz[0];
  NX = images[0]->sz[1];
  NY = images[0]->sz[2];
  
  int x,y;
  
  float_image_t* novo   = float_image_new(3,NX,NY);
  for(x = 0; x < NX; x++){
    for(y = 0; y < NY; y++){
      
      double dx,dy,dz;
      dx = float_image_get_sample(images[0],0,x,y);
      dy = float_image_get_sample(images[1],0,x,y);
      dz = float_image_get_sample(images[2],0,x,y);
      
            
        /*Correct what PovRay messes up. The formula is 
	Np = 0.5 + asin(Nr)/Pi
	where Np is the normal given by PovRay and Nr the real normal value
	so, we compute the inverse
	Nr = sin(Pi*(Np - 0.5))
	*/
      
      dx = sin(M_PI*(dx - 0.5));
      dy = sin(M_PI*(dy - 0.5));
      dz = sin(M_PI*(dz - 0.5));
      
      
   //   fprintf(stderr,"dx %lf dy %lf dz %lf\n",dx,dy,dz);
    /*   
      dx = 2.0*dx - 1.0;
      dy = 2.0*dy - 1.0;
      dz = 2.0*dz - 1.0;
      
    */
      
      double m = sqrt(dx*dx + dy*dy + dz*dz);
      
   //   fprintf(stderr,"dx %lf dy %lf dz %lf M %lf\n",dx,dy,dz,m);
      
      dx = dx/m;
      dy = dy/m;
      dz = dz/m;

      
    //  fprintf(stderr,"DX %lf DY %lf DZ %lf\n",dx,dy,dz);
      
    //  char c;
    //  scanf("%c",&c);
      
      if(m == 0) dz = dy = dx = 0;
      
      float_image_set_sample(novo,0,x,y,dx);
      float_image_set_sample(novo,1,x,y,dy);
      float_image_set_sample(novo,2,x,y,dz);

    }
  }
  
  return novo;
}

int main(int argc, char** argv){
  
  if(argc < 5){
    fprintf(stderr,"Program usage:\nrgb3_to_fni <X> <Y> <z>  <outfile>\n");
    return 1;
  }
  
  char* in_filename[3];
  in_filename[0] = argv[1];
  in_filename[1] = argv[2];
  in_filename[2] = argv[3];
  char* out_filename = argv[4];
  
  float_image_t** images = (float_image_t**)malloc(sizeof(float_image_t*)*3);
  int i;
  for(i = 0;i < 3; i++){
    //images[i] = float_pnm_image_read(in_filename[i], 1.0,0.0, TRUE,TRUE,FALSE);
    images[i] = Read_PAM(in_filename[i]);
  }
  
  float_image_t* novo = JoinNormals(images);
  
  FILE* out_file = open_write(out_filename,TRUE);
  float_image_write(out_file,novo);
  fclose(out_file);
  
  return 0;  
}

