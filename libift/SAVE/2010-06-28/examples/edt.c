/*
  Copyright (C) <2003> <Alexandre Xavier Falcão>
  
  This program is free software; you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation; either version 2 of the License, or
  (at your option) any later version.

  This program is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.
  
  You should have received a copy of the GNU General Public License
  along with this program; if not, write to the Free Software
  Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA

  please see full copyright in COPYING file.
  -------------------------------------------------------------------------
  written by A.X. Falcão <afalcao@ic.unicamp.br>, March 27th 2003

  This program computes the exact Euclidean distance transform and
  discrete Voronoi diagram of points given in the seeds.txt file. One
  can optimize it by accumulating the absolute displacements dx and dy
  along the path, and by using a look-up table for all possible
  squared values inside the image. In other words, one can use
  path-cost function:

  (sum |x_{i+1}-x_i|)^2 + (sum |y_{i+1}-y_i|)^2, 

  for a path <(x_1,y_1),(x_2,y_2),...,(x_n,y_n)>, where (x_i,y_i) is
  the coordinate of pixel i.
*/

#include "gift.h"

int main(int argc, char **argv)
{
  timer tic,toc;
  FILE *ftime=fopen("time.txt","w");
  AnnImg *aimg=NULL;
  AdjRel *A=NULL;
  Image *lambda,*edt;
  int nx,ny,x,y,i,n,dist,p;
  FILE *fp;
  Pixel *S;

  /* The following block must the remarked when using non-linux machines 

  void *trash = malloc(1);                 
  struct mallinfo info;   
  int MemDinInicial, MemDinFinal;
  free(trash); 
  info = mallinfo();
  MemDinInicial = info.uordblks;
  */

  fp     = fopen("seeds.txt","r");
  fscanf(fp,"%d %d %d",&nx,&ny,&n);
  lambda = LambdaInit(nx,ny);
  S      = (Pixel *)calloc(n,sizeof(Pixel));
  i=1;
  while ((fscanf(fp,"%d %d",&x,&y))!=EOF){
    lambda->val[x+lambda->tbrow[y]]=x+lambda->tbrow[y];
    
    S[i-1].x = x; S[i-1].y = y; 
    i++;
  }
  fclose(fp);
  A = Circular(20.0);

  gettimeofday(&tic,NULL);  
  aimg = IFTFIFO(lambda,lambda,NULL,A,Fedt);
  gettimeofday(&toc,NULL);
  fprintf(ftime,"EDT in %f milliseconds\n",(toc.tv_sec-tic.tv_sec)*1000.0 + (toc.tv_usec-tic.tv_usec)*0.001);

  /* Check if it is correct and label DVD */
      
   for (y=0; y < ny; y++) 
    for (x=0; x < nx; x++) {      
      p    = x + aimg->img->tbrow[y];
      for (i=0; i < n; i++) {
	dist = (x-S[i].x)*(x-S[i].x) + (y-S[i].y)*(y-S[i].y);
	if (dist < aimg->cost->val[p]){
	  fprintf(stderr,"Error: (%d,%d) is closer to (%d,%d)\n",x,y,S[i].x,S[i].y);
	  exit(-1);
	}
      }
    }

  edt = SQRT(aimg->cost);
  WriteImage(edt,"edt.pgm");
  WriteImage(aimg->label,"dvd.pgm");

  DestroyImage(&edt);
  DestroyAdjRel(&A);
  DestroyImage(&lambda);
  DestroyAnnImg(&aimg);
  free(S);

  /* The following block must the remarked when using non-linux machines 

  info = mallinfo();
  MemDinFinal = info.uordblks;
  if (MemDinInicial!=MemDinFinal)
    printf("\n\nDinamic memory was not completely deallocated (%d, %d)\n",
	   MemDinInicial,MemDinFinal);   
  */
  fclose(ftime);

  return(0);

}
