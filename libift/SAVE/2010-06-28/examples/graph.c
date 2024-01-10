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

  This program illustrates some examples of IFTs for graphs defined in
  the graph1.txt and graph2.txt, under different tiebreaking policies.

*/

#include "gift.h"
#include "libps/ps.h"
#include "libps/js.h"

void WriteEPS(AnnImg *aimg, AdjRel *A, char *figname) 
{
  Pixel u,v;
  int i, ifile, pass;
  double unit = 32.0; /* Pixel spacing in pt */
  double mm = 72/25.4; /* One mm in pt. */
  double hx = 0.18*unit; /* Label x displacement in pt. */
  double hy = 0.25*unit; /* Label y displacement in pt. */

  for (ifile = 0; ifile < 2; ifile++)
    { 
      char *filename = txtcat(figname, (ifile ==0 ? "-a.eps" : "-b.eps"));
      FILE *fp = fopen(filename, "w"); 
      ps_begin_figure(fp, 0.0, unit*(aimg->img->ncols), 0.0, unit*(aimg->img->nrows));
      ps_set_label_font(fp, "TimesRoman", 18.0);
      ps_set_pen(fp, 0,0,0, 0.15, 0,0);
      for (pass = 0; pass < 2; pass++)
	{
	  for (u.y=0; u.y < aimg->img->nrows; u.y++)
    	    {
    	      for (u.x=0; u.x < aimg->img->ncols; u.x++) 
    		{
    		  int p = u.x + aimg->img->tbrow[u.y];
    		  double ux = unit*(u.x + 0.25);
    		  double uy = unit*(aimg->img->nrows-1 - u.y + 0.25);
    		  bool proot = (aimg->pred->val[p] == NIL);
    		  double pradius = unit*(proot ? 0.150 : 0.100);  /* radius in pt */
                  if (pass == 0)
                    {
		      /* Draw arrows: */
		      for (i=1; i < A->n; i++) 
			{
			  v.x = u.x + A->dx[i];
			  v.y = u.y + A->dy[i];
			  if (ValidPixel(aimg->img,v.x,v.y))
			    {
			      int q  = v.x + aimg->img->tbrow[v.y]; /* Index of neighbor. */
			      /* Decide whether arc should be drawn: */
			      bool drawit;
			      if (ifile == 0)
				{ drawit = true; }
			      else
				{ drawit = (q == aimg->pred->val[p]); }
			      if (drawit)
				{
				  bool qroot = (aimg->pred->val[q] == NIL);
				  double qradius = unit*(qroot ? 0.150 : 0.100);  /* radius in pt */
				  double vx = unit*(v.x+0.25);
				  double vy = unit*(aimg->img->nrows-1 - v.y + 0.25);
				  /* A vector along the arc: */
				  double dlm = hypot(vx - ux, vy - uy);
				  double dlx = (vx - ux)/dlm;
				  double dly = (vy - uy)/dlm;
				  /* A vector orthogonal to the arc: */
				  double dtx = -dly;
				  double dty = dlx;
				  /* Arrowhead tip: */
				  double ax = vx - qradius*dlx;
				  double ay = vy - qradius*dly;
				  /* Arrowhead corners: */
				  double bx = ax - unit*(0.15*dlx + 0.08*dtx);
				  double by = ay - unit*(0.15*dly + 0.08*dty);
				  double cx = ax - unit*(0.15*dlx - 0.08*dtx);
				  double cy = ay - unit*(0.15*dly - 0.08*dty);
				  /* Draw arc from p to q: */
				  ps_draw_segment(fp, ux,uy, vx,vy);
				  /* Paint arrowhead: */
				  ps_fill_triangle(fp, ax, ay, bx, by, cx, cy,  0,0,0);
				} 
			    }
			}
		    }
                  else
		    { 
		      int value = (ifile == 0 ? aimg->img->val[p] : aimg->cost->val[p]); /* Image value to show. */
		      char *cval = fmtint(value,1);
		      int dhx, dhy;
		      /* Paint dot and image value: */
                      ps_set_pen(fp, 1,1,1, 0.15, 0,0);
		      for (dhy = -1; dhy <= 1; dhy++)
                        for (dhx = -1; dhx <= 1; dhx++)
                          ps_put_label(fp, cval, ux+hx+1.0*dhx, uy+hy+1.0*dhy, 0,0);
		      ps_set_pen(fp, 0,0,0, 0.15, 0,0);
		      ps_put_label(fp, cval, ux+hx, uy+hy, 0,0);    
		      ps_fill_dot(fp, ux,uy,pradius/mm,0,0,0);
                      free(cval);
		    }
    		}
    	    }
	}
      ps_end_figure (fp);
      fclose(fp);
      free(filename);
    }
}

int main(int argc, char **argv)
{
  timer tic,toc;
  FILE *ftime=fopen("time.txt","w");
  AnnImg *aimg=NULL;
  AdjRel *A=NULL;
  Image *lambda,*img,*handicap;
  FILE *fp;
  int nx,ny,i,x,y,n;

  /* The following block must the remarked when using non-linux machines 

  void *trash = malloc(1);                 
  struct mallinfo info;   
  int MemDinInicial, MemDinFinal;
  free(trash); 
  info = mallinfo();
  MemDinInicial = info.uordblks;
  */

  /*---- Example of an IFT ------------------------------------------*/

  fp     = fopen("graph1.txt","r");
  fscanf(fp,"%d %d %d",&nx,&ny,&n);
  lambda   = LambdaInit(nx,ny);
  handicap = CreateImage(nx,ny);
  img      = CreateImage(nx,ny);
  for (i=0; i < n; i++){
    fscanf(fp,"%d %d",&x,&y);
    lambda->val[x+lambda->tbrow[y]]=i+1;
  }
  for (y=0; y < ny; y++)
    for (x=0; x < nx; x++){
      fscanf(fp,"%d",&img->val[x+img->tbrow[y]]);
    }
  fclose(fp);


  A = Circular(1.0);
  gettimeofday(&tic,NULL);  
  aimg = IFTFIFO(img,lambda,handicap,A,Fsuprec);
  gettimeofday(&toc,NULL);
  fprintf(ftime,"Watershed in %f milliseconds\n",(toc.tv_sec-tic.tv_sec)*1000.0 + (toc.tv_usec-tic.tv_usec)*0.001);

  WriteEPS(aimg,A,"graph1");

  DestroyImage(&handicap);
  DestroyImage(&lambda);
  DestroyImage(&img);
  DestroyAnnImg(&aimg);
  DestroyAdjRel(&A);

  /*----- Example of tiebreak FIFO with 4-neighborhood -------*/

  fp     = fopen("graph2.txt","r");
  fscanf(fp,"%d %d %d",&nx,&ny,&n);
  lambda   = LambdaInit(nx,ny);
  handicap = CreateImage(nx,ny);
  img      = CreateImage(nx,ny);
  for (i=0; i < n; i++){
    fscanf(fp,"%d %d",&x,&y);
    lambda->val[x+lambda->tbrow[y]]=i+1;
  }
  for (y=0; y < ny; y++)
    for (x=0; x < nx; x++){
      fscanf(fp,"%d",&img->val[x+img->tbrow[y]]);
    }
  fclose(fp);

  A = Circular(1.0);
  gettimeofday(&tic,NULL);  
  aimg = IFTFIFO(img,lambda,handicap,A,Fsuprec);
  gettimeofday(&toc,NULL);
  fprintf(ftime,"Watershed in %f milliseconds\n",(toc.tv_sec-tic.tv_sec)*1000.0 + (toc.tv_usec-tic.tv_usec)*0.001);

  WriteEPS(aimg,A,"graph2-1");

  /*----- Example of tiebreak LIFO with 4-neighborhood -------*/

  DestroyAnnImg(&aimg);

  gettimeofday(&tic,NULL);  
  aimg = IFTLIFO(img,lambda,handicap,A,Fsuprec);
  gettimeofday(&toc,NULL);
  fprintf(ftime,"Watershed in %f milliseconds\n",(toc.tv_sec-tic.tv_sec)*1000.0 + (toc.tv_usec-tic.tv_usec)*0.001);

  WriteEPS(aimg,A,"graph2-2");

  /*----- Example of tiebreak FIFO with 8-neighborhood -------*/

  DestroyAnnImg(&aimg);
  DestroyAdjRel(&A);
  A = Circular(1.5);
  
  gettimeofday(&tic,NULL);  
  aimg = IFTFIFO(img,lambda,handicap,A,Fsuprec);
  gettimeofday(&toc,NULL);
  fprintf(ftime,"Watershed in %f milliseconds\n",(toc.tv_sec-tic.tv_sec)*1000.0 + (toc.tv_usec-tic.tv_usec)*0.001);

  WriteEPS(aimg,A,"graph2-3");

  DestroyImage(&handicap);
  DestroyImage(&lambda);
  DestroyImage(&img);
  DestroyAnnImg(&aimg);
  DestroyAdjRel(&A);

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
