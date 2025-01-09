#define _GNU_SOURCE
#include <vetorn.h>
#include <normais.h>
#include <tabela.h>
#include <affirm.h>
#include <float_image.h>
#include <float_image_masked.h>
#include <float_image_io_pnm.h>
#include <jsfile.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <assert.h>
#include <float.h>
#include <r2.h>
#include <r3.h>
#include <rn.h>
#include <r3x3.h>
#include <math_utils.h>

#define MINMAG 0.01

struct LinhaTabela{
  r3_t normal;      // normal - normal da superfície.
  double *intdir;   // intdir[num_luzes] - assinatura normalizada
  double intmag;    // intmag norma euclidiana da assinatura
};

struct Tabela{
  int num_linhas;
  int num_luzes;
  r3_t view_dir; //stores the camera direction as seen from the center of the gauge
  LinhaTabela* tabela;
};


r3_t calcula_normal_esfera(r2_t *uv)
  { double u = uv->c[0];
    double v = uv->c[1];
    /* Check if inside unit {u,v} disk: */
    double uv2 = u*u + v*v;
    if (uv2 >= 0.999999)
      { /* Outside - return null normal: */
        return (r3_t) {{ 0.0, 0.0, 0.0 }};
      }
    else
      { /* Inside - compute {w} coordinate: */
        double w = sqrt(1.000001 - uv2);
        /* The normal is the points's position: */
        return (r3_t) {{ u, v, w }};
      }
  }

void gera_pontos_no_gabarito_old(int resolucao, r2_t centro, double raio, r2_t **pontoP, r3_t **normalP, int *num_pontosP)
{
  assert(sizeof(int) >= 4*sizeof(char));

  fprintf(stderr,  "Gerando pontos no gabarito\n");
  fprintf(stderr,  "  Resolução: %d\n", resolucao);
  fprintf(stderr,  "  Centro: %f %f\n", centro.c[0], centro.c[1]);
  fprintf(stderr,  "  Raio: %f\n", raio);
  fprintf(stderr,  "\n");

  r2_t *ponto = NULL;
  r3_t *normal = NULL; 
  int np = 0;
  int passo;
  for (passo = 0; passo < 2; passo++) {
    /* Gera pontos no disco unitário, em {M} círculos concentricos. */
    int kp = 0;
    int kr, kt;
    for (kr = 0; kr < resolucao; kr++){
      /* Escolhe o raio do círculo número {kr}: */
      /* double r = ((double) kr)/resolucao;  Igualmente espaçados em X,Y */
      /* Igualmente espaçados na esfera. */
      /* Escolhe o número de pontos {N} nesse círculo. Note que {N=1} quando {kr=0}: */
      int N = (int)floor(M_PI*kr) +1;
      if (ponto == NULL) {
        np += N; 
      } else {
      	for (kt = 0; kt < N; kt++){
      	  /* Escolhe ângulo {t} no círculo, com deslocamento alternado: */
	  double r = (kr)/(double)resolucao;
	  double t = (2*M_PI*kt)/N;
      	      	  /* Calcula normal {dx,dy,dz} da esfera nesse ponto: */
      	  double dx = r*cos(t);
      	  double dy = r*sin(t);
      	  double dz = sqrt(1 - (dx*dx) - (dy*dy));
      	  /* Calcula coordenadas {dx,dy} do ponto no disco unitário: */
      	  double qX = centro.c[0] + raio*dx;
      	  double qY = centro.c[1] + raio*dy;
      	  ponto[kp] = (r2_t){{qX, qY}}; 
      	  normal[kp] = (r3_t){{dx, dy, dz}}; 
      	  fprintf(stderr, "%d %6.3lf %6.3lf\n",kp,qX,qY);
      	  kp++;
      	} 
      }
    }
    if (ponto == NULL) {
      ponto = (r2_t *)malloc(np*sizeof(r2_t)); 
      normal = (r3_t *)malloc(np*sizeof(r3_t)); 
    } else {
      assert(kp == np); 
    }
  }
  fprintf(stderr,  "  Num pontos: %d\n", np);
  (*pontoP) = ponto;
  (*normalP) = normal;
  (*num_pontosP) = np;
}

r3_t compute_view_dir_from_stretch(r2_t gauge_stretch, double radius){
	r3_t view_dir;
	double m = r2_norm(&gauge_stretch);
	if(m != 0 ){
		double cR = radius/(radius + m);
  		double sR = sqrt(1 - cR*cR);
		view_dir = (r3_t){{ (gauge_stretch.c[0]*sR)/m, -(gauge_stretch.c[1]*sR)/m , cR}};
		(void)r3_dir(&view_dir,&view_dir);//just in case...
		
	}else{
		view_dir = (r3_t){{ 0,0,1}};
	}
	
	return view_dir;

}



void gera_normais_na_esfera_uniforme(int resolucao, r3_t** normalP, int* num_pontosP, r3_t view_dir){
	double step = (M_PI*0.5)/((double) resolucao);
	int k;
	int np = 0;
	r3_t* normal = NULL;
	int passo;
	int kp = 0;
	for(passo = 0; passo < 2; passo++){
		for(k = 0; k < resolucao; k++){
			double theta = (M_PI*0.5)*(k/(double)resolucao);
			double r = sin(theta);
			int N = (int)ceil((2.0*M_PI*r)/step);
			int i;
			if(normal == NULL){
				np+= N;
			}else{
				for(i = 0; i < (N-1); i++){
					double phi = (2.0*M_PI*i)/N;
					//double phi = (M_PI*(2*i + (k%2)))/N;
					double x,y,z;
					x = r*cos(phi);
					y = r*sin(phi);
					z = cos(theta);
					r3_t normal_bruta = (r3_t){{ x,y,z }};
					normal[kp] = normal_bruta;
					r3_print(stderr,&normal_bruta);
					fprintf(stderr,"\n");
					kp++;
				}
			}
		}
		if(normal == NULL){
			normal = (r3_t*)malloc(sizeof(r3_t)*np);
		}
	}
	*normalP = normal;
	*num_pontosP = np;
}

void gera_normais_no_gabarito(int resolucao,double thetaMax, r3_t** normalP, int* num_pontosP, r3_t view_dir){
  fprintf(stderr,  "Gerando normais no gabarito\n");
  fprintf(stderr,  "  Resolução: %d\n", resolucao);
  fprintf(stderr,  "  ThetaMax: %lf\n", thetaMax);
  fprintf(stderr,  "\n");

  r3_t *normal = NULL; 
  int np = 0;
  int passo;
  
   r3x3_t roda_normal = compute_normal_correction_matrix(view_dir);
   for (passo = 0; passo < 2; passo++) {
    	/* Gera pontos no disco unitário, em {M} círculos concentricos. */
    	int kp = 0;
    	int kr, kt;
    	for (kr = 0; kr < resolucao; kr++){
      		/* Escolhe o raio do círculo número {kr}: */
      		/* double r = ((double) kr)/resolucao;  Igualmente espaçados em X,Y */
// 		double fr = (((double)kr)/(resolucao -1))*(1.0 - (0.25/resolucao) );
//       		double r = sin(thetaMax*fr); /* Igualmente espaçados na esfera. */
      		double r = sin(thetaMax*((double) kr)/(resolucao)); /* Igualmente espaçados na esfera. */
      		/* Escolhe o número de pontos {N} nesse círculo. Note que {N=1} quando {kr=0}: */
      		int N = (int)floor(M_PI*r*resolucao) + 1;
      		if (normal == NULL) {
	        	np += N; 
      		}else {
      			for (kt = 0; kt < N; kt++){
      	  			/* Escolhe ângulo {t} no círculo, com deslocamento alternado: */
      	  			double t = (M_PI*(2*kt + (kr%2)))/N;
      	   			/* Calcula a posição {uv} relativa ao círculo unitário */
      	  			r2_t uv;
	  			uv.c[0] = r*cos(t);
      	  			uv.c[1] = r*sin(t);
      	  			r3_t normal_bruta = calcula_normal_esfera(&uv);
      	  			
				r3x3_map_row(&(normal_bruta),&roda_normal,&(normal[kp]));
				double dot = r3_dot(&view_dir,&normal[kp] );
				if(dot < -1.0e-6) {
					fprintf(stderr,"entry %d\n",kp);
					fprintf(stderr,"normal_bruta = "); r3_print(stderr,&normal_bruta); fprintf(stderr,"\n");
					fprintf(stderr,"normal_ajustada = "); r3_print(stderr,&normal[kp]); fprintf(stderr,"\n");
					assert(FALSE);
				}
      	       			//  fprintf(stderr, "%d %6.3lf %6.3lf\n",kp,qX,qY);
				kp++;
      			} 
		}
    	}
    	if (normal == NULL) {
      		normal = (r3_t *)malloc(np*sizeof(r3_t)); 
    	} else {
      		assert(kp == np); 
    	}
   }

  fprintf(stderr,  "  Num normais: %d\n", np);
  (*normalP) = normal;
  (*num_pontosP) = np;
	
}

void gera_pontos_no_gabarito_eliptico(int resolucao,double thetaMax, r2_t centro, double raio, r2_t estica, r2_t **pontoP, r3_t **normalP, int *num_pontosP,r3_t view_dir){
	
  assert(sizeof(int) >= 4*sizeof(char));

  fprintf(stderr,  "Gerando pontos no gabarito eliptico\n");
  fprintf(stderr,  "  Resolução: %d\n", resolucao);
  fprintf(stderr,  "  ThetaMax: %lf\n", thetaMax);
  fprintf(stderr,  "  Centro: %f %f\n", centro.c[0], centro.c[1]);
  fprintf(stderr,  "  Raio: %f\n", raio);
  fprintf(stderr,  "  Estica: %f %f\n", estica.c[0], estica.c[1]);
  fprintf(stderr,  "\n");

  r2_t *ponto = NULL;
  r3_t *normal = NULL; 
  int np = 0;
  
  r2_t dir_estica;
  double en = r2_norm(&estica);
  fprintf(stderr,"EN = %24.16e\n",en);
  if((estica.c[0] == 0.0) && (estica.c[1] == 0.0)){
	dir_estica = (r2_t){{1.0,0.0}};
  }
  else{
	r2_dir(&estica,&dir_estica);
  }

  gera_normais_no_gabarito(resolucao,thetaMax,&normal,&np,view_dir);
  ponto = (r2_t *)malloc(np*sizeof(r2_t)); 
  int k;
  for(k = 0; k < np; k++){
	r3_t uvw = normal[k];
	r2_t uv;
	uv.c[0] = uvw.c[0];
      	uv.c[1] = uvw.c[1];
	/* Calcula o fator de esticamento (b) */
	double b = r2_dot(&uv,&dir_estica);
	//calcula a posição q do ponto no gabarito
      	r2_t q;
	q.c[0] = (raio*uv.c[0] + b*estica.c[0]) + centro.c[0];
	q.c[1] = (raio*uv.c[1] + b*estica.c[1]) + centro.c[1];
      	ponto[k] = q; 
  }

  fprintf(stderr,  " Num pontos: %d\n", np);
  (*pontoP) = ponto;
  (*normalP) = normal;
  (*num_pontosP) = np;

		
}


void gera_pontos_com_mascara(float_image_t* mask,r2_t centro, double raio, r2_t estica, r2_t **pontoP, r3_t **normalP, int *num_pontosP,r3_t view_dir){
	
  assert(sizeof(int) >= 4*sizeof(char));

  fprintf(stderr,  "Gerando pontos no gabarito com mascara\n");
  fprintf(stderr,  "  Centro: %f %f\n", centro.c[0], centro.c[1]);
  fprintf(stderr,  "  Raio: %f\n", raio);
  fprintf(stderr,  "  Estica: %f %f\n", estica.c[0], estica.c[1]);
  fprintf(stderr,  "\n");

 
  r2_t dir_estica;
  double compr_estica = r2_norm(&estica);
  fprintf(stderr,"COMPR_ESTICA = %24.16e\n",compr_estica);
  if((estica.c[0] == 0.0) && (estica.c[1] == 0.0)){
	dir_estica = (r2_t){{1.0,0.0}};
  }
  else{
	r2_dir(&estica,&dir_estica);
  }

  r2_t dir_perp = (r2_t){{-dir_estica.c[1],dir_estica.c[0]}};

  r3x3_t roda_normal = compute_normal_correction_matrix(view_dir);
 /*conta quantos pontos a imagem possui dentro da máscara*/
  int x,y;
  
  demand(mask->sz[0] == 1,"Mask File with more than 1 channel");
  int NX = (int)mask->sz[1];
  int NY = (int)mask->sz[2];
   r2_t * ponto = NULL;
  r3_t * normal = NULL;
  int passo;
  int np = 0;
  for(passo = 0; passo < 2; passo++){
	if(passo == 1){
		ponto = (r2_t *)malloc(np*sizeof(r2_t)); 
		normal = (r3_t *)malloc(np*sizeof(r3_t));
	}
	int kp = 0;	
	bool_t debug  = FALSE;
	for(y = 0; y < NY; y++){
		for(x = 0; x < NX; x++){
			if(float_image_get_sample(mask,0,x,y) != 0){
				r2_t p = (r2_t){{x+0.5,y+0.5}};
				r2_t q; 
				r2_sub(&p,&(centro),&q);
				double u = r2_dot(&dir_estica,&q)/(raio + compr_estica);
				double v = r2_dot(&dir_perp,&q)/(raio);
				r2_t nxy;
				r2_mix(u,&(dir_estica),v,&(dir_perp),&(nxy));
				double nxy2 = r2_norm_sqr(&nxy);
				assert(nxy2 <= (1.00001+(1.42/raio)) );
				if( nxy2 <  1.00 ){
					double nz = sqrt( 1.0 - fmin(nxy2, 1.0));
					r3_t normal_bruta = (r3_t){{ nxy.c[0],nxy.c[1], nz}};
					(void)r3_dir(&normal_bruta,&normal_bruta);
					if(debug){
						fprintf(stderr,"x=%d y=%d ",x,y);
						fprintf(stderr,"NXY= %9.6f %9.6f ",nxy.c[0],nxy.c[1]);
						fprintf(stderr,"ABS= %9.6f ",sqrt(nxy2));
						r3_print(stderr,&normal_bruta);
						fprintf(stderr,"\n");
					
						debug = FALSE;
					}
					r3_t normal_rodada;
					r3x3_map_row(&(normal_bruta),&roda_normal,&(normal_rodada));
					if(passo == 0){ np++;}
					else{
						ponto[kp] = p;
						normal[kp] = normal_rodada;
						kp++;
					}
				}
			}
		}
	}
  }

  fprintf(stderr,  "  Num pontos: %d\n", np);
  (*pontoP) = ponto;
  (*normalP) = normal;
  (*num_pontosP) = np;

		
}

void gera_pontos_no_gabarito(int resolucao, r2_t centro, double raio, r2_t **pontoP, r3_t **normalP, int *num_pontosP)
{
  assert(sizeof(int) >= 4*sizeof(char));

  fprintf(stderr,  "Gerando pontos no gabarito\n");
  fprintf(stderr,  "  Resolução: %d\n", resolucao);
  fprintf(stderr,  "  Centro: %f %f\n", centro.c[0], centro.c[1]);
  fprintf(stderr,  "  Raio: %f\n", raio);
  fprintf(stderr,  "\n");

  r2_t *ponto = NULL;
  r3_t *normal = NULL; 
  int np = 0;
  int passo;
  for (passo = 0; passo < 2; passo++) {
    /* Gera pontos no disco unitário, em {M} círculos concentricos. */
    int kp = 0;
    int kr, kt;
    for (kr = 0; kr < resolucao; kr++){
      /* Escolhe o raio do círculo número {kr}: */
      /* double r = ((double) kr)/resolucao;  Igualmente espaçados em X,Y */
      double r = sin(M_PI/2*((double) kr)/resolucao); /* Igualmente espaçados na esfera. */
      /* Escolhe o número de pontos {N} nesse círculo. Note que {N=1} quando {kr=0}: */
      int N = (int)floor(M_PI*r*resolucao) + 1;
      if (ponto == NULL) {
        np += N; 
      } else {
      	for (kt = 0; kt < N; kt++){
      	  /* Escolhe ângulo {t} no círculo, com deslocamento alternado: */
      	  double t = (M_PI*(2*kt + (kr%2)))/N;
      	  /* Calcula normal {dx,dy,dz} da esfera nesse ponto: */
      	  double dx = r*cos(t);
      	  double dy = r*sin(t);
      	  double dz = sqrt(1 - (dx*dx) - (dy*dy));
      	  /* Calcula coordenadas {dx,dy} do ponto no disco unitário: */
      	  double qX = centro.c[0] + raio*dx;
      	  double qY = centro.c[1] + raio*dy;
      	  ponto[kp] = (r2_t){{qX, qY}}; 
      	  normal[kp] = (r3_t){{dx, dy, dz}}; 
      	//  fprintf(stderr, "%d %6.3lf %6.3lf\n",kp,qX,qY);
      	  kp++;
      	} 
      }
    }
    if (ponto == NULL) {
      ponto = (r2_t *)malloc(np*sizeof(r2_t)); 
      normal = (r3_t *)malloc(np*sizeof(r3_t)); 
    } else {
      assert(kp == np); 
    }
  }
  fprintf(stderr,  "  Num pontos: %d\n", np);
  (*pontoP) = ponto;
  (*normalP) = normal;
  (*num_pontosP) = np;
}

int ponto_dentro_do_gabarito(double radius, r2_t stretch,r2_t dp);
int ponto_dentro_do_gabarito(double radius, r2_t stretch,r2_t dp){
	r2_t e; //direção do stretch
  	double en = r2_norm(&stretch);
  	if(en == 0.0){
		e = (r2_t){{1.0,0.0}};
  	}
  	else{
		r2_dir(&stretch,&e);
  	}
	
	r2_t f = (r2_t){{-e.c[1],e.c[0]}}; //perpendicular ao {e}
	double dpe = r2_dot(&dp,&e)/radius;
	double dpf = r2_dot(&dp,&f)/(radius + en);
	 
	return ((dpe*dpe + dpf*dpf) <= 1);
	
  
}

float_image_t  *gera_mascara_do_gabarito(int nx, int ny, r2_t centro, double raio,r2_t stretch)
{
  float_image_t  *M = float_image_new(1, nx, ny);
  int x, y;
  for (x = 0; x < nx; x++) {
    for (y = 0; y < ny; y++) {
      /* Determina o centro {px,py} do pixel {x,y}: */
      double px = x + 0.5; 
      double py = y + 0.5;
      /* Determina o ponto {px,py} do pixel que está mais longe do centro do gabarito: */
      if (px < centro.c[0]) { px -= 0.5; } else { px += 0.5; }
      if (py < centro.c[1]) { py -= 0.5; } else { py += 0.5; }
      /* Decide se o pont está dentro do gabarito: */
      double dx = px - centro.c[0];
      double dy = py - centro.c[1];
      //int dentro = (dx*dx + dy*dy <= raio*raio);
      int dentro = ponto_dentro_do_gabarito(raio,stretch,(r2_t){{dx,dy}});
      /* Guarda na máscara: */
      float v = (dentro ? 1.0 : 0.0);
      float_image_set_sample(M, 0, x, y, v);
    }
  }
  return M;
}

Tabela* aloca_tabela_vazia(int num_luzes, int num_linhas, r3_t view_dir){
	Tabela* tab = (Tabela*)malloc(sizeof(Tabela));
  	tab->tabela = (LinhaTabela*)malloc(sizeof(LinhaTabela)*num_linhas);
	tab->num_linhas = num_linhas;
  	tab->num_luzes = num_luzes;
  	tab->view_dir = view_dir;
	int i;
	for (i = 0; i < num_linhas; i++){
		tab->tabela[i].intdir = (double*)malloc(sizeof(double)*num_luzes);
	}
	return tab;
}

Tabela* cria_tabela
  ( float_image_t  *G[],
    float_image_t  *M, 
    double albedo, 
    int num_luzes, 
    int canal, 
    r2_t ponto[], 
    r3_t normal[], 
    int num_pontos,
    r3_t view_dir,
    bool_t interpolate_pixels
  )
{

  fprintf(stderr,  "Criando tabela fotométrica\n");
  fprintf(stderr,  "  Num luzes: %d\n", num_luzes);
  fprintf(stderr,  "  Canal: %d\n", canal);
  fprintf(stderr,  "  Num pontos dados: %d\n", num_pontos);
  fprintf(stderr,  "  View Dir = "); r3_print(stderr,&view_dir); fprintf(stderr,"\n");

  Tabela* tab = (Tabela*)malloc(sizeof(Tabela));
  tab->tabela = (LinhaTabela*)malloc(sizeof(LinhaTabela)*num_pontos);

  tab->num_luzes = num_luzes;
  tab->view_dir = view_dir;
  /* Guarda vetores de observação úteis na tabela: */
  int k, luz;
  int count_minus =0;
  int count_norm =0;
  int count_plus = 0;
  int count_wzero = 0;
  int count_small = 0;
  int num_ok = 0;  /* Conta pontos OK. */
  int num_rej = 0; /* Conta pontos rejeitados. */
  for (k = 0; k < num_pontos; k++){
    double GO[num_luzes]; /* Vetor de observações. */
    double GW[num_luzes]; /* Vetor de pesos. */
    double Gmag2 = 0;
    r2_t *q = &(ponto[k]);
    double x = q->c[0];
    double y = q->c[1];
    r3_t *u = &(normal[k]); 
    //normal must be a unitary vector
    r3_dir(u,u);
    assert(!isnan(u->c[0]));
    int nbugs = 0; /* Numero de erros neste ponto. */
    for(luz = 0; luz < num_luzes; luz++){
      float_image_masked_t GM = (float_image_masked_t){ .img = G[luz], .msk = M };
      float val, wht;
      if(interpolate_pixels){
      	float_image_masked_interpolate(&GM, canal, x, y, 2, &val, &wht);
      }else{
	int ix = (int)floor(x);
	int iy = (int)floor(y);
	val = float_image_get_sample(G[luz],canal,ix, iy);
	wht = float_image_get_sample(M,0,ix, iy);
      }
     // float_image_masked_interpolate_exclusive(&GM, canal, x, y, 1, &val, &wht);
      GO[luz] = val; GW[luz] = wht;
      if (GW[luz] == 0.0) { count_wzero++; nbugs++; }
      if (GO[luz] < 0.0) { count_minus++; nbugs++; }
      if (GO[luz] > 1.0) { count_plus++; nbugs++; }
      if (u->c[2] < 0.0) { count_norm++; nbugs++; }
      Gmag2 += GO[luz]*GO[luz];
    }
    double Gmag = sqrt(Gmag2);
    if (Gmag < MINMAG) { count_small++; nbugs++; }
    if ((nbugs > 0) && (count_wzero+count_minus+count_plus+count_small == nbugs)) {
      /* Primeiro erro nesta tabela. */
      fprintf(stderr, "pixel problematico no gabarito");
      fprintf(stderr, "  luz = %2d  q = (%7.2f %7.2f)\n", luz, x, y);
      for(luz = 0; luz < num_luzes; luz++){
        fprintf(stderr, "  GO[%2d] = %+9.6f wht = %12.6e\n", luz, GO[luz], GW[luz]);
      }
    }
    
    if (nbugs == 0) {
      /* Aloca vetor de observação normalizado: */
      double *go = (double*)malloc(sizeof(double)*num_luzes);
      for(luz = 0; luz < num_luzes; luz++){ go[luz] = GO[luz]/Gmag; }
      tab->tabela[num_ok].intdir = go;
      tab->tabela[num_ok].intmag = Gmag / albedo;
      tab->tabela[num_ok].normal = (*u);
      num_ok++;
    } else {
      num_rej++;
    }
  }
  tab->num_linhas = num_ok;
  fprintf(stderr,  "  Num pontos aceitos: %d\n", num_ok);
  fprintf(stderr,  "  Num pontos rejeitados: %d\n", num_rej);
  fprintf(stderr,  "  Num pontos com intensidade abaixo de zero: %d\n",count_minus);
  fprintf(stderr,  "  Num pontos com intensidade acima de um: %d\n",count_plus);
  fprintf(stderr,  "  Num pontos com normal Z negativa: %d\n",count_norm);
  fprintf(stderr,  "\n");
  return tab;
}


Tabela* cria_tabela_virtual
  ( r3_t  directions[],
    double radius,r2_t centro,
    double albedo, 
    int num_luzes, 
    int canal, 
    r2_t ponto[], 
    r3_t normal[], 
    int num_pontos,
    r3_t view_dir
  )
{

  fprintf(stderr,  "Criando tabela fotométrica virtual\n");
  fprintf(stderr,  "  Num luzes: %d\n", num_luzes);
  fprintf(stderr,  "  Canal: %d\n", canal);
  fprintf(stderr,  "  Num pontos dados: %d\n", num_pontos);
  fprintf(stderr,  "  View Dir = "); r3_print(stderr,&view_dir); fprintf(stderr,"\n");
 
  Tabela* tab = (Tabela*)malloc(sizeof(Tabela));
  tab->tabela = (LinhaTabela*)malloc(sizeof(LinhaTabela)*num_pontos);

  tab->num_luzes = num_luzes;
  tab->view_dir = view_dir;
  /* Guarda vetores de observação úteis na tabela: */
  int k, luz;
  int count_minus =0;
  int count_plus = 0;
  int count_wzero = 0;
  int count_small = 0;
  int num_ok = 0;  /* Conta pontos OK. */
  int num_rej = 0; /* Conta pontos rejeitados. */
  /*entorta direcoes da luz*/
  r3_t true_directions[num_luzes];
  r3x3_t roda_normal = compute_normal_correction_matrix(view_dir);
  for(luz = 0; luz < num_luzes;luz++){
	r3x3_map_row(&(directions[luz]),&roda_normal,&(true_directions[luz]));
  }
  

  for (k = 0; k < num_pontos; k++){
    double GO[num_luzes]; /* Vetor de observações. */
    double GW[num_luzes]; /* Vetor de pesos. */
    double Gmag2 = 0;
    r2_t *q = &(ponto[k]);
    double x = q->c[0];
    double y = q->c[1];
    r3_t *u = &(normal[k]); 
    int nbugs = 0; /* Numero de erros neste ponto. */
    for(luz = 0; luz < num_luzes; luz++){
      float val, wht;
      //float_image_masked_interpolate(&GM, canal, x, y, 2, &val, &wht);
      //val = virtual_gab_intensity(x,y,radius,cor.c[0],cor.c[1], cor.c[2], centro.c[0], centro.c[1]);
      val = (float)lambertian_shading(true_directions[luz],albedo,*u);
      wht = 1.0;
     // float_image_masked_interpolate_exclusive(&GM, canal, x, y, 1, &val, &wht);
      GO[luz] = val; GW[luz] = wht;
      if (GW[luz] == 0.0) { count_wzero++; nbugs++; }
      if (GO[luz] < 0.0) { count_minus++; nbugs++; }
      if (GO[luz] > 1.0) { count_plus++; nbugs++; }
      Gmag2 += GO[luz]*GO[luz];
    }
    double Gmag = sqrt(Gmag2);
    if (Gmag < MINMAG) { count_small++; nbugs++; }
    if ((nbugs > 0) && (count_wzero+count_minus+count_plus+count_small == nbugs)) {
      /* Primeiro erro nesta tabela. */
      fprintf(stderr, "pixel problematico no gabarito");
      fprintf(stderr, "  luz = %2d  q = (%7.2f %7.2f)\n", luz, x, y);
      for(luz = 0; luz < num_luzes; luz++){
        fprintf(stderr, "  GO[%2d] = %+9.6f wht = %12.6e\n", luz, GO[luz], GW[luz]);
      }
    }
    
    if (nbugs == 0) {
      /* Aloca vetor de observação normalizado: */
      double *go = (double*)malloc(sizeof(double)*num_luzes);
      for(luz = 0; luz < num_luzes; luz++){ go[luz] = GO[luz]/Gmag; }
      tab->tabela[num_ok].intdir = go;
      tab->tabela[num_ok].intmag = Gmag / albedo;
      tab->tabela[num_ok].normal = (*u);
      num_ok++;
    } else {
      num_rej++;
    }
  }
  tab->num_linhas = num_ok;
  fprintf(stderr,  "  Num pontos aceitos: %d\n", num_ok);
  fprintf(stderr,  "  Num pontos rejeitados: %d\n", num_rej);
  fprintf(stderr,  "  Num pontos intensidade abaixo de zero: %d\n",count_minus);
  fprintf(stderr,  "  Num pontos intensidade acima de um: %d\n",count_plus);
  fprintf(stderr,  "\n");
  return tab;
}


Tabela* criaSubTabela(Tabela* tab,int subsetSize,int* subSets, int** index){
	Tabela* nova = (Tabela*)malloc(sizeof(Tabela));
        nova->num_luzes = subsetSize;
	//first count how lines would be used
	int num_linhas = tab->num_linhas;
	int i,valid_lines;
	valid_lines = 0;
	int is_valid_line[num_linhas];
	fprintf(stderr,  "Criando subTabela fotométrica a partir da Tabela principal\n");
	
	for(i = 0; i < num_linhas;i++){
		double GO[subsetSize];
		double Gmag = 0;
		const double* goTab = get_intdir(tab, i);
		double gmagTab = get_intmag(tab, i);
		int j;
		for(j = 0; j < subsetSize; j++){
			int ind = subSets[j];
			GO[j] = goTab[ind]*gmagTab;
			Gmag+= GO[j]*GO[j];
		}
		Gmag = sqrt(Gmag);
		is_valid_line[i] = 0;
		if (Gmag > MINMAG) { valid_lines++; is_valid_line[i] = 1;}
	}
	fprintf(stderr,  "Linhas selecionadas: %d de %d\n",valid_lines,num_linhas);
	fprintf(stderr,  "Preenchendo estrutura da tabela\n");
	//now we fill our table
	*index = (int*)malloc(sizeof(int)*valid_lines);
	nova->num_linhas = valid_lines;
	nova->tabela = (LinhaTabela*)malloc(sizeof(LinhaTabela)*valid_lines);
        nova->view_dir = tab->view_dir;
	int k = 0;
	for(i = 0; i < num_linhas; i++){
		if(is_valid_line[i] == 1){
			double* go = (double*)malloc(sizeof(double)*subsetSize);
			double gmag = 0;
			//get GO and gmag from Table
			const double* goTab = get_intdir(tab, i);
			double gmagTab = get_intmag(tab, i);
			int j;
			for(j = 0; j < subsetSize; j++){
				int ind = subSets[j];
				go[j] = goTab[ind]*gmagTab;
				gmag+= go[j]*go[j];
			}
			gmag = sqrt(gmag);
			//normalize vector
			for(j = 0; j < subsetSize; j++){ go[j] = go[j]/gmag; }
			//store data 
			nova->tabela[k].intdir = go;
			nova->tabela[k].intmag = gmag;
			nova->tabela[k].normal = get_normal(tab, i);
			(*index)[k] = i;
			k++;
		}
	}
  	fprintf(stderr,  "SubTabela fotométrica gerada\n");
	return nova;
	
}




double distPoints(double x1, double x2, double y1, double y2);
double get_intdir_selecionado(Tabela* tab,int num_linha,int arquivo,int canal);
void normaliza_assinatura(double **intdir, double *intmag, int num_luzes);
double ObservationMagnitude(double** intdir,int num_luzes,int canal);

double distPoints(double x1, double x2, double y1, double y2){
  double X = (x1 - x2)*(x1 - x2);
  double Y = (y1 - y2)*(y1 - y2);
  return sqrt(X+Y);
}

/* FUNÇÕES DE INTERFACE DA TABELA*/

const double *get_intdir(Tabela* tab, int linha)
{
  return tab->tabela[linha].intdir;
}

double get_intmag(Tabela* tab, int linha)
{
  return tab->tabela[linha].intmag;
}


void set_intdir(Tabela* tab, int linha,double g[]){
	int i;
	for(i = 0; i < tab->num_luzes; i++){
		tab->tabela[linha].intdir[i] = g[i];
	}
}

void set_intmag(Tabela* tab, int linha,double Gmag){
	tab->tabela[linha].intmag = Gmag;
}

r3_t get_view_dir(Tabela* tab){
	return tab->view_dir;
}

void print_linha(Tabela* tab, int linha)
{
  int i;
  fprintf(stderr,  "Linha %d Assinatura", linha);
  const double *go =tab->tabela[linha].intdir;
  for (i = 0; i < tab->num_luzes; i++) { fprintf(stderr,  " %f", go[i]); }
  fprintf(stderr,  " Magnitude %f \n", tab->tabela[linha].intmag);
	
}

int get_num_linhas(Tabela* tab)
{
  return (tab == NULL ? 0 : tab->num_linhas);
}

void set_num_linhas(Tabela* tab,int num_linhas)
{
	tab->num_linhas = num_linhas; 
}


int get_num_luzes(Tabela* tab)
{
  return tab->num_luzes;
}

r3_t get_normal(Tabela* tab, int linha)
{
  return tab->tabela[linha].normal;
}

void set_normal(Tabela* tab, int linha, r3_t normal){
	tab->tabela[linha].normal = normal;	
}

double ObservationMagnitude(double** intdir,int num_luzes,int canal){
	int i;
	double soma = 0.000001;
    	for(i = 0; i < num_luzes; i++){
      		double vi = intdir[i][canal];
      		soma += vi*vi;
    	}
    	return sqrt(soma);
    	
	
}

void normaliza_assinatura(double **intdir, double *intmag, int num_luzes){
  int canal;
  for (canal = 0; canal < 3; canal++){
    double soma = 0.000001;
    int i;
    for(i = 0; i < num_luzes; i++){
      double vi = intdir[i][canal];
      soma += vi*vi;
    }
    intmag[canal] = sqrt(soma);
    for(i = 0; i < num_luzes; i++){
      intdir[i][canal] /= intmag[canal];
    }
  }
}

int localiza_linha_por_normal(Tabela* tab, r3_t *sn)
{
  int tam = get_num_linhas(tab);
  int imin = -1;
  int i;
  double distmin = HUGE_VAL;
  for(i = 0; i < tam; i++) {
    double dist;
    r3_t gn = get_normal(tab,i);
    dist = r3_dist(&gn, sn);
    if (dist < distmin) { distmin = dist; imin = i;	}
		
  }
  return imin;
}

void ShowTableData(char* prefix, Tabela* tab,r2_t ponto[], r3_t normal[], int num_pontos){
	char* filename = NULL;
	char *filename = jsprintf("%s_TableData.txt",prefix);
	FILE* arq;
	arq = fopen(filename,"wt");
	if(arq == NULL){
		fprintf(stderr,"showTableData %s FAILED !\n",filename);
		return;
	}
	fprintf(arq,"# Numero linhas\n");
	fprintf(arq,"%d\n",tab->num_linhas);
	fprintf(arq,"# Numero luzes\n");
	fprintf(arq,"%d\n",tab->num_luzes);
	fprintf(arq,"#(x,y) (normal) (intmag) (intdir)\n");
	int i;
	int linha = 0;
	for(i = 0; i < num_pontos; i++){
		int j;
		int test = 1;
		for(j = 0; j < 3; j++){
			if( normal[i].c[j] != tab->tabela[linha].normal.c[j] ){
				test = 0;
			}
		}
		if(test){
			fprintf(arq,"%9.6f %9.6f     ",ponto[i].c[0], ponto[i].c[1]);
			for(j = 0; j < 3; j++){
				fprintf(arq,"%9.6f ",tab->tabela[linha].normal.c[j]);
			}
			fprintf(arq,"    %9.6f    ",tab->tabela[linha].intmag);
			for(j = 0; j < tab->num_luzes;j++){
				fprintf(arq,"%9.6f ",tab->tabela[linha].intdir[j]);
			}
			fprintf(arq,"\n");
			linha++;
		}
	}
	fclose(arq);

}

void SaveTable(char* filename,Tabela* tab, bool_t check){
	FILE* arq;
	arq = fopen(filename,"wt");
	if(arq == NULL){
		fprintf(stderr,"SaveTable %s FAILED !\n",filename);
		return;
	}
	fprintf(arq,"# Numero linhas\n");
	fprintf(arq,"%d\n",tab->num_linhas);
	fprintf(arq,"# Numero luzes\n");
	fprintf(arq,"%d\n",tab->num_luzes);
        fprintf(arq,"# View dir\n");
	r3_gen_print (arq, &tab->view_dir, "%+10.7f", ""," ", "\n");
	fprintf(arq,"#(x,y) (normal) (intmag) (intdir)\n");
	int i;
	bool_t bug = FALSE;
	for(i = 0; i < tab->num_linhas; i++){
	  const double *intdir = tab->tabela[i].intdir;
	  r3_t *normal = &(tab->tabela[i].normal);
	  double intmag = tab->tabela[i].intmag;

	  int j;
	  for(j = 0; j < 3; j++){
	    fprintf(arq,"%e ",normal->c[j]);
	  }
	  fprintf(arq,"    %e    ",intmag);
	  for(j = 0; j < tab->num_luzes;j++){
	    fprintf(arq,"%e ",intdir[j]);
	  }
	  fprintf(arq,"\n");
	  if (check){
	    // Check normal:
	    double len = r3_norm(normal);
	    if (fabs(len - 1.0) > 0.000001){
	      fprintf(stderr, "normal %d is not normalized\n", i); bug = TRUE;
	    }
	    double dot = r3_dot(&tab->view_dir,normal);
	    if(dot < -1.0e-6) {
	      fprintf(stderr, "normal %d is not visible\n", i); bug = TRUE;
	    }
	    // Check signature:
            for(j = 0; j < tab->num_luzes;j++){
	      if (intdir[j] < 0.0){
	      fprintf(stderr, "signature %d light %d is negative\n", i, j); bug = TRUE;
	      }
	    }
            double mag = rn_norm(tab->num_luzes, (double *)intdir);
            if (fabs(mag - 1.0) > 0.000001){
	      fprintf(stderr, "signature %d is not normalized\n", i); bug = TRUE;
	    }
	    if (intmag < 0.000001){
	      fprintf(stderr, "magniture %d = %9.6f is not positive\n", i,intmag); bug = TRUE;
	    }
            fflush(arq);
	    assert(! bug);
	  }
	}
	fclose(arq);
}


void loadTableData(char* filename, Tabela** table,r2_t** pontos_P, r3_t** normals_P ){
	
	auto void advanceUntilEqual(FILE* arq,char ch);
	void advanceUntilEqual(FILE* arq,char ch){
		int c;
		c = fgetc(arq);
		while( c != EOF){
                  if((char)c == ch) break;
			c = fgetc(arq);
		}
	}
	
	FILE* arq;
	arq = fopen(filename,"rt");
	if(arq == NULL){
		*table = NULL;
		*pontos_P = NULL;
		*normals_P = NULL;
		return ;
	}
	Tabela* tab;
	r2_t* pontos;
	r3_t* normals;
	tab = (Tabela*)malloc(sizeof(Tabela));
	advanceUntilEqual(arq,'\n');
	fscanf(arq,"%d",&(tab->num_linhas));	
 	advanceUntilEqual(arq,'#');
	advanceUntilEqual(arq,'\n');
	fscanf(arq,"%d",&(tab->num_luzes));
	advanceUntilEqual(arq,'#');
	advanceUntilEqual(arq,'\n');
	pontos = (r2_t*)malloc(sizeof(r2_t)*(tab->num_linhas));
	normals = (r3_t*)malloc(sizeof(r3_t)*(tab->num_linhas));
	tab->tabela = (LinhaTabela*)malloc(sizeof(LinhaTabela)*(tab->num_linhas));
	
	int i;
	for(i = 0; i < tab->num_linhas; i++){
		r2_t pt;
		r3_t norm;
		double mag;		

		fscanf(arq,"%lf %lf",&pt.c[0],&pt.c[1]);
		pontos[i] = pt;
		fscanf(arq,"%lf %lf %lf",&norm.c[0],&norm.c[1],&norm.c[2]);
		normals[i] = norm;
		fscanf(arq,"%lf",&mag);
		tab->tabela[i].intdir = (double*)malloc(sizeof(double)*(tab->num_luzes));
		tab->tabela[i].normal = norm;
		tab->tabela[i].intmag = mag;
		int j;
		for(j = 0; j < tab->num_luzes; j++){
			fscanf(arq,"%lf",&(tab->tabela[i].intdir[j]));
		}
	}
	fclose(arq);
	*table = tab;
	*pontos_P = pontos;
	*normals_P = normals;
	return ;
	
}

void LoadTable(char* filename, Tabela** table){
	
	auto void advanceUntilEqual(FILE* arq,char ch);
	void advanceUntilEqual(FILE* arq,char ch){
		int c = fgetc(arq);
		while( c != EOF){
                  if((char)c == ch) break;
			c = fgetc(arq);
		}
	}
	
	FILE* arq;
	arq = fopen(filename,"rt");
	if(arq == NULL){
		*table = NULL;
		return ;
	}
	Tabela* tab;
	tab = (Tabela*)malloc(sizeof(Tabela));
	advanceUntilEqual(arq,'\n');
	fscanf(arq,"%d",&(tab->num_linhas));	
 	advanceUntilEqual(arq,'#');
	advanceUntilEqual(arq,'\n');
	fscanf(arq,"%d",&(tab->num_luzes));
        advanceUntilEqual(arq,'#');
	advanceUntilEqual(arq,'\n');
	int res = fscanf(arq,"%lf %lf %lf",&(tab->view_dir.c[0]),&(tab->view_dir.c[1]),&(tab->view_dir.c[2]) );
	assert(res == 3);
	advanceUntilEqual(arq,'#');
	advanceUntilEqual(arq,'\n');
	tab->tabela = (LinhaTabela*)malloc(sizeof(LinhaTabela)*(tab->num_linhas));
	
	int i;
	for(i = 0; i < tab->num_linhas; i++){
		r3_t norm;
		double mag;		
		fscanf(arq,"%lf %lf %lf",&norm.c[0],&norm.c[1],&norm.c[2]);
		fscanf(arq,"%lf",&mag);
		tab->tabela[i].intdir = (double*)malloc(sizeof(double)*(tab->num_luzes));
		tab->tabela[i].normal = norm;
		tab->tabela[i].intmag = mag;
		int j;
		for(j = 0; j < tab->num_luzes; j++){
			fscanf(arq,"%lf",&(tab->tabela[i].intdir[j]));
		}
	}
	fclose(arq);
	*table = tab;
	return ;
	
}

void PrintTableStats(Tabela* tab){
	fprintf(stderr,"Lights: %d\nLines: %d\n View Dir (%lf,%lf,%lf)\n",tab->num_luzes, tab->num_linhas,tab->view_dir.c[0],tab->view_dir.c[1],tab->view_dir.c[2]);
}

void LiberaTabela(Tabela* tab){
	int i ;
	for( i = 0; i < tab->num_linhas; i++){
		free(tab->tabela[i].intdir);
	}
	free(tab->tabela);
	free(tab);
	
}

Tabela* FixTableData(Tabela* tab){
  int n = get_num_linhas(tab);
  int m = get_num_luzes(tab);
  int l;
  bool_t valid_lines[n];
  int count_good_lines = 0;
  int not_norm = 0;
  int not_vis = 0;
  int neg_sig = 0;
  int not_norm_sig = 0;
  int neg_mag = 0;
  for(l = 0; l < n; l++){
    const double *intdir = tab->tabela[l].intdir;
    r3_t *normal = &(tab->tabela[l].normal);
    double intmag = tab->tabela[l].intmag;
    //at first we consider this line correct
    bool_t bug = FALSE;

    //case 1 - not normalized normal
    double len = r3_norm(normal);
    if (fabs(len - 1.0) > 0.000001){
    //  fprintf(stderr, "normal %d is not normalized\n", l);
      not_norm++;
      bug = TRUE;
    }

    //case 2 - not visible normal from view direction
    double dot = r3_dot(&tab->view_dir,normal);
    if(dot < -1.0e-6) {
 //     fprintf(stderr, "normal %d is not visible\n", l);
      not_vis++;
      bug = TRUE;
    }
    
    //case 3 -  Check signature anomalies:
    int j;
    for(j = 0; j < m;j++){
      if (intdir[j] < 0.0){
	//fprintf(stderr, "signature %d light %d is negative\n", l, j);
	neg_sig++;
	bug = TRUE;
      }
    }

    double mag = rn_norm(m, (double *)intdir);
    if (fabs(mag - 1.0) > 0.000001){
   //   fprintf(stderr, "signature %d is not normalized\n", l);
     not_norm_sig++;
      bug = TRUE;
    }
    if (intmag < 0.000001){
     // fprintf(stderr, "magniture %d = %9.6f is not positive\n", l,intmag);
      neg_mag++;
      bug = TRUE;
    }

    valid_lines[l] = !bug;
    if(!bug) count_good_lines++;
  }

  //now we allocate the fixed table
  Tabela* fixed_tab = aloca_tabela_vazia(m,count_good_lines,tab->view_dir);
  int count_l = 0;
  for(l = 0; l < n; l++){
    if(valid_lines[l]){
      double *intdir = tab->tabela[l].intdir;
      r3_t normal = (tab->tabela[l].normal);
      double intmag = tab->tabela[l].intmag;
      
      set_intdir(fixed_tab,count_l,intdir);
      set_normal(fixed_tab,count_l,normal);
      set_intmag(fixed_tab,count_l,intmag);
      count_l++;
    }
  }

  fprintf(stderr,"Fixed Table - %d lines, %d discarded \n",count_l, n-count_l);
  fprintf(stderr,"Not Normalized Normals: %d\n",not_norm);
  fprintf(stderr,"Not Visible Normals: %d\n",not_vis);
  fprintf(stderr,"Negative Signatures: %d\n",neg_sig);	
  fprintf(stderr,"Not Normalized Signatures: %d\n",not_norm_sig);
  fprintf(stderr,"Negative Magnitudes: %d\n",neg_mag );
  return fixed_tab;
  
}

Tabela* cria_tabela_from_model(int num_lights,int resolution,r3_t view_dir,double thetaMax, approx_model_t* lm,void** l_data){
  
  r2_t* points;
  r3_t* normals;
  int num_lines;
  
  r2_t gauge_center = (r2_t){{0,0}};
  r2_t gauge_stretch = (r2_t){{0,0}};
  double radius = 100;
  
  gera_pontos_no_gabarito_eliptico(resolution,thetaMax,gauge_center,radius,gauge_stretch,&points,&normals,&num_lines,view_dir);
  Tabela* tab = aloca_tabela_vazia(num_lights,num_lines,view_dir);
  int i;
  for( i = 0; i < num_lines; i++){
    set_normal(tab,i,normals[i]);
    double GO[num_lights];
    int l;
    for(l = 0; l < num_lights; l++){
      GO[l] = lm->evaluate(normals[i].c,l_data[l]);
      if (GO[l] < 0) GO[l] = 0;
    }
    double go[num_lights];
    double Gmag = rn_dir(num_lights,GO,go);
    set_intmag(tab,i,Gmag);
    set_intdir(tab,i,go);
  }
   Tabela* fixed_tab = FixTableData(tab);
   LiberaTabela(tab);
  free(normals);
  free(points);
  
   return fixed_tab;

}

void check_dup(argparser_t *pp, bool_t dup, char *key);
  /* If {dup} is true, prints an error "duplicate keyword {key}" and aborts. */
  
void check_dup(argparser_t *pp, bool_t dup, char *key)
  { if (dup)
      { char *msg = NULL;
        char *msg = jsprintf("duplicate keyword \"%s\"", key);
        argparser_error(pp, msg);
        free(msg);
      }
  }

gauge_data_t *parse_gauge_args(argparser_t *pp, bool_t input)
  {
    gauge_data_t *ga = (gauge_data_t*)malloc(sizeof(gauge_data_t));
    ellipse_crs_t *E = &(ga->E);
    
    /* Initialize fields to NAN: */
    E->ctr = (r2_t){{ NAN, NAN }};
    E->rad = NAN;
    E->str = (r2_t){{ NAN, NAN }};
    
    ga->image = NULL;
    ga->view = (r3_t){{ NAN, NAN, NAN }};
    int c; 
    for (c = 0; c < MAX_NC; c++) { ga->albedo[c] = NAN; }
    
    ga->magnify = -1;
    ga->margin = -1;
    
    ga->trim = NAN;
    ga->gamma = NAN;
    
    /* Parse fields of a gauge spec: */
    bool_t done;
    do 
      { 
        done = FALSE;
        if (argparser_keyword_present_next(pp, "image"))
          { check_dup(pp, ga->image != NULL, "image");
            ga->image = argparser_get_next(pp);
          }
        else if (argparser_keyword_present_next(pp, "albedo"))
          {
            check_dup(pp, ! isnan(ga->albedo[0]), "albedo");
            ga->albedo[0] = argparser_get_next_double(pp, 0.000, 1000.0);
            int c = 1;
            while ((c < MAX_NC) && argparser_next_is_number(pp))
              { ga->albedo[c] = argparser_get_next_double(pp, 0.000, 1000.0);
                c++;
              }
            /* Replicate last value through all remaining channels: */
            while (c < MAX_NC) { ga->albedo[c] = ga->albedo[c-1]; c++; }
          }
        else if (input)
          { /* Fields allowed only in input gauge spec: */
            if (argparser_keyword_present_next(pp, "center"))
              { check_dup(pp, ! isnan(E->ctr.c[0]), "center");
                E->ctr.c[0] = argparser_get_next_double(pp, -100000, 100000.0);
                E->ctr.c[1] = argparser_get_next_double(pp, -100000, 100000.0);
              }
            else if (argparser_keyword_present_next(pp, "radius"))
              { check_dup(pp, ! isnan(E->rad), "radius");
                E->rad = argparser_get_next_double(pp, 0.0, 100000.0);
              }
            else if (argparser_keyword_present_next(pp, "stretch"))
              { check_dup(pp, ! isnan(E->str.c[0]), "stretch");
                E->str.c[0] = argparser_get_next_double(pp, -100000.0, 100000.0);
                E->str.c[1] = argparser_get_next_double(pp, -100000.0, 100000.0);
              }
            else if (argparser_keyword_present_next(pp, "view"))
              { check_dup(pp, ! isnan(ga->view.c[0]), "view");
                ga->view.c[0] = argparser_get_next_double(pp, -100000.0, +100000.0);
                ga->view.c[1] = argparser_get_next_double(pp, -100000.0, +100000.0);
                ga->view.c[2] = argparser_get_next_double(pp, -100000.0, +100000.0);
                /* Normalize to unit length: */
                (void)r3_dir(&(ga->view), &(ga->view));
              }
            else if( argparser_keyword_present_next(pp, "trim")){
	        check_dup(pp, ! isnan(ga->trim), "trim");
		ga->trim =  argparser_get_next_double(pp, 0.0, +100000.0);
	    }
	    else if( argparser_keyword_present_next(pp, "gamma")){
	        check_dup(pp, ! isnan(ga->gamma), "gamma");
		ga->gamma =  argparser_get_next_double(pp, 0.01, +100.0);
	    }
            else
              { done = TRUE; }
          }
        else
          { /* Fields allowed only in output gauge spec: */
            if (argparser_keyword_present_next(pp, "magnify"))
              {
                check_dup(pp, ga->magnify > 0, "magnify");
                ga->magnify = (int)argparser_get_next_int(pp, 1, 64);
              }
            else if (argparser_keyword_present_next(pp, "margin"))
              {
                check_dup(pp, ga->margin >= 0, "margin");
                ga->margin = (int)argparser_get_next_int(pp, 0, 1024);
              }
            else
              { done = TRUE; }
          }
      }
    while (! done);
      
    /* Check for required fields and provide defaults: */
    if (ga->image == NULL) { argparser_error(pp, "missing \"image\" for gauge."); }
    if (isnan(ga->albedo[0])) 
      { for (c = 0; c < MAX_NC; c++) { ga->albedo[c] = 1.0; } }
    
    if (input)
      { /* Fields allowed/required only in input gauge spec: */
        if (isnan(E->rad)) { argparser_error(pp, "missing \"radius\" for gauge."); }
        if (isnan(E->ctr.c[0])) { argparser_error(pp, "missing \"center\" for gauge."); }
        if (isnan(E->str.c[0])) { E->str = (r2_t) {{ 0.0, 0.0}}; }
        if (isnan(ga->view.c[0])) { ga->view = (r3_t) {{ 0.0, 0.0, 1.0}}; }
        if (isnan(ga->trim)) { ga->trim = 0 ; }
        if (isnan(ga->gamma)) { ga->gamma = 1 ; }
      }
    else
      { /* Fields allowed/required only in output gauge spec: */
        if (ga->magnify < 0) { ga->magnify = 1; }
        if (ga->margin < 0) { ga->margin = 0; }
      }
      
    return ga;
  }
  

r3_t gauge_normal_at_point(r2_t *q, ellipse_crs_t *E)
  {
    r2_t w = ellipse_crs_relative_coords(E, q);
    double r2 = r2_norm_sqr(&w);
    demand(r2 <= 1.0, "point is not on gauge");
    r3_t nrm = (r3_t){{ w.c[0], w.c[1], sqrt(1 - r2) }};
    return nrm;
  }

    
r3_t gauge_normal_in_pixel(int x, int y, ellipse_crs_t *E)
  {
    /* Find the mean normal by sampling: */
    int NS = 5;  /* Pixel subsampling points along each axis for normal averaging. */
    int kx, ky;
    r3_t nrm = (r3_t){{ 0, 0, 0 }}; /* Sum of normal vectors. */
    for (kx = 0; kx < NS; kx++)
      { for (ky = 0; ky < NS; ky++)
          { /* Pick a subsampling point in pixel: */
            r2_t q = (r2_t){{ x + (kx + 0.5)/NS, y + (ky + 0.5)/NS }};
            assert(ellipse_crs_inside(E, &q));
            r3_t nq = gauge_normal_at_point(&q, E);
            assert(! isnan(nq.c[0]));
            r3_add(&nq, &nrm, &nrm);
          }
      }
    r3_scale(1.0/(NS*NS), &nrm, &nrm);
    double m = r3_dir(&nrm, &nrm);
    assert(m > 0);
    return nrm;
  }


int pixel_position_in_gauge_image(int x, int y, ellipse_crs_t *E)
  {
    double debug = FALSE;
    if (debug) { fprintf(stderr, "position of pixel [%d %d]", x, y); }
    
    /* Find the center {p} of pixel {x,y}: */
    r2_t p = (r2_t){{ x + 0.5, y + 0.5 }};
    /* If the four corners are inside, the pixel is inside: */
    bool_t inside = TRUE;
    int kx, ky;
    for (kx = -1; (kx <= 1) && inside; kx += 2)
      { for (ky = -1; (ky <= 1) && inside; ky += 2)
          { /* Get corner, slightly outwards of pixel: */
            r2_t q = p;
            q.c[0] += 0.50001*kx;
            q.c[1] += 0.50001*ky;
            if (debug) { r2_gen_print(stderr, &q, "%12.7f", "  q = (", " ", ")\n"); }
            inside &= ellipse_crs_inside(E, &q);
          }
      }
    if (inside) { return +1; }
    
    /* Else, if {p} is outside the fattened ellipse, the pixel is outside: */
    ellipse_crs_t EF = (ellipse_crs_t){ .ctr = E->ctr, .rad = E->rad + 0.707107, .str = E->str };
    bool_t outside = ! ellipse_crs_inside(&EF, &p);
    if (outside) { return -1; }
      
    /* Else give up: */
    return 0;
  }



    
double gauge_coverage_of_pixel(int x, int y, ellipse_crs_t *E)
  {
    double debug = FALSE;
    int pos = pixel_position_in_gauge_image(x, y, E);
    if (debug) { fprintf(stderr, "pixel [%d %d]  position = %+d ", x, y, pos); }
    double cov;
    if (pos < 0)
      { /* Fully outside */
        if (debug) { fprintf(stderr, " (outside)\n"); }
        cov = 0.0;
      }
    else if (pos > 0)
      { /* Fully inside */
        if (debug) { fprintf(stderr, " (inside)\n"); }
        cov = 1.0;
      }
    else
      { /* Straddling the border; find {cov} by sampling: */
        if (debug) { fprintf(stderr, " (straddles)"); }
        int NS = 15;  /* Pixel subsampling points along each axis for coverage. */
        int kx, ky;
        int nin = 0; /* Counts subsampling points inside the pixel. */
        for (kx = 0; kx < NS; kx++)
          { for (ky = 0; ky < NS; ky++)
              { /* Pick a subsampling point in pixel: */
                r2_t q = (r2_t){{ x + (kx + 0.5)/NS, y + (ky + 0.5)/NS }};
                /* Check whether {q} is inside the sphere's projection: */
                bool_t inside = ellipse_crs_inside(E, &q);
                if (inside) { nin++; }
              }
          }
        if (debug) { fprintf(stderr, " nin = %d", nin); }
        /* Make sure that {cov} is sufficiently away from 0 and 1 to avoid confusion: */
        double cov_min = 1.0/255.0;
        double cov_max = 254.0/255.0;
        cov = ((double)nin + 1)/(NS*NS + 2);
        if (cov < cov_min) { cov = cov_min; }
        if (cov > cov_max) { cov = cov_max; }
      }
    if (debug) { fprintf(stderr, " cov = %8.6f\n", cov); }
    return cov;
  }

void extract_data_from_gauge_image
  ( float_image_t *img,
    float_image_t *xtr,
    float_image_t* mask,
    ellipse_crs_t *E,
    double cutRadius,
    r3_t *view,
    r3_t **XP,
    double ***FP,
    int *NPP
  )
  {
    /* Get the image dimensions: */
    int NC = (int)img->sz[0];
    int NX = (int)img->sz[1];
    int NY = (int)img->sz[2];
    
    int c;
    
    /* Max width of gauge projection: */
    int WPmax = (int)(2*ceil(E->rad + r2_norm(&(E->str)))) + 1;
    
    /* Max number of pixels: */
    int NPmax = WPmax*WPmax;
    
    /* Allocate the data arrays for the max pixel count: */
    r3_t *X = (r3_t *)notnull(malloc(NPmax*sizeof(r3_t)), "no mem");
    double **F = (double **)notnull(malloc(NC*sizeof(double*)), "no mem");
    for (c = 0; c < NC; c++)
      { F[c] = (double *)notnull(malloc(NPmax*sizeof(double)), "no mem"); }
    int NP = 0;
    
    /* Get the normal correction matrix: */
    r3x3_t VM = compute_normal_correction_matrix(*view);
    
    int x, y;
    for (x = 0; x < NX; x++) 
      { for (y = 0; y < NY; y++) 
          { /* Compute fractional coverage {cov}: */
            double cov = gauge_coverage_of_pixel(x, y, E);

           
	    if(mask != NULL){
	      cov = cov * float_image_get_sample(mask,0,x,y);
	    }
	    /*This is for test cdist*/
 	    double dist = sqrt( ((x - E->ctr.c[0])*(x - E->ctr.c[0])) + ((y - E->ctr.c[1]))*(y-E->ctr.c[1]));
	    if(cutRadius > 0){
	      if( (cov > 0) && ( dist > (E->rad -cutRadius))){
// 		fprintf(stderr,"Pixel excluded - border %lf < %lf (%d,%d) - %lf\n",E->rad - dist, cutRadius,x,y,dist);
		cov = 0;
	      }
	    }
 	    if (xtr != NULL) { float_image_set_sample(xtr, 0,x,y, (float)cov); }
            if (cov == 1.0)
              { /* Compute mean view-relative surface normal in pixel: */
                r3_t nrm_rel = gauge_normal_in_pixel(x, y, E);
                /* Adjust normal for viewing direction: */
                r3_t nrm_abs;
                r3x3_map_row(&nrm_rel, &VM, &nrm_abs);
                /* Store normal and intensity in the vectors: */
                X[NP] = nrm_abs;
                for (c = 0; c < NC; c++) { F[c][NP] = float_image_get_sample(img,c,x,y); }
                NP++;
              }
          }
      }
              
    /* Resize vectors: */
    X = (r3_t *)realloc(X, NP*sizeof(r3_t));
    for (c = 0; c < NC; c++) { F[c] = (double *)realloc(F[c], NP*sizeof(double)); }
    /* return what should be returned */
    (*XP) = X;
    (*FP) = F;
    (*NPP) = NP;
  }


    
float_image_t *read_gauge_image(char *filename, double gamma)
  { char *ext = strrchr(filename, '.');
    demand(ext != NULL, "image name has no extension");
    float_image_t *img;
    if ((strcmp(ext, ".ppm") == 0) || (strcmp(ext, ".pgm") == 0))
      { img = float_image_read_pnm_named(filename,FALSE, gamma, 0.0, TRUE,TRUE,FALSE);
      }
    else if (strcmp(ext, ".fni") == 0)
      { FILE *rd = open_read(filename, TRUE);
        img = float_image_read(rd);
        if (rd != stdin) { fclose(rd); }
      }
    else
     { demand(FALSE, "unknown image name extension"); }
    return img;
  }
  
  
  
  
  
  


void synthetic_gauge_color_at_point
  ( r2_t *p,
    ellipse_crs_t *E, 
    r3x3_t *VM,
    int NC,
    void** l_data,
    evaluate_function* shading,
    double value[]  /* (OUT) */
  )
  {
    int c;
    for (c = 0; c < NC; c++)
      { /* Subtract the center and undo the stretching: */
        r2_t q = ellipse_crs_relative_coords(E, p);
        double r2 = r2_norm_sqr(&q);
        //light_field_t *Lc = &(L[c]);
        if(r2 <= 1.0)
          { /* Compute the view-relative normal to the sphere: */
            r3_t normal = (r3_t){{ q.c[0], q.c[1], sqrt(1 - r2) }}; /* View-relative normal */
            /* Compute the absolute normal: */
            r3x3_map_row(&normal, VM, &normal);
            //value[c] = compute_shade(D[c], &normal, Lc);
	    double x[3];
	    x[0] = normal.c[0];
	    x[1] = normal.c[1];
	    x[2] = normal.c[2];
	    //value[c] = shading(&normal,l_data[c]);
	    value[c] = shading(x,l_data[c]);
          }
        else
          { value[c] = 0; }
        if(!isfinite(value[c])){
	  fprintf(stderr,"Failure to generate virtual gab (channel %d) at:\n",c);
	  r3_t normal = (r3_t){{ q.c[0], q.c[1], sqrt(1 - r2) }}; /* View-relative normal */
            /* Compute the absolute normal: */
          r3x3_map_row(&normal, VM, &normal);
	  fprintf(stderr,"[%lf,%lf] - (%lf,%lf,%lf)\n",q.c[0], q.c[1],normal.c[0],normal.c[1],normal.c[2]);
	  fprintf(stderr,"VAL - %lf\n",value[c]);
	  double x[3];
	  x[0] = normal.c[0];
	  x[1] = normal.c[1];
	  x[2] = normal.c[2];
	  value[c] = shading(x,l_data[c]);
	}
        assert(isfinite(value[c]));
      }
  }



void synthetic_gauge_color_in_pixel
  ( int x,
    int y,
    ellipse_crs_t *E, 
    r3x3_t *VM,
    int NC,
    void** l_data,
    evaluate_function* shading,
    double value[]  /* (OUT) */
  )
  {
    int NS = 5;  /* Pixel subsampling points along each axis for normal averaging. */
    int c;
    for (c = 0; c < NC; c++) { value[c] = 0; }
    double qvalue[NC];
    int kx, ky;
    for (kx = 0; kx < NS; kx++)
      { for (ky = 0; ky < NS; ky++)
          { /* Pick a subsampling point in pixel: */
            r2_t q = (r2_t){{ x + (kx + 0.5)/NS, y + (ky + 0.5)/NS }};
            synthetic_gauge_color_at_point(&q, E, VM, NC, l_data,shading, qvalue);
            for (c = 0; c < NC; c++) { value[c] += qvalue[c]; }
          }
      }
    double nq = NS*NS;
    for (c = 0; c < NC; c++) { value[c] /= nq; }
  }



void paint_synthetic_gauge_image
  ( float_image_t *img,
    ellipse_crs_t *E, 
    r3_t *view,
    void** l_data,
    evaluate_function* shading
  )
  {
    /* Get the normal correction matrix: */
    r3x3_t VM = compute_normal_correction_matrix(*view);
    
    int NC = (int)img->sz[0];
    int NX = (int)img->sz[1];
    int NY = (int)img->sz[2];

    int nneg[NC];     /* Number of negative samples in each channel. */
    int nbig[NC];     /* Number of samples greater than 1.0 in each channel. */
    double value[NC]; /* Temp for pixel color. */
    int c;
    for (c = 0; c < NC; c++) { nneg[c] = 0; }
    int x, y;
    for (x = 0; x < NX; x++)
      { for (y = 0; y < NY; y++)
          { synthetic_gauge_color_in_pixel(x, y, E, &VM, NC,l_data,shading, value);
            for (c = 0; c < NC; c++)
              { if (value[c] < 0) { nneg[c]++; value[c] = 0; }
                if (value[c] > 1) { nbig[c]++; value[c] = 1; }
                float_image_set_sample(img, c, x, y, (float)(value[c]));
              }
          } 
      }
    for (c = 0; c < NC; c++) 
      {
        if (nneg[c] > 0)
          { fprintf(stderr, "** warning: %d negative samples in channel %d of synthetic image.\n", nneg[c], c); }
        if (nbig[c] > 0)
          { fprintf(stderr, "** warning: %d samples over 1.0 in channel %d of synthetic image.\n", nbig[c], c); }
      }
  }
  
