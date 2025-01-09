/* Last edited on 2023-01-14 01:00:37 by stolfi */

#define _GNU_SOURCE
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <float.h>
#include <assert.h>
#include <affirm.h>
#include <jsfile.h>
#include <float_image.h>
#include <float_image_io_pnm.h>
#include <float_image_heights.h>
#include <normais.h>
#include <estlogprob.h>
#include <rn.h>



#define TINY (1.0e-200)
/* Um valor minúsculo usado para evitar divisão por zero. */

#define SQRT_2_PI (2.5066282746310005024)



struct NoIndice{
  int indice;
  double valor;
};

typedef struct NoIndice noind;

/* Explicar estes pesos !!! */
double wRed = 0.4;
double wGrn = 0.4;
double wBlu = 0.2;

double corGab[3] = {1.0,1.0,1.0};

vetorn* cria_vetor(int tam);
noind* copy_vetorn(vetorn* v);
void libera_noind(noind* vet);
int compare (const void * a, const void * b);
double dist_alpha2(vetorn* s, vetorn* g);

double getRed(void){
  return wRed;
}

double getGrn(void){
  return wGrn;
}

double getBlu(void){
  return wBlu;
}

void set_Gab_Color(double R, double G, double B){
  corGab[0] = R;
  corGab[1] = G;
  corGab[2] = B;
}

void set_Peso(double Red, double Green, double Blue){
  wRed = Red;
  wGrn = Green;
  wBlu = Blue;
	
}

double media_ponderada(double R, double G, double B){
  double temp;
  temp = wRed*R;
  temp = temp + (wGrn*G);
  temp = temp + (wBlu*B);
  //temp = fabs(temp);// fabs - funcao do math.h que retorna o valor absoluto de um  double
  return temp;
}

double nova_media(double R, double G, double B, double* w){
  double divisor;
  double epslon = 0.01;
  double temp;
  divisor = (R*R) + (B*B) + (G*G) + epslon;
  *w = 1 /divisor;
	
  temp = wRed*R;
  temp = temp + (wGrn*G);
  temp = temp + (wBlu*B);

  return temp;
	
	
}

void extrai_assinatura(const double OBS[], double ass[], double *mag, int num_luzes)
{
  int i;
  double soma = TINY;
  for(i = 0; i < num_luzes; i++){ soma = soma + OBS[i]*OBS[i]; }
  (*mag) = sqrt(soma);
  for(i = 0; i < num_luzes; i++) { ass[i] = OBS[i]/(*mag); }
}

double dist_euclid(const double u[], const double v[], int n){
  int i;
  double soma = 0.0;
  for(i = 0; i< n;i++ ){
    double aux = (u[i] - v[i]);
    aux = (aux*aux);
    soma = soma + aux;
  }
	
  return sqrt(soma);
}



double dist_point_box(double p[],r2_t box[],int n){
	//compute the distance from point and a n-dimensional box 
	int i;
	int dist = 0;
	for(i = 0; i < n; i++){
		double blo = box[i].c[0];
		double bhi = box[i].c[1];
		if(p[i] < blo){
			dist+= ((blo - p[i])*(blo - p[i]));
		}else if(p[i] > bhi){
			dist+= ((bhi - p[i])*(bhi - p[i]));
		}
	}
	return  sqrt(dist);
	
}

double dist(vetorn* u, vetorn* v){
  return dist_euclid(u->valor,v->valor,u->tam);
}

int compare (const void * a, const void * b)
{
  double valora;
  double valorb;
  noind*  noa;
  noind*  nob;
  noa = (noind*) a;
  nob = (noind*) b;

  valora = noa->valor;
  valorb = nob->valor;
	
  //	fprintf(stderr,  "A:%f B:%f \n",valora,valorb);
  if(valora > valorb){
    return -1;
  }
  if(valora == valorb ) return 0;
  return 1;
}

int CompareDouble(const void* a, const void* b){
  double va,vb;
  va = *((double*)a);
  vb = *((double*)b);
  return va < vb;
}

bool_t debug_dist = FALSE;
/* Se esta variável for TRUE, as funções {EstLogPrSG_{NN}} se debugam. */

int localiza_simples
  ( Tabela* tab, 
    const double SO[], 
    estima_log_prob_S_G_t *estLogPrSG, 
    estima_albedo_t *estAlbedo, 
    double sigma, 
    double omg0, 
    double omg1, 
    double *logPrSG, 
    double *albedo, 
    normal_debug_opts_t *dbopt
  )
{
  int num_luzes = get_num_luzes(tab);
  int num_linhas = get_num_linhas(tab);

  /* Normaliza o vetor de observações {SO} obtendo a assinatura {so}: */
  double so[num_luzes];
  double Smag;
  extrai_assinatura(SO, so, &Smag, num_luzes);

  /* Melhor linha (candidata a resposta): */
  int melhor_linha = -1;
  double maior_logPrSG = -HUGE_VAL;
  
  /* Para debugging: */
  double menor_logPrSG = +HUGE_VAL;
  int pior_linha = -1;
  int segunda_melhor_linha = -1;

  /* Procura a entrada de {tab} que maximiza {estLogPrSG}: */
  int linha;
  debug_dist = (dbopt != NULL);
  for(linha = 0; linha < num_linhas; linha++) {
    const double *go = get_intdir(tab,linha);
    double Gmag = get_intmag(tab, linha);
    double logPrSG_linha = estLogPrSG(so, Smag, go, Gmag, num_luzes, sigma, omg0, omg1);
    
    if (logPrSG_linha > maior_logPrSG) {
      segunda_melhor_linha = melhor_linha;
      melhor_linha = linha;
      maior_logPrSG = logPrSG_linha;
    }
    if (logPrSG_linha < menor_logPrSG) {
      pior_linha = linha;
      menor_logPrSG = logPrSG_linha;
    }

    if (dbopt != NULL) {
      fprintf(dbopt->arq_debug, "SMP --- ");
      r3_t SN = get_normal(tab,linha);
      r3_t RN = dbopt->normal_ref;
      escreve_log_prob_e_alphas
        (dbopt->arq_debug, dbopt->hp, dbopt->vp, &RN, &SN, logPrSG_linha, so, Smag, go, Gmag, num_luzes);
    }
  }
  assert(melhor_linha != -1);
  assert(maior_logPrSG != -HUGE_VAL);
   
  /* Gera plots de {SO[i] x GO[i]} e {Pr(S)} em função de alpha, se solicitado: */
  if ((dbopt != NULL) && (dbopt->gera_plots_q)) { 
    /* Escolhe algumas linhas interessantes {lin[0..m-1]} da tabela: */
    int m = 5;
    int lin[m];
    r3_t zenite = (r3_t){{ 0,0,1 }};
    lin[0] = melhor_linha; /* Linha escolhida pela função de busca. */
    lin[1] = segunda_melhor_linha; /* Segunda melhor opção da função de busca. */
    lin[2] = localiza_linha_por_normal(tab, &(dbopt->normal_ref)); /* Linha da normal de referencia. */
    lin[3] = localiza_linha_por_normal(tab, &zenite); /* Linha com normal {(0,0,1)}. */
    lin[4] = pior_linha; /* Linha considerada pior pela função de busca. */
    escreve_linhas_interessantes(dbopt, tab, SO, estLogPrSG, sigma, omg0, omg1, lin, m);
    gera_plot_Si_Gi(dbopt, tab, SO, lin, m);
    gera_plot_LogPrSGalb(dbopt, tab, SO, sigma, omg0, omg1, lin, m);
     
  }

  /* Gera mapa de logprob sobre o gabarito, se solicitado: */
  if ((dbopt != NULL) && (dbopt->mapa_gabarito)) { 
    gera_mapa_prob_gabarito(dbopt, tab, SO, estLogPrSG, sigma, omg0, omg1);
  }

  /* Devolve o log da probabilidade: */
  (*logPrSG) = maior_logPrSG;
  
  /* Calcula e devolve o albedo mais provável: */
  { const double *go = get_intdir(tab, melhor_linha);
    double Gmag = get_intmag(tab, melhor_linha);
    (*albedo) = estAlbedo(so, Smag, go, Gmag, num_luzes, sigma, omg0, omg1);
  }

  /* Devolve o índice da linha encontrada. */
  return melhor_linha;
}

int localiza_alternativa
  ( Tabela* tab, 
    const double SO[], 
    estima_log_prob_S_G_t *estLogPrSG, 
    estima_albedo_t *estAlbedo, 
    double sigma, 
    double omg0, 
    double omg1, 
    double *logPrSG, 
    double *albedo, 
    normal_debug_opts_t *dbopt
  )
{
  int linha_res;
  int logPrSG_res;
  double albedo_res;

  /* Primeiro, procura a entrada com {EstLogPrSG_00} (supondo que não há sombras): */
  double logPrSG_e = -HUGE_VAL;
  double albedo_e = 0.0;
  int linha_e = localiza_simples (
    tab, SO, 
    EstLogPrSG_00, estAlbedo, sigma, omg0, omg1, 
    &logPrSG_e, &albedo_e, 
    dbopt
  );
  assert(linha_e != -1);
  assert(logPrSG_e != -HUGE_VAL);

  /* Calcula o erro esperado se não houvesse sombras: */ 
  int num_linhas = get_num_linhas(tab);
  double logPrSG_corte = sqrt(M_PI/num_linhas); /* !!! Determinar a fórmula !!! */

  if(logPrSG_e <= logPrSG_corte) { 
      /* Aceitável: */ 
      logPrSG_res = logPrSG_e;
      linha_res = linha_e;
      albedo_res = albedo_e;
  } else { 
      /* Repete busca com estimador {estLogPrSG}: */
      double logPrSG_w = -HUGE_VAL;
      double albedo_w = 0.0;
      int linha_w = localiza_simples (
        tab, SO, 
        estLogPrSG, estAlbedo, sigma, omg0, omg1, 
        &logPrSG_w, &albedo_w, 
        dbopt
      );
      assert(linha_w != -1);
      assert(logPrSG_w != -HUGE_VAL);
      linha_res = linha_w;
      logPrSG_res = logPrSG_w;
      albedo_res = albedo_w;
  }

  /* Devolve a resposta: */
  (*logPrSG) = logPrSG_res;
  (*albedo) = albedo_res;
  return linha_res;
}

void escreve_linhas_interessantes (
    normal_debug_opts_t *dbopt,
    Tabela* tab, 
    const double SO[], 
    estima_log_prob_S_G_t *estLogPrSG, 
    double sigma, 
    double omg0, 
    double omg1,
    int lin[],
    int m
  ) 
  {
    char *nome_arq = NULL;
    char *nome_arq = jsprintf("%s_%d_%04d_%04d_linhas.txt", dbopt->prefixo, dbopt->c, dbopt->hp, dbopt->vp);
    FILE *arq = open_write(nome_arq, TRUE);
    int n = get_num_luzes(tab);

    double so[n];
    double Smag;
    extrai_assinatura(SO, so, &Smag, n);

    int j, i;
    for (j = 0; j <= m; j++) {
      fprintf(arq, "%3d %6d", j, (j < 0 ? -1 : lin[j]));
      if (j == m) {
        fprintf(arq, "  %+9.6f %+9.6f %+9.6f", 0.0,0.0,0.0);
        /* Auto-gabarito: */
        double logPrSG_self = estLogPrSG(so, Smag, so, Smag, n, sigma, omg0, omg1);
        fprintf(arq, "  %+12.6lf", logPrSG_self);
        for (i = 0; i < n; i++) { fprintf(arq, " %5.3f", SO[i]); }
      } else {
        int linha = lin[j];
        const double *go = get_intdir(tab,linha);
        double Gmag = get_intmag(tab,linha);
        r3_t gn = get_normal(tab, linha);
        fprintf(arq, "  %+9.6f %+9.6f %+9.6f", gn.c[0], gn.c[1], gn.c[2]);
        double logPrSG_linha = estLogPrSG(so, Smag, go, Gmag, n, sigma, omg0, omg1);
        fprintf(arq, "  %+12.6lf", logPrSG_linha);
        for (i = 0; i < n; i++) { fprintf(arq, " %5.3f", Gmag*go[i]); }
      }
      fprintf(arq, "\n");
    }
    fclose(arq);
    free(nome_arq);
    fprintf(stderr, "pronto.\n");
  }


void gera_mapa_prob_gabarito (
    normal_debug_opts_t *dbopt,
    Tabela* tab, 
    const double SO[], 
    estima_log_prob_S_G_t *estLogPrSG, 
    double sigma, 
    double omg0, 
    double omg1
  )
  {
    int nx = 100;
    int ny = 100;
    float_image_t  *img = float_image_new(1, nx, ny);
    float_image_fill(img, -INF); 

    int n = get_num_luzes(tab);
    double so[n];
    double Smag;
    extrai_assinatura(SO, so, &Smag, n);
    
    int num_linhas = get_num_linhas(tab);
    int linha;
    for (linha = 0; linha < num_linhas; linha ++)
      { const double *go = get_intdir(tab,linha);
        double Gmag = get_intmag(tab,linha);
        r3_t gn = get_normal(tab, linha);
        double logPrSG = estLogPrSG(so, Smag, go, Gmag, n, sigma, omg0, omg1);
        int x = (int)floor(nx*(gn.c[0] + 1)/2);
        int y = (int)floor(ny*(gn.c[1] + 1)/2);
        float *p = float_image_get_sample_address(img, 0, x, y);
	if (logPrSG > (*p)) { (*p) = logPrSG; }
      }
  
    /* Find the range {[vmin _ vmax]} of all finite values in image: */
    float vmin = 0.0f, vmax = +1.0e-38f;
    float_image_update_sample_range(img, 0, &vmin, &vmax);
    float_image_rescale_samples(img, 0, vmin, vmax, 0.0, 1.0);
    char *nome_arq = NULL;
    char *nome_arq = jsprintf("%s_%d_%04d_%04d_pr_gab.pgm", dbopt->prefixo, dbopt->c, dbopt->hp, dbopt->vp);
    float_image_write_pnm_named(nome_arq, img,TRUE, 1.000, 0.000,FALSE,TRUE,TRUE);
    free(nome_arq);
    float_image_free(img);
    fprintf(stderr, "pronto.\n");
  }

void calcula_alphas(const double so[], const double go[], int num_luzes, double epsilon, double alpha[]) 
{
  int luz;
  for (luz = 0; luz < num_luzes; luz++) {
    alpha[luz] = calcula_alpha(so[luz],go[luz],epsilon);
  }
}

int CompareAlphaObs (const void * a, const void * b){
  alpha_obs_t *va = (alpha_obs_t *)a;
  alpha_obs_t *vb = (alpha_obs_t *)b;
  return va->alpha < vb->alpha;
}

void calcula_alphas_ordenados
  ( const double so[], 
    double Smag, 
    const double go[], 
    double Gmag, 
    int num_luzes, 
    double epsilon,  
    alpha_obs_t ao[]
  ) 
{
  int i;
  for(i = 0; i < num_luzes; i++) {
    ao[i].S_p = so[i]*Smag;
    ao[i].G_q = go[i]*Gmag;
    ao[i].alpha = calcula_alpha(ao[i].S_p, ao[i].G_q, epsilon);
    ao[i].luz = i;
  }
  qsort (ao, num_luzes, sizeof(alpha_obs_t), CompareAlphaObs);
}

void escreve_candidato_Sp_Gq(FILE *arq, int hp, int vp, r3_t *snp, r3_t *normal_ref, int num_luzes, double SO[], double GO[])
{
  fprintf(arq,"%04d %04d", hp, vp);
  fprintf(arq,"  %+9.6lf %+9.6lf %+9.6lf", snp->c[0], snp->c[1], snp->c[2]);
  fprintf(arq,"  %+9.6lf %+9.6lf %+9.6lf", normal_ref->c[0], normal_ref->c[1], normal_ref->c[2]);
  int j;
  for(j = 0; j < num_luzes; j++){
    fprintf(arq,"  %8.6lf %8.6lf", SO[j], GO[j]);
  }
  fprintf(arq,"\n");
  fflush(arq);
}

void escreve_log_prob_e_alphas
( FILE *output, 
  int hp, int vp, 
  r3_t *rnp,
  r3_t *snp, 
  double logPrSG_est, 
  const double so[], 
  double Smag, 
  const double go[], 
  double Gmag, 
  int num_luzes
  )
{
  fprintf(output,"%04d %04d",hp,vp);
  fprintf(output,"  %+8.5lf %+8.5lf %+8.5lf", rnp->c[0], rnp->c[1], rnp->c[2]);
  fprintf(output,"  %+8.5lf %+8.5lf %+8.5lf", snp->c[0], snp->c[1], snp->c[2]);
  double distNormals = dist_euclid(snp->c, rnp->c, 3);
  fprintf(output,"  %8.6lf %+12.6lf ", distNormals, logPrSG_est);
  // double logRatio = (2*log(distNormals+0.000001) - logPrSG_est)/log(10.0);
  // fprintf(output,"  %+10.6lf", logRatio);
  alpha_obs_t ao[num_luzes];
  calcula_alphas_ordenados(so, Smag, go, Gmag, num_luzes, OBS_NOISE, ao) ;
  int i;
  for(i = 0; i < num_luzes; i++){
    fprintf(output,"  %12.6lf",ao[i].alpha);
    fprintf(output," %5.3lf %5.3lf %02d", ao[i].S_p, ao[i].G_q, ao[i].luz);
  }
  fprintf(output,"\n");
}

double calcula_alpha(double S,double G,double epsilon)
{
  return (sqrt(S*S+epsilon*epsilon) - epsilon)/sqrt(G*G + epsilon*epsilon);
  //return sqrt(S*S+(epsilon*epsilon*0.5))/sqrt(G*G + epsilon*epsilon);
}






/* ESTIMATIVAS DO ALBEDO **************************************************/

double EstAlbedo_00(const double so[], double Smag, const double go[], double Gmag, int n, double sigma, double omg0, double omg1)
{
  double alb = (Smag + TINY)/(Gmag + TINY);
  return fmax(0.0, fmin(1.0, alb));
}


estima_albedo_t *escolhe_estAlbedo(int num_funcao)
{
  switch(num_funcao) {
    case 0: return &EstAlbedo_00;
    default: demand(FALSE, "numero invalido da funcao de albedo"); return NULL;
  }
}

double calcula_albedo(Tabela* tab, int resposta, const double SO[])
{
  int num_luzes = get_num_luzes(tab);
  int i;
  double soma = TINY;
  for(i = 0; i < num_luzes; i++ ){ soma += SO[i]*SO[i]; }
  double Smag = sqrt(soma);
  double Gmag = get_intmag(tab, resposta);
  return Smag/Gmag;
}

r3_t calcula_media_das_normais
  ( double peso_0, r3_t norm_0,
    double peso_1, r3_t norm_1,
    double peso_2, r3_t norm_2
  )
{
  /* Calcula a média ponderada das normais: */
  r3_t norm_media;
  int a;
  for (a = 0; a < 3; a++)
    { norm_media.c[a] = peso_0*norm_0.c[a] + peso_1*norm_1.c[a] + peso_2*norm_2.c[a]; }
  /* Normaliza {norm_media} para comprimento unitário: */
  double mag = r3_norm(&norm_media);
  if (mag != 0) { r3_scale(1/mag, &norm_media, &norm_media); }
  return norm_media;
}

double calcula_media_dos_pesos
  ( double dist_0, r3_t norm_0,
    double dist_1, r3_t norm_1,
    double dist_2, r3_t norm_2,
    r3_t norm_media
  )
{
  double eps = 0.01;
  double peso_0 = eps/hypot(dist_0, eps);
  double peso_1 = eps/hypot(dist_1, eps);
  double peso_2 = eps/hypot(dist_2, eps);
  //double peso = sqrt(wRed*peso_0*peso_0 + wGrn*peso_1*peso_1 + wBlu*peso_2*peso_2);
  /* O peso é menor quando a superfície é quase vertical: */

  double peso = eps/sqrt((wRed*peso_0*peso_0) +(wGrn*peso_1*peso_1) + (wBlu*peso_2*peso_2) + eps*eps);
  //double peso = eps/sqrt((wRed*dist_0*dist_0) +(wGrn*dist_1*dist_1) + (wBlu*dist_2*dist_2) + eps*eps);

  peso = peso * fmax(0.0, norm_media.c[2]);
  return peso;
}




double virtual_gab_phi1(double h,double r, double v){
	return h/r;
	//derivada em a1 da funcao phi
}

double virtual_gab_phi2(double h, double r, double v){
	return v/r;
}

double virtual_gab_phi3(double h, double r, double v){
	double aux;
	aux = r*r;
	aux = aux - (v*v);
	aux = aux - (h*h);
	if(aux < 0){
		aux = 0;
		//printf(" Menor %f %f.",h,v);
	}
	aux = sqrt(aux);
	aux = aux/r;
	//;printf("PHI3: %f",aux);
	return aux;
}

double virtual_gab_phi(double h, double r, double v,double a1, double a2, double a3){
	double resp;
	resp = a1*virtual_gab_phi1(h,r,v);
	resp = resp + (a2*virtual_gab_phi2(h,r,v));
	resp = resp + (a3*virtual_gab_phi3(h,r,v));
	if(resp < 0){
		resp =0;
	}
	return resp;

}


double virtual_gab_intensity(double x, double y,double radius,double a1,double a2, double a3, double x_center, double y_center){
	double h,v;
	v = y_center - y;
	h = x_center - x;
	
	return  virtual_gab_phi(h,radius,v,-a1,-a2,a3);
	
}

double lambertian_shading(r3_t dir_luz,double albedo, r3_t normal){
	double s = r3_dot(&dir_luz,&normal);
	return (s > 0 ? s*albedo: 0);
}

void converte_derivadas_para_normais(float_image_t* IDX, float_image_t* IDY, float_image_t* IN){
	assert(IN->sz[0] == 3);
	int nx = IN->sz[1];
	int ny = IN->sz[2];
	assert(IDX->sz[1] == nx);
	assert(IDX->sz[2] == ny);
	assert(IDY->sz[1] == nx);
	assert(IDY->sz[2] == ny);
	int x,y;
	for(x = 0; x < nx; x++){
		for(y = 0; y < ny; y++){
			r3_t vec;
			vec.c[0] = -float_image_get_sample(IDX,0,x,y); 
			vec.c[1] = -float_image_get_sample(IDY,0,x,y);
			vec.c[2] = 1.0 ;
			r3_dir(&vec,&vec);
			float_image_set_sample(IN,0,x,y,vec.c[0]);
			float_image_set_sample(IN,1,x,y,vec.c[1]);
			float_image_set_sample(IN,2,x,y,vec.c[2]);
		}
	}
} 

void interpola_normais(float_image_t* nrm,float_image_t* IW,char* prefDebug){
	
	
	
	float_image_t* IDX = calcula_derivadas(nrm, 'x', IW);
	float_image_t* IDY = calcula_derivadas(nrm, 'y', IW);
	interpola_derivadas_recursivo(IDX, IDY, IW,0, prefDebug);
	converte_derivadas_para_normais(IDX,IDY,nrm);
	float_image_free(IW);
}

void interpola_normais_unsafe(float_image_t* nrm,float_image_t* IW,char* prefDebug){
	
	
	
	float_image_t* IDX = calcula_derivadas_unsafe(nrm, 'x', IW);
	float_image_t* IDY = calcula_derivadas_unsafe(nrm, 'y', IW);
	interpola_derivadas_recursivo(IDX, IDY, IW,0, prefDebug);
	converte_derivadas_para_normais(IDX,IDY,nrm);
	float_image_free(IW);
}


/*RANSAC-de-burro*/

r3_t compute_normal_by_RANSAC(double SO[],Tabela* tab,double sigma,double probOutlier,double albGab,double*albedo_res,double* ransac_weight){
  int nLights = get_num_luzes(tab);
  int nComb =  nLights*(nLights-1)*(nLights-2)/6.0;
  int valids[nComb];
  r3_t normals[nComb];
  int tablines[nComb];
  double SmagK[nComb];
  double GmagK[nComb];
  int i,j,k;
  
  int count = 0;
  
  for(i = 0; i < nLights; i++){
    for(j = (i+1); j < nLights; j++){
      for(k = (j+1); k < nLights;k++){
	
	double c_SO[3];
	c_SO[0] =  SO[i];
	c_SO[1] =  SO[j];
	c_SO[2] =  SO[k];
	
	double c_so[3];
	
	double c_Smag = rn_dir(3,c_SO,c_so);
	SmagK[count] = c_Smag;
	
	int l;
	int num_linhas = get_num_linhas(tab);
	int melhor_linha = 0;
	double melhor_logprob = INF;
	double c_go_best[3];
	double c_Gmag_best;
	
	//return (r3_t){{0,0,1}};
	
	for(l = 0 ; l < num_linhas; l++){
	  const double* go = get_intdir(tab,l);
	  double Gmag = get_intmag(tab,l);
	  
	  double c_GO[3];
	  c_GO[0] = go[i]*Gmag;
	  c_GO[1] = go[j]*Gmag;
	  c_GO[2] = go[k]*Gmag;
	  
	  double c_go[nLights];
	  double c_Gmag = rn_dir(3,c_GO,c_go);
	    
	  double dist  = rn_dist(3,c_so,c_go);
	  if(dist < melhor_logprob){
	    melhor_linha = l;
	    melhor_logprob = dist;
	    rn_copy(3,c_go,c_go_best);
	    c_Gmag_best = c_Gmag;
	  }
	}
	
	normals[count] = get_normal(tab,melhor_linha);
	tablines[count] = melhor_linha;
	GmagK[count] = c_Gmag_best;
	valids[count] = 1;
	
	
	count++;
	assert(count <= nComb);
      }
    }
  }
  
    double Var[nComb];
   // double P[nComb];
    double weights[nComb];
    double previousW[nComb];
//    double A[nComb];
    //double albedo;
    
    r3_t avg_norm;
  
    //init weights
    for(i = 0; i < nComb; i++){
      weights[i] = 0;
      if(valids[i]){
	weights[i] = 1.0/nComb;
	//P[i] = 1.0/(double)nComb;
      }
    }

    
    double diff = 10;
    int iter = 0;
    double gaussian_norm_const = 1.0/(sigma*sigma*2*M_PI);
    
    while((diff  > 10e-06) && (iter < 200)){
     
      r3_zero(&avg_norm);
      /*Save old weights*/
      rn_copy(nComb,weights,previousW);
      /*Compute the averaged normal*/
      for(i = 0; i < nComb;i++){
	if(valids[i]){
	  r3_mix_in(weights[i],&(normals[i]),&(avg_norm));
	}
      }
      double avg_len = r3_dir(&(avg_norm),&(avg_norm));
      /*Crash functions */
      assert(avg_len != 0);
      assert(!isnan(avg_norm.c[0]));
      assert(!isnan(avg_norm.c[1]));
      assert(!isnan(avg_norm.c[2]));
      /*recompute the weights*/
      double sumW = 0;
      for(i = 0; i < nComb;i++){
	double dist = r3_dist(&avg_norm,&(normals[i]));
	double probk = gaussian_norm_const*exp(dist*dist/(2*sigma*sigma));/*probability to be inlier*/
	//probk = exp((-(dist*dist)/4.0 + 1.0)*exp(1));
	double prob_ink = ((1 - probOutlier)*probk)/((probOutlier/2.0*M_PI) +  ((1-probOutlier)*probk));
	double prob_ouk = 1 - prob_ink;
	
	Var[i] = (prob_ink*(sigma*sigma)) + (prob_ouk*0.5);
	weights[i] = 1/Var[i];
	sumW+= weights[i];
	//weights[i] = (1 - probOutlier)/(probOutlier +  ((1-probOutlier)*probk));
      }
      /*normalize the weights*/
      
      for(i = 0; i < nComb;i++){
	weights[i]/=sumW;
      }
      
      diff = rn_dist(nComb,weights,previousW);
      iter++;
    }
    
    //usar maior dos pesos ou a média ?
    double finalW = 0;
    double sumW = 0;
    double albedo = 0;
    for(i = 0; i < nComb; i++){
	if(valids[i]){
	  
	   /* Sets {r := r + s * a}. */
	  double albk = SmagK[i]/(GmagK[i]*albGab);
	  albedo+= albk*weights[i];
	  finalW+= weights[i]*weights[i]*Var[i];
	  sumW+=weights[i];
	}
    }
    
    
    *ransac_weight = 1/(finalW/(sumW*sumW));
    //r3_dir(&avg_norm,&avg_norm);
//     r3_scale(albedo,&avg_norm,&avg_norm);
//     n_res[0] = avg_norm.c[0];
//     n_res[1] = avg_norm.c[1];
//     n_res[2] = avg_norm.c[2];
//       
  
  r3_t n_res_snp = avg_norm;
  assert(!isnan(avg_norm.c[0]));
  assert(!isnan(avg_norm.c[1]));
  assert(!isnan(avg_norm.c[2]));
  if(isnan(albedo)){
    albedo = 0;
  }
//   n_res_snp.c[0] = n_res[0];
//   n_res_snp.c[1] = n_res[1];
//   n_res_snp.c[2] = n_res[2];
//   
  *albedo_res = albedo;
  return n_res_snp;
  
  
}



r3_t compute_normal_by_fast_RANSAC(double SO[],int nLights,fast_hash_t** FH_set,int** seq_num, int num_sets,double sigma,double probOutlier,double albGab,double*albedo_res,double* ransac_weight,int* valid_lights){
//   int nLights = FH_set[0]->nLights;
  int nComb =  num_sets;
  int valids[nComb];
  r3_t normals[nComb];
  double SmagK[nComb];
  double GmagK[nComb];
  double alber[nComb];
  int i,j,k;
  
  int count = 0;
  int count_valids = 0;
 
//   for(i = 0; i < nLights; i++){
//     for(j = (i+1); j < nLights; j++){
//       for(k = (j+1); k < nLights;k++){
// 
  for(count = 0; count < nComb; count++){
	double c_SO[3];
	i = seq_num[count][0];
	j = seq_num[count][1];
	k = seq_num[count][2];
	
	c_SO[0] =  SO[i];
	c_SO[1] =  SO[j];
	c_SO[2] =  SO[k];
	    
	double c_so[3];
	    
	double c_Smag = rn_dir(3,c_SO,c_so);
	SmagK[count] = c_Smag;
	    
	double alb;
	double logProb;
	
	    
	normals[count] = fast_hash_compute_normal(FH_set[count],c_SO,sigma,0,0,&alb,&logProb);
	alber[count] = alb;
	GmagK[count] = SmagK[count]/alb;
	double len = r3_norm(&(normals[count]));
	if( len < 0.99) {
	  logProb  = 0;
	}
	if(logProb != 0){
	  valids[count] = 1;
	  count_valids++;
	}else{
	  valids[count] = 0;
	}
	
	/*This should be put back in case of restore previous function*/
	//count++;
	assert(count <= nComb);
	
  }
	
//       }
//     }
//   }
//     
    double Var[nComb];
   // double P[nComb];
    double weights[nComb];
    double previousW[nComb];
//    double A[nComb];
    //double albedo;
  
    *valid_lights = count_valids;
    
    r3_t avg_norm;
  if(count_valids == 0){
    *albedo_res = 0;
    *ransac_weight=  0;
    
    return (r3_t){{0,0,0}};
    
  }
  
  
  
  if(count_valids == 1) {
    for(i = 0; i < nComb; i++){
      if(valids[i] == 1) break;
    }
    double Vfinal = probOutlier*1 + (1-probOutlier)*sigma;
    
   *albedo_res = alber[i];
   *ransac_weight = 1/Vfinal;
   return normals[i];
  }
  
  if(count_valids == 2){
    r3_t normal_valid[2];
    double albs[2];
    int count_norm = 0;
    
    for(i = 0; i < nComb; i++){
      if(valids[i] == 1) {
	normal_valid[count_norm] = normals[i];
	albs[count_norm] = alber[i];
	count_norm++;
      }
    }
    
    r3_t avg_norm;
    r3_add(&(normal_valid[0]),&(normal_valid[1]),&avg_norm);
    r3_dir(&avg_norm,&avg_norm);
    double albedo = (albs[0]+albs[1])/2.0;
    double Vfinal = sigma + ( 1 - r3_dot(&(normal_valid[0]),&(normal_valid[1])));
    *albedo_res = albedo;
    *ransac_weight = 1/Vfinal;
     return avg_norm;
    
  }
    //init weights
    for(i = 0; i < nComb; i++){
      weights[i] = 0;
      if(valids[i]){
	weights[i] = 1.0/nComb;
	//P[i] = 1.0/(double)nComb;
      }
    }

    
    double diff = 10;
    int iter = 0;
    double gaussian_norm_const = 1.0/(sigma*sigma*2*M_PI);
    
    while((diff  > 10e-06) && (iter < 200)){
     
      r3_zero(&avg_norm);
      /*Save old weights*/
      rn_copy(nComb,weights,previousW);
      /*Compute the averaged normal*/
      for(i = 0; i < nComb;i++){
	if(valids[i]){
	  r3_mix_in(weights[i],&(normals[i]),&(avg_norm));
	}
      }
      double avg_len = r3_dir(&(avg_norm),&(avg_norm));
      /*Crash functions */
      assert(avg_len != 0);
      assert(!isnan(avg_norm.c[0]));
      assert(!isnan(avg_norm.c[1]));
      assert(!isnan(avg_norm.c[2]));
      /*recompute the weights*/
      double sumW = 0;
      for(i = 0; i < nComb;i++){
	double dist = r3_dist(&avg_norm,&(normals[i]));
	double probk = gaussian_norm_const*exp(dist*dist/(2*sigma*sigma));/*probability to be inlier*/
	//probk = exp((-(dist*dist)/4.0 + 1.0)*exp(1));
	double prob_ink = ((1 - probOutlier)*probk)/((probOutlier/2.0*M_PI) +  ((1-probOutlier)*probk));
	double prob_ouk = 1 - prob_ink;
	
	Var[i] = (prob_ink*(sigma*sigma)) + (prob_ouk*0.5);
	weights[i] = 1/Var[i];
	sumW+= weights[i];
	//weights[i] = (1 - probOutlier)/(probOutlier +  ((1-probOutlier)*probk));
      }
      /*normalize the weights*/
      
      for(i = 0; i < nComb;i++){
	weights[i]/=sumW;
      }
      
      diff = rn_dist(nComb,weights,previousW);
      iter++;
    }
    
    //usar maior dos pesos ou a média ?
    double finalW = 0;
    double sumW = 0;
    double albedo = 0;
    for(i = 0; i < nComb; i++){
	if(valids[i]){
	  
	   /* Sets {r := r + s * a}. */
	  double albk = SmagK[i]/(GmagK[i]*albGab);
	  albedo+= albk*weights[i];
	  finalW+= weights[i]*weights[i]*Var[i];
	  sumW+=weights[i];
	}
    }
    
    
    *ransac_weight = 1/(finalW/(sumW*sumW));
    //r3_dir(&avg_norm,&avg_norm);
//     r3_scale(albedo,&avg_norm,&avg_norm);
//     n_res[0] = avg_norm.c[0];
//     n_res[1] = avg_norm.c[1];
//     n_res[2] = avg_norm.c[2];
//       
  
  r3_t n_res_snp = avg_norm;
  assert(!isnan(avg_norm.c[0]));
  assert(!isnan(avg_norm.c[1]));
  assert(!isnan(avg_norm.c[2]));
  if(isnan(albedo)){
    albedo = 0;
  }
//   n_res_snp.c[0] = n_res[0];
//   n_res_snp.c[1] = n_res[1];
//   n_res_snp.c[2] = n_res[2];
//   
  *albedo_res = albedo;
  return n_res_snp;
  
  
}
