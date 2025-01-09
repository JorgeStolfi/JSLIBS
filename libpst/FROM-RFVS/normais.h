#ifndef normais_H
#define normais_H

#include <vetorn.h>
#include <tabela.h>
#include <float_image.h>
#include <fast_hash.h>

#define OBS_NOISE 0.001
/* Erro esperado em cada componente de um vetor de observações
  (supondo intensidades em [0_1]). */

typedef struct Vetor vetor;

typedef struct alpha_obs_t { double alpha; double S_p; double G_q; int luz; } alpha_obs_t;
/* Estrutura para ordenação de alphas com vetores de observações. 
  {S_p} é a intensidade de um pixel da cena, {G_q} a intensidade
  do pixel associado no gabarito, {luz} é o índice da luz. */

typedef struct normal_debug_opts_t {
  int c;              /* Canal da imagem sendo processado. */
  int hp;             /* Coluna do pixel {p} da cena (estilo PNM). */
  int vp;             /* Linha do pixel {p} da cena (estilo PNM). */
  FILE *arq_debug;
  char *prefixo;
  bool_t mapa_gabarito;    /* Se TRUE, gera uma imagem do gabarito mostrando probabilidades. */
  r3_t normal_ref;         /* Normal correta, ou {(0,0,0)} se desconhecida. */
  bool_t gera_plots_q;     /* Se TRUE, gera arquivos com varios plots para {q} escolhidos. */
} normal_debug_opts_t;


typedef double estima_log_prob_S_G_t(const double so[], double Smag, const double go[], double Gmag, int n, double sigma, double omg0, double omg1); 
/* Tipo de uma função que estima a probabilidade do vetor de observações da cena 
  {SO[0..n-1} supondo que o vetor de de observações do gabarito é {GO[0..n-1]}.
  As assinaturas {SO} e {GO} são as assinaturas normalizadas {so,go}
  multiplicadas por {Smag,Gmag}. */

/* DETERMINAÇÃO DA NORMAL
  
  Cada função {busca_XXX} abaixo localiza a entrada na tabela {tab}
  cujo vetor de observações {GO[0..n-1]} maximiza a verossimilhanca
  {Pr(SO|GO)}, estimada pela função {estLogPrSG}; onde {SO[0..n-1} é
  o vetor de observações de um pixel, e {n} é o número de luzes.
  
  A função devolve o índice dessa entrada como resultado, coloca em
  {*logPrSG} o valor retornado por {estLogPrSG}, e coloca em {*albedo}
  o albedo mais provável supondo que essa é a entrada correta.

  Se {dbopt} não é NULL, a função pode gerar vários arquivos de debug.
  Nesse caso ela usa o canal {c} e os índices do pixel {hp,vp} (para
  identificação e nomes de arquivos apenas) e a normal de referência
  {*normal_ref}.  */


typedef double estima_albedo_t(const double so[], double Smag, const double go[], double Gmag, int n, double sigma, double omg0, double omg1); 
/* Tipo de uma função que determina o albedo mais provável de um pixel
  da cena, dado seu vetor de observações {SO[0..n-1} e o vetor de de
  observações {GO[0..n-1]} do pixel do gabarito com a mesma normal.
  As assinaturas {SO} e {GO} são as assinaturas normalizadas {so,go}
  multiplicadas por {Smag,Gmag}. */

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
  );
/* Faz busca exaustiva (compara com todas as entradas da tabela).  Usa
  a distância-alpha de tipo {which_eval} com parâmetro {sigma} para
  determinar a similaridade das entradas. */

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
  );
/* Faz primeiro uma busca exaustiva usando o estimador simples {EstLogPrSG_00}
  (baseado na distância euclidiana de assinaturas normalizadas). Se o resultado
  for ruim ({*logPrSG} muito baixo, considerada a granularidade da tabela),
  refaz a busca usando o estimador dado {estLogPrSG} */

void extrai_assinatura(const double OBS[], double ass[], double *mag, int num_luzes);
double dist_euclid(const double *u, const double *v, int n);
double dist_point_box(double p[],r2_t box[],int n);
void calcula_alphas(const double so[], const double go[], int num_luzes, double epsilon, double alpha[]) ;
double calcula_alpha(double S,double G,double epsilon);
  
void calcula_alphas_ordenados
  ( const double so[], 
    double Smag, 
    const double go[], 
    double Gmag, 
    int num_luzes, 
    double epsilon,  
    alpha_obs_t ao[]
  );
/* Preenche o vetor {ao[0..num_luzes-1]} com os vetores de observação
  {SO} e {GO}, obtidos multiplicando as assinaturas
  {so[0..num_luzes-1]} e {go[0..num_luzes-1]} por {Smag} e {Gmag},
  respectivamente. Também calcula o campo {alpha} de cada entrada, e
  ordena as entradas por {alpha} descrescente.  O parâmetro {epsilon}
  deve ser o desvio padrão dos ruídos de quantização e medida em cada
  componente dos vetores de observação. */

double calcula_albedo(Tabela* tab, int linha, const double SO[]);
/* Calcula o albedo de um pixel {p} da cena com vetor de observação {SO}, 
  supondo que a entrada correspondente{linha} de tab} é um ponto do gabarito
  com a mesma normal que {p} e albedo 1.0, E QUE NÃO HÁ SOMBRAS NEM
  HIGHLIGHTS.  Era {cor_original} */

r3_t calcula_media_das_normais
  ( double peso_0, r3_t norm_0,
    double peso_1, r3_t norm_1,
    double peso_2, r3_t norm_2
  );
  /* Calcula a normal de um pixel {p}, dadas as normais
    {norm0,norm1,norm2} determinada em cada canal, e os respectivos
    pesos de confiabilidade {peso_0,peso_1,peso_2}. */
  

double calcula_media_dos_pesos
  ( double peso_0, r3_t norm_0,
    double peso_1, r3_t norm_1,
    double peso_2, r3_t norm_2,
    r3_t norm_media
  );
  /* Calcula o peso de confiabilidade de um pixel {p}, dadas as
    normais {norm0,norm1,norm2} determinada em cada canal, os
    respectivos pesos de confiabilidade {peso_0,peso_1,peso_2}, e a normal
    calculada por {calcula_media_das_normais}. */

void escreve_candidato_Sp_Gq
  ( FILE *arq, 
    int hp, 
    int vp, 
    r3_t *snp, 
    r3_t *normal_ref, 
    int num_luzes, 
    double SO[], 
    double GO[]
  );
/* Escreve no arquivo {arq} os dados de um par de pixels {p,q}, um da cena e outro do gabarito.
  As coordenadas {hp,vp} são gravadas simplesmente para identificação. */

void escreve_log_prob_e_alphas
  ( FILE *output, 
    int hp, 
    int vp, 
    r3_t *rnp,   /* Normal correta do ponto (hp,vp), ou (0,0,0) se desconhecida. */
    r3_t *snp,   /* Normal determinada pelo algoritmo para o ponto (hp,vp). */
    double logPrSG_est, 
    const double so[], 
    double Smag, 
    const double go[], 
    double Gmag, 
    int num_luzes
  );
/* Escreve no arquivo {output} os dados de um par de pixels {p,q}, um da cena e outro do gabarito.
  Calcula os alphas e ordena as luzes em ordem de alpha decrescente.
  O parâmetro {logPrSG_est} é {log(Pr(S|G))} estimada pelo algoritmo que determinou {snp}. 
  As coordenadas {hp,vp} são gravadas simplesmente para identificação. */

double media_ponderada(double R, double G, double B);
double nova_media(double R, double G, double B, double* w);
void set_Peso(double Red, double Green, double Blue);
double getRed(void);
double getGrn(void);
double getBlu(void);
void set_Gab_Color(double R, double G, double B);
double dist(vetorn* u, vetorn* v);
int CompareDouble(const void* a, const void* b);

double calcula_albedo(Tabela* tab, int linha, const double SO[]);
/* Calcula o albedo de um pixel {p} da cena com vetor de observação {SO}, 
  supondo que a entrada correspondente{linha} de tab} é um ponto do gabarito
  com a mesma normal que {p} e albedo 1.0, E QUE NÃO HÁ SOMBRAS NEM
  HIGHLIGHTS.  Era {cor_original} */

int CompareAlphaObs (const void * a, const void * b);
/* Compara os campos {alpha} de dois {alpha_obs_t}. */


  

/* ESTIMADORES DE ALBEDO */
  
estima_albedo_t *escolhe_estAlbedo(int num_funcao);

double EstAlbedo_00(const double so[], double Smag, const double go[], double Gmag, int n, double sigma, double omg0, double omg1);
/* Primário: {albedo = Smag/Gmag}, evitando divisão por 0 e clipado em [0 _ 1]. */



/* FUNÇÕES DE DEBUG INTERNAS */

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
  );



void gera_mapa_prob_gabarito (
    normal_debug_opts_t *dbopt,
    Tabela* tab, 
    const double SO[], 
    estima_log_prob_S_G_t *estLogPrSG, 
    double sigma, 
    double omg0, 
    double omg1
  );

/*Functions related to virtual gauges*/
//DEPRECATED
double virtual_gab_phi1(double h, double r, double v);
double virtual_gab_phi2(double h, double r, double v);
double virtual_gab_phi3(double h, double r, double v);
double virtual_gab_phi(double h, double r, double v,double a1, double a2, double a3);
double virtual_gab_intensity(double x, double y,double radius,double a1,double a2, double a3, double x_center, double y_center);
// New shading function
double lambertian_shading(r3_t dir_luz,double albedo, r3_t normal);

void interpola_normais(float_image_t* nrm,float_image_t* IW,char* prefDebug);
void interpola_normais_unsafe(float_image_t* nrm,float_image_t* IW,char* prefDebug);
void converte_derivadas_para_normais(float_image_t* IDX, float_image_t* IDY, float_image_t* IN);
/*
Interpola normais {nrm} preenchendo todas as normais que são {0,0,0}. 
*/

r3_t compute_normal_by_RANSAC(double SO[],Tabela* tab,double sigma,double probOutlier,double albGab,double*albedo_res,double* ransac_weight);
/*Dumb RANSAC normal computation, returns the norma, the albedo and the weight*/

r3_t compute_normal_by_fast_RANSAC(double SO[],int nLights,fast_hash_t** FH_set,int** seq_num, int num_sets,double sigma,double probOutlier,double albGab,double*albedo_res,double* ransac_weight,int* valid_lights);
/*RANSAC using FAST hash normal computation, it supports variable number of sets, the parameter {seq_num} is a array {num_sets}x3
  which contains the correspondent lights of each set
 */

#endif
