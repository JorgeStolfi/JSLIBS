#ifndef ESTLOGPROB_H
#define ESTLOGPROB_H

#include <stdio.h>
#include <normais.h>
#include <tabela.h>
#define K_09 1.5

double ProbSiGialb_canonical(double Si, double Gi, double alb, double sigma);
  /* Calcula {Pr(S[i]|G[i],alb)} pelo modelo contínuo adotado, com parâmetros {sigma} supondo que Si é canonico. */
double ProbSiGialb_shadow(double Si, double Gi, double alb, double sigma);
  /* Calcula {Pr(S[i]|G[i],alb)} pelo modelo contínuo adotado, com parâmetros {sigma} supondo que Si é sombreado. */
double ProbSiGialb_highlight(double Si, double Gi, double alb, double sigma);
  /* Calcula {Pr(S[i]|G[i],alb)} pelo modelo contínuo adotado, com parâmetros {sigma} supondo que Si tem highlights. */

double ProbSiGialb(double Si, double Gi, double alb, double sigma, double omg0, double omg1);
  /* Calcula {Pr(S[i]|G[i],alb)} pelo modelo contínuo adotado, com parâmetros {sigma,omg0,omg1}. */

double LogProbSiGialb(double Si, double Gi, double alb, double sigma, double omg0, double omg1);
  /* Calcula {log(Pr(S[i]|G[i],alb))} pelo modelo contínuo adotado, com parâmetros {sigma,omg0,omg1}. */

double LogPrSGalb(const double SO[], const double GO[], int n, double alb, double sigma, double omg0, double omg1);
/* Calcula {log(Pr(SO|GO,alb))} para o valor dado de {alb}, somando {LogProbSiGialb(SO[i],GO[i],alb,...)}
  para todo {i}.  Supõe que os erros nas componentes {SO[i]} são independentes. */  
    
double LogPrSG_CbP(const double SO[],const double GO[], int n, double sigma, double omg0, double omg1,double K, FILE *arq);
/* Calcula uma aproximação do log da integral definida de {Pr(SO|GO,alb)} para {alb} entre 0 e 1.  Supõe que os erros nas componentes
  {SO[i]} são independentes.  Usa uma aproximção constante-por-partes dos fatores, com intervalo central de largura {K*sigma}.
  Se {arq} é diferente de NULL, escreve nele a aproximação usada.  */  

double LogPrSG_UnS(const double SO[],const double GO[], int n, double sigma, double omg0, double omg1,int  nsteps);
/* Calcula uma aproximação do log da integral definida de {Pr(SO|GO,alb)} para {alb} entre 0 e 1.  Supõe que os erros nas componentes
  {SO[i]} são independentes.  Usa uma somatória com {nsteps} parcelas. 
 */  

/* LISTAS DE INTERVALOS DE ALBEDO
  
  Para as funções a seguir, cada par {S[i],G[i]} define um intervalo
  de albedos com extremos {alo[i]} e {ahi[i]}.  Os extremos de todos
  esses intervalos são armazenados em {vx[0..2*n-1]}, misturados, e
  identificados pelos inteiros {ix[0..2*n-1]}, com a seguinte
  convenção: se {ix[j]} é positivo, {vx[j]} é o extremo superior
  {ahi[ix[j]-1]}; se {ix[j]} é negativo, {vx[j]} é o extremo inferior
  {alo[1-ix[j]]}.  Por exemplo, se {ix[0..5] = {-2,+2,-1,-3,+1,+3}},
  os extremos {vx[0..5]} são {alo[1],ahi[1],alo[0],alo[2],ahi[0],ahi[2]}.

  Os intervalos são truncados de modo a estar estritamente no interior
  do interval {[0 _ 1]}. */

void CriaExtremosDeAlbedo_08(const double SO[], const double GO[], int n, double R, int ix[], double vx[]);
  /* Cria intervalos de albedo com {alo[i] = (S[i] - R)/G[i]} e {ahi[i] = (S[i] + R)/G[i]}.
    O parâmetro {R} é a faixa de incerteza de {Si} para um albedo e {Gi} fixos.  */

void CriaExtremosDeAlbedo_09(const double SO[], const double GO[], int n, double R, int ix[], double vx[]);
  /* Cria os intervalos de albedo necessários para o estimador {EstLogProb_09}.
    O parâmetro {R} é a faixa de incerteza de {Si} para um albedo e {Gi} fixos. */

void OrdenaExtremosDeAlbedo(int ix[], double vx[], int nex);
  /* Ordena a lista de extremos de intervalos de albedos {vx[0..nex-1]} e {ix[0..nex-1]}, 
    em ordem crescente de {vx}. */

/* ESTIMADORES DE VEROSSIMILHANÇA



  Estas funções calculam o valor aproximado de {log(Pr(SO|GO))} por vários métodos.
  O parâmetro {sigma} é o desvio padrão dos erros de medida `normais' em cada {SO[i]}.
  Os parâmetros {omg0,omg1} determinam a probabilidade do valor de {SO[i]} ser `anômalo'
  (outlier), respectivamente menor e maior que o valor esperado pela normal e albedo.  */





estima_log_prob_S_G_t *escolhe_estLogPrSG(int num_funcao);

double EstLogPrSG_00(const double so[], double Smag, const double go[], double Gmag, int n, double sigma, double omg0, double omg1);
/* Modelo simples, supõe que não há highlights, ignora {omg0,omg1}. */

double EstLogPrSG_01(const double so[], double Smag, const double go[], double Gmag, int n, double sigma, double omg0, double omg1);
/* Modelo mais simples ainda,ignora {go,Gmag,omg0,omg1}.Decide prob. apenas baseado no {Smag/sigma} */

double EstLogPrSG_07(const double so[], double Smag, const double go[], double Gmag, int n, double sigma, double omg0, double omg1);
/* não exige alphas ordenados.  ACEITA HIGHLIGHTS. Usa Bayes. Integra {Pr(S|G,alb) dalb} sobre
  todo o intervalo {[0,1]}. Demorada. */

double EstLogPrSG_08(const double so[], double Smag, const double go[], double Gmag, int n, double sigma, double omg0, double omg1);
/* não exige alphas ordenados.  ACEITA HIGHLIGHTS. Usa Bayes. Demorada. Similar a {EstLogPrSG_07}, mas procura limitar a 
  integral a um sub-intervalo onde há maior número de fatores com picos. */

double EstLogPrSG_09(const double so[], double Smag, const double go[], double Gmag, int n, double sigma, double omg0, double omg1);
/* não exige alphas ordenados.  ACEITA HIGHLIGHTS. Usa Bayes. Demorada. Similar a {EstLogPrSG_07}, mas usa
  uma aproximação constante-por-partes dos fatores. */

void gera_plot_Si_Gi (
    normal_debug_opts_t *dbopt,
    Tabela* tab, 
    const double SO[], 
    int lin[],
    int m
  );
  /*Write into a file specified by {dbopt} a plot of the pairs {SO[i],GO[k][i]} for k in [0..m] and i
   in [0...n-1], where GO[k] is the row of {tab} with index lin[k], and n is the number of lights (within {tab}); except when k=m in which case GO[k] = SO.
   There is one line in the file for each i with m+2 columns, {SO[i], GO[0][i]... GO[m-1][i], SO[i]}.
  */

void gera_plot_LogPrSGalb (
    normal_debug_opts_t *dbopt,
    Tabela* tab, 
    const double SO[], 
    double sigma, 
    double omg0, 
    double omg1,
    int lin[],
    int m
  );
  
  /*Write into a file specified by {dbopt} a plot of LogPrSGalb(SO,GO[k],alb) for k in [0..m] and alb
   varying between zero and one, where GO[k] is the row of {tab} with index lin[k]; except when k=m in which case GO[k] = SO.
   There is one line in the file for each i with 2 columns, alb and LogPrSGalb.
  */



#endif