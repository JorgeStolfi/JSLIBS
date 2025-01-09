#define _GNU_SOURCE
#include  <stdio.h>
#include <estlogprob.h>
#include <bool.h>
#include <affirm.h>
#include <assert.h>
#include <math.h>
#include <jsfile.h>


#define TINY (1.0e-200)

estima_log_prob_S_G_t *escolhe_estLogPrSG(int num_funcao)
{
  switch(num_funcao) {
    case 0: return &EstLogPrSG_00;
    case 1: return &EstLogPrSG_01;
   // case 2: return &EstLogPrSG_02;
    case 7: return &EstLogPrSG_07;
    case 8: return &EstLogPrSG_08;
    case 9: return &EstLogPrSG_09;
  default: demand(FALSE, "numero invalido da funcao de probabilidade"); return NULL;
  }
}


double ProbSiGialb_canonical(double Si, double Gi, double alb, double sigma)
  { /*
  Supõe que {Pr(Si | Gi,alb) = (1/sigma)*exp(-r*r/2)/sqrt(2*pi)}
      onde
        {Mi = alb*Gi},
        {r = (Si - Mi)/sigma},
      !!! Deveriamos usar uma aproximação que evite calcular {exp} e {log}.
      !!! Compensar quando {Mi} é menor que {3*sigma} ou maor que {1-3*sigma} para integral ser 1. */

    double Mi = alb*Gi;  /* Valor esperado de {Si} */
    double r = (Si - Mi)/sigma;
    double G = (fabs(r) > 5.0 ? 0.0 : (1.0/sigma)*exp(-r*r/2.0)/sqrt(2.0*M_PI));
    return G;
  }
  
double ProbSiGialb_shadow(double Si, double Gi, double alb, double sigma)
  { /* Supõe que {Pr(Si | Gi,alb) = (1/Mi)*(0.5 - 0.5*erf(r))}
      onde
        {Mi = alb*Gi},
        {r = (Si - Mi)/sigma},
      !!! Deveriamos usar uma aproximação que evite calcular {exp} e {log}.
      !!! Compensar quando {Mi} é menor que {3*sigma} ou maor que {1-3*sigma} para integral ser 1. */

    double R = 2.0*sigma;
    double Mi = alb*Gi;  /* Valor esperado de {Si} */

    double r0 = (Si - Mi - R)/sigma;
    double ERF0 = (r0 < -5.0 ? 0.0 : (r0 > +5.0 ? 1.0 : 0.5*(1 + erf(r0))));
    double E = (1/(Mi + R))*(1 - ERF0);
    
    return E ;
  }
  
  
double ProbSiGialb_highlight(double Si, double Gi, double alb, double sigma)
  { /* Supõe que 
        {Pr(Si | Gi,alb) = (1/(1-Mi))*(0.5 + 0.5*erf(r)}
      onde
        {Mi = alb*Gi},
        {r = (Si - Mi)/sigma},
      !!! Deveriamos usar uma aproximação que evite calcular {exp} e {log}.
      !!! Compensar quando {Mi} é menor que {3*sigma} ou maor que {1-3*sigma} para integral ser 1. */

    double R = 2.0*sigma;
    double Mi = alb*Gi;  /* Valor esperado de {Si} */

   
    double r1 = (Si - Mi + R)/sigma;
    double ERF1 = (r1 < -5.0 ? 0.0 : (r1 > +5.0 ? 1.0 : (0.5 + 0.5*erf(r1))));
    double H = (1/(1-Mi+R))*ERF1;

    return H;
  } 

double ProbSiGialb(double Si, double Gi, double alb, double sigma, double omg0, double omg1)
  { /* Supõe que 
        {Pr(Si | Gi,alb) = omg0*E + omg1*H + (1-omg0-omg1)*G}
      onde
        {E = (1/Mi)*(0.5 - 0.5*erf(r))}      # FDP supondo sombra.
        {H = (1/(1-Mi))*(0.5 + 0.5*erf(r))}  # FDP supondo highlight.
        {G = (1/sigma)*exp(-r*r/2)/sqrt(2*pi)}  # FDP supondo iluminação normal.
        {Mi = alb*Gi},
        {r = (Si - Mi)/sigma},
      !!! Deveriamos usar uma aproximação que evite calcular {exp} e {log}.
      !!! Compensar quando {Mi} é menor que {3*sigma} ou maor que {1-3*sigma} para integral ser 1. */

    double E = ProbSiGialb_shadow(Si,Gi,alb,sigma);
    double H = ProbSiGialb_highlight(Si,Gi,alb,sigma);
    double G = ProbSiGialb_canonical(Si,Gi,alb,sigma);
    double Pr = omg0*E + omg1*H + (1-omg0-omg1)*G;
    return Pr;
  }

double LogProbSiGialb(double Si, double Gi, double alb, double sigma, double omg0, double omg1)
  { return log(ProbSiGialb(Si, Gi, alb, sigma, omg0, omg1)); }

double LogPrSGalb(const double SO[], const double GO[], int n, double alb, double sigma, double omg0, double omg1)
  { double logPrSG = 0.0;
    int i;
    for(i = 0; i < n; i++) {
      //  if (i != 1) { continue; } /* !!!!! TIRAR ISTO !!!!!!! */
      double lpri = LogProbSiGialb(SO[i], GO[i], alb, sigma, omg0, omg1);
      logPrSG += lpri;
    }
    return logPrSG;
  }

double LogPrSG_CbP(const double SO[], const double GO[], int n, double sigma, double omg0, double omg1,double K, FILE *arq)
  {     
     
    double R = K*sigma; /* Largura do intervalo central no eixo {S}. */
    
    /* Obtém a lista ordenada dos extremos de intervalos: */
    int nex = 2*n;  /* Número de extremos de intervalos de albedo. */
    int ix[nex];    /* Identificadores de extremos em ordem crescente. */
    double vx[nex]; /* Valores dos extremos em ordem crescente. */
    CriaExtremosDeAlbedo_09(SO, GO, n, R, ix, vx);
    OrdenaExtremosDeAlbedo(ix, vx, nex);

    /* Calcula os três valores de cada fator: */
    double PrLoAlb[n]; /* Valor de Pr(S[i] | G[i],alb)} para {alb < alo[i]}. */ 
    double PrMdAlb[n]; /* Valor de Pr(S[i] | G[i],alb)} para {alb} entre {alo[i],ahi[i]}. */ 
    double PrHiAlb[n]; /* Valor de Pr(S[i] | G[i],alb)} para {alb > ahi[i]}. */ 
    int i;
    for (i = 0; i < n; i++) {

      double Si = SO[i];
      double Gi = GO[i];
    
      /* Limites da fórmula constante-por-partes de {Pr(Si|Gi,alb)}: */
      double alo = (Si < 2*R ? 0.0 : (Si > Gi + R ? 1.0 : (Si - R)/Gi ));
      double ahi = ((Si > Gi - R) || (Si > 1 - 2*R) ? 1.0 : (Si + R)/Gi);
      double amd = 0.5*(alo + ahi);
      
      double lam = 0.57735026918962576450; /* sqrt(1/3) para quadratura de Gauss. */
      double grd = lam * (ahi - amd); 

      PrLoAlb[i] = ProbSiGialb(Si,Gi,0.5*(0 + alo),sigma,omg0,omg1);
      PrHiAlb[i] = ProbSiGialb(Si,Gi,0.5*(1 + ahi),sigma,omg0,omg1);
      double P0 = ProbSiGialb(Si,Gi,amd - grd,sigma,omg0,omg1);
      double P1 = ProbSiGialb(Si,Gi,amd + grd,sigma,omg0,omg1);
      PrMdAlb[i] = 0.5*(P0 + P1);
      if (isnan(PrLoAlb[i]) || isnan(PrHiAlb[i]) || isnan(PrMdAlb[i])) { fprintf(stderr, "#@#! %d\n", i); }
    }
    
    /* Calcula o valor de {Pr(S|G,alb)} para {alb = 0}: */
    double Prob = 1.0;
    for (i = 0; i < n; i++) { Prob *= PrLoAlb[i]; }
    if(isnan(Prob)){
      return 0;
    }
    assert(Prob >= 1.0e-290);
      
    /* Percorre a lista {ix,vx}, atualizando o produto {Prob} e acumulando a integral. */
    double Intg = 0.0; /* Integral de {Pr(S|G,alb)dalb} até o albedo corrente. */
    int j;
    if (arq != NULL) { fprintf(arq, "%3d %8.6f %24.16e\n", -1, 0.0, Prob); }
    for (j = 0; j < nex; j++)
      { /* Neste momento o albedo corrente é {vx[j] - epsilon}. */
        if (j > 0) { assert(vx[j] >= vx[j-1]); }
        if (arq != NULL) { fprintf(arq, "%3d %8.6f %24.16e\n", j, vx[j], Prob); }
        /* Acumula a integral: */
        double dalb = vx[j] - (j == 0 ? 0.0 : vx[j-1]); /* Largura do trecho que terminou. */
        Intg += Prob*dalb;
        /* Atualiza o valor de {Prob}. */
        /* Determina o índice {i} do intervalo ao qual pertence o extremo {ix[j],vx[j]}:  */
        int i = (ix[j] > 0 ? +ix[j] - 1 : -ix[j] - 1);
	// if (i != 1) { continue; } /* !!!!! TIRAR ISTO !!!!!!! */
        if (ix[j] < 0) 
          { /* Extremo inferior de novo intervalo: */
            Prob /= PrLoAlb[i];
            Prob *= PrMdAlb[i];
          }
        else
    	  { /* Extremo superior de um intervalo: */
            Prob /= PrMdAlb[i];
            Prob *= PrHiAlb[i];
    	  } 
	/*   if (isnan(Prob)) { fprintf(stderr, " %d\n", i); }
	     if (TRUE) {  !!!!!! Prob < 1.0e-290 
          fprintf(stderr, "extremo vx[%d] = %24.16e S[%d] = %24.16e  G[%d] = %24.16e\n", j, vx[j], i, SO[i], i, GO[i]);
          fprintf(stderr, "PrLoAlb = %24.16e\n", PrLoAlb[i]);
          fprintf(stderr, "PrMdAlb = %24.16e\n", PrMdAlb[i]);
          fprintf(stderr, "PrHiAlb = %24.16e\n", PrHiAlb[i]);
          fprintf(stderr, "Prob = %24.16e\n", Prob);
	  assert(isfinite(Prob) && (Prob >= 1.0e-290));
 	}*/
        if (arq != NULL) { fprintf(arq, "%3d %8.6f %24.16e\n", j, vx[j], Prob); }
      }
    if (arq != NULL) { fprintf(arq, "%3d %8.6f %24.16e\n", nex, 1.0, Prob); }
    double dalb = 1 - vx[nex-1]; /* Largura do último trecho da integração. */
    Intg += Prob*dalb;
    
    return log(Intg);
  }

double LogPrSG_UnS(const double SO[],const double GO[], int n, double sigma, double omg0, double omg1,int nsteps){
  
  auto double IntPrSGalb(double amin, double amax);
  /* Calcula {\int_{amin}^{amax} Pr(S|G,alb)Pr(alb)dalb}, supondo {alb}
    uniformemente distribuida em [0_1]
    (isto é, {Pr(alb) = 1}), e que {Pr(S|G,alb)} é bem-compurtada no
    intervalo {[amin _ amax]}. */

  auto double PrSGalb(double alb);
  /* Calcula {Pr(S|G,alb)}, para um intervalo pequeno. */
  
  double PrSG = 0.0; /* Valor calculado de {A*Pr(G|S)}. */
  int i;
  for(i = 0; i < nsteps; i++) {	
    double amin = ((double)i)/nsteps;
    double amax = ((double)i+1)/nsteps;
    double intPr = IntPrSGalb(amin, amax);
    PrSG += intPr;
  }
  return log(PrSG);
  

  double IntPrSGalb(double amin, double amax)
    { /* Integral de Gauss com dois pontos. */
      double da = (amax - amin)*(0.5/sqrt(3));
      double f0 = PrSGalb(amin + da);
      double f1 = PrSGalb(amax - da);
      return 0.5*(f0 + f1)*(amax - amin);
    }

  double PrSGalb(double alb)
    { return exp(LogPrSGalb(SO, GO,n, alb, sigma, omg0, omg1)); }
    
}


void CriaExtremosDeAlbedo_08(const double SO[], const double GO[], int n, double R, int ix[], double vx[])
  {
    int nex = 2*n; /* Número de extremos. */
    int i;

    /* double QUASE_UM =   0.999999999; */
    /* double QUASE_ZERO = 0.000000001; */
    
    /* Coloca todos os extremos em {ix[0..2*n-1], vx[0..2*n-1]}: */
    int j = 0;
    for (i = 0; i < n; i++)
      { double Si = SO[i];
        double Gi = GO[i];
        ix[j] = -i-1;  /* Reprsenta {alo[i]} */
	vx[j] = fmin(1.0, fmax(0.0, (Si - R)/(Gi + TINY)));  /* {alo[i]} */
        j++;
        ix[j] = +i+1; /* Representa {ahi[i]} */
	vx[j] = fmin(1.0, fmax(0.0, (Si + R)/(Gi + TINY))); /* {ahi[i]} */
        j++;
      }
    assert(j == nex);
  }

void CriaExtremosDeAlbedo_09(const double SO[], const double GO[], int n, double R, int ix[], double vx[])
  {
    int nex = 2*n; /* Número de extremos. */
    int i;

    /* double QUASE_UM =   0.999999999; */
    /* double QUASE_ZERO = 0.000000001; */
    
    /* Coloca todos os extremos em {ix[0..2*n-1], vx[0..2*n-1]}: */
    int j = 0;
    for (i = 0; i < n; i++)
      { double Si = SO[i];
        double Gi = GO[i];
    
        /* Limites da fórmula constante-por-partes de {Pr(Si|Gi,alb)}: */
        double aLO = (Si < 2*R ? 0.0 : (Si > Gi + R ? 1.0 : (Si - R)/Gi ));
	double aHI = ((Si > Gi - R) || (Si > 1 - 2*R) ? 1.0 : (Si + R)/Gi);
    
       /* Guarda no vetor: */
        ix[j] = -i-1;  /* Reprsenta {alo[i]} */
	vx[j] = aLO;  /* {alo[i]} */
        j++;
        ix[j] = +i+1; /* Representa {ahi[i]} */
	vx[j] = aHI; /* {ahi[i]} */
        j++;
      }
    assert(j == nex);
  }

void OrdenaExtremosDeAlbedo(int ix[], double vx[], int nex)
  {
    /* Insertion sort dos extremos: */
    int i;
    for (i = 0; i < nex; i++)
      { /* Insere {ix[i]} entre {ix[0..i-1]}: */
        int j = i;
        double vxi = vx[i];
        int ixi = ix[i];
        while ((j > 0) && (vxi < vx[j-1])) { vx[j] = vx[j-1]; ix[j] = ix[j-1];  j--; }
        ix[j] = ixi;
        vx[j] = vxi;
      }
    return;
  }

  
/* ESTIMADORES DE VEROSSIMILHANÇA */

double EstLogPrSG_00(const double so[], double Smag, const double go[], double Gmag, int n, double sigma, double omg0, double omg1)
{
//   double so_go = 0.0;
//   int i;
//   for (i = 0; i < n;i++){ so_go += so[i]*go[i]; }
//   double L1 = -n*log(sigma*SQRT_2_PI);
//   double t = Smag/sigma;
//   double L2 = -t*t*(1 - so_go);
//   return L1 + L2;'

//  double alb = (Smag + TINY)/(Gmag + TINY);
//  double sum  = 0.0;
//  int i;
//  for (i = 0; i < n;i++){ 
// 	double term  = (so[i]*Smag) - (alb*go[i]*Gmag);
// 	sum += term*term;
//  }
//  return -sum/(2.0*sigma*sigma);
   
    double dist_eu  = dist_euclid(so,go,n);
    double csi = (dist_eu*Smag/sigma);
    return -0.5*(csi*csi) + 2.0*log(Smag);
  
}

double EstLogPrSG_01(const double so[], double Smag, const double go[], double Gmag, int n, double sigma, double omg0, double omg1)
{
   return  2.0*log(Smag/sigma);
  
}


double EstLogPrSG_07(const double so[], double Smag, const double go[], double Gmag, int n, double sigma, double omg0, double omg1)
{
  /* Calcula {Pr(S|G) = \int_0^1 Pr(S|G,alb)Pr(alb)dalb} onde {alb} é o albedo hipotético 
    do pixel da cena.  O problema é que {Pr(S|G,alb)} é cheia de picos estreitos onde
    está toda a probabilidade.  Chutamos que esses picos podem ter 
    largura mínima {sigma/sqrt(n)}.  Portanto vamos dividir o intervalo de integração
    em {m} intervalos dessa largura, e usar a formula
    {C_PrG_S = \sum_{i=0}^{m-1} \int_{i/m}^{(i+1)/m} Pr(S|G,alb)Pr(alb)dalb} */

  int nsteps = (int)ceil(sqrt(n)/sigma);
  assert(nsteps >= 1);
  double SO[n];
  double GO[n];
  int i;
  for(i = 0; i < n;i++){ SO[i] = so[i]*Smag; GO[i] = go[i]*Gmag; }
  return LogPrSG_UnS(SO,GO,n,sigma,omg0,omg1,nsteps);
    
}



double EstLogPrSG_08(const double so[], double Smag, const double go[], double Gmag, int n, double sigma, double omg0, double omg1)
{
  /* Estima {Pr(S|G)} pela integral {\int_{aMIN}^{aMAX}
    Pr(S|G,alb)Pr(alb)dalb}, onde {alb} é o albedo hipotético do pixel
    da cena e {[aMIN _ aMAX]} é um intervalo de albedos onde se
    acredita que a integral seja significativa.

    Escolhemos {aMIN} e {aMAX} procurando primeiro o albedo {amed}
    onde há o maior número de pares {SO[i],GO[i]} ativos.  Um par
    {SO[i],GO[i]} é ativo para um albedo {alb} se {\abs{SO[i] -
    alb*GO[i]} <= K*sigma}, onde {K} é um parametro a ajustar
    (digamos, {K=3}).  Uma vez determinado {amed}, escolhemos {aMIN =
    amed - C} e {aMAX = amed + C}, onde {C} é outro
    parâmetro a ajustar. */

  double K = 4.0; /* Chute. */
  
  /* Obtém a lista ordenada dos extremos de intervalos: */
  double SO[n];
  double GO[n];
  int k ;
  for(k = 0; k < n;k++){ GO[k] = go[k]*Gmag; SO[k]= so[k]*Smag;}
  
  int nex = 2*n; /* Número de extremos de intervalos de albedo. */
  int ix[nex];    /* Identificadores de extremos em ordem crescente. */
  double vx[nex]; /* Valores dos extremos em ordem crescente. */
  CriaExtremosDeAlbedo_09(SO, GO, n, K*sigma, ix, vx);
  OrdenaExtremosDeAlbedo(ix, vx, nex);

  /* Procura o albedo {amed}: */
  int n_ativos = 0;           /* Número de intervalos que incluem o albedo corrente. */
  double sumGi2 = 0.0;        /* Soma de {G[i]^2} para esses intervalos. */

  int max_ativos = 0;         /* Número máximo de pares ativos para qualquer albedo. */
  double melhor_amed = -1.0;  /* Albedo para o qual esse máximo foi atingido. */
  double melhor_sumGi2 = 0.0; /* Valor correspondente de {sumGi2}. */
  
  int j;
  for (j = 0; j < nex; j++)
    { /* Neste momento o albedo corrente é {vx[j] - epsilon}. */
      if (j > 0) { assert(vx[j] >= vx[j-1]); }
      /* Determina o índice {i} do intervalo ao qual pertence o extremo {ix[j],vx[j]}:  */
      int i = (ix[j] > 0 ? +ix[j] - 1 : -ix[j] - 1);
      double Gi = Gmag*go[i];
      if (ix[j] < 0) 
        { /* Extremo inferior de novo intervalo: */
          n_ativos++; 
          sumGi2 += Gi*Gi;
        }
      else
	{ /* Extremo superior de um intervalo: */
          if (n_ativos > max_ativos)
	    { /* O trecho de albedo de onde estamos saindo é o melhor que já encontramos: */
              max_ativos = n_ativos;
	      assert(j > 0);
	      melhor_amed = (vx[j-1] + vx[j])/2;
	      melhor_sumGi2 = sumGi2;
	    }
          n_ativos--;
          sumGi2 -= Gi*Gi;
	}
    }
  assert(n_ativos == 0);
  bool_t debug_dist = (max_ativos >= 9);

  /* Estima o desvio padrão do produto das gaussianas ativas em {melhor_amed}: */
  double sigma_prod = sigma/sqrt(melhor_sumGi2 + TINY);

  if (debug_dist)
    {/* fprintf(stderr, "melhor albedo = %5.3f ativos = %2d sigma_prod = %7.5f", melhor_amed, max_ativos, sigma_prod);*/ }
  
  /* Escolhe o intervalo de integração: */
  double C = 3.0*sigma_prod; /* Outro chute. */
  double aMIN = fmin(1.0, fmax(0.0, melhor_amed - C));
  double aMAX = fmin(1.0, fmax(0.0, melhor_amed + C));
  
  /* Calcula o número de passos: */

  int nsteps = (int)ceil((aMAX-aMIN)/sigma_prod);
  assert(nsteps >= 1);
  
  if (debug_dist)
    { /*fprintf(stderr, "  intg = [ %5.3f _ %5.3f ] nsteps = %3d\n", aMIN, aMAX, nsteps);*/ }

  auto double IntPrSGalb(double amin, double amax);
  /* Calcula {\int_{amin}^{amax} Pr(S|G,alb)Pr(alb)dalb}, supondo {alb}
    uniformemente distribuida em [0_1]
    (isto é, {Pr(alb) = 1}), e que {Pr(S|G,alb)} é bem-compurtada no
    intervalo {[amin _ amax]}. */

  auto double PrSGalb(double alb);
  /* Calcula {Pr(S|G,alb)}, para um intervalo pequeno. */
  
  double PrSG = 0.0; /* Valor calculado de {A*Pr(G|S)}. */
  int i;
  for(i = 0; i < nsteps; i++) {	
    double tmin = ((double)i)/nsteps;
    double amin = (1-tmin)*aMIN + tmin*aMAX;
    double tmax = ((double)i+1)/nsteps;
    double amax = (1-tmax)*aMIN + tmax*aMAX;
    double intPr = IntPrSGalb(amin, amax);
    PrSG += intPr;
  }
  return log(PrSG);
  

  double IntPrSGalb(double amin, double amax)
    { /* Integral de Gauss com dois pontos. */
      double da = (amax - amin)*(0.5/sqrt(3));
      double f0 = PrSGalb(amin + da);
      double f1 = PrSGalb(amax - da);
      return 0.5*(f0 + f1)*(amax - amin);
    }

  double PrSGalb(double alb)
    { return exp(LogPrSGalb(SO, GO, n, alb, sigma, omg0, omg1)); }
  
}


/* K utilizado na EstLogProb_09*/

double EstLogPrSG_09(const double so[], double Smag, const double go[], double Gmag, int n, double sigma, double omg0, double omg1)
{
  /* Calcula {Pr(S|G) = \int_0^1 Pr(S|G,alb)Pr(alb)dalb} onde {alb} é o albedo hipotético 
    do pixel da cena.  Nesta versão, usamos uma aproximação constante-por-partes
    de cada fator {Pr(S[i]|G[i],alb)}. O intervalo {[0_1]} 
    é dividido em três partes, sendo a parte do meio um intervalo
    {[alo[i] _ ahi[i]]}, definido como na função {OrdenaExtremosDeAlbedo}
    com {R = K_09*sigma} para algum {K_09} 

    Uma vez que cada fator é constante-por-partes, o mesmo vale para 
    o produto {Pr(S|G,alb)}.  Basta então percorrer a lista ordenada 
    dos extremos dos intervalos {alo[i],ahi[i]}, calculando o valor
    em cada intervalo.  Custo: ordenação de {2*n} valores, mais
    constante vezes {n}.  */
  
  double SO[n];
  double GO[n];
  int i;
  for(i = 0; i < n;i++){ SO[i] = so[i]*Smag; GO[i] = go[i]*Gmag; }
  return LogPrSG_CbP(SO, GO, n, sigma, omg0, omg1, K_09, NULL);
}


void gera_plot_Si_Gi (
    normal_debug_opts_t *dbopt,
    Tabela* tab, 
    const double SO[], 
    int lin[],
    int m
  ) 
  {
    char *nome_arq = NULL;
    char *nome_arq = jsprintf("%s_%d_%04d_%04d_Si_Gi.txt", dbopt->prefixo, dbopt->c, dbopt->hp, dbopt->vp);
    FILE *arq = open_write(nome_arq, TRUE);
    int n = get_num_luzes(tab);
    int i;
    for (i = 0; i < n; i++) {
      fprintf(arq, "%3d %9.6f", i, SO[i]);
      int j;
      for (j = 0; j <= m; j++) {
        if (j == m) {
          /* Auto-gabarito: */
          fprintf(arq, " %9.6f", SO[i]);
        } else {
          int linha = lin[j];
          const double *go = get_intdir(tab,linha);
          double Gmag = get_intmag(tab,linha);
          fprintf(arq, " %9.6f", Gmag*go[i]);
	}
      }
      fprintf(arq, "\n");
    }
    fclose(arq);
    free(nome_arq);
    fprintf(stderr, "pronto.\n");
  }


void gera_plot_LogPrSGalb (
    normal_debug_opts_t *dbopt,
    Tabela* tab, 
    const double SO[], 
    double sigma, 
    double omg0, 
    double omg1,
    int lin[],
    int m
  )
  {
    char *nome_arq = NULL;
    char *nome_arq = jsprintf("%s_%d_%04d_%04d_Pr_S_a.txt", dbopt->prefixo, dbopt->c, dbopt->hp, dbopt->vp);
    FILE *arq = open_write(nome_arq, TRUE);
    int n = get_num_luzes(tab);
    
    double so[n];
    double Smag;
    extrai_assinatura(SO, so, &Smag, n);
    
    int nsteps = (int)ceil(3.0*sqrt(n)/sigma);
    int k;
    for (k = 0; k <= nsteps; k++) {
      double albedo = ((double)k)/((double)nsteps);
      fprintf(arq, "%5d %9.6f ", k, albedo);
      int j;
      for (j = 0; j <= m; j++) {
        double logPrSG;
	double GO[n];
	int kk;
        if (j == m) {
	  /* Auto-gabarito: */
	  for(kk = 0; kk < n; kk++){ GO[kk] = SO[kk];}
	} else {
          int linha = lin[j];
          const double *go = get_intdir(tab,linha);
          double Gmag = get_intmag(tab,linha);
	  for(kk = 0; kk < n; kk++){ GO[kk] = Gmag*go[kk];}
	}
	logPrSG =  LogPrSGalb(SO, GO, n, albedo, sigma, omg0, omg1);
        fprintf(arq, " %+9.6lf", logPrSG);
      }
      fprintf(arq, "\n");
    }
    fclose(arq);
    free(nome_arq);
    fprintf(stderr, "pronto.\n");
  }



