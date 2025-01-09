/* Resolução de sistemas esparsos de equações. */
/* Last edited on 2011-03-08 17:40:08 by stolfi */ 

#ifndef sistema_H
#define sistema_H

#define MAXCOEFS 5

typedef struct equacao_t
  { long int ind[5];  /* {ind[k]} é o índice de uma incógnita que entra na equação. */
    double coef[5];   /* Coeficiente dessa incógnita. */
    double indep;     /* Termo constante da equação. */
  } equacao_t;
  /* Um registro do tipo {equacao_t} representa uma equação linear 
     envolvendo até {MAXCOEFS} incógnitas, da forma
      {SOMA {coef[k]*Z[ind[k]] : k = 0..MAXCOEFS-1} == indep},
     onde {Z[0..N-1]} é o vetor com todas as incógnitas.
     Entradas vazias são indicadas por {ind[k] == -1}. */

typedef struct sistema_t
  { long int N;
    equacao_t* eq;
  } sistema_t;
  /* Um {sistema_t} representa {N} equações lineares {eq[0..N-1]} sobre
    {N} incógnitas {Z[0..N-1}. */

sistema_t* cria_sistema(long int N);
  /* Cria um sistema linear com {N} equações para {N} incógnitas. 
     Não inicializa nenhum coeficiente, índice ou termo independente. */

void resolve_sistema
  ( sistema_t* S, 
    double *Z, 
    long int max_iter, 
    double tol,
    int para, 
    int szero
  );
  /* Resolve o sistema {S}. A solução é armazenada no vetor {Z[0..N-1}}
     onde {N = S->N}. 
     
     Usa o método iterativo de Gauss-Jordan, e portanto supõe que a
     incógnita {Z[i]} aparece na equação {S->eq[i]}, para {i =
     0..S->N-1}, com coeficiente ``suficientemente grande''. 
     O chute inicial é o valor de {Z} na entrada do procedimento. Executa
     no máximo {max_iter} iterações, mas termina quando duas iterações
     consecutivas não mudam nenhuma incógnita mais do que a tolerância
     {tol}.

     Se {para == 1}, usa forma ``paralela'' da iteração de Gauss, senão 
     usa a forma sequencial (Gauss-Jordan). Se {szero == 1},
     ajusta a solução de modo a ter soma zero, após cada iteração. */
 
double calcula_dif_total(double* Zvelho, double* Z, long int N);
  /* Retorna o máximo de {abs(Zvelho[i]-Z[i])}, para {i = 0..N-1}. */

#endif
