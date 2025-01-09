/* Resolu��o de sistemas esparsos de equa��es. */
/* Last edited on 2011-03-08 17:40:08 by stolfi */ 

#ifndef sistema_H
#define sistema_H

#define MAXCOEFS 5

typedef struct equacao_t
  { long int ind[5];  /* {ind[k]} � o �ndice de uma inc�gnita que entra na equa��o. */
    double coef[5];   /* Coeficiente dessa inc�gnita. */
    double indep;     /* Termo constante da equa��o. */
  } equacao_t;
  /* Um registro do tipo {equacao_t} representa uma equa��o linear 
     envolvendo at� {MAXCOEFS} inc�gnitas, da forma
      {SOMA {coef[k]*Z[ind[k]] : k = 0..MAXCOEFS-1} == indep},
     onde {Z[0..N-1]} � o vetor com todas as inc�gnitas.
     Entradas vazias s�o indicadas por {ind[k] == -1}. */

typedef struct sistema_t
  { long int N;
    equacao_t* eq;
  } sistema_t;
  /* Um {sistema_t} representa {N} equa��es lineares {eq[0..N-1]} sobre
    {N} inc�gnitas {Z[0..N-1}. */

sistema_t* cria_sistema(long int N);
  /* Cria um sistema linear com {N} equa��es para {N} inc�gnitas. 
     N�o inicializa nenhum coeficiente, �ndice ou termo independente. */

void resolve_sistema
  ( sistema_t* S, 
    double *Z, 
    long int max_iter, 
    double tol,
    int para, 
    int szero
  );
  /* Resolve o sistema {S}. A solu��o � armazenada no vetor {Z[0..N-1}}
     onde {N = S->N}. 
     
     Usa o m�todo iterativo de Gauss-Jordan, e portanto sup�e que a
     inc�gnita {Z[i]} aparece na equa��o {S->eq[i]}, para {i =
     0..S->N-1}, com coeficiente ``suficientemente grande''. 
     O chute inicial � o valor de {Z} na entrada do procedimento. Executa
     no m�ximo {max_iter} itera��es, mas termina quando duas itera��es
     consecutivas n�o mudam nenhuma inc�gnita mais do que a toler�ncia
     {tol}.

     Se {para == 1}, usa forma ``paralela'' da itera��o de Gauss, sen�o 
     usa a forma sequencial (Gauss-Jordan). Se {szero == 1},
     ajusta a solu��o de modo a ter soma zero, ap�s cada itera��o. */
 
double calcula_dif_total(double* Zvelho, double* Z, long int N);
  /* Retorna o m�ximo de {abs(Zvelho[i]-Z[i])}, para {i = 0..N-1}. */

#endif
