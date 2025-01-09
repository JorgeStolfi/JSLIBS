/* Veja sistema.h. */
/* Last edited on 2025-01-07 08:12:13 by stolfi */ 

#include <sistema.h>

#include <assert.h>
#include <stdio.h>
#include <stdlib.h>

sistema_t* cria_sistema(long int N)
  {
    sistema_t* S = (sistema_t*)malloc(sizeof(sistema_t));
    S->N = N;
    S->eq = (equacao_t*)malloc(sizeof(equacao_t)*N);
    return S;
  }

void resolve_sistema
  ( sistema_t* S, 
    double *Z, 
    long int max_iter, 
    double tol,
    int para, 
    int szero
  )
  {
    double dif;
    long int N = S->N;
    long int i;

    /* Solu�o calculada na itera�o anterior: */
    double *Zvelho = (double*)malloc(sizeof(double)*N);
    
    dif = 100.0;
    
    long int iter;
    for (iter = 0; iter < max_iter; iter++)
      {
        /* Mais uma passada, recalcula todas as inc�nitas {Z[i]}: */
        for (i = 0; i < N; i++)
          {
            /* Guarda solu�o corrente em {Zvelho}: */
            Zvelho[i] = Z[i];
            /* Pega a equa�o {i}: */
            equacao_t *eqi = &(S->eq[i]);
            /* Calcula {Z[i]} usando a equa�o {i}: */
            double soma = eqi->indep; /* Lado direito da equa�o. */
            double coef_i = 0.0; /* Coeficiente de {Z[i]} na equa�o. */
            int k;
            for(k = 0; k < MAX_COEFFS; k++)
              {
                /* Pega mais uma inc�nita {Z[j] que entra na equa�o {i}: */
                int j = eqi->ind[k];
                double coef_j = eqi->coef[k];
                if (j == i) 
                  { /* Este termo usa {Z[i]}, guarde o coeficiente: */
                    coef_i = coef_j;
                  }
                else if ((j >= 0) && (j < N) && (coef_j != 0))
                  {
                    /* A inc�nita {Z[j]} �distinta de {Z[i]}. */
                    /* Pega o valor apropriado (velho ou novo) de {Z[j]}: */
                    double Zj = (para && (j < i) ? Zvelho[j] : Z[j]);
                    /* Subtrai o termo do lado direito da equa�o: */
                    soma = soma - Zj * coef_j;
                  }
              }
            
            /* �bom que a equa�o {eqi} dependa de {Z[i]}: */
            assert(coef_i != 0.0); 
            
            /* Resolve a equa�o: */
            Z[i] = soma / coef_i;
          }
          
        if (szero)
          { /* Normaliza a m�ia em zero: */
            double soma = 0;
            for (i = 0; i < N; i++) { soma += Z[i]; }
            double media = soma/N;
            for (i = 0; i < N; i++) { Z[i] -= media; }
          }
        iter++;
        dif = calcula_dif_total(Zvelho, Z, N);
        if (dif <= tol) { return; }
      }
  }

double calcula_dif_total(double* Zvelho, double* Z, long int N)
  {
    int i;
    double maior_dif = 0;
    long int maior_i = 0;
    for(i = 0; i< N; i++)
      {
        double di = Z[i] - Zvelho[i];
        if(di < 0) { di = - di; }
        if(di > maior_dif) { maior_i = i; maior_dif = di; }
      }
    return maior_dif;
  }


