#define _GNU_SOURCE
#include <stdio.h>
#include <stdint.h>
#include <stdlib.h>
#include <affirm.h>
#include <tf_kalman.h>

kalman_parameters_t tf_kalman_new_parameters (int32_t order)
{
   kalman_parameters_t K = (kalman_parameters_t)malloc(sizeof(struct _kalman_parameters_t));

   K->order = order;
   K->avg = 0.0;
   K->dev = 1.0;
   K->coeffs = (double *)malloc(order * sizeof(double));
   int32_t i;
   for (i = 0; i < K->order; i++) { K->coeffs[i] = 0.0; }
   return K;
}

kalman_parameters_t tf_kalman_read_parameters (FILE *f)
{
   int32_t order;

   fscanf(f, "%d",  &order);
   kalman_parameters_t K = tf_kalman_new_parameters (order);

   fscanf(f, "%lf", &(K->avg));
   fscanf(f, "%lf", &(K->dev));
   demand(K->dev >= 0.0, "invalid Kalman deviation");

   int32_t i;
   for (i = 0; i < K->order; i++) {
     fscanf(f, "%lf", &(K->coeffs[i]));  
   }

   return K;
}

void tf_kalman_write_parameters (kalman_parameters_t K, FILE *f)
{
   
   fprintf(f, "%d",  K->order);
   fprintf(f, " %24.16le",  K->avg);
   fprintf(f, " %24.16le",  K->dev);

   int32_t i;
   for (i = 0; i < K->order; i++) {
     fprintf(f, " %24.16le", K->coeffs[i]);  
   }
   fprintf(f, "\n");
   fflush(f);
}

void tf_kalman_show_parameters (kalman_parameters_t K, FILE *f)
{
   fprintf(f, "Kalman order = %d",  K->order);
   fprintf(f, "  model Z[k] = (");
   int32_t i;
   for (i = 0; i < K->order; i++) {
     fprintf(f, " %+lf*Z[k-%d]", K->coeffs[i], i+1);  
   }
   fprintf(f, " %+lf )",  K->avg);
   fprintf(f, " +/- %lf\n",  K->dev);
   fflush(f);
}

void tf_kalman_predict_parameter (kalman_parameters_t K, double data[], double *avg, double *dev)
{
   int32_t i;
   double sum = 0.0;   
   for (i = 0; i < K->order; i++) {
      sum += K->coeffs[i] * data[i];
   }
   *avg = sum + K->avg;
   *dev = K->dev;
}

void tf_kalman_free_parameters (kalman_parameters_t K) 
{
  free(K->coeffs);
  free(K);
}
















