#ifndef tf_kalman_H 
#define tf_kalman_H 

#define _GNU_SOURCE
#include <stdio.h>

typedef struct _kalman_parameters_t {
    int order;
    double avg;
    double dev;
    double *coeffs;
} *kalman_parameters_t;

kalman_parameters_t tf_kalman_new_parameters (int order);

kalman_parameters_t tf_kalman_read_parameters (FILE *f);

void tf_kalman_write_parameters (kalman_parameters_t K, FILE *f);

void tf_kalman_show_parameters (kalman_parameters_t K, FILE *f);

void tf_kalman_predict_parameter (kalman_parameters_t K, double data[], double *avg, double *dev);
/* Predicts the value of a parameter given the kalman model {K} and the
   previous values {data[0..K->order-1]}. Assumes that {data[0]} is the
   most recent value. Returns the expected value in {*avg}, and the 
   standard deviation (uncertainty) of the prediction in {*dev}. */

void tf_kalman_free_parameters (kalman_parameters_t K);


#endif
