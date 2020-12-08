#ifndef btc_bubble_nl_opt_add_integer_parameter_bias_H
#define btc_bubble_nl_opt_add_integer_parameter_bias_H

/* Favor the guessed vaues of integer BTC price bubble parameters. */
/* Last edited on 2015-04-20 00:47:27 by stolfilocal */

double btc_bubble_nl_opt_add_integer_parameter_bias(double Q, int npi, int pi_a[], int pi_b[], double alpha);
  /* Adds to the goal function value {Q} a small bias term equal to 
    {alpha} times the square of the distance between the integer vector {pi_a[0..npi]} and
    the reference vctor {pi_b[0..npi-1]}. */
    

#endif
