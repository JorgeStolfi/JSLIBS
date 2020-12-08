#ifndef btc_bubble_t_H
#define btc_bubble_t_H

/* Parameters of a bubble model. */
/* Last edited on 2015-04-29 23:49:54 by stolfilocal */

typedef struct btc_bubble_t /* Parameters of one bubble in the model. */
  { double coef;     /* Fitted coefficient. */
    int id_ini_sg;   /* Day index of start of relevant period, {0..nd-1}. */
    double rt_up;    /* Daily rate of increase during rally. */
    int id_fin_up;   /* Day index of end of rally phase {0..nd-1}. */
    int wd_plat;     /* Width of plateau, in days (non-negative). */
    double rt_dn;    /* Daily rate of decrease during decay. */
    int id_fin_sg;   /* Day index of end of relevant period, {0..nd-1}. */
    char* tag;       /* Name tag of bubble. */
    char *color;     /* RGB color for plotting, as six hex digits. */
  } btc_bubble_t;

#endif
