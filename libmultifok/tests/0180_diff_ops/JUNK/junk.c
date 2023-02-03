/* Last edited on 2023-02-01 18:37:39 by stolfi */

    argparser_get_keyword(pp, "-winSize");
    o->winSize = (int32_t)argparser_get_next_int(pp, 3, 99);  
    if (o->winSize % 2 != 1) { argparser_error(pp, "window size must be odd"); }

void mfdo_read_term_index_table
  ( char *fname, 
    int32_t NB,
    char *belName[],
    int32_t *NP_P,
    multifok_term_prod_t **prix_P,
    int32_t *NT_P,
    char ***termName_P
  );
  /* Calls {multifok_term_read_index_table} on the file "{fname}". */

void mfdo_read_term_index_table
  ( char *fname, 
    int32_t NB,
    char *belName[],
    int32_t *NP_P,
    multifok_term_prod_t **prix_P,
    int32_t *NT_P,
    char ***termName_P
  )
  { bool_t verbose = TRUE;
    FILE *rd = open_read(fname, verbose);
    multifok_term_read_index_table(rd, NB, belName, NP_P, prix_P, NT_P, termName_P, verbose);
    fclose(rd);
  }

    fprintf(stderr, "ignored %d pixels for low sharpness\n", NL_blur); 

int32_t NL_blur = 0;    /* Number of pixels discarded because of low {sharp}. */
    
            double sharp_min = 0.2;      /* Ignore pixels with sharpness below this. */

 if (sharp < sharp_min)
              { NL_blur++;
                mkval = 0.25;
              }
            else
