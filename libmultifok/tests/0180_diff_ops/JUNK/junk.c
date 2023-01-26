/* Last edited on 2023-01-25 06:06:29 by stolfi */

    argparser_get_keyword(pp, "-winSize");
    o->winSize = (int32_t)argparser_get_next_int(pp, 3, 99);  
    if (o->winSize % 2 != 1) { argparser_error(pp, "window size must be odd"); }

