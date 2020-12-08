    /* Choose a suitable {maxval} for it: */
    uint32_t maxmaxval = uint16_image_MAX_SAMPLE;
    sample_conv_choose_maxval(chns, imaxval, maxmaxval, &omaxval);
    assert(maxmaxval <= 65535);  /* So that the conversion does not overflow. */

    bool_t debug_conv = FALSE;
    
    for (int32_t row = 0; row < rows; row++)
      { 
        uint16_t *pnm_P = img->smp[row];
        for (int32_t col = 0; col < cols; col++)
          { for (int32_t chn = 0; chn < chns; chn++)
              {

                /* Convert to PNM sample {osmp}: */
                uint32_t den = imaxval[chn];
                uint32_t osmp;
                if (den == omaxval) 
                  { /* No conversion needed: */
                    osmp = ismp;
                  }
                else
                  { /* Linear scale with rounding: */
                    osmp = (ismp*omaxval + (den >> 1))/den;
                  }
                if (debug_conv) { fprintf(stderr, "ismp = %u  osmp = %u\n", (uint32_t)ismp, (uint32_t)osmp); }
                assert(osmp <= img->maxval);
                pnm_P++;
              }
          }
      }    
