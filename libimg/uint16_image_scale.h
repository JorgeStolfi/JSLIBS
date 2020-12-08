  
#define uint16_image_io_png_read_SAMPLE_CONV_INFO \
  "If the file max values {imaxval[0..NC-1]} are all equal, the max sample" \
  " value {img.maxval} will be that value, and the samples will" \
  " be as stored in the file.  Otherwise the procedure will choose" \
  " suitable {img.maxval}, not less than any of them and" \
  " scale the input samples in each channel {i} by the" \
  " fraction {s=img.maxval/imaxval[i]} so that 0 remains 0, and" \
  " {imaxval[i]} becomes {img.maxval}.  This conversion" \
  " may involve some rounding that causes the differeces" \
  " between successive input sample values to vary" \
  " between {floor(s)} and {ceil(s)}."
