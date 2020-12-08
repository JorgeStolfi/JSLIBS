#! /bin/bash
# Last edited on 2012-01-16 02:57:21 by stolfilocal

# Mask in/RAW/white-msk.png was created with gimp (400x400)

# Raw files from photographs (various sizes):

convert in/RAW/gras2-big.png -crop '420x300+33+12' -resize 'x100'     in/gras2-img.ppm
convert in/RAW/white-msk.png                       -resize '!140x100' in/gras2-msk.pgm

convert in/RAW/grass-big.png -crop '360x300+20+50' -resize 'x100'     in/grass-img.ppm
convert in/RAW/white-msk.png                       -resize '!120x100' in/grass-msk.pgm

# Raw files created with gimp (400x400, centered between pixels):

convert in/RAW/janus-big.png -resize 'x100'  in/janus-img.ppm
convert in/RAW/white-msk.png -resize 'x100'  in/janus-msk.pgm

# convert in/JUNK/jdiag-big.png -resize 'x100'  in/jdiag-img.ppm
# convert in/RAW/white-msk.png  -resize 'x100'  in/jdiag-msk.pgm
# 
# convert in/JUNK/odots-big.png -resize 'x100'  in/odots-img.ppm
# convert in/RAW/white-msk.png  -resize 'x100'  in/odots-msk.pgm
# 
# convert in/JUNK/sqres-big.png -resize 'x100'  in/sqres-img.ppm
# convert in/RAW/white-msk.png  -resize 'x100'  in/sqres-msk.pgm
# 
# convert in/JUNK/trigs-big.png -resize 'x100'  in/trigs-img.ppm
# convert in/RAW/white-msk.png  -resize 'x100'  in/trigs-msk.pgm

# Raw files created with make_bullseye (512x512, centered between pixels):

convert in/RAW/bueye-big.png -crop '200x200+156+156' -resize 'x100'  in/bueye-img.ppm
convert in/RAW/white-msk.png                         -resize 'x100'  in/bueye-msk.pgm
    
# Raw files created with make_test_slope_maps (513x513, centered at pixel)

convert in/RAW/bdots-big.png -crop '505x505+4+4' -resize 'x404' -crop '400x400+2+2' -resize 'x100'  in/bdots-img.ppm
convert in/RAW/white-msk.png                                                        -resize 'x100'  in/bdots-msk.pgm

convert in/RAW/bpent-big.png -crop '505x505+4+4' -resize 'x404' -crop '400x400+2+2' -resize 'x100'  in/bpent-img.ppm
convert in/RAW/white-msk.png                                                        -resize 'x100'  in/bpent-msk.pgm

convert in/RAW/btria-big.png -crop '505x505+4+4' -resize 'x404' -crop '400x400+2+2' -resize 'x100'  in/btria-img.ppm
convert in/RAW/white-msk.png                                                        -resize 'x100'  in/btria-msk.pgm

convert in/RAW/ptpyr-big.png -crop '505x505+4+4' -resize 'x404' -crop '400x400+2+2' -resize 'x100'  in/ptpyr-img.ppm
convert in/RAW/white-msk.png                                                        -resize 'x100'  in/ptpyr-msk.pgm

convert in/RAW/wavec-big.png -crop '505x505+4+4' -resize 'x404' -crop '400x400+2+2' -resize 'x100'  in/wavec-img.ppm
convert in/RAW/white-msk.png                                                        -resize 'x100'  in/wavec-msk.pgm
