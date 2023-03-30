#! /bin/bash
# Last edited on 2008-10-18 10:32:05 by stolfi

# To be used internally by {plot_spectrum_exact}.
# Reads the temporary file and excludes the {EQL} entries.

gawk '($2+0 != 0){ print; }'
