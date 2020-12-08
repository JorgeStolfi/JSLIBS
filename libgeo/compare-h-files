#! /bin/csh -f
# Last edited on 2002-12-29 23:34:54 by stolfi

set filea = "$1"; shift;
set fileb = "$1"; shift;

foreach f ( ${filea} ${fileb} )
  cat $f \
    | sed \
        -e 's:r[0-9m]x[0-9n]:rNxN:g' \
        -e 's:r[0-9n]:rN:g' \
    > /tmp/$f
end
diff -b /tmp/${filea} /tmp/${fileb}
rm /tmp/${filea} /tmp/${fileb}
