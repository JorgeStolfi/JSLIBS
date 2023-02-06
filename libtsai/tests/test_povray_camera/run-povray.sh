#! /bin/bash

# Runs POV-ray on the "camera_dat.inc" produced by test_povray_camera,
# and displays the result.

name="main"
outfile=out/${name}.png

/bin/rm -f ${outfile}
nice povray +K0.0000 \
            +FN +Q9 +MB1 \
            +W640 +H480 \
            +AM1 +A0.0 +R2 \
            +D +SP32 +EP4 \
            +Lpovray \
            +Lpovray/tt-fonts \
            +Ipovray/${name}.pov \
            +O${outfile} \
          3>&2 | tee out/${name}.log

display ${outfile} &
