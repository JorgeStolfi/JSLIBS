#! /bin/bash
# Last edited on 2024-10-26 07:06:06 by stolfi

src1="projects/image-collections/texture-bank/ppm-400x400"

imgs1=( \
    cloth-84 decor-82 falls-87 fiber-01 gecko-02 grave-02 gravl-87 husks-01 india-81 leafy-09 \
    mosic-96 pebbl-87 rebar-03 sqirt-01 tampi-86 tarmc-81 trees-07 wavys-33 wheat-81 \
  )
for img in ${imgs1[@]}; do 
  img_src="${src1}/${img}.ppm"
  img_dst="./${img/-/}.png"
  convert ${img_src} -colorspace Gray -normalize -resize '50%' ${img_dst}
done

src2="../../../../libimg/tests/330_test_images/out"
imgs2=( grittie )
for img in ${imgs2[@]} ; do
  img_src="${src2}/test-1024x1024-1-${img}.png"
  img_dst="./${img}.png"
  cp -a ${img_src} ${img_dst}
done
  
src3="../../../../libimg/tests/425_noise/out"
imgs3=( img-0.000000-0.000000-cmF-sqF:noise01 )
for img in ${imgs3[@]} ; do
  img_src="${src3}/png-1024x1024/${img/:*/}.png"
  img_dst="./${img/*:/}.png"
  cp -a ${img_src} ${img_dst}
done
