# Last edited on 2020-10-11 15:26:05 by jstolfi

# Creates test PNG files from master images

MAXVALS_SMA := 00001 00003 00015 00255
MAXVALS_ALL := 00001 00003 00015 00255 65535
MAXVALS_BIG := 00255 65535

PREFIX := test

VARIANTS_GRY_A0 := ${foreach mxv,${MAXVALS_ALL},${addsuffix -${mxv},GRY-A0}}
VARIANTS_GRY_A1 := ${foreach mxv,${MAXVALS_BIG},${addsuffix -${mxv},GRY-A1}}

VARIANTS_RGB_A0 := ${foreach mxv,${MAXVALS_BIG},${addsuffix -${mxv},RGB-A0}}
VARIANTS_RGB_A1 := ${foreach mxv,${MAXVALS_BIG},${addsuffix -${mxv},RGB-A1}}

VARIANTS_MAP_A0 := ${foreach mxv,${MAXVALS_SMA},${addsuffix -${mxv},MAP-A0}}
VARIANTS_MAP_A1 := ${foreach mxv,${MAXVALS_SMA},${addsuffix -${mxv},MAP-A1}}

VARIANTS := \
  ${VARIANTS_GRY_A0} ${VARIANTS_GRY_A1} \
  ${VARIANTS_RGB_A0} ${VARIANTS_RGB_A1} \
  ${VARIANTS_MAP_A0} ${VARIANTS_MAP_A1} 
 
PNG_OUT_FILES :=  ${addsuffix .png,${addprefix ${PREFIX}-,${VARIANTS}}}

.PHONY:: all clean testfiles 

all: check

check:  clean testfiles

SHOW := display -title '%f'
        
testfiles: ${PNG_OUT_FILES}

clean::
	/bin/rm -f ${PNG_OUT_FILES}
  
# ====================================================================
# Data files:

MASTER_A1 := ../master/${PREFIX}-master-alpha.png
MASTER_A0 := ../master/${PREFIX}-master-plain.png

# ----------------------------------------------------------------------
# GRY colorspace, no alpha
        
${PREFIX}-GRY-A0-00001.png: ${MASTER_A0}
	convert ${MASTER_A0} -depth 16 -colorspace Gray -resize '320x240' +repage -depth 1 -define png:bit-depth=1  $@
	identify $@ ; display -title '%f' $@

${PREFIX}-GRY-A0-00003.png: ${MASTER_A0}
	convert ${MASTER_A0} -depth 16 -colorspace Gray -resize '320x240' +repage -depth 2 -define png:bit-depth=2  $@
	identify $@ ; display -title '%f' $@

${PREFIX}-GRY-A0-00015.png: ${MASTER_A0}
	convert ${MASTER_A0} -depth 16 -colorspace Gray -resize '320x240' +repage -depth 4 -define png:bit-depth=4  $@
	identify $@ ; display -title '%f' $@

${PREFIX}-GRY-A0-00255.png: ${MASTER_A0}
	convert ${MASTER_A0} -depth 16 -colorspace Gray -resize '320x240' +repage -depth 8 -define png:bit-depth=8  $@
	identify $@ ; display -title '%f' $@

${PREFIX}-GRY-A0-65535.png: ${MASTER_A0}
	convert ${MASTER_A0} -depth 16 -colorspace Gray -resize '320x240' +repage -depth 16 -define png:bit-depth=16  $@
	identify $@ ; display -title '%f' $@
        
# ----------------------------------------------------------------------
# GRAY colorspace, with alpha
        
${PREFIX}-GRY-A1-00255.png: ${MASTER_A1}
	convert ${MASTER_A1} -depth 16 -colorspace Gray -resize '320x240' +repage -depth 8 -define png:bit-depth=8  $@
	identify $@ ; display -title '%f' $@

${PREFIX}-GRY-A1-65535.png: ${MASTER_A1}
	convert ${MASTER_A1} -depth 16 -colorspace Gray -resize '320x240' +repage -depth 16 -define png:bit-depth=16  $@
	identify $@ ; display -title '%f' $@

# ----------------------------------------------------------------------
# RGB colorspace, no alpha

${PREFIX}-RGB-A0-00255.png: ${MASTER_A0}
	convert ${MASTER_A0} -depth 16 -resize '320x240' +repage -depth 8 -define png:bit-depth=8  $@
	identify $@ ; display -title '%f' $@

${PREFIX}-RGB-A0-65535.png: ${MASTER_A0}
	convert ${MASTER_A0} -depth 16 -resize '320x240' +repage -depth 16 -define png:bit-depth=16  $@
	identify $@ ; display -title '%f' $@
        
# ----------------------------------------------------------------------
# RGB colorspace, with alpha

${PREFIX}-RGB-A1-00255.png: ${MASTER_A1}
	convert ${MASTER_A1} -depth 16 -resize '320x240' +repage -depth 8 -define png:bit-depth=8  $@
	identify $@ ; display -title '%f' $@

${PREFIX}-RGB-A1-65535.png: ${MASTER_A1}
	convert ${MASTER_A1} -depth 16 -resize '320x240' +repage -depth 16 -define png:bit-depth=16  $@
	identify $@ ; display -title '%f' $@
        
# ----------------------------------------------------------------------
# RGB colorspace, colormapped, no alpha

${PREFIX}-MAP-A0-00001.png: ${MASTER_A0}
	convert ${MASTER_A0} -depth 16 -resize '320x240' +repage -depth 8 -colors 2 \
           -define png:bit-depth=1 -define png:color-type=3  $@
	identify $@ ; display -title '%f' $@

${PREFIX}-MAP-A0-00003.png: ${MASTER_A0}
	convert ${MASTER_A0} -depth 16 -resize '320x240' +repage -depth 8 -colors 4 \
           -define png:bit-depth=2 -define png:color-type=3  $@
	identify $@ ; display -title '%f' $@

${PREFIX}-MAP-A0-00015.png: ${MASTER_A0}
	convert ${MASTER_A0} -depth 16 -resize '320x240' +repage -depth 8 -colors 14 \
           -define png:bit-depth=4 -define png:color-type=3  $@
	identify $@ ; display -title '%f' $@

${PREFIX}-MAP-A0-00255.png: ${MASTER_A0}
	convert ${MASTER_A0} -depth 16 -resize '320x240' +repage -depth 8 -colors 250 \
          -define png:bit-depth=8 -define png:color-type=3  $@
	identify $@ ; display -title '%f' $@
        
# ----------------------------------------------------------------------
# RGB colorspace, colormapped, with alpha

${PREFIX}-MAP-A1-00001.png: ${MASTER_A1}
	convert ${MASTER_A1} -depth 16 -resize '320x240' +repage -depth 8 -colors 2 \
           -alpha Set -fuzz '50%' -transparent 'rgb(150,150,150)' \
           -define png:include-chunk=tRNS -define png:bit-depth=1 -define png:color-type=3  $@
	identify $@ ; display -title '%f' $@

${PREFIX}-MAP-A1-00003.png: ${MASTER_A1}
	convert ${MASTER_A1} -depth 16 -resize '320x240' +repage -depth 8 -colors 4 \
           -alpha Set -fuzz '50%' -transparent 'rgb(150,150,150)' \
           -define png:include-chunk=tRNS -define png:bit-depth=2 -define png:color-type=3  $@
	identify $@ ; display -title '%f' $@

${PREFIX}-MAP-A1-00015.png: ${MASTER_A1}
	convert ${MASTER_A1} -depth 16 -resize '320x240' +repage -depth 8 -colors 14  \
           -alpha Set -fuzz '50%' -transparent 'rgb(150,150,150)' \
           -define png:include-chunk=tRNS -define png:bit-depth=4 -define png:color-type=3  $@
	identify $@ ; display -title '%f' $@

${PREFIX}-MAP-A1-00255.png: ${MASTER_A1}
	convert ${MASTER_A1} -depth 16 -resize '320x240' +repage -depth 8 -colors 250 \
          -alpha Set -fuzz '50%' -transparent 'rgb(150,150,150)' \
          -define png:include-chunk=tRNS -define png:bit-depth=8 -define png:color-type=3  $@
	identify $@ ; display -title '%f' $@


