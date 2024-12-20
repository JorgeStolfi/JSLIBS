#! /usr/bin/python3
# Last edited on 2024-12-16 09:40:13 by stolfi

import os, sys, subprocess
from subprocess import run
  
err = sys.stderr

def main():

  run_prog = False

  imageSize = 2
  
  pix_smp = 2
  dir_smp = 10
  imgType = 'noise01'
  
  NX = imageSize*160
  NY = imageSize*120

  pixDeb = choose_pixDeb(NX, NY)
  stackDir = stack_dir_name(NX, NY, pix_smp, dir_smp, imgType)
  if run_prog:
    clean_stack(stackDir)
    make_stack(NX, NY, pix_smp, dir_smp, imgType, stackDir, pixDeb)
  show_results(stackDir,pixDeb)
  os.system(f"( cd out/ && rmempty.sh )");
  return 0
  
def stack_dir_name(NX, NY, pix_smp, dir_smp, imgType):
  # Returns the name of the top-level directory for all output files
  # of this run. Includes the "out/" prefix.
  
  sizeTag = f"{NX:04d}x{NY:04d}"
  samplingTag = f"hs{pix_smp:02d}-kr{dir_smp:02d}"
  stackDir = f"out/img-{imgType}-{sizeTag}-{samplingTag}"
  return stackDir
  # ------------------------------------------------------------
  
def clean_stack(stackDir):
  # Returns the name of the top-level directory for all output files
  # of this run. Includes the "out/" prefix.
  
  os.system(f"rm -rf '{stackDir}'")
  os.system(f"mkdir -pv {stackDir}")
  return stackDir
  # ------------------------------------------------------------

def choose_pixDeb(NX, NY):
  # Returns a list {(ix,iy)} with the indices of the pixel 
  # for which the ray-tracing, sampling and debugging data is to 
  # be writen.
  
  pixDeb = ( NX//2, NY//2 )
  return pixDeb
  # ------------------------------------------------------------

def make_stack(NX, NY, pix_smp, dir_smp, imgType, stackDir, pixDeb):
  # Runs the C program generating the rendered frame images  and 
  # sampling information files. Also compiles the program.
  
  prog = 'test_mfok_raytrace'
  
  err.write("compiling program ...\n")
  status = os.system(f"make {prog}");
  assert status == 0, f"** 'make {prog}' returned status = {status}"

  if not os.path.exists(prog):
    err.write(f"** 'make {prog}' failed\n")
  else:
    err.write(f"running {prog} ...\n")
    status = run(
        [ f"{prog}",
          "-imageSize",   f"{NX}", f"{NY}",
          "-imageType",   f"{imgType}",
          "-pixSampling", f"{pix_smp}",
          "-dirSampling", f"{dir_smp}",
          "-debugPixel",  f"{pixDeb[0]}", f"{pixDeb[1]}",
          "-stackDir",    f"{stackDir}"
        ]
      )
    assert status.returncode == 0, f"** '{prog}' failed - returned status = {status}"
  return
  # ------------------------------------------------------------
 
def show_results(stackDir, pixDeb):  
  # Runs 'mfok_plot_rays.sh {file}' script on the pixel sampling information {file}
  # of each frames.  Then runs 'display' with all rendered frame images.
  #
  # Assumes that each frame image file is is "{stackDir}/{frameName}"
  # where {frameName} is each string returned by {list_frames}.
  #
  # Assumes that the pixel information files are called
  # "{stackDir}/{frameName}/pixel-rays-{XXXX}-{YYYY}.txt"
  # where {XXXX} and {YYYY} are the pixel indices in {pixDeb}
  # formatted as "%04d".
  
  err.write("collecting the images and plotting ray data ...\n")
  images = []
  for frameName in list_frames(stackDir):
    err.write(f"  frame {frameName} ...\n");
    frame_dir = f"{stackDir}/{frameName}"
    plot_rays(stackDir, frameName, pixDeb);

    imageName = "img-blur"
    images.append(f"{frameName}/{imageName}.png")

  # Display the images>
  args = " ".join(images)
  err.write(f"displaying the images from {stackDir} ...\n")
  err.write("args = «" + args + "»\n")
  status = os.system(f"( cd {stackDir} && display -title '%d' -filter Box -resize '300%' {args} )")
  
  plot_pixels(stackDir, pixDeb)
  return
  # ------------------------------------------------------------
    
def plot_pixels(stackDir, pixDeb):
  # Plots the computed pixel properties of pixel {pixDeb}
  # across all frames in directory {stackDir},
  # using the script {plotter} below.
 
  plotter = 'mfok_plot_pixels.sh'
  
  err.write("plotting pixel data ...\n")
  pixelFile = f"pixel-data-{pixDeb[0]:04d}-{pixDeb[1]:04d}.txt"
  if not os.path.isfile(f"{stackDir}/{pixelFile}"):
    err.write(f"** file {pixelFile} not found\n")
  else:
    status = os.system(f"( cd {stackDir} && ../../{plotter} {pixelFile} )")
    assert status == 0, f"** '{plotter} ...' returned status = {status}"

def plot_rays(stackDir, frameName, pixDeb):   
  # Display the sampling points and ray hits fro pixel {pixDeb}
  # of the frame {frameName} in directory {stackDir},
  # using the script {plotter} below.
 
  plotter = 'mfok_plot_rays.sh'
  
  rayName = f"pixel-rays-{pixDeb[0]:04d}-{pixDeb[1]:04d}.txt"
  rayFile = f"{frameName}/{rayName}"
  if not os.path.isfile(f"{stackDir}/{rayFile}"):
    err.write(f"** file {rayFile} not found\n")
  else:
    status = os.system(f"( cd {stackDir} && ../../{plotter} {rayFile} )")
    assert status == 0, f"** '{plotter} ...' returned status = {status}"
  
def list_frames(stackDir):
  # Returns a list of the names of all frame folders in directory
  # {stackDir}.  The folder names are assumed to begin with "zf".
  
  LS = os.listdir(stackDir); LS.sort()
  FS = [ x for x in LS if x[0:2] == "zf" ]
  return FS
  
main()
