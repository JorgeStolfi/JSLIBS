#! /usr/bin/python3
# Last edited on 2025-02-04 18:07:14 by stolfi

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
  stackFolder = stack_dir_name(NX, NY, pix_smp, dir_smp, imgType)
  if run_prog:
    clean_stack(stackFolder)
    make_stack(NX, NY, pix_smp, dir_smp, imgType, stackFolder, pixDeb)
  show_results(stackFolder,pixDeb)
  os.system(f"( cd out/ && rmempty.sh )");
  return 0
  
def stack_dir_name(NX, NY, pix_smp, dir_smp, imgType):
  # Returns the name of the top-level directory for all output files
  # of this run. Includes the "out/" prefix.
  
  sizeTag = f"{NX:04d}x{NY:04d}"
  samplingTag = f"hs{pix_smp:02d}-kr{dir_smp:02d}"
  stackFolder = f"out/img-{imgType}-{sizeTag}-{samplingTag}"
  return stackFolder
  # ------------------------------------------------------------
  
def clean_stack(stackFolder):
  # Removes the folder {stackFolder} and recreats it empty.
  
  os.system(f"rm -rf '{stackFolder}'")
  os.system(f"mkdir -pv {stackFolder}")
  return stackFolder
  # ------------------------------------------------------------

def choose_pixDeb(NX, NY):
  # Returns a list {(ix,iy)} with the indices of the pixel 
  # for which the ray-tracing, sampling and debugging data is to 
  # be writen.
  
  pixDeb = ( NX//2, NY//2 )
  return pixDeb
  # ------------------------------------------------------------

def make_stack(NX, NY, pix_smp, dir_smp, imgType, stackFolder, pixDeb):
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
          "-stackFolder", f"{stackFolder}"
        ]
      )
    assert status.returncode == 0, f"** '{prog}' failed - returned status = {status}"
  return
  # ------------------------------------------------------------
 
def show_results(stackFolder, pixDeb):  
  # Runs 'mfok_plot_rays.sh {file}' script on the pixel sampling information {file}
  # of each frames.  Then runs 'display' with all rendered frame images.
  #
  # Assumes that each frame image file is is "{stackFolder}/{frameName}"
  # where {frameName} is each string returned by {list_frames}.
  #
  # Assumes that the pixel information files are called
  # "{stackFolder}/{frameName}/pixel-rays-{XXXX}-{YYYY}.txt"
  # where {XXXX} and {YYYY} are the pixel indices in {pixDeb}
  # formatted as "%04d".
  
  err.write("collecting the images and plotting ray data ...\n")
  images = []
  for frameName in list_frames(stackFolder):
    err.write(f"  frame {frameName} ...\n");
    frame_dir = f"{stackFolder}/{frameName}"
    plot_rays(stackFolder, frameName, pixDeb);

    imageName = "img-blur"
    images.append(f"{frameName}/{imageName}.png")

  # Display the images>
  args = " ".join(images)
  err.write(f"displaying the images from {stackFolder} ...\n")
  err.write("args = «" + args + "»\n")
  status = os.system(f"( cd {stackFolder} && display -title '%d' -filter Box -resize '300%' {args} )")
  
  plot_pixels(stackFolder, pixDeb)
  return
  # ------------------------------------------------------------
    
def plot_pixels(stackFolder, pixDeb):
  # Plots the computed pixel properties of pixel {pixDeb}
  # across all frames in directory {stackFolder},
  # using the script {plotter} below.
 
  plotter = 'mfok_plot_pixels.sh'
  
  err.write("plotting pixel data ...\n")
  pixelFile = f"pixel-data-{pixDeb[0]:04d}-{pixDeb[1]:04d}.txt"
  if not os.path.isfile(f"{stackFolder}/{pixelFile}"):
    err.write(f"** file {pixelFile} not found\n")
  else:
    status = os.system(f"( cd {stackFolder} && ../../{plotter} {pixelFile} )")
    assert status == 0, f"** '{plotter} ...' returned status = {status}"

def plot_rays(stackFolder, frameName, pixDeb):   
  # Display the sampling points and ray hits fro pixel {pixDeb}
  # of the frame {frameName} in directory {stackFolder},
  # using the script {plotter} below.
 
  plotter = 'mfok_plot_rays.sh'
  
  rayName = f"pixel-rays-{pixDeb[0]:04d}-{pixDeb[1]:04d}.txt"
  rayFile = f"{frameName}/{rayName}"
  if not os.path.isfile(f"{stackFolder}/{rayFile}"):
    err.write(f"** file {rayFile} not found\n")
  else:
    status = os.system(f"( cd {stackFolder} && ../../{plotter} {rayFile} )")
    assert status == 0, f"** '{plotter} ...' returned status = {status}"
  
def list_frames(stackFolder):
  # Returns a list of the names of all frame folders in directory
  # {stackFolder}.  The folder names are assumed to begin with "zf".
  
  LS = os.listdir(stackFolder); LS.sort()
  FS = [ x for x in LS if x[0:2] == "zf" ]
  return FS
  
main()
