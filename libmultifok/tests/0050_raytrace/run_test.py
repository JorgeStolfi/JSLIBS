#! /usr/bin/python3
# Last edited on 2025-04-19 08:14:20 by stolfi

import os, sys, subprocess
  
err = sys.stderr

def main():

  runProg = True

  # Output image size:
  imageRelSize = 2
  NX = imageRelSize*160
  NY = imageRelSize*120
  
  pixSmp = 2
  dirSmp = 10
  imageDesc = ('melon16', 'noise01', 'spots06', 512, 512)

  debugPixel = choose_debug_pixel(NX, NY)
  stackFolder = stack_dir_name(NX, NY, pixSmp, dirSmp, inputDesc[0])
  bash(f"mkdir -p {stackFolder}");
  if runProg:
    clean_stack(stackFolder)
    make_stack(NX, NY, pixSmp, dirSmp, inputDesc, stackFolder, debugPixel)
  show_results(stackFolder,debugPixel)
  bash(f"( cd out/ && rmempty.sh )");
  return 0
  
def stack_dir_name(NX, NY, pixSmp, dirSmp, inputImageName):
  # Returns the name of the top-level directory for all output files
  # of this run. Includes the "out/" prefix.
  
  sizeTag = f"{NX:04d}x{NY:04d}"
  samplingTag = f"hs{pixSmp:02d}-kr{dirSmp:02d}"
  stackFolder = f"out/img-{inputImageName}-{sizeTag}-{samplingTag}"
  return stackFolder
  # ------------------------------------------------------------
  
def clean_stack(stackFolder):
  # Removes the folder {stackFolder} and recreats it empty.
  
  bash(f"rm -rf '{stackFolder}'")
  bash(f"mkdir -pv {stackFolder}")
  return stackFolder
  # ------------------------------------------------------------

def choose_debug_pixel(NX, NY):
  # Returns a list {(ix,iy)} with the indices of the pixel 
  # for which the ray-tracing, sampling and debugging data is to 
  # be writen.
  
  debugPixel = ( NX//2, NY//2 )
  return debugPixel
  # ------------------------------------------------------------

def make_stack(NX, NY, pixSmp, dirSmp, inputDesc, stackFolder, debugPixel):
  # Runs the C program generating the rendered frame images  and 
  # sampling information files. Also compiles the program.
  #
  # The {inputDesc} should be a list {(inputImageName, imageBase1, imageBase2)}
  # where "{imageBase1}.png" and "{imageBase2}.png" are two image files
  # in "../00-DATA/textures" that must be multiplied together to produce 
  # the background image to be raytraced, which will be stored into
  # "in/{inputImageName}.png"
  
  prog = 'test_mfok_raytrace'
  
  err.write("compiling program ...\n")
  bash(f"make {prog}");
  
  inputImageName = make_input_image(inputDesc);
  
  if not os.path.exists(prog):
    err.write(f"** 'make {prog}' failed\n")
  else:
    err.write(f"running {prog} ...\n")
    run_process(
        [ f"{prog}",
          "-imageSize",   f"{NX}", f"{NY}",
          "-inputImage",  f"{inputImageName}",
          "-pixSampling", f"{pixSmp}",
          "-dirSampling", f"{dirSmp}",
          "-debugPixel",  f"{debugPixel[0]}", f"{debugPixel[1]}",
          "-stackFolder", f"{stackFolder}"
        ]
      )
  return
  # ------------------------------------------------------------
 
def make_input_image(inputDesc):
  # Creates the input image file as described by {inputDesc}.
  # Returns the input image file name, without the folder ("in/") or
  # extension (".png").
  resDir = "in"
  resName = inputDesc[0]
  baseDir = "../00-DATA/texture-bank/pgm-512x512"
  base1Name = inputDesc[1]
  base2Name = inputDesc[2]
  resNX = inputDesc[3]
  resNY = inputDesc[4]

  # The input image "{resDir}/{resName}.png" is obtained by multiplying the 
  # textures "{baseDir}/{base1Name}.png" and "{baseDir}/{base2Name}.png".
  resFile = f"{resDir}/{resName}.png"
  base1File = f"{baseDir}/{base1Name}.png"
  base2File = f"{baseDir}/{base2Name}.png"
  bash(f"${baseDir}/multiply_gray_textures.sh {base1File} {base2File} {resNX} {resNY} {resFile}"
  assert file_OK(resFile), "{resFile} not created"
  return resName
  # ----------------------------------------------------------------------

def show_results(stackFolder, debugPixel):  
  # Runs 'mfok_plot_rays.sh {file}' script on the pixel sampling information {file}
  # of each frames.  Then runs 'display' with all rendered frame images.
  #
  # Assumes that each frame image file is is "{stackFolder}/{frameName}"
  # where {frameName} is each string returned by {list_frames}.
  #
  # Assumes that the pixel information files are called
  # "{stackFolder}/{frameName}/pixel-rays-{XXXX}-{YYYY}.txt"
  # where {XXXX} and {YYYY} are the pixel indices in {debugPixel}
  # formatted as "%04d".
  
  err.write("collecting the images and plotting ray data ...\n")
  images = []
  for frameName in list_frames(stackFolder):
    err.write(f"  frame {frameName} ...\n");
    frame_dir = f"{stackFolder}/{frameName}"
    plot_rays(stackFolder, frameName, debugPixel);

    imageName = "img-blur"
    images.append(f"{frameName}/{imageName}.png")

  # Display the images>
  args = " ".join(images)
  err.write(f"displaying the images from {stackFolder} ...\n")
  err.write("args = «" + args + "»\n")
  bash(f"( cd {stackFolder} && display -title '%d' -filter Box -resize '300%' {args} )")
  
  plot_pixels(stackFolder, debugPixel)
  return
  # ------------------------------------------------------------
    
def plot_pixels(stackFolder, debugPixel):
  # Plots the computed pixel properties of pixel {debugPixel}
  # across all frames in directory {stackFolder},
  # using the script {plotter} below.
 
  plotter = 'mfok_plot_pixels.sh'
  
  err.write("plotting pixel data ...\n")
  pixelFile = f"pixel-data-{debugPixel[0]:04d}-{debugPixel[1]:04d}.txt"
  if not os.path.isfile(f"{stackFolder}/{pixelFile}"):
    err.write(f"** file {pixelFile} not found\n")
  else:
    bash(f"( cd {stackFolder} && ../../{plotter} {pixelFile} )")

def plot_rays(stackFolder, frameName, debugPixel):   
  # Display the sampling points and ray hits fro pixel {debugPixel}
  # of the frame {frameName} in directory {stackFolder},
  # using the script {plotter} below.
 
  plotter = 'mfok_plot_rays.sh'
  
  rayName = f"pixel-rays-{debugPixel[0]:04d}-{debugPixel[1]:04d}.txt"
  rayFile = f"{frameName}/{rayName}"
  if not os.path.isfile(f"{stackFolder}/{rayFile}"):
    err.write(f"** file {rayFile} not found\n")
  else:
    bash(f"( cd {stackFolder} && ../../{plotter} {rayFile} )")
  
def list_frames(stackFolder):
  # Returns a list of the names of all frame folders in directory
  # {stackFolder}.  The folder names are assumed to begin with "zf".
  
  LS = os.listdir(stackFolder); LS.sort()
  FS = [ x for x in LS if x[0:2] == "zf" ]
  return FS

def file_OK(path):
  return path != None and os.path.exists(path) and os.path.getsize(path) > 0
  # ......................................................................

def Er(msg):
  sys.stderr.write(msg)
  # ......................................................................

def run_command(command):
  result = subprocess.run(command, text = True)
  if result.returncode != 0:
    print(result.stderr)
    print(result.stdout)
    assert False, f"** {command[0]} failed - returned status = {result}"
  return
  # ......................................................................
  
def bash(cmd):
  # Execute the string {cmd} with "/bin/bash".
  result = subprocess.run([ cmd ], shell = True, executable = "/bin/bash")
  if result.returncode != 0:
    print(result.stderr)
    print(result.stdout)
    assert False, f"** {cmd} failed - returned status = {result}"
  # ......................................................................

main()
