# Last edited on 2020-12-07 22:20:18 by jstolfi

repository JSLIBS
Jorge Stolfi's C libraries

This repository contains a set of C libraries that I have written over
almost 30 years. They were originally intended for "private" use, by
myself and maybe my graduate students. I am placing them on GitHub in
case someone else finds them useful.

--Jorge Stolfi, IC-UNICAMP 2020-12-07


libjs
  Miscellaneous hacks: text concatenation, random
  numbers, bool data type, assertions with custom messages,
  simple parsing tools, structured file parsing, 
  command-line parsing, etc..
  
libflt        
  Floating-point tools.
  
libia
  Interval arithmetic (IA) with guaranteed rounding.
  
libaa         
  Affine Arithmetic (AA) with guaranteed rounding.
  
libfgraph  
  Plotting functions using IA or AA.

libbbopt      
  Branch-and-Bound optimization using IA or AA.
  
libbezier     
  Multi-dimensional Bézier patches.
  
libdelaunay   
  Delaunay triangulation.
  
libdygrid     
  Multi-dimensional dyadic grids.
  
libgeo        
  Cartesian and Projective geometry in 2--6 dimensions.
  
libgmo 
  Dot-line-triangle geometric models.
  
libimg
  Image processing tools.
  
libintg
  Integration of ordinary differential equations.
  
  Miscellaneous tools for PBM/PGM/PPM images.
  
libminn
  Minimization of functions of N variables.
  
libminu
  Minimization of univariate functions.
  
libppv
  Representation and tools for multi-dimensional signals.
  
libps
  Creating Postscript drawings.
  
libpsplot
  [FUTURE] Creating Postscript drawings.
  
libquad
  The quad-edge data structure (for orientable manifolds only).
  
libsexp
  Lisp-like symbolic expressions.
  
libstreetmap
  Representation and plotting of street maps.

libapprox     

LIBRARIES

  Each library is in a directory of its own.

    * Splines:

      libbezier       
        Tools for N-dimensional tensor polynomial splines.

    * Miscellaneous utilities:

      libjs
        text concatenation, random numbers, bool data type, assertions
        with custom messages, simple parsing tools, structured file
        parsing, command-line parsing, etc..

    * Postscript plotting:

      libps
        Creating Postscript drawings.

    * Basic geometry:

      libgeo        
        Cartesian and Projective geometry in 2--6 dimensions.

    * Quad-edge and Voronoi/Delanay:

      libquad
        The quad-edge data structure (for orientable manifolds only).

      libdelaunay   
        Delaunay triangulation.

    * Validated numeric computation (interval and affine arithmetic):

      libflt        
        Floating-point tools, used by "libia" and "libaa".

      libia
        Interval arithmetic (IA) with guaranteed rounding.

      libaa         
        Affine Arithmetic (AA) with guaranteed rounding.

      libfgraph  
        Plotting functions using IA or AA.

      libzf
        Univariate zero finding using validated arithmetic (IA or AA).

      libbbopt      
        Branch-and-Bound optimization using IA or AA.



MAKEFILES

  Each library has its own Makefile, and a sub-directory "tests".
  containing one or more programs that check the library and/or
  provide examples of its use. These makefiles use the generic
  templates {GENERIC-LIB.make} and {GENERIC-LIB-TEST.make},
  respectively, which in turn uses the script {extract-ho-deps} to
  obtain the module dependency graph. which should have been made
  avaliable with in the distribution.

  The Makefiles assume an environment variable ${PLATFORM} that
  indicates the machine and operating system. This is to allow the
  same libraries to be compiled for different target platforms.
  
  The main "Makefile" targets are:

    clean
      Deletes all derived files.

    depend 
      Recreates the files "Deps.make" in the directory
      and all sub-directories. Must be executed manually
      whenever any "#include" directives have changed.
      
    install
      Rebuilds all libraries and installs the public 
      files in the global repositories.

    uninstall
      Deletes all public files from the global repositories.
      
    check 
      Runs the validation and example programs in the "tests"
      directories.

  Specifically, the makefiles will install the `public' header files
   in the global repositories
     
     ${STOLFIHOME}/include - header files (extensions ".h" and ".ho").
     ${STOLFIHOME}/lib - lisp library code.
     ${STOLFIHOME}/bin - shell scripts.
     ${STOLFIHOME}/lib/${PLATFORM} - library files (extension ".a").
     
  When a library imports files from another, it gets them from the
  global repository.
  
INSTALLATION

  Instructions for /bin/bash users:

    (*) I your home directory is shared among machines with different
        and incompatible OS/architecture combinations (e.g., Linux and
        Solaris, Intel and SPARC), set the variable ${PLATFORM} to a
        string that identifies the hardware and OS for which the
        libraries are to be compiled. For example
        
          export PLATFORM="Intel-Linux"
          
        If you have only one OS/architecture in 
        your network, you may set
        
          export PLATFORM="."
       
    (*) Set the variable ${STOLFIHOME} to the name of a directory
        where to place the libraries and associated files, for example

          export STOLFIHOME="${HOME}/stolfi-stuff"
          
    (*) Set the variable ${SRCDIR} to the name of the sub-directory 
        of ${STOLFIHOME} where the library sources will live:

          export SRCDIR="${STOLFIHOME}/programs/c"
          mkdir -p ${SRCDIR}
          
    (*) Move the tarfile to that directory:
    
          mv ~/Downloads/stolfi-JSLIBS-2009-02-10.tgz ${SRCDIR}/
    
    (*) Unpack the tar file with, for example

          cd ${SRCDIR}
          tar -xvzf stolfi-JSLIBS-2009-02-10.tgz

    (*) Compile the library:
    
          make build
          
    (*) Install the library
    
          make install
          
          
