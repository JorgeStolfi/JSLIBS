# JSLIBS - Jorge Stolfi's C libraries

% Last edited on 2021-04-12 20:12:23 by jstolfi

This repository contains a set of C libraries that I have written over
almost 30 years. They were originally intended for "private" use, by
myself and maybe my graduate students. I am placing them on GitHub in
case someone else finds them useful.

  --Jorge Stolfi, IC-UNICAMP 2020-12-07
  
Creator: Jorge Stolfi
Supervisor: Jorge Stolfi
Intended users: Jorge Stolfi, general public

##LIBRARIES

Each library is in a package of its own.

### Data structures:

  * libspmat
    Functions for generic sparse matrix definition and manipulation,
    with user-specified element types.

  * libdgraph
    Directed graphs as sparse matrices of booleans,
    based on {libspmat}.

  * libjsarray

    N-dimensional arrays whose elements are floats, intervals,
    or other user-specfied scalar type.  See also {vec.h} in {libjs}.

  * libsexp
    Functions for reading and printing Lisp-like
    S-expressions.

### Function approximation and optimization:

  * liblsq 
    Least squares fitting, also with outlier rejection.

  * libapprox
    Object-oriented functions for least-squares functional approximation.
    Developed with Anamaria Gomide for her Ph. D. Thesis.
    See also {liblsq}.
    
  * libpspulse
    Polynomial spline pulses and tents.

  * libminn
    General non-linear n-variate minimization routines.

  * libminu
    Minimization of univariate non-linear functions.

  * libbbopt
    Tools for finding maximum or minimum of functions
    by the bounding box method, using  IA or AA estimators.
  
  * libbezier
    Bézier representation of univariate polynomials.

  * libclassif
    Procedures for vector clustering and classification,
    using various methods including the Falcão's Optimal Path Forest.
    Mainly to compare the accuracy etc of those methods.
    
### Numerical calculus:

  * libintg
    J. Stolfi's implementation of the Euler and Runge-Kutta
    methods for integration of 

### Computational geometry and graphics:

  * libgeo
    Basic geometry routines, including linear algebra in {R^2},{R^3}, {R^4},
    {R^6} and {R^n}, matrix operations, Gaussian elimination.  
    Also some tools for homogeneous coordinates in {T^2}, {T^3} and {T^n},
    for integer tuples in {Z^2}, {Z^3}, and {Z^n}, for ellipse representations,
    for quadratic function optimization by the edge-divided simplex method.
    Also tools for parsing vector arguments from files or the command line.
    See also {libjsarray}.
  
  * libdelaunay
    Planar Voronoi/Delaunay construction with quad-edge data structure
    
  * libeps
    Creation of Encapsulated Postscript files,
    with emphasis on graphics as opposed to text.

  * libpsextra
    Functions to plot Bézier maps of {R^2}.
    >> TO BE DECOMISSIONED
    
  * libquad
    Representation and tools of the quad-edge and oct-edge
    data structures. Uses low-order bits of address to indicate flip and rot.

  * libgem
    Implementation of the GEM data structure.
    Partly developed by Arnaldo Montagner and Lucas Moutinho Bueno
    for their thesis projects.

  * libgmap
    Representation and manipulation of n-dimensional maps
    (complexes of topological polytopes), represented internally by the
    gems of their barycentric subdivisions.
    
  * libstmesh
    A data structure for riangle meshes with quantized vertices
    and semi-topological structure.  Related to the UTFPR 3D slicing project 
    by Minetto, Volpatto, Habib, and Stolfi (2015).
    
  * libstpoly
    Representation of polygonal figures by unstructured lists 
    of segments (the "STP" format, a 2D analog of the popular STL format)

  * libmkgr
    Functions to create grids of marks, e.g. calibration grids,
    quadrille paper, etc.

### Images, videos, sounds, tomograms:

  * libimg
    Representation and processing of multichannel two-dimensional images
    as arrays of {float} values. Image operations include interpolation,
    deformation, and Fourier-Hartley transform, Also conversion of such
    images to/from various image file formats (PPM, PNG, JPG). Also
    functions to represent and manipulate colors as vectors of 3
    {float}s in various color spaces.
    
  * libift
    Functions to build the Image Forest Transform (IFT) of an image.

  * libcamfirewire
    Controlling camera through the Firewire interface.
    Developed mainly by Rafael Saracchini for his PhD project.

  * libjsaudio
    J. Stolfi's tools for handling audio files.

  * libppv
    Portable multi-dimensional arrays of integer samples
    for general signal  processing (images, videos, sound, spectra, etc.).a

  * libvoxm
    Tools for voxel-based 3D solid modeling.

### Computer vision and pattern matching:

  * libmsmatch
    Routines for multiscale sequence matching.
    Derived from Helena Leitão's Ph. D. thesis project
    with contributions by Rafael Saracchini.

  * libmultifok
    Functions to do multifocus stereo.
    >> IN DEVELOPMENT
    
  * libpst
    Functions for photometric stereo.
    Developed with Rafael Saracchiní for his Ph. D. thesis project.
    
  * libtsai
    Tsai's camera calibration routines.
    Used in Rodrigo Minetto's PhD project.

### Bioinformatics:

  * libdnaenc
    Encoding of DNA strings as sequences of points in 3D.

  * libdnaview
    Visualizing DNA strings as curves in 3D.
 
### Interval and affine arithmetic:

  * libia
    Standard interval arithmetic (IA) operations. They provide
    guaranteed error bounds for the results, accounting for roundoff
    errors.  Also axes-aligned boxes (multidimensional intervals). 
    >> {ia_box.h} SUPERSEDED BY {interval.h}

  * libaa
    Basic Affine Arithmetic (AA) operations.  They provide
    guaranteed error bounds for the results, accounting for roundoff
    errors, that are usually tighter than those of interval arithmetic (IA).
    This implementation requires manual allocation of AA forms..

  * libflt
    Operations on floating-point values that parallel those of {libia}
    and {libaa}.   Also parsing of algebraic formulas
    into a stack-machine pseudocode.

  * libaacmp
    Compiling formulas into affine arithmetic function calls.
    
  * libaafuncs 
    Some 1-argument and 2-argument test functions for the
    interval arithmetic library {libia} and the affine arithmetic
    library {libaa}, as well as the "isomorhic" floating point
    arithmetic library {libflt}.
    
  * libaagraph
    Procedures to plot univariate functions using libraries for
    interval arithmetic {libia}, affine arithmetic {libaa}.
    and "isomorphic" floating point {libflt}
    
### Natural language processing:

  * libdicio
    Representing word lists as acyclic automata,
    as in the Dicio project (lucchesi/stolfi).
    This is a C rewriting of the original Modula-3 libraries.

### Drawing buidings, furniture, maps:

  * libarchdraw
    Simple (but not at all user-friendly) 
    tools to generate architectural ground plan sketches. 
    
  * libsheetcut
    Functions for planning the cutting of rectangular
    plates out of rectangular stock sheets.
    
  * libstreetmap
    Functions for reading, plotting, anda manipulating street maps.
    
### Miscellaneous

  * libjs
    Miscellaneous hacks: text concatenation, random
    numbers, bool data type, sign data type, generic 
    self-sized extensible vector data types,
    integer ranges, assertions with custom messages,
    simple file parsing tools, structured file parsing, 
    command-line parsing, indexing in multimensional arrays,
    heaps of integers, merging and sorting arrays of integers, 
    subsets of {0..31}, linear search and interpolation
    in tables, program timing, etc..

  * libjsextra
    Various interfaces which should 
    be in {jslib.h} or the like, but which aren't
    working yet.
    >> UNDER DEVELOPMENT

  * libbtc
    Data analysis related to bitcoin, such as the
    sum-of-bubbles model and log-quadratic fitting.

  * libcryptoy
    Trivial criptography tools.

### NeuroMat project:

  * libneuro
    EEG procesinng tools for the NeuroMat project.

  * libnmsim
    Data structures and basic functions
    for modeling of neurons and neuronal networks.

  * libnmsim_e
    Tools to simulate a neuronal net model as described in {libnmsim},
    at the {elem} level (individual neurons and synapses).

## MAKEFILES

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
  
## INSTALLATION

  Instructions for /bin/bash users:

    (*) I your home directory is shared among machines with different
        and incompatible OS/architecture combinations (e.g., Linux and
        Solaris, Intel and SPARC), set the variable ${PLATFORM} to a
        string that identifies the hardware and OS for which the
        * libraries are to be compiled. For example
        
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
          
          
