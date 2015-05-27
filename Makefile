#================================================================================#
# This makefile copied from B. Jonkman by J. Michalakes on 29-Jan-2013,          #
# adapted from Crunch (M. Buhl on 25-Jan-2013).                                  #
#                                                                                #
# This makefile has been tested on Windows 7 with gfortran.                      #
# This makefile works with mingw32-make.exe.                                     #
#                                                                                #
# It was designed to be used with:                                               #
#     BeamDyn        (v1.00.04,      ???????????)                       #
#     NWTC Subroutine Library (CSC version    29-Jan-2013)                       #
#                                                                                #
# Older versions of ModuleName and NWTC Library may not work with this makefile. #
#================================================================================#

################################################
# Added by Eliot for BeamDyn interface
BEAMDYN_LIB = libBeamDyn.${LIBEXT}
#LIBEXT = dylib
#LAPACK_DIR = /usr/local/Cellar/lapack/3.5.0/lib
LIBEXT = so
LAPACK_DIR = ${HOME}/lapack-3.5.0
CPP_DRIVER = cppbeam
LIB_INST_DIR = ${FOAM_USER_LIBBIN}
INC_INST_DIR = $(WM_PROJECT_USER_DIR)/platforms/$(WM_OPTIONS)/include
################################################

   # 32-bit or 64-bit?
#BITS = 32
BITS = 64


   # Location of source files for ModuleName and the NWTC Library.  You will probably need to change these for your system.
OS = Linux

# these settings are overridden below
ifeq ($(OS),Windows_NT)
   LIB_DIR  = C:/Users/bjonkman/Documents/DATA/DesignCodes/miscellaneous/nwtc_subs/SVNdirectory/trunk/source
   REGISTRY = Registry
else
   LIB_DIR    = ${HOME}/NWTC_Library/trunk/source
   REGISTRY   = ${HOME}/NWTC_Library/FAST_Registry/registry.exe
   NETLIB_DIR = ${LIB_DIR}/dependencies/lapack
   LAPACK_LINK  = -L${LAPACK_DIR} -llapack -lblas
endif

MODNAME_DIR = .

   # Name of compiler to use and flags to use.

FC     = gfortran
CC	   = gcc
OPT    = 
#FFLAGS = -O2 -m$(BITS) -fbacktrace -ffree-line-length-none -x f95-cpp-input -ffree-line-length-none
#FFLAGS = -O2 -g -fbacktrace -ffree-line-length-none -x f95-cpp-input -fdefault-real-8 -fdefault-double-8
#LDFLAGS = -O2 -g -fbacktrace
#FFLAGS = -g -O2 -m32 -fbacktrace -ffree-line-length-none -x f95-cpp-input -fdefault-real-8 -fdefault-double-8 #-fcheck=bounds

# from Qi:
#FFLAGS = $(OPT) -O2 -m32 -fdefault-real-8 -fdefault-double-8 -fbacktrace -ffree-line-length-none -x f95-cpp-input -Wsurprising 
#LDFLAGS = $(OPT) -O2 -m32 -fbacktrace

#FFLAGS = $(OPT) -O2 -fdefault-real-8 -fdefault-double-8 -fbacktrace -ffree-line-length-none -x f95-cpp-input -Wsurprising -DDOUBLE_PRECISION
FFLAGS = $(OPT) -O2 -fdefault-real-8 -fdefault-double-8 -fbacktrace -ffree-line-length-none -x f95-cpp-input -Wsurprising -DDOUBLE_PRECISION -DFPE_TRAP_ENABLED
LDFLAGS = $(OPT) -O2 -fbacktrace

CFLAGS = $(OPT) -O2

   # Precision.

# Use "SingPrec" for single precision and "DoubPrec" for double precision.  You may also need to change an option switch to make constants DP.
PREC = SingPrec
#PREC = DoubPrec # set DOUBLE_PRECISION flag instead

   #==========================================================#
   # You should not need to change anything beyond this point #
   #==========================================================#

   # System-specific settings.

ifeq ($(OS),Windows_NT)
      # Windows
   DEL_CMD   = del
   EXE_EXT   = _gwin$(BITS).exe
   INTER_DIR = Obj_win$(BITS)
   MD_CMD    = @mkdir
   OBJ_EXT   = .obj
   PATH_SEP  = \\
   SYS_FILE  = SysGnuWin
else
      # Linux
   DEL_CMD   = rm -f
   EXE_EXT   = _glin$(BITS)
   INTER_DIR = Obj_lin$(BITS)
   MD_CMD    = @mkdir -p
   OBJ_EXT   = .o
   PATH_SEP  = /
   SYS_FILE  = SysGnuLinux
endif

   # Destination and RootName for executable

OUTPUT_NAME = BeamDyn_AM2
DEST_DIR    = .

   # Library files.

LIB_SOURCES =             \
	$(PREC).f90        \
        NWTC_Base.f90      \
	$(SYS_FILE).f90    \
	NWTC_Library_Types.f90   \
	NWTC_IO.f90        \
        NWTC_Num.f90       \
	ModMesh_Types.f90  \
	ModMesh.f90        \
        ModMesh_Mapping.f90  \
	NWTC_Library.f90   \

NETLIB_SOURCES=             \
	NWTC_LAPACK.f90 

MODNAME_SOURCES   =             \
	BeamDyn_SP.f90             \
	BeamDyn_Types.f90 
	#Driver_BeamDyn_AM2.f90                \

#FITPACK_SOURCES = \
#	fpback.f \
#	fpbspl.f \
#	fpchec.f \
#	fpcurf.f \
#	fpdisc.f \
#	fpgivs.f \
#	fpknot.f \
#	fprati.f \
#	fprota.f \
#	curfit.f

vpath %.f90 $(LIB_DIR) $(NETLIB_DIR) $(MODNAME_DIR)
vpath %.mod $(INTER_DIR)
vpath %.obj $(INTER_DIR)

ALL_SOURCES = $(LIB_SOURCES) $(NETLIB_SOURCES) $(MODNAME_SOURCES)
tmp_objs1   = $(ALL_SOURCES:.f90=.obj)
tmp_objs2   = $(tmp_objs1:.F90=.obj)       #note the upper case here (from IceFloe)
ALL_OBJS    = $(tmp_objs2:.f=.obj)

FITPACK_OBJS = $(FITPACK_SOURCES:.f=.obj)

#NETLIB_OBJS    = $(NETLIB_SOURCES:.f90=.obj)
#LIB_OBJS          = $(LIB_SOURCES:.f90=.obj)
#MODNAME_OBJS  = $(MODNAME_SOURCES:.f90=.obj)

   # Rule to do everything.

all:     default
default: BeamDyn_Types.f90 $(INTER_DIR) $(DEST_DIR)/$(OUTPUT_NAME)$(EXE_EXT)

   # General rule for making the files.

# -B is needed for MinGW version of Gfortran

#%.obj: %.f
#	$(FC) -I $(INTER_DIR) -g -c $< -o $(INTER_DIR)/$@ -J $(INTER_DIR)

#spline.obj: spline.f90
#	$(FC) -I $(INTER_DIR) $(FFLAGS) -g -c -fPIC $< -o $(INTER_DIR)/$@ -J $(INTER_DIR)
#%.obj: %.f90
#	$(FC) -I $(INTER_DIR) $(FFLAGS) -g -c $< -o $(INTER_DIR)/$@ -J $(INTER_DIR)

%.obj: %.f90
	$(FC) -I $(INTER_DIR) $(FFLAGS) -g -c -fPIC $< -o $(INTER_DIR)/$@ -J $(INTER_DIR)

%.obj: %.C
	$(CC) -I $(INTER_DIR) $(CFLAGS) -g -c $< -o $(INTER_DIR)/$@ -J $(INTER_DIR)

   #  Dependency rules.

#NWTC Library dependency rules:
NWTC_Base.obj:              $(PREC).obj
$(SYS_FILE).obj:            NWTC_Base.obj
NWTC_Library_Types.obj:     $(SYS_FILE).obj
NWTC_IO.obj:                NWTC_Library_Types.obj
NWTC_Num.obj:               NWTC_IO.obj
ModMesh_Types.obj:          NWTC_Num.obj
ModMesh.obj:                ModMesh_Types.obj
ModMesh_Mapping.obj:        ModMesh.obj NWTC_LAPACK.obj
NWTC_Library.obj:           ModMesh.obj  ModMesh_Mapping.obj

NWTC_LAPACK.obj:            NWTC_Base.obj

#BeamDyn_Types.f90:       Registry_BeamDyn.txt
BeamDyn_Types.obj:       NWTC_Library.obj
BeamDyn_SP.obj: BeamDyn_Types.obj DynamicSolution_AM2.f90 GenerateDynamicElement_AM2.f90 *.f90

Driver_BeamDyn_AM2.obj: BeamDyn_SP.obj NWTC_Library.obj 

#$(OUTPUT_NAME)$(EXE_EXT): Driver_Dynamic_AM2.obj

   # Make sure the destination directory for the intermediate files exist.

$(INTER_DIR):
	$(MD_CMD) $(INTER_DIR)


   # Run the registry if the input file changes.

BeamDyn_Types.f90:
	#$(REGISTRY) Registry_BeamDyn.txt
	@echo "Skipping registry"
	cp -v $@_orig $@


   # For compiling the driver/glue code.

$(DEST_DIR)/$(OUTPUT_NAME)$(EXE_EXT): Driver_BeamDyn_AM2.obj $(ALL_OBJS) | $(INTER_DIR)
	$(FC) $(LDFLAGS) -I $(INTER_DIR) -o $(DEST_DIR)/$(OUTPUT_NAME)$(EXE_EXT) \
	$(INTER_DIR)/Driver_BeamDyn_AM2.obj $(foreach src, $(ALL_OBJS), $(addprefix $(INTER_DIR)/,$(src))) $(MAP_lib) $(LAPACK_LINK)


   # For interface with OpenFOAM (EWQ 04/01/2015)

#fitpack: $(FITPACK_OBJS)
#	$(FC) -shared -o $(FITPACK_LIB) \
#	$(foreach src, $(FITPACK_OBJS), $(addprefix $(INTER_DIR)/,$(src)))

lib: spline.obj $(ALL_OBJS) Library_BeamDyn_AM2.obj | $(INTER_DIR)
	$(FC) -shared -o $(BEAMDYN_LIB) \
	$(INTER_DIR)/spline.obj \
	$(INTER_DIR)/Library_BeamDyn_AM2.obj \
	$(foreach src, $(ALL_OBJS), $(addprefix $(INTER_DIR)/,$(src)))

cpp: $(CPP_DRIVER)

$(CPP_DRIVER): install $(CPP_DRIVER).obj
	$(CC) $(LDFLAGS) -I $(INTER_DIR) -o $(CPP_DRIVER) $(INTER_DIR)/$(CPP_DRIVER).obj $(BEAMDYN_LIB) -lstdc++ $(LAPACK_LINK)

   # Cleanup afterwards.

clean:
	$(DEL_CMD) $(INTER_DIR)$(PATH_SEP)*.mod $(INTER_DIR)$(PATH_SEP)*.obj $(OUTPUT_NAME)$(EXE_EXT) BeamDyn_Types.f90 *.dat *.out

install: lib
	@echo "Make sure that you've sourced the appropriate OpenFOAM rc file so that this installs correctly."
	cp -v beamDynInterface.H $(INC_INST_DIR)/
	cp -v $(BEAMDYN_LIB) $(LIB_INST_DIR)/
	cp -v $(LAPACK_DIR)/liblapack.* $(LIB_INST_DIR)/
	#cp -v $(LAPACK_DIR)/libblas.* $(LIB_INST_DIR)/

