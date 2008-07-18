# Copyright 2004-2006 Spencer Olson

MKOCTFILE = @MKOCTFILE@

# Destination for octave files.
octave_dir = $(libexecdir)

# Make sure that our .oct files do not have extensions
EXEEXT=

bin_PROGRAMS =          \
  LPAxicon.oct          \
  LPBegin.oct           \
  LPCircAperture.oct    \
  LPCircScreen.oct      \
  LPConvert.oct         \
  LPForward.oct         \
  LPForvard.oct         \
  LPFresnel.oct         \
  LPGain.oct            \
  LPGaussAperture.oct   \
  LPGaussScreen.oct     \
  LPInterpol.oct        \
  LPLens.oct            \
  LPLensForvard.oct     \
  LPLensFresnel.oct     \
  LPNormal.oct          \
  LPPipFFT.oct          \
  LPRectAperture.oct    \
  LPRectScreen.oct      \
  LPStrehl.oct          \
  LPSuperGaussAperture.oct \
  LPSuperGaussScreen.oct \
  LPTilt.oct            \
  LPToroidalLens.oct    \
  LPZernike.oct

LPOctave_EXTRAHEADERS = oct-Field.h

LPAxicon_oct_SOURCES  = LPAxicon.C $(LPOctave_EXTRAHEADERS)
LPBegin_oct_SOURCES  = LPBegin.C $(LPOctave_EXTRAHEADERS)
LPCircAperture_oct_SOURCES  = LPCircAperture.C $(LPOctave_EXTRAHEADERS)
LPCircScreen_oct_SOURCES  = LPCircScreen.C $(LPOctave_EXTRAHEADERS)
LPConvert_oct_SOURCES  = LPConvert.C $(LPOctave_EXTRAHEADERS)
LPForward_oct_SOURCES  = LPForward.C $(LPOctave_EXTRAHEADERS)
LPForvard_oct_SOURCES  = LPForvard.C $(LPOctave_EXTRAHEADERS)
LPFresnel_oct_SOURCES  = LPFresnel.C $(LPOctave_EXTRAHEADERS)
LPGain_oct_SOURCES  = LPGain.C $(LPOctave_EXTRAHEADERS)
LPGaussAperture_oct_SOURCES  = LPGaussAperture.C $(LPOctave_EXTRAHEADERS)
LPInterpol_oct_SOURCES     = LPInterpol.C $(LPOctave_EXTRAHEADERS)
LPGaussScreen_oct_SOURCES  = LPGaussScreen.C $(LPOctave_EXTRAHEADERS)
LPLens_oct_SOURCES  = LPLens.C $(LPOctave_EXTRAHEADERS)
LPLensForvard_oct_SOURCES  = LPLensForvard.C $(LPOctave_EXTRAHEADERS)
LPLensFresnel_oct_SOURCES  = LPLensFresnel.C $(LPOctave_EXTRAHEADERS)
LPNormal_oct_SOURCES  = LPNormal.C $(LPOctave_EXTRAHEADERS)
LPPipFFT_oct_SOURCES  = LPPipFFT.C $(LPOctave_EXTRAHEADERS)
LPRectAperture_oct_SOURCES  = LPRectAperture.C $(LPOctave_EXTRAHEADERS)
LPRectScreen_oct_SOURCES  = LPRectScreen.C $(LPOctave_EXTRAHEADERS)
LPStrehl_oct_SOURCES  = LPStrehl.C $(LPOctave_EXTRAHEADERS)
LPSuperGaussAperture_oct_SOURCES  = LPSuperGaussAperture.C $(LPOctave_EXTRAHEADERS)
LPSuperGaussScreen_oct_SOURCES  = LPSuperGaussScreen.C $(LPOctave_EXTRAHEADERS)
LPTilt_oct_SOURCES  = LPTilt.C $(LPOctave_EXTRAHEADERS)
LPToroidalLens_oct_SOURCES  = LPToroidalLens.C $(LPOctave_EXTRAHEADERS)
LPZernike_oct_SOURCES  = LPZernike.C $(LPOctave_EXTRAHEADERS)

# I haven't finished the C++ library for these functions.
# LPGaussHermite
# LPPhaseUnwrap             # Probably just use the unwrap function in octave.
# LPReflectMultilayer
# LPSteps
# LPStrehl
# LPTransmitMultilayer

# We're going to put it in the libexec dir.
bindir = $(octave_dir)

# Some of these compare to functions in the non-gpl matlab version of
# lightpipes.  These are all implemented in octave.  
mfilesdir = $(octave_dir)
dist_mfiles_DATA =  \
    LPBeamMix.m \
    LPCompatible.m \
    LPIntensity.m \
    LPPolarizer.m \
    LPRandomIntensity.m \
    LPRandomPhase.m


# These will not be reimplemented for octave; use octave methods, they are too
# easy.  They will get installed anyway, soas to provide a warning when
# someone attempts to use them.
deadmfilesdir = $(octave_dir)
dist_deadmfiles_DATA = \
    LPIntAttenuator.m \
    LPMultHelp.m \
    LPMultIntensity.m \
    LPMultPhase.m \
    LPPhase.m \
    LPSubIntensity.m \
    LPSubPhase.m

AM_CPPFLAGS = $(shell ${MKOCTFILE} -p INCFLAGS) $(LIGHTPIPES_INCLUDES)

LIGHTPIPES_INCLUDES   = -I../src/include
LIGHTPIPES_LIBS       = -L../src/lib -llightpipes

#SH_LDFLAGS_OCT=$(foreach opt, $(shell ${MKOCTFILE} -p SH_LDFLAGS), -Xlinker $(opt))
#DL_LDFLAGS_OCT=$(foreach opt, $(shell ${MKOCTFILE} -p DL_LDFLAGS), -Xlinker $(opt))
SH_LDFLAGS_OCT=$(shell ${MKOCTFILE} -p SH_LDFLAGS)
DL_LDFLAGS_OCT=$(shell ${MKOCTFILE} -p DL_LDFLAGS)

PWD=$(shell pwd)

#%.o:%.C
#	$(CXX) $(DEFS) $(DEFAULT_INCLUDES) $(INCLUDES) $(AM_CPPFLAGS) \
#	$(CPPFLAGS) $(AM_CFLAGS) $(CFLAGS) $(LIGHTPIPES_INCLUDES) -c $< -o $(PWD)/$@
# can't seem to get CXXCOMPILE variable created by automake.

LPAxicon.oct          : $(LPAxicon_oct_OBJECTS)
LPBegin.oct           : $(LPBegin_oct_OBJECTS)
LPCircAperture.oct    : $(LPCircAperture_oct_OBJECTS)
LPCircScreen.oct      : $(LPCircScreen_oct_OBJECTS)
LPConvert.oct         : $(LPConvert_oct_OBJECTS)
LPForward.oct         : $(LPForward_oct_OBJECTS)
LPForvard.oct         : $(LPForvard_oct_OBJECTS)
LPFresnel.oct         : $(LPFresnel_oct_OBJECTS)
LPGain.oct            : $(LPGain_oct_OBJECTS)
LPGaussAperture.oct   : $(LPGaussAperture_oct_OBJECTS)
LPInterpol.oct        : $(LPInterpol_oct_OBJECTS)
LPGaussScreen.oct     : $(LPGaussScreen_oct_OBJECTS)
LPLens.oct            : $(LPLens_oct_OBJECTS)
LPLensForvard.oct     : $(LPLensForvard_oct_OBJECTS)
LPLensFresnel.oct     : $(LPLensFresnel_oct_OBJECTS)
LPNormal.oct          : $(LPNormal_oct_OBJECTS)
LPPipFFT.oct          : $(LPPipFFT_oct_OBJECTS)
LPRectAperture.oct    : $(LPRectAperture_oct_OBJECTS)
LPRectScreen.oct      : $(LPRectScreen_oct_OBJECTS)
LPStrehl.oct          : $(LPStrehl_oct_OBJECTS)
LPSuperGaussAperture.oct : $(LPSuperGaussAperture_oct_OBJECTS)
LPSuperGaussScreen.oct : $(LPSuperGaussScreen_oct_OBJECTS)
LPTilt.oct            : $(LPTilt_oct_OBJECTS)
LPToroidalLens.oct    : $(LPToroidalLens_oct_OBJECTS)
LPZernike.oct         : $(LPZernike_oct_OBJECTS)
%.oct:%.o
	@echo "[$(MKOCTFILE) $@]"
	@export CC="$(CC)" \
	export CXX="$(CXX)" \
	export SH_LD="$(CXX)" \
	export DL_LD="$(CXX)" \
	export F77="$(F77)" \
	export CFLAGS="$(CFLAGS)"; \
	export CPPFLAGS="$(LIGHTPIPES_INCLUDES) $(CPPFLAGS)"; \
	export CXXFLAGS="$(CXXFLAGS)"; \
	export FPICFLAG="$(FPICFLAG)"; \
        export SH_LDFLAGS="$(SH_LDFLAGS_OCT) $(LDFLAGS)"; \
        export DL_LDFLAGS="$(DL_LDFLAGS_OCT) $(LDFLAGS)"; \
        $(MKOCTFILE) $< $(LIBS) $(LIGHTPIPES_LIBS) -o $@
