ACLOCAL_AMFLAGS = -I m4

pkgconfigdir = $(libdir)/pkgconfig
pkgconfig_DATA = lightpipes.pc

if OCTAVE_ENABLED
    OCTAVE_DIR = octave
endif

SUBDIRS = src $(OCTAVE_DIR)
DIST_SUBDIRS = src octave
EXTRA_DIST = manual scripts

tests :
	cd examples; make

apidoc :
	cd src/lib; make apidoc
