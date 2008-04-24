ACLOCAL_AMFLAGS = -I m4

pkgconfigdir = $(libdir)/pkgconfig
pkgconfig_DATA = lightpipes.pc

if OCTAVE_ENABLED
    OCTAVE_DIR = octave
endif

SUBDIRS = src $(OCTAVE_DIR)
DIST_SUBDIRS = src octave

tests :
	cd examples; make





include doxygen.am

MOSTLYCLEANFILES = $(DX_CLEANFILES)

if DX_COND_doc
all-local: doxygen-doc
APIDOCS = doc/api
endif

EXTRA_DIST = manual scripts $(DX_CONFIG) $(APIDOCS)

