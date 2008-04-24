# Copyright 2006 Spencer Olson

lib_LIBRARIES = liblightpipes.a

liblightpipes_a_SOURCES = \
        Field.C Field.h \
        fftn.c fftn.h

# Make it easy to install the headers
liblightpipes_a_includedir = $(includedir)/lightpipes
nobase_liblightpipes_a_include_HEADERS = \
        Field.h \
        $(PACKAGE)-config.h

LP_INC=../include/lightpipes
$(shell ${mkdir_p} ${LP_INC})
$(shell ${INSTALL} -p ${srcdir}/Field.h ${PACKAGE}-config.h ${LP_INC}/ )

PWD=$(shell pwd)
#OUTESC=$(shell pwd | sed -e "s/\//\\\\\//g")
OUTESC=$(shell ${top_srcdir}/fullpath.sh ${top_builddir} | sed -e "s/\//\\\\\//g")
