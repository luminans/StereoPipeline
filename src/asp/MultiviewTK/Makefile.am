# __BEGIN_LICENSE__
# Copyright (C) 2006-2010 United States Government as represented by
# the Administrator of the National Aeronautics and Space Administration.
# All Rights Reserved.
# __END_LICENSE__


#########################################################################
# sources
#########################################################################

if MAKE_MODULE_MULTIVIEWTK

include_HEADERS =

if ENABLE_EXCEPTIONS

MULTIVIEWTK_LOCAL_LIBS = @MODULE_MULTIVIEWTK_LIBS@

multiview_SOURCES = multiview.cc
multiview_LDADD   = $(MULTIVIEWTK_LOCAL_LIBS)

bin_PROGRAMS = multiview

endif

endif

#########################################################################
# general
#########################################################################

AM_CPPFLAGS = @ASP_CPPFLAGS@
AM_LDFLAGS = @ASP_LDFLAGS@

includedir = $(prefix)/include/asp/MultiviewTK

include $(top_srcdir)/config/rules.mak
