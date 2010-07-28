# __BEGIN_LICENSE__
# Copyright (C) 2006-2010 United States Government as represented by
# the Administrator of the National Aeronautics and Space Administration.
# All Rights Reserved.
# __END_LICENSE__


########################################################################
# sources
########################################################################

if MAKE_MODULE_MULTIVIEWTK

TestPatchMultiView_SOURCES = TestPatchMultiView.cxx
TestAlbedoOptimizers_SOURCES = TestAlbedoOptimizers.cxx

TESTS = TestPatchMultiView TestAlbedoOptimizers

endif

########################################################################
# general
########################################################################

AM_CPPFLAGS = @ASP_CPPFLAGS@
AM_LDFLAGS  = @ASP_LDFLAGS@ @PKG_MULTIVIEWTK_LIBS@

check_PROGRAMS = $(TESTS)

include $(top_srcdir)/config/rules.mak
include $(top_srcdir)/config/tests.am
