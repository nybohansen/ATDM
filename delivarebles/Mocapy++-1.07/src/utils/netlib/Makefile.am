SLATECDIR = slatec
CEPHESDIR = cephes

noinst_LTLIBRARIES = libintegrate.la libfunctions.la

libintegrate_la_SOURCES = $(SLATECDIR)/dqags.f $(SLATECDIR)/dqagse.f $(SLATECDIR)/d1mach.f $(SLATECDIR)/dqelg.f $(SLATECDIR)/dqk21.f $(SLATECDIR)/dqpsrt.f $(SLATECDIR)/fdump.f $(SLATECDIR)/i1mach.f $(SLATECDIR)/j4save.f $(SLATECDIR)/xercnt.f $(SLATECDIR)/xerhlt.f $(SLATECDIR)/xermsg.f $(SLATECDIR)/xerprn.f $(SLATECDIR)/xersve.f $(SLATECDIR)/xgetua.f integrate.cpp integrate.h
libintegrate_la_LDFLAGS = -static

libfunctions_la_SOURCES = $(CEPHESDIR)/i0.c $(CEPHESDIR)/i1.c $(CEPHESDIR)/chbevl.c $(CEPHESDIR)/gamma.c $(CEPHESDIR)/const.c $(CEPHESDIR)/mtherr.c $(CEPHESDIR)/polevl.c $(CEPHESDIR)/isnan.c $(CEPHESDIR)/mconf.h  $(CEPHESDIR)/zeta.c functions.h 
libfunctions_la_LDFLAGS = -static

EXTRA_DIST = $(SLATECDIR)/readme $(CEPHESDIR)/readme $(CEPHESDIR)/cephes.doc
