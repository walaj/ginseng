bin_PROGRAMS = ginseng

ginseng_CPPFLAGS = \
	-I$(top_srcdir)/SeqLib \
	-I$(top_srcdir)/SeqLib/htslib \
	-I$(top_srcdir)/apophenia-pkg/include \
	-I/xchip/gistic/Jeremiah/software/gsl-2.2.1/include -Wno-sign-compare

ginseng_LDADD = $(top_builddir)/SeqLib/src/libseqlib.a $(top_builddir)/SeqLib/htslib/libhts.a \
		$(top_builddir)/apophenia-pkg/lib/libapophenia.a -lgsl -lgslcblas -lsqlite3

ginseng_LDFLAGS = -L/xchip/gistic/Jeremiah/software/gsl-2.2.1/lib -fopenmp

ginseng_SOURCES = FishModel.cpp gzstream.C ginseng.cpp swap.cpp Histogram.cpp Matrix.cpp \
		prng_engine.hpp Fractions.cpp FishHookInterval.cpp sim.cpp fishhook.cpp \
		Model.cpp BEDIntervals.cpp Events.cpp
