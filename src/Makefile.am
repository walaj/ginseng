bin_PROGRAMS = ginseng

ginseng_CPPFLAGS = \
	-I$(top_srcdir)/SeqLib \
	-I$(top_srcdir) \
	-I$(top_srcdir)/SeqLib/htslib \
	-Wno-sign-compare

ginseng_LDADD = $(top_builddir)/SeqLib/src/libseqlib.a $(top_builddir)/SeqLib/htslib/libhts.a

ginseng_SOURCES = FishModel.cpp gzstream.C ginseng.cpp swap.cpp Histogram.cpp Matrix.cpp \
		prng_engine.hpp Fractions.cpp FishHookInterval.cpp sim.cpp fishhook.cpp \
		Model.cpp BEDIntervals.cpp Events.cpp
