bin_PROGRAMS = swap

swap_CPPFLAGS = \
	-I$(top_srcdir)/SeqLib \
	-I$(top_srcdir)/SeqLib/htslib 

swap_LDADD = $(top_builddir)/SeqLib/src/libseqlib.a $(top_builddir)/SeqLib/htslib/libhts.a

swap_LDFLAGS = -pthread -std=c++11

swap_SOURCES = gzstream.C ginseng.cpp swap.cpp Histogram.cpp Matrix.cpp prng_engine.hpp Fractions.cpp FishHookInterval.cpp sim.cpp fishhook.cpp

