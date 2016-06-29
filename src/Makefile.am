bin_PROGRAMS = swap

swap_CPPFLAGS = \
	-I$(top_srcdir)/SnowTools/src \
	-I$(top_srcdir)/SnowTools/htslib \
	-I$(top_srcdir)/SnowTools/multifast-v1.4.2 

swap_LDADD = $(top_builddir)/SnowTools/src/libsnowtools.a $(top_builddir)/SnowTools/htslib/libhts.a

swap_LDFLAGS = -pthread -std=c++11

swap_SOURCES = swap.cpp Matrix.cpp prng_engine.hpp

