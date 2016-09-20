
##AM_CPPFLAGS = -I$(top_srcdir)/include -O3 -fno-tree-pre -ftree-vectorize -fopt-info -march=native -mfpmath=sse -fopt-info-vec-optimized
AM_CPPFLAGS = -I$(top_srcdir)/include -O0
AM_LDFLAGS = -llapack -lblas
bin_PROGRAMS = promising blah
promising_SOURCES = coreroutines.cpp CompleteGraphScorer.cpp PValueModuleScorer.cpp main.cpp FastScorer.cpp
blah_SOURCES = graph_kernels.cpp
