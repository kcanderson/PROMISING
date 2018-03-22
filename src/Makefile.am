##AM_CPPFLAGS = -I$(top_srcdir)/include -O3 -fno-tree-pre -ftree-vectorize -fopt-info -march=native -mfpmath=sse -fopt-info-vec-optimized
AM_CPPFLAGS = -I$(top_srcdir)/include -std=c++0x
AM_LDFLAGS = -llapack -lblas
bin_PROGRAMS = promising reglaplacian pullentriesfrommat
promising_SOURCES = coreroutines.cpp CompleteGraphScorer.cpp PValueModuleScorer.cpp main.cpp FastScorer.cpp
reglaplacian_SOURCES = graph_kernels.cpp
pullentriesfrommat_SOURCES = coreroutines.cpp pullentriesfrommat.cpp
