#
# Makefile for Cray XC30.
#

#For XC30
#CC = /home/saitoutk/bin/bin/ccache cc
CC = cc

#FOR GNU
#OPTIMIZE = -std=c99
#OPTIMIZE += -O3 -finline-functions -ffast-math -funroll-loops -frerun-cse-after-loop -fexpensive-optimizations -fomit-frame-pointer -mfpmath=sse -fprefetch-loop-arrays
#OPTIMIZE += -Winline -Wundef -Wmissing-noreturn -Wstrict-prototypes
#OPTIMIZE += -g

#FOR Intel
#OPTIMIZE = -std=c99  -O3
#OPTIMIZE = -std=c99 -g
#OPTIMIZE = -std=c99 -O3 
#OPTIMIZE += -ipo 
#OPTIMIZE += -Ofast -fp-model fast=2

#FOR Cray
OPTIMIZE = -O3 
OPTIMIZE += -hfp3

#OPTIMIZE += -ra
## OPTIMIZE += -h list=a
#OPTIMIZE += -h bounds
# #OPTIMIZE += -fprofile-arcs -ftest-coverage  # for gcov
# #OPTIMIZE += -pg 
#OPTIMIZE += -K trap=divz,fp,inv
#OPTIMIZE += -G0 


GRAPE="XC50"
ifeq ($(UNAME),masakawats)
INCLUDEPATH += -I/home/masakawats/lib/avx_phantom
#INCLUDEPATH += -I/home/matsuihd/lib/gsl/include
INCLUDEPATH += -I${GSL_DIR}/include
INCLUDEPATH += -I/home/masakawats/lib/celib/celib_hirai2
INCLUDEPATH += -I/home/masakawats/lib/CloudyCooling
OPTIONS += -DHAVE_AVX_PHANTOM_GRAPE
LDFLAGS += -L/home/masakawats/lib/avx_phantom -lpg5
#LDFLAGS += -L/home/matsuihd/lib/gsl/lib -lgsl -lgslcblas
LDFLAGS += -L${GSL_DIR}/lib -lgsl -lgslcblas
LDFLAGS += -L/home/masakawats/lib/celib/celib_hirai2 -lCELib
LDFLAGS += -L/home/masakawats/lib/CloudyCooling -lCloudyCooling
endif
#ifeq ($(UNAME),hiraiyt)
#INCLUDEPATH += -I/home/hiraiyt/ASURA_ver3.5.0/lib/avx_phantom
#INCLUDEPATH += -I/home/saitoh/lib/gsl/include
#INCLUDEPATH += -I${GSL_DIR}/include
#INCLUDEPATH += -I/home/hiraiyt/ASURA_ver3.5.0/lib/CELib
#INCLUDEPATH += -I/home/hiraiyt/celib/src_20170720
#INCLUDEPATH += -I/home/hiraiyt/celib/src_20170919
#INCLUDEPATH += -I/home/hiraiyt/celib/src_20170919_c
#INCLUDEPATH += -I/home/hiraiyt/ASURA_ver3.5.0/lib/CloudyCooling
#OPTIONS += -DHAVE_AVX_PHANTOM_GRAPE
#LDFLAGS += -L/home/hiraiyt/ASURA_ver3.5.0/lib/avx_phantom -lpg5
#LDFLAGS += -L/home/saitoh/lib/gsl/lib -lgsl -lgslcblas
#LDFLAGS += -L${GSL_DIR}/lib -lgsl -lgslcblas
#LDFLAGS += -L/home/hiraiyt/ASURA_ver3.5.0/lib/CELib -lCELib
#LDFLAGS += -L/home/hiraiyt/celib/src_20170720 -lCELib
#LDFLAGS += -L/home/hiraiyt/celib/src_20170919 -lCELib
#LDFLAGS += -L/home/hiraiyt/celib/src_20170919_c -lCELib
#LDFLAGS += -L/home/hiraiyt/ASURA_ver3.5.0/lib/CloudyCooling -lCloudyCooling
#endif
#ifeq ($(UNAME),hiraiyt)
#INCLUDEPATH += -I/home/saitoutk/lib/saitoh_phantom
#INCLUDEPATH += -I/home/saitoutk/lib/gsl/include
#INCLUDEPATH += -I/home/saitoutk/lib/CloudyCooling
#INCLUDEPATH += -I/home/saitoutk/lib/CELib
#OPTIONS += -DHAVE_PHANTOM_GRAPE
#LDFLAGS += -L/home/saitoutk/lib/saitoh_phantom -lpg5pot
#LDFLAGS += -L/home/saitoutk/lib/gsl/lib -lgsl -lgslcblas
#LDFLAGS += -L/home/saitoutk/lib/CELib -lCELib
#LDFLAGS += -L/home/saitoutk/lib/CloudyCooling -lCloudyCooling
#endif
