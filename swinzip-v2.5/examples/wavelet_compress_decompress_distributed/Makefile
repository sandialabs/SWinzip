include ../../src/CURRENT_defines

INC = /usr/local
OPT = 3
CPPFLAGS = -g -pedantic -Wunused -Wformat -Wimplicit -std=c++11
OPTIMFLAGS = -ftemplate-depth-100 -funroll-loops -O3 -finline-functions -Wno-inline -DNDEBUG -DBOOST_UBLAS_CHECK_FALSE
DFLAGS  = 
LDFLAGS =

TARGET1 = wcompress
TARGET2 = wdecompress
SRCS1 = main_compress.cpp
OBJS1 = $(SRCS1:.cpp=.o)
SRCS2 = main_decompress.cpp
OBJS2 = $(SRCS2:.cpp=.o)
.SUFFIXES : .cpp

LIBINCDIR = ../../include
LIBDIR    = ../../lib
INCDIRS = -I . -I $(LIBINCDIR) -I$(COMPDIR)/include

LIBS = -L$(COMPDIR)/lib -L/usr/local/lib -L$(LIBDIR) -lAlpert -lUtils -L/usr/lib64/ -llapack -lgfortran

ifdef OPENBLAS
  LIBS+= -L$(OPENBLAS) -lopenblas
else
  LIBS+= -lblas
endif

ifdef OPENMP
  LIBS+= -fopenmp
endif

all: $(TARGET1) $(TARGET2)

$(TARGET1): $(OBJS1) $(LDEP)
	$(MPICXX) $(CPPFLAGS) $(OPTIMFLAGS) $(INCDIRS) -o $(TARGET1) $(OBJS1) $(LIBS)

$(TARGET2): $(OBJS2) $(LDEP)
	$(MPICXX) $(CPPFLAGS) $(OPTIMFLAGS) $(INCDIRS) -o $(TARGET2) $(OBJS2) $(LIBS)

clean:
	rm -f *~ $(OBJS1) $(TARGET1) $(OBJS2) $(TARGET2)

.cpp.o:
	$(MPICXX) $(CPPFLAGS) $(OPTIMFLAGS) $(INCDIRS) -c $*.cpp

