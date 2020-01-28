# Trilinos install folder
include ../../src/CURRENT_defines

INC = /usr/local
OPT = 3
CPPFLAGS = -g -pedantic -Wunused -Wformat -Wimplicit -std=c++11
OPTIMFLAGS = -ftemplate-depth-100 -funroll-loops -O3 -finline-functions -Wno-inline -DNDEBUG -DBOOST_UBLAS_CHECK_FALSE
DFLAGS  = 
LDFLAGS =

TARGET = test_case
SRCS = main.cpp 
OBJS = $(SRCS:.cpp=.o)
.SUFFIXES : .cpp

LIBINCDIR = ../../include
LIBDIR    = ../../lib
INCDIRS = -I . -I $(LIBINCDIR) -I $(TRILINOSINSTALL)/include -I$(COMPDIR)/include $(BOOST_NUMERIC_BINDINGS)

LIBS = -L$(COMPDIR)/lib -L/usr/local/lib -L$(LIBDIR) -lSampling -lAlpert -lUtils -lcstomp -llapack 
ifdef OPENBLAS
  LIBS+= -L$(OPENBLAS) -lopenblas
else
  LIBS+= -lblas
endif

ifdef OPENMP
  LIBS+= -fopenmp
endif

all: $(TARGET)

$(TARGET): $(OBJS) $(LDEP)
	$(CXX) $(CPPFLAGS) $(OPTIMFLAGS) $(INCDIRS) -o $(TARGET) $(OBJS) $(LIBS)

clean:
	rm -f *~ $(OBJS) $(TARGET)

.cpp.o:
	$(CXX) $(CPPFLAGS) $(OPTIMFLAGS) $(INCDIRS) -c $*.cpp

