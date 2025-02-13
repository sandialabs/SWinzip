# Many macros are defined in the below included site-specific configuration file
include ../CURRENT_defines
#
#INC = /Users/ndfabia/boost_install
CPPFLAGS = -g -pedantic -Wunused -Wformat -fPIC
ifdef OPENMP
  CPPFLAGS+= -DOPENMP
endif

OPTIMFLAGS = -ftemplate-depth-100 -funroll-loops -O3 -finline-functions -Wno-inline -DNDEBUG -DBOOST_UBLAS_CHECK_FALSE -DWITH_LAPACK=OFF
ifdef OPENMP
  OPTIMFLAGS+= -fopenmp
endif

DFLAGS  = 
LDFLAGS =
# To use other versions of BLAS, LAPACK or SLATEC, just change the macros below
# As a disclaimer: the Makefiles have only been tested with the versions of these
# libraries that come with the UQTk distribution.
#
AR = ar cr
INC = ../..

TARGET = libAlpert.a
CPPFILES = Alpert_Matrix.cpp Alpert_Transform.cpp
INCFILES = Alpert_Matrix.hpp Alpert_Transform.hpp
SRCS = $(CPPFILES)
OBJS = $(CPPFILES:.cpp=.o)    $(CPPFILES:.cpp=.so) 
INCDIRS = -I $(TRILINOSINSTALL)/include -I $(INC)/include $(BOOST_NUMERIC_BINDINGS) -I ../Utilts/
LDEP = Alpert_Matrix.cpp Alpert_Transform.cpp Alpert_Matrix.hpp Alpert_Transform.hpp

all: $(TARGET)

$(TARGET): $(OBJS) $(LDEP)
	$(AR) $(TARGET) $(OBJS)
	$(RANLIB) $(TARGET)
	rm -f ../../lib/$(TARGET)
	cp $(TARGET) ../../lib/$(TARGET)
	for f in *.hpp;\
	 do rm -f ../../include/$$f;\
         cp $$f ../../include/$$f;\
        done;

clean:
	rm -f *~ $(OBJS) $(MODS) $(TARGET) *.o

#.cpp.o:
#	$(CXX) $(CPPFLAGS) $(OPTIMFLAGS) $(DEFS) $(INCDIRS) -c $(SRCS)


%.so : %.cpp
	$(CXX) $(CPPFLAGS) $(OPTIMFLAGS) $(DEFS) $(INCDIRS) -shared -c $< -o $@
%.o : %.cpp
	$(CXX) $(CPPFLAGS) $(OPTIMFLAGS) $(DEFS) $(INCDIRS) -c $< -o $@
