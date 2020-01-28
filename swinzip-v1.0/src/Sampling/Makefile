# Many macros are defined in the below included site-specific configuration file
include ../CURRENT_defines

SHELL = /bin/sh
INC = /usr/local
CPPFLAGS = -g -pedantic -Wunused -Wformat -Wimplicit -std=c++11  -fPIC
OPTIMFLAGS = -ftemplate-depth-100 -funroll-loops -O3 -finline-functions -Wno-inline -DNDEBUG -DBOOST_UBLAS_CHECK_FALSE
RANLIB = ranlib
AR = ar cr

TARGET = libSampling.a
CPPFILES = Sampling_Matrix.cpp
INCFILES = Sampling_Matrix.hpp
SRCS = $(CPPFILES)
OBJS = $(CPPFILES:.cpp=.o) 
INCDIRS = -I $(TRILINOSINSTALL)/include -I $(INC)/include
LDEP = Sampling_Matrix.cpp Sampling_Matrix.hpp

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

.cpp.o:
	$(CXX) $(CPPFLAGS) $(OPTIMFLAGS) $(DEFS) $(INCDIRS) -c $(SRCS)
