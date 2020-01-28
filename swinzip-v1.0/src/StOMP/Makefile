# Many macros are defined in the below included site-specific configuration file
INC = /usr/local
CPPFLAGS = -g -pedantic -Wunused -Wformat -Wimplicit 
OPTIMFLAGS = -O3 -Wno-inline -DNDEBUG -fPIC -std=c++11
DFLAGS  = 
LDFLAGS =
# To use other versions of BLAS, LAPACK or SLATEC, just change the macros below
# As a disclaimer: the Makefiles have only been tested with the versions of these
# libraries that come with the UQTk distribution.
#
RANLIB = ranlib
AR = ar cr

TARGET = libcstomp.a
SRCS = $(CPPFILES)
OBJS = $(CPPFILES:.cpp=.o) 
INCDIRS = -I $(INC)/include

all: $(TARGET)

$(TARGET): ccstomp.o
	$(AR) $(TARGET) $(LIB_OBJS) ccstomp.o
	$(RANLIB) $(TARGET)
	rm -f ../../lib/$(TARGET)
	cp $(TARGET) ../../lib/$(TARGET)
	for f in *.hpp;\
         do rm -f ../../include/$$f;\
         cp $$f ../../include/$$f;\
        done;

ccstomp.o:ccstomp.cpp
	$(CXX) -fpic $(OPTIMFLAGS) -c ccstomp.cpp

clean:
	rm -f *~ $(OBJS) $(MODS) $(TARGET) *.o

.cpp.o:
	$(CXX) $(CPPFLAGS) $(OPTIMFLAGS) $(DEFS) $(INCDIRS) -c $(SRCS)
