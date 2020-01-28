# Many macros are defined in the below included site-specific configuration file
include ../CURRENT_defines

CPPFLAGS = -g -pedantic -Wunused -Wformat -std=c++11 -fPIC
OPTIMFLAGS = -ftemplate-depth-100 -funroll-loops -O3 -finline-functions -Wno-inline -DNDEBUG -DBOOST_UBLAS_CHECK_FALSE
DFLAGS  = 
LDFLAGS =
# To use other versions of BLAS, LAPACK or SLATEC, just change the macros below
# As a disclaimer: the Makefiles have only been tested with the versions of these
# libraries that come with the UQTk distribution.
#
AR = ar cr

TARGET = libUtils.a
CPPFILES = Utils.cpp
INCFILES = Utils.hpp
SRCS = $(CPPFILES)
OBJS = $(CPPFILES:.cpp=.o) 
INCDIRS = -I $(TRILINOSINSTALL)/include -I $(INC)/include $(BOOST_NUMERIC_BINDINGS) 
LDEP = Utils.cpp Utils.hpp

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



%.o : %.cpp
	$(CXX) $(CPPFLAGS) $(OPTIMFLAGS) $(DEFS) $(INCDIRS) -c $< -o $@
