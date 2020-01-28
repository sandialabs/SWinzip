# Trilinos install folder
#include ../src/MNS_defines
include ../CURRENT_defines
# Openmpi folder
INC = /usr/local
OPT = 3
CPPFLAGS = -g -pedantic -Wunused -Wformat -std=c++11  -fPIC
CFLAGS = -g -Wunused -Wformat
#OPTIMFLAGS = -ftemplate-depth-100 -funroll-loops -O3 -finline-functions -Wno-inline -DNDEBUG -DBOOST_UBLAS_CHECK_FALSE
DFLAGS  = 
LDFLAGS =

AR = ar cr


TARGET = libexo_rw.a
CSRCS = exo_rw.c
COBJS = $(CSRCS:.c=.o)
.SUFFIXES : .cpp .c

LIBINCDIR = ../../include
LIBDIR    = ../../lib
INCDIRS = -I . -I $(LIBINCDIR) -I $(TRILINOSINSTALL)/include -I$(COMPDIR)/include -I$(NETCDF_DIR)/include


all: $(TARGET)

$(TARGET): $(COBJS)
	$(AR) $(TARGET) $(COBJS)
	$(RANLIB) $(TARGET)
	rm -f ../../lib/$(TARGET)
	cp $(TARGET) ../../lib/$(TARGET)
	for f in *.h;\
         do rm -f ../../include/$$f;\
         cp $$f ../../include/$$f;\
        done;

clean:
	rm -f *~ $(OBJS) $(COBJS) $(TARGET)

.cpp.o:
	$(CXX) $(CPPFLAGS) $(OPTIMFLAGS) $(INCDIRS) -c $*.cpp

.c.o:
	$(CC) $(CFLAGS) $(OPTIMFLAGS) $(INCDIRS) -c $*.c

clean:
	rm -f *~ $(OBJS) $(MODS) $(TARGET) *.o
