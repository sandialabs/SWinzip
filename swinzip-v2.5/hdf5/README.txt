This folder contains the C code of the Alpert wavelets library. It is equivalent
to the one in swinzip-v2.0/src/Alpert

This C code can be compiled as a library but the intent here is to use it to
make a HDF5 filter for data compression using Alpert wavelets. This compression
currently works on regular grid data only, which is the norm in HDF5.

In order to compile the code with HDF5 and use it, copy the files:
init_alpert.c
init_alpert.h
alpert_transform.c
alpert_transform.h
alpert_utils.c
alpert_utils.h
H5Zalpert.c

into the hdf5-x.xx.x/src source folder and add to hdf5-x.xx.x/src/Makefile.in
the following line:
H5Zalpert.c alpert_utils.c init_alpert.c alpert_transform.c \
in the block:
am__libhdf5_la_SOURCES_DIST =

and the following line:
H5Zalpert.lo alpert_utils.lo init_alpert.lo alpert_transform.lo \
in the block:
am_libhdf5_la_OBJECTS =

and the following line:
H5Zalpert.c alpert_utils.c init_alpert.c alpert_transform.c \
in the block:
libhdf5_la_SOURCES =

Then compile and build HDF5.

There is an example in the code swinzip-v2.0/hdf5/example/test_Alpert_3D.c
to show how the Alpert filter can be used. It can compiled and built using:

gcc -O2 -I /home/mnsallo/hdf5/include -o test_Alpert_3D test_Alpert_3D.c -lm -lblas -llapack -L /home/mnsallo/hdf5/lib -lhdf5

Running the executable will output the original data followed by the decompressed data.
For comparison, run the matlab file run_compare_hdf5.m found in the same folder.

