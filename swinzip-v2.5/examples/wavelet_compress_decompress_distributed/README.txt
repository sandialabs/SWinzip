This example provides a pipeline of parallel compression of data files and
decompression using optimal threshold-based Alpert wavelets transforms.

These codes can be tweaked to accomodate any given data files and formats.

The data provided here assumes that the first line contains the number of lines
and columns in the file

Usage:
-----

Compression:
mpirun -np [number of processes] ./wcompress [threshold value] [path to data files+filename root] [number of physical dimensions] [wavelet order]

Examples:
mpirun -np 64 ./wcompress 0.025e0 ../../data/Greenland/gis 3 3
mpirun -np 62 ./wcompress 0.025e0 ../../data/freejet/jet 3 4
mpirun -np 40 ./wcompress 0.025e0 ../../data/channel/ksgs 3 2

Deompression:
mpirun -np [number of processes] ./wdecompress [mesh filename root] [compressed data filename root] [reconstructed data filename root]

Examples:
mpirun -np 64 ./wdecompress mesh_ fw_ out_r
mpirun -np 62 ./wdecompress mesh_ fw_ out_r
mpirun -np 40 ./wdecompress mesh_ fw_ out_r


The matlab .m files are used to plot the results in the .jpg images