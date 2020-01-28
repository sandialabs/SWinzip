This example provides a pipeline of compression and decompression using
compressed sensing and Alpert wavelets.

In this version of SWinzip, this pipeline relies on openMP multicore computing.
Future versions will also accommodate multi-node computing.

please refer to /src/CURRENT_defines for inforamtion on how to set the
multicore settings.

Usage:
-----

./test_case --seed 'seed_value' --reduction 'reduction_value' [mesh input_file] [data input_file]

If the data input file is provided it should have the same number of pointsas
the mesh, if not a default data field will be assigned

- The 'seed' is used for random number generation.
- The 'reduction' refers to the compression ratio required
- The first line in [mesh input_file] should contain the number of pointsand the
number of dimensions. The rest of the lines contain the points coordinates i.e.
the number of columns should be the same as the second number in the first line.
- [data input_file] should consist of one column where the number of lines is
equal to the number of points in the mesh i.e. the first number in the first
line of [mesh input_file]

Both [mesh input_file] and [data input_file] should be in ASCII format.
