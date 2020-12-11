SWinzip v2.5 (Sensing Wavelets in-situ zip)
-------------------------------------------

SWinzip is a Matlab and C++ library for scientific lossy data compression and
reconstruction using compressed sensing and tree-wavelets transforms. These
methods are known for their large compression and usefulness in data analytics
such as features extraction. Compressed sensing and wavelets methods rely
heavily on sparse and dense linear algebra operations implemented through the
Boost codes in our library. SWinzip accommodates data represented on both
regular grids (e.g. image data) and point-clouds (e.g. unstructured meshes).

In SWinzip v2.5, a HDF5 Alpert wavelet compression filter is added. It is based
on te matrix-free Alpert wavelet transform. Please check
swinzip-v2.5/hdf5/README.txt on how to use it.

Background information that is generally available is the following papers:

•	Salloum, M., Johnson, K.L., Bishop, J.E., Aytac, J.M., Dagel, D., van
Bloemen Waander, B.G. "Adaptive wavelet compression of large additive
manufacturing experimental and simulation datasets" Computational Mechanics,
63(3): 491–510, 2019.

•	Salloum, M., Fabian, N.D., Hensinger, D.M., Lee, J., Allendorf, E.M.,
Bhagatwala, A., Blaylock, M.L., Chen, J.H., Templeton, J.A., Tezaur, I. "
Optimal Compressed Sensing and Reconstruction of Unstructured Mesh Datasets "
Data Science and Engineering, 3(1): 1–23, 2018.

•	D. L. Donoho, Y. Tsaig, I. Drori, and J.-L. Starck. Sparse solution of
underdetermined systems of linear equations by stagewise orthogonal matching
pursuit. IEEE Transactions on Information Theory, 58(2):1094–1121, 2012.

•	B. Alpert, G. Beylkin, R. Coifman, and V. Rokhlin. Wavelet-like bases
for the fast solution of second-kind integral equations. SIAM Journal on
Scientific Computing, 14(1):159–184, 1993.

Author(s):  Maher Salloum, Sandia National Laboratories, mnsallo@sandia.gov
            David Hensinger, Sandia National Laboratories
            Nathan Fabian, Sandia National Laboratories
            Jina Lee, Sandia National Laboratories
            Kyle Karlson, Sandia National Laboratories
                           
The library is primarily intended for computers running the Unix, Linux,
operating system with the GNU Compiler Collection and/or Matlab. Building
and using the library requires Matlab (www.mathworks.com), the Boost library
(www.boost.org) and the Boost BLAS and LAPACK bindings to be installed on the
same operating system. The C++ library performance is limited so far to serial
computations however it can be linked to OpenBlas (http://www.openblas.net/) for
multicore computing. 


The results accuracy will depend on the CPU architecture used. The tests have
shown that it is best to run the library on Xeon CPUs with 4 or more cores.

3rd party open-source software used in the library:
---------------------------------------------------

- SparseLab: (Apache 2.0 license) available at: https://sparselab.stanford.edu/ 
The Stagewise Orthogonal Matching Pursuit (StOMP) Matlab algorithm
(SolveStOMP.m) was imported and improved.

- Boost 1.58: (Boost Software License 1.0) available at: http://www.boost.org/
Boost header files where included in C++ code.

- sparse-msrf: (BSD license) available at:
http://www.sandia.gov/~jairay/kfcs/kfcs.html#Releases
The sampling matrices construction Matlab code is imported

- Alpert Transform Toolbox: (BSD license) available at:
https://github.com/gpeyre/matlab-toolboxes/tree/master/toolbox_alpert
The Alpert wavelet matrix building Matlab algorithm was imported and improved.

- Python 3.0


Contents:
--------

/data: Contains data files that are used in the Matlab and C++ examples. A
README.txt file in the folder describes each dataset. 

C++ source code:
  src/Alpert
  src/Sampling
  src/StOMP
  src/Utils

C++ examples:
  /examples/compress_reconstruct:
A data compression/decompression pipeline using compressed sensing and wavelet
compression with fixed rate. Please refer to
/examples/compress_reconstruct/README.txt for more info

/examples/compress_wavelets:
Data compression/decompression using wavelet compression with fixed rate. Please
refer to /examples/compress_wavelets/README.txt for more info

/examples/wavelet_compress_decompress_distributed:
Optimal parallel data compression/decompression using wavelet compression with
given threshold Please refer to
/examples/wavelet_compress_decompress_distributed/README.txt for more info

    
Matlab source code:
  /matlab/src

Matlab examples:
  /matlab/examples
    
Include and Library folders: (do not delete!)
  /include
  /lib

Boost Lapack bindings
  /bindings
This folder has to be copied to: [boost_folder]/numeric/.

HDF5 filter
  /hdf5
This folder contains the C source code of Alpert compression as well as the HDF5
filter code H5Zalpert.C. It also contains a folder with an example code using
the filter.

Compiling the library and example:
----------------------

For the library:
- Go to /src
- adjust the required settings in 'CURRENT_defines'
(the C++ code will compile with gcc 4.8.5 but might compile with other versions.
It is only tested with Boost 1.58)
- run 'make'

For the example:
- Go to /examples/compress_reconstruct
- run 'make'

