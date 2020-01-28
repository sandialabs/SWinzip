SWinzip v1.0 (Sensing Wavelets in-situ zip)
--------------------

SWinzip is a Matlab and C++ library for scientific lossy data compression and
reconstruction using compressed sensing and tree-wavelets transforms. These
methods are known for their large compression and usefulness in data analytics
such as features extraction. Compressed sensing and wavelets methods rely
heavily on sparse and dense linear algebra operations implemented through the
Boost codes in our library. SWinzip accommodates data represented on both
regular grids (e.g. image data) and point-clouds (e.g. unstructured meshes).
Background information that is generally available is the following papers:

•	Salloum, M., Fabian, N., Hensinger, D.M. Templeton, J.A. " Compressed
Sensing and Reconstruction of Unstructured Mesh Datasets ", 2015, SAND
2015-4995C Available at http://arxiv.org/pdf/1508.06314

•	D. L. Donoho, Y. Tsaig, I. Drori, and
J.-L. Starck. Sparse solution of underdetermined systems of linear equations by
stagewise orthogonal matching pursuit. IEEE Transactions on Information Theory,
58(2):1094–1121, 2012.

•	B. Alpert, G. Beylkin, R. Coifman, and V. Rokhlin. Wavelet-like bases
for the fast solution of second-kind integral equations. SIAM Journal on
Scientific Computing, 14(1):159–184, 1993.

Author(s):  Maher Salloum, Sandia National Laboratories, mnsallo@sandia.gov
            David Hensinger, Sandia National Laboratories
            Nathan Fabian, Sandia National Laboratories
            Jina Lee, Sandia National Laboratories
                           
The library is primarily intended for computers running the Unix, Linux,
operating system with the GNU Compiler Collection and/or Matlab. Building
and using the library requires Matlab (www.mathworks.com), the Boost library
(www.boost.org) and the Boost BLAS and LAPACK bindings to be installed on the
same operating system. The C++ library performance is limited so far to serial
computations however it can be linked to OpenBlas (http://www.openblas.net/) for
multicore computing. Future releases will be based on more sophisticated
parallel computing libraries.

The results accuracy will depend on the CPU architecture used. The tests have
shown that it is best to run the library on Xeon CPUs with 4 or more cores.

3rd party open-source software used in the library:
---------------------------------------------------

- SparseLab: (Apache 2.0 license) available at: https://sparselab.stanford.edu/ 
The Stagewise Orthogonal Matching Pursuit (StOMP) Matlab algorithm
(SolveStOMP.m) was imported and improved.

- Boost: (Boost Software License 1.0) available at: http://www.boost.org/
Boost header files where included in C++ code.

- sparse-msrf: (BSD license) available at:
http://www.sandia.gov/~jairay/kfcs/kfcs.html#Releases
The sampling matrices construction Matlab code is imported

- Alpert Transform Toolbox: (BSD license) available at:
https://github.com/gpeyre/matlab-toolboxes/tree/master/toolbox_alpert
The Alpert wavelet matrix building Matlab algorithm was imported and improved.

Contents:
--------

/data: Contains data files that are used in the Matlab and C++ examples. A
README.txt file in the folder describes each dataset. 

C++ source code:
  src/Alpert
  src/Sampling
  src/StOMP
  src/Utils

C++ example:
  /examples/compress_reconstruct: A data compression/decompression pipeline,
please refer to /examples/compress_reconstruct/README.txt for more info
    
Matlab source code:
  /matlab/src

Matlab examples:
  /matlab/examples
    
Include and Library folders: (do not delete!)
  /include
  /libraries


Compiling the library and example:
----------------------

For the library:
- Go to /src
- adjust the required settings in 'CURRENT_defines'
- run 'make'

For the example:
- Go to /examples/compress_reconstruct
- run 'make'

Notes:
- The Boost library has to be properly installed along its blas and lapack bindings.
- OpenMP and openblas also have to be properly installed in order to use the
 multicore capabilities.
- The optimal number of threads is 8. It is set by the following Linux command:
export OMP_NUM_THREADS=8 

To do in future releases:
-------------------------

- Implement the Optimized StOMP in C++
- Implement the gradient computation for unstructured meshes in Matlab and C++
- Implement other data formats e.g. HDFS
- Implement global parallel data compression


Questions and comments?
-----------------------

Please email Maher Salloum: mnsallo@sandia.gov

