# SWinzip
Compression of data represented on unstructured meshes and point-clouds using
Alpert wavelets and Compressed Sensing

SWinzip v1.0 (Sensing Wavelets in-situ zip): is a Matlab and C++ library for
scientific lossy data compression and reconstruction using compressed sensing
and tree-wavelets transforms. These methods are known for their large
compression and usefulness in data analytics such as features extraction.
Compressed sensing and wavelets methods rely heavily on sparse and dense linear
algebra operations implemented through the Boost codes in our library. SWinzip
accommodates data represented on both regular grids (e.g. image data) and
point-clouds (e.g. unstructured meshes).




SWinzip v2.0 (Sensing Wavelets in-situ zip): In SWinzip v2.0, there is more
focus on wavelet compression and the following features are added:

- The matrix-free (direct) wavelet transform is implemented in both C++ and
Matlab. A threshold-based compression can now be used. Examples for distributed
compression are also added.
- Python implementation of Alpert wavelet transforms
- Field data comparison using wavelets in Matlab with global and local metrics. 

