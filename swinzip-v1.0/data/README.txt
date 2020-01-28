This is a description of the datasets used (or could be used) in the examples:

1 - mesh2D_irreg.txt
A 2D unstructured mesh consisting of 33,062 mesh points in ASCII format, no data
is present here. The first line contains the number of points and
dimensionality. The remaining lines form two columns for the x and y coordinates

2 - mesh2D_irreg_[coarse,medium,fine].txt
Similar meshes with different levels of refinement. The data format is the same
as mesh2D_irreg.txt

3 - mesh3D_irreg.txt
A 3D unstructured mesh consisting of 53,209 mesh points in ASCII format, no data
is present here. The first line contains the number of points and
dimensionality. The remaining lines form two columns for the x, y and z
coordinates

4 - mesh[2,3]D_irreg.mat
Same as mesh[2,3]D_irreg.txt but in matlab format

5 - rcci_*.dat
2D structured data (here the mesh is obvious) obtained from a RCCI turbulent
combustion simulation provided by Jackie Chen and Ankit Bhagatwala from Sandia
National Labs.

rcci_U_*.dat:     Velocity magnitude
rcci_Y_CO2_*.dat: CO2 concentration

rcci_*_0002: Data before the ignition time (smooth with very few features)
rcci_*_0003: Data after the ignition time (exhibits discontinuities with
             interesting features such as combustion fronts)

6. data_*.txt: Sample data fields that could be used in the C++ example.