# TensorView for Matlab
TensorView for Matlab is a tool to visualize chemical shift/shielding tensors in a molecular context. TensorView for Matlab can read arbitrary .pdb and .xyz files for molecular visualization. Any 3D tensor can be used for visualisation though the chemical shift/shielding tensor is used as an input in script.  
   
   TensorView for Matlab is licenced with creative commons CC BY. https://creativecommons.org/licenses/
Free to share and adapt. Give appropriate credits to authors:

Credits:

github.com/LeoSvenningsson/TensorViewforMatlab

mathworks.com/matlabcentral/fileexchange/55231-molecoule3d

onlinelibrary.wiley.com/doi/full/10.1002/mrc.4793

   Version: 1.12

Authors: 
Dr. Leo Svenningsson (leo.svenningsson@chalmers.se), 
Dr. André Ludwig (aludwig@alumni.ethz.ch),
Prof. Leonard Mueller (leonard.mueller@ucr.edu)

   TensorView for Matlab is a collaboration with works derrived from molecule3D.m (André Ludwig: mathworks.com/matlabcentral/fileexchange/55231-molecule3d) and TensorView (Prof. Leonard Mueller: https://onlinelibrary.wiley.com/doi/full/10.1002/mrc.4793 and https://sites.google.com/ucr.edu/Mueller/home/tensorview)


Versions for windows and mac reside in their respective folder which only requires the latest matab compiler. The Scource folder contains the matlab code, but also include a matlab script which can be run with matlab. As of now, there are more camera options for the script version using matlab, which is why it is included separately.

### Changelog 1.12

Tensoview now simplifies symmetric tensors in the Unique-x PAS form. Example: [1 2 2]. 
Unique-x ordering lack lotations within the tensor symmetry, therefore only this system is solved nummerically with a least square method.

A few variable name changes within the code.

### Changelog 1.11
Fixed bug with some forms of atom labels. 

Added atomic number read for .xyz

Updated contact information

Fixed bug with color slider

Tensoview now simplifies symmetric tensors in the Unique-z PAS form. Example: [2 2 1]. 

Recompiled for MacOS "Big Sur"
