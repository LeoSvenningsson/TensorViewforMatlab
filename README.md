# TensorView for Matlab

TensorView for Matlab is a GUI and script based tool to visualize tensors in a molecular context. TensorView for Matlab includes: Reading arbitrary .pdb and .xyz files for molecular visualization. Ovaloid and ellipsoid tensor visualisation. 3D model file conversion to .glb and .wrl. Euler angle and relative angle decoder. All angles are cross-validated with the input tensor.
   
   TensorView for Matlab is licenced with creative commons CC BY. https://creativecommons.org/licenses/
Free to share and adapt. Give appropriate credits to authors:

Please reference article:
https://doi.org/10.1016/j.ssnmr.2022.101849

Credits:

github.com/LeoSvenningsson/TensorViewforMatlab

mathworks.com/matlabcentral/fileexchange/55231-molecoule3d

onlinelibrary.wiley.com/doi/full/10.1002/mrc.4793

https://se.mathworks.com/matlabcentral/fileexchange/109264-matlab2glb

Version: 1.152 (June 1st 2024)

Authors: 
Dr. Leo Svenningsson () 
Prof. Leonard Mueller (leonard.mueller@ucr.edu)
Contributors:
Dr. André Ludwig (aludwig@alumni.ethz.ch, molecule3d)
Dmitri Sastin (https://se.mathworks.com/matlabcentral/fileexchange/109264-matlab2glb, matlab2glb)

------ 
   TensorView for Matlab is a collaboration with works derrived from molecule3D.m (André Ludwig: mathworks.com/matlabcentral/fileexchange/55231-molecule3d), matlab2glb ( Dmitri Sastin: https://se.mathworks.com/matlabcentral/fileexchange/109264-matlab2glb), and TensorView (Prof. Leonard Mueller: https://onlinelibrary.wiley.com/doi/full/10.1002/mrc.4793 and https://sites.google.com/ucr.edu/Mueller/home/tensorview)


Versions for windows and mac reside in their respective folder which only requires the latest matab compiler. The Scource folder contains the matlab code, but also include a matlab script which can be run with matlab. As of now, there are more camera options for opening the GUI through matlab. Matlab2022 or the free "matlab2022 runtime" or later versions are required.

In matlab, set your working folder to the "source" folder and klick on TensorView.mlapp (in the same folder) to start the app.

.glb files are easely opened with the free Blender software or (windows) 3D paint.
With the inclusion of .wrl support, a short tutorial for Blender is included below. 
1: press file and import and choose the .wrl format.<br/>
2: Blender has a few defoult workspaces, select the Layout workspace from the top tab if not already selected.<br/>
3: select all objects either in the viewer or from the side panel.<br/>
4: click one additional time on one of the objects. in the bottom right, select the materials properties tab. Click on base color and select vertex colors. Now the color has been added to the selected object. <br/>
5: Click control+l and choose materials. This will link the materials between all of the selected objects.<br/>
6: press on the globes in the viewport top right corner to see the materials and shading in the viewport.<br/>
7: in the top left corner, select object and choose "shade smooth" to make it look pretty. <br/>

### Wanted list
Wigner rotations

### Known issues
the windows and mac apps are only compiled to patch 1.151
mac app is currently only compiled for apple silicon
### Changelog 1.152
A bug was found where an erronius criteria if acos(U(3,3)) == 1 was found in three places, and changed to if U(3,3) == 1. These are used for simplifications of the specific angle beta=0.
It is unlikly that these angles simplifications where ever erroniously activated in since U(3,3) would have needed to be a very specific irrational number. 

### Changelog 1.151
Recompiled Tensorview for matlab without any changes. This somehow fixed a bug where the precompiled versions would not display figures.

### Changelog 1.15 (Major update)
New features: Export 3D models as .glb format. Module for relative angles of two tensors. Module for reference frame plotting. Feature to load and save lists of tensors as text files (.dat). Feature to load lists of reference frames. Feature to save lists of relative angles from two tensors. Removed all iterative calculations of euler angles for all analytical solutions though the function for sqmineuler still remains in the function folder if someone wants to use it in their own scripts. The app also now includes an extra validation step that will give notification in the matlab command promt if for any reason the obtained MFtoEuler or relative euler angles should not rectreat the correct rotation matrix/tensor. The new cross validation feature is an insurance that the angles are calculated correctly and will give a notification in the command window if the  cross validation could not exactly reproduce the input tensor. This is mainly to find all of the edge cases, usually from floating point jittering. 


### Changelog 1.14
Tensorview for Matlab now supports .wrl export of the 3d model. The "Exportsurfwrl.m" script will export any surface figure. Try fig = figure; surf(X,Y,Z) for a surface defined by X, Y and Z.  Exportsurfwrl(fig,filename,path). The .wrl format is best used in combination with "Blender" https://www.blender.org/. The .wrl format supports primitives, however this is not yet implemented, which means that the 3d files are larger than its most compressed form... for now.

### Changelog 1.13
Semicomplete conections lists can now be used.
Updated the .pdb interpreter. New version should have better compatability with different .pdb files.
Bugfix on atom names longer then one letter.

### Changelog 1.12

Tensoview now simplifies symmetric tensors in the Unique-x PAS form. Example: [1 2 2]. 
Unique-x ordering lack rotations within the tensor symmetry, therefore, only this system is solved nummerically with a least square method.

A few variable name changes within the code.

### Changelog 1.11
Fixed bug with some forms of atom labels. 

Added atomic number read for .xyz

Updated contact information

Fixed bug with color slider

Tensoview now simplifies symmetric tensors in the Unique-z PAS form. Example: [2 2 1]. 

Recompiled for MacOS "Big Sur"
