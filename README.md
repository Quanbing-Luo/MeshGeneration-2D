# MeshGeneration-2D
This repository includes the source code, program, and other materials for the following  paper: 
*  Quanbing Luo; Hua Li, A mesh generator for two-dimensional triangulated geometries written in C++ and interfaced with MATLAB, Submitted to Journal. 
<!-- , [Engineering with Computers](https://doi.org/10.1007/s00366-020-01262-x), 2021 (Published Online) -->

---

There are some notes before using the MATLAB mesh generator. 

* Please open (and run) **startup.m** to add **"Functions"** folder to the MATLAB Search Path at Startup before using the MATLAB mesh generator. 

* Please read the examples in the **"Examples"** folder before using the MATLAB mesh generator to generate your own meshes.   

* The MATLAB mesh generator have successfully tested at **MATLAB R2023b** on Windows 11, the mesh generator might not run successfully at other environments. 

* The core source code are written in C++, your compiler must support **C++23** if you need to recompile the C++ source code (**meshgeneration.cpp**) for some reasons. 


In the "Functions" folder, the `meshgeneration` function of the MATLAB mesh generator is written in  **meshgeneration.m** file. The real core source code is written in C++ (**meshgeneration.cpp**), a console application (**meshgeneration.exe**) was obtained by compile the source code. Meanwhile, **geometryread.m**, **geometrywrite.m**, **meshradiiread.m**, and **meshradiiread.m** are provided to read and write binary geometry and mesh files.   

---
The meshes could be generated with `meshgeneration` function at different **name-value arguments**. For examples, 
```
meshgeneration('geometry.dat','meshradii.dat');
meshgeneration('geometry.dat','meshradii.dat',Hmax=0.03); 
meshgeneration('geometry.dat','meshradii.dat',Hvertex={4,0.02},Hgrad=1.1);
meshgeneration('geometry.dat','meshradii.dat',Hedge={4,0.01},Hgrad=1.1);
```
In the above code, `'geometry.dat'`  is the binary data file of triangulated geometry, a MATLAB `geometrywrite` function is provided simultaneously as output function of the file. 
`'meshradii.dat'` is the  binary data file of generated mesh, a MATLAB `meshradiiread` function  is provided as input function of the file. 

---
The **name-value arguments** of `meshgeneration` function were designed similar to the name-value arguments of `generateMesh` function in MATLAB Partial Differential Equation Toolbox. 
The four name-value arguments used for `meshgeneration` function are described as follows. 


*	**Hmax**:  The `Hmax` is maximum diameter of node bubbles used for mesh generation, whose length is almost at the level of required mesh edge length around the node. If the value of `Hmax` is not provided, the mesh generator would be set as a very large value. Then, the mesh size of generated mesh would be controlled by narrow regions and mesh growth rate (gradient).   
	
*	**Hvertex**: The `Hvertex` is maximum diameter of node bubbles used for mesh generation at  selected nodes/vertices. 
	
*	**Hedge**:  The `Hedge` is maximum diameter of node bubbles used for mesh generation at selected boundary edges.  
	
*	**Hgrad**: The 	`Hgrad` is mesh growth rate (gradient), mesh growth rate was specified as a number greater than 1 and less than or equal to 2. The default value of `Hgrad` is set as 1.1 (rather than 1.5 in the `generateMesh` function) for higher-quality mesh. Here, the value of `Hgrad` should not be set as 1 as it might bring some side effects. If even-distributed mesh elements are required globally, please change the value of `Hmax` preferentially rather than decrease the value of `Hgrad`. 	
	
	
**Note:** For the name-value arguments of `Hvertex` and `Hedge`, the value are specified as a cell array containing an even number of elements. Odd-indexed elements are positive integers or vectors of positive integers specifying vertex/edge IDs. Even-indexed elements are positive numbers specifying the target sizes for the corresponding vertex/edges. The syntax is the same to the  name-value arguments of `generateMesh` function in MATLAB Partial Differential Equation Toolbox. 


