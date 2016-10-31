## Instructions on how to run the code in fem-elasticity repository:

**IMPORTANT NOTE:** 
Mesh files have been provided in separate folders in this repository. Make sure
all mesh folders are removed from MATLAB path before running code. Different mesh
files have the same names in different mesh folders. Due to MATLAB's strange path system, 
wrong mesh files may be imported by the code if all mesh folders are on MATLAB path 
before the code is run! The example code (e.g. runL2projection.m, runPoisson.m, etc) 
show how to add and remove the necessary mesh folder for those problems.

**For all any problem:**

  1. Clone the entire fem-elastisity repository
  2. Place the fem-elastisity folder **only** on MATLAB path. 
  3. Make sure the **mesh** and **test_units** folders are removed from Path (see above 
   **Important Note**)

**Runnning L2 Projection problem:**

  1. Run _runL2projection.m_  
  2. MATLAB should generate a plot that looks like _jpg_L2proj.jpg_

  ![alt text][L2]

   [L2]: https://github.com/ArashMehraban/fem-elasticity/blob/master/jpg_L2proj.jpg "L2 projection"

 **Note:** The _u_ (given_u in _runL2projection.m_) is the solutions of L2 projection problem in 2D
       and 3D. So, if _u_ is changed, then _userf_L2_2d.m_ or _userf_L2_3d.m_ must be modified accordingly.


**Runnning Poisson problem with Dirichlet boundary condition only:**

   1. Run _runPoisson.m_  
   2. MATLAB should generate a plot that looks like _jpg_poisson.jpg_
  
  ![alt text][Poisson]

   [Poisson]: https://github.com/ArashMehraban/fem-elasticity/blob/master/jpg_poisson.jpg "Poisson with Dirichlet Boundary conditions only"

 Note 1: The _u_'s in _runPoisson.m_ are the solutions of Poisson PDE in 2D and 3D that 
         have been constructed using _construct_poisson_2d.m_ and 
         _construct_poisson_3d.m_ files. The returns by these script files
         must be used _userf_poisson_2d.m_ and _userf_poisson_3d.m_ files. So, if _u_ 
         is changed, then the g variable(s) in  _userf_poisson_2d.m_ and 
         _userf_poisson_3d.m_ must be modified accordingly.

Note 2:  _userf_poisson_2d.m_ and _userf_poisson_3d.m_ are the physics of the problem 
         that the user provides. User must also provide the constituents for the
         **consistent tangent** of each problem. That is perfomed by modifying 
         _user**d**f_poisson_2d.m_ and _user**d**f_poisson_3d.m_
         Note that, the code in _userdf_poisson_2d.m_ and _userdf_poisson_3d.m_ are 
         almost identical to _userf_poisson_2d.m_ and _userf_poisson_3d.m_ as Poisson 
         is a linear problem.


**Runnning Plane Strain problem with Dirichlet boundary condition only:**

   1. Run _runPlaneStrain.m_  
   2. MATLAB should generate a plot that looks like _jpg_poisson.jpg_
  
  ![alt text][Plane Strain]

   [Plane Strain]: https://github.com/ArashMehraban/fem-elasticity/blob/master/jpg_planeStrain.jpg "Plane Straine with Dirichlet Boundary conditions only"

 Note 1: The _u_'s in _runPlaneStrain.m_ are the solutions of Plane Strain PDE in 2D that 
         have been constructed using _construct_plane_strain.m_ file. The returns by this 
         script file must be used in _userf_planeStrain.m_ file. So, if _u_ is changed, 
         then the g variable(s) in  _userf_planeStrain.m_ must be modified accordingly.

 Note 2: _userf_planeStrain.m_  contains the physics of the problem that the user provides. 
         User must also provide the constituents for the **consistent tangent** of the
         problem. That is perfomed by modifying _user**d**f_planeStrain.m_. 
         Note that, the code in _userf_planeStrain.m_ is almost identical to _userdf_planeStrain.m_  
         as Plane Strane is a linear problem.

**Runnning 3D Elastisity problem with Dirichlet boundary condition only:**

   1. Run _run3dElas.m_  
   2. MATLAB should generate a plot that looks like _jpg_3dElas.jpg_
  
  ![alt text][3D Elastisity]

   [3D Elastisity]: https://github.com/ArashMehraban/fem-elasticity/blob/master/jpg_3dElas.jpg "3D Elastisity with Dirichlet Boundary conditions only"

 Note 1: The _u_'s in _run3dElas.m_ are the solutions of 3D Elastisity PDE in 3D that 
         have been constructed using _construct_linear_3d_elas.m_ file. The returns by this 
         script file must be used in _userf_3d_elas.m_ file. So, if _u_ is changed, 
         then the g variable(s) in  _userf_3d_elas.m_ must be modified accordingly.

 Note 2: _userf_3d_elas.m_  contains the physics of the problem that the user provides. 
         User must also provide the constituents for the **consistent tangent** of the
         problem. That is perfomed by modifying _user**d**f_3d_elas.m_. 
         Note that, the code in _userdf_3d_elas.m_ is almost identical to _userf_3d_elas.m_  
         as 3D Elastisity is a linear problem.


