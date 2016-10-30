#### Instructions on how to run the code in fem-elasticity repository:

**IMPORTANT NOTE:** 
Mesh files have been provided in separate folders in this repository. Make sure
all mesh folders are removed from MATLAB path before running code. The code will 
add and remove the necessary mesh folder automatically. Different mesh files have
the same names in different mesh folders. Due to MATLAB's strange path system, 
wrong mesh files may be imported by the code if all mesh folders are on MATLAB path 
before the code is run!

**Runnning L2 Projection problem:**

  1. Save all functions (.m files) starting with **get_** along with **runL2projection.m**
      in a folder in MATLAB path
  2. Save "_mesh_no_ns_ss_" folder in the same folder as above, but do NOT add it to
     MATLAB path (see above IMPORTANT NOTE). The contents of mesh folder must be 
     in a separate folder from .m files.
  3. Run _runL2projection.m_ 
  4. MATLAB should generate a plot that looks like _jpg_L2proj.jpg_

  ![alt text][L2]

   [L2]: https://github.com/ArashMehraban/fem-elasticity/blob/master/jpg_L2proj.jpg "L2 projection"

 **Note:** The _u_'s in _runL2projection.m_ are the solutions of L2 projection problem in 2D
       and 3D. So, if _u_ is changed, then _get_userf.m_ must be modified accordingly.

  How to modify _get_userf_ with a different _u_?
 
  + Make sure the part of the code in _get_userf.m_ corresponding to L2 projection problem 
     is uncommented while everything else is commented, EXCEPT where in the code says 
     "do not comment/remove". The corresponding part to L2 in _get_userf_ is placed in between:

   %==== L2 Starts====%

   %==== L2 ends======% 

  + Once the corresponding part of the code to L2 projection in _get_userf_ is uncommented,
       look for variable _**g**_ under 2D and/or 3D sections for L2 Projection in _get_userf.m_ 
       and place the new _u_ for your 2D and/or 3D problem.

======

**Runnning Poisson problem with Dirichlet boundary condition only:**

   1. Save all functions (.m files) starting with **get_** along with **runPoisson.m**
      in a folder in MATLAB path
   2. Save "_mesh_Dirich_only_" folder in the same folder as above but do NOT add it to
     MATLAB path (see above IMPORTANT NOTE). The contents of mesh folder must be 
     separate from .m files.
   3. Run _runPoisson.m_ 
   4. MATLAB should generate a plot that looks like _jpg_runPoisson_Dirich.jpg_
  
  ![alt text][Poisson]

   [Poisson]: https://github.com/ArashMehraban/fem-elasticity/blob/master/jpg_runPoisson_Dirich.jpg "Poisson with Dirichlet Boundary conditions only"

 Note: The _u_'s in runPoisson.m are the solutions of Poisson PDE in 2d and 3D that 
       have been constructed using _construct_neg_laplacian2d.m_ and 
       _construct_neg_laplacian3d.m_ script files. The returns by these script files
       have been provided in _get_user_f.m_ file. So, if _u_ is changed, then _get_userf.m_ 
       must be modified accordingly.

  How to modify _get_userf_ with a different _u_?
  Answer coming soon!

======

**Runnning Plane Strain problem:**
   Under Construction


