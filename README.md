Instructions on how to run the code in fem-elasticity repository:

IMPORTANT NOTE: 
Mesh files have been provided in separate folders in this repository. Make sure
you do NOT add them to MATLAB path manually. The code will add and remove the necessary
mesh folder automatically. Due to MATLAB 's strange path system, wrong mesh files may
be imported by the code. Make sure all mesh folders have been removed from path before 
running code.

I) Runnning L2 Projection problem:

  1) Save all functions (.m files) starting with "get_" along with runL2projection.m
      in a folder in MATLAB path

  2) Save "mesh_no_ns_ss" folder in the same folder as above, but do NOT add it to
     MATLAB path (see above IMPORTANT NOTE). The contents of mesh folder must be 
     separate from .m files.

  3) Run runL2projection.m 

  4) MATLAB should generate a plot that looks like jpg_L2proj.jpg 

 Note: The u's in runL2projection.m are the solutions of L2 projection problem in 2D
       and 3D. So, if u is changed, then get_userf.m must be modified accordingly.

  How to change get_userf with a different u?
 
  i) Make sure the part of the code in get_userf.m corresponding to L2 projection problem 
     is uncommented while everything else, EXCEPT where in the code says do not comment/remove,
     is commented. The coressponding part is in between

   %==== L2 Starts====%

   %==== L2 ends======% 

   ii) Look for variable g under 2D or 3D sections for L2 Projection in get_userf.m 
       and place the new u for your 2D or 3D problem.

++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

II) Runnning Poisson problem with Dirichlet boundary condition only:

   1) Save all functions (.m files) starting with "get_" along with runPoisson.m
      in a folder in MATLAB path

   2) Save "mesh_Dirich_only" folder in the same folder as above ibut do NOT add it to
     MATLAB path (see above IMPORTANT NOTE). The contents of mesh folder must be 
     separate from .m files.

   3) Run runPoisson.m 

   4) MATLAB should generate a plot that looks like jpg_runPoisson_Dirich.jpg 

 Note: The u's in runPoisson.m are the solutions of Poisson PDE in 2d and 3D that 
       have been constructed using construct_neg_laplacian2d.m and 
       construct_neg_laplacian3d.m script files. The rerurns by these script files
       have been provided in get_user_f.m file. So, if u is changed, then get_userf.m 
       must be modified accordingly.

  How to change get_userf with a different u?
  Answer coming soon!

++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

III) Runnning Plane Strain problem:
   Under Construction


