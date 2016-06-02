Instructions on how to run the code in fem-elasticity repository:

I) Runnning L2 Projection problem:(Under Construction)

  coming soon!



II) Runnning Poisson problem with Dirichlet boundary condition only:

   1) Save all functions (.m files) starting with "get_" along with runPoisson.m
      in a folder in MATLAB path

   2) Save "mesh_Dirich_only" folder in the same folder as above in MATLAB path. 
      Keep this in its own folder. The contents of mesh folder must be separate 
      from .m files.

   3) Run runPoisson.m 

   4) MATLAB should generate a plot that looks like jpg_runPoisson_Dirich.jpg 

 Note: The u's in runPoisson.m are the solutions of Poisson PDE in 2d and 3D that 
       have been constructed using construct_neg_laplacian2d.m and 
       construct_neg_laplacian3d.m ascript files. The rerurns by these script files
       have been provided in get_user_f.m. So, if u is changed, then get_userf.m 
       must be modified accordingly.

  How to change get_userf with a different u?
  Answer coming soon!



III) Runnning Plane Strain problem:
   Under Construction


