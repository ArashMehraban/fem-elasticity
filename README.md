Instructions on how to run the code in fem-elasticity repository:

I) Runnning L2 Projection problem:(Under Construction)
   1) Make sure all code (.m files) are in one folder in MATLAB path
   2) Make sure all the filenames in variable files in runL2Projection.m 
      are in the same folder as the code
   3) Unzip the folder exo for all exodus files that need to be used with 
      runL2Projection.m. 
   4) Place all .exo files in the same folder that contains the code 
   5) run runL2Projection.m 
 Note: runL2Projection.m works with the u provided in runL2Projection.m. If 
       u is changed, the function get_userf must be modified accordingly. 
    
 How to change get_userf for with a different u?
 Answer soon!


II) Runnning Poisson problem:
   1) Make sure all code (.m files) are in one folder in MATLAB path
   2) Make sure all the filenames in variable files in runPoisson.m 
      are in the same folder as the code
   3) Unzip the folder exo for all exodus files that need to be used with 
      runL2Projection.m. 
   4) Place all .exo files in the same folder that contains the code 
   5) run runPoisson.m 
 Note: runPoisson.m works with the u provided in runPoisson.m. If 
       u is changed, the function get_userf must be modified accordingly.
       construct_neg_laplacian_2D or construct_neg_laplacian_3D can be used
       for manufactured solutions.

  How to change get_userf for with a different u?
  Answer soon!


III) Runnning Plane Strain problem:
   Under Construction


