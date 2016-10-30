####Instructions on how to run the code in fem-elasticity repository:


**IMPORTANT NOTE:** 
Mesh files have been provided in separate folders in this repository. Make sure
all mesh folders are removed from MATLAB path before running code. The code will 
add and remove the necessary mesh folder automatically. Different mesh files have
the same names in different mesh folders. Due to MATLAB's strange path system, 
wrong mesh files may be imported by the code if all mesh folders are on MATLAB path 
before the code is run!

**Runnning L2 Projection problem:**
1. Run _runL2projection.m_
2. MATLAB should generate a plot that looks like _jpg_L2proj.jpg_

  ![alt text][L2]

   [L2]: https://github.com/ArashMehraban/fem-elasticity/blob/master/jpg_L2proj.jpg "L2 projection"

 **Note:** The _u_'s in _runL2projection.m_ are the solutions of L2 projection problem in 2D
       and 3D. So, if _u_ is changed, then userf.m_ must be modified accordingly.

