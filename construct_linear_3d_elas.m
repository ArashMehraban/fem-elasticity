function [g1,g2,g3] = construct_linear_3d_elas()
%CONSTRUCT_LINEAR_3D_ELAS manufactures the Linear 3D Elastisity PDE for a 
% give right hand side u.
%
%  div.(sigma) + g = 0  where  g and sigma are vectors. So it is written as : 
%
%   \partial_sigma_11    \partial_sigma_12      \partial_sigma_13
%    ---------------- +  ------------------  +  ------------------  + g1 = 0
%      \partial_x           \partial_y             \partial_z
%
%   \partial_sigma_22    \partial_sigma_21      \partial_sigma_23
%    ---------------- +  ------------------  + ------------------ +  g2 = 0
%      \partial_y           \partial_x             \partial_z
%
%   \partial_sigma_33    \partial_sigma_31       \partial_sigma_32
%    ---------------- +  ------------------  +  ------------------ + g3 = 0
%      \partial_z           \partial_x              \partial_z
%
% Therefore,
%
%           \partial_sigma_11    \partial_sigma_12      \partial_sigma_13
%  g1 = -( ---------------- +  ------------------  +  ------------------ )
%               \partial_x           \partial_y             \partial_z
%
%          \partial_sigma_22    \partial_sigma_21      \partial_sigma_23
%  g2 = -( ---------------- +  ------------------  +   ----------------- )
%               \partial_y           \partial_x             \partial_z
%
%          \partial_sigma_33    \partial_sigma_31       \partial_sigma_32 
%  g3 = -( ---------------- +  ------------------  +   ----------------- )
%            \partial_z           \partial_x              \partial_z
%
%  
%
% g1, g2 and g3 should be used in userf_3d_elas 


   
   % Assume u is the solution: (u is displpacement)
   u{1}=@(x,y,z)exp(2*x).*sin(3*y).*cos(4*z);
   u{2}=@(x,y,z)exp(3*y).*sin(4*z).*cos(2*x); 
   u{3}=@(x,y,z)exp(4*z).*sin(2*x).*cos(3*y);
   
   % Manufacture the 3D Elastisity PDE:
   syms x y z nu E 
   
   C =(E/((1+nu)*(1-2*nu)))*[1-nu,    nu,   nu,            0,  0,            0; 
                               nu,  1-nu,   nu,            0,  0,            0;
                               nu,    nu, 1-nu,            0,  0,            0;
                               0,      0,    0, 0.5*(1-2*nu),  0,            0;
                               0,      0,    0,            0,  0.5*(1-2*nu), 0;
                               0,      0,    0,            0,  0,            0.5*(1-2*nu)];   
   strain = [diff(u{1},x) ; diff(u{2},y) ; diff(u{3},z) ;
             0.5*(diff(u{1},y) + diff(u{2},x));
             0.5*(diff(u{1},z) + diff(u{3},x));
             0.5*(diff(u{2},z) + diff(u{3},y))];
   
   %  [sigma11 = 1]
   %  [sigma22 = 2]
   %  [sigma33 = 3]
   %  [sigma12 = 4]
   %  [sigma13 = 5]
   %  [sigma23 = 6]
   sigma = C * strain;
   
   g1_4_pde = -(diff(sigma(1),x) + diff(sigma(4),y) + diff(sigma(5),z)); 
   g2_4_pde = -(diff(sigma(2),y) + diff(sigma(4),x) + diff(sigma(6),z));
   g3_4_pde = -(diff(sigma(3),z) + diff(sigma(5),x) + diff(sigma(6),y));
    
   g1=matlabFunction(g1_4_pde);
   g2=matlabFunction(g2_4_pde);
   g3=matlabFunction(g3_4_pde);
end