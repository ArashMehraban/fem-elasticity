function [g1,g2] = construct_plane_strain()
%CONSTRUCT_PLAIN_STRAIN manufactures the Plain Strain PDE for a given 
%right hand side u:
%
%   div.(sigma) + g = 0  where  g and sigma are vectors. So it is written as : 
%
%   \partial_sigma_11    \partial_sigma_12
%    ---------------- +  ------------------  + g1 = 0
%      \partial_x           \partial_y  
%
%   \partial_sigma_22    \partial_sigma_21
%    ---------------- +  ------------------  + g2 = 0
%      \partial_y           \partial_x    
%
% Therefore,
%
%           \partial_sigma_11    \partial_sigma_12
%  g1 = - ( ---------------- +  ------------------- )
%              \partial_x           \partial_y  
%
%
%            \partial_sigma_22    \partial_sigma_21
%  g2 = - (  ---------------- +  ------------------ )
%              \partial_y           \partial_x    
%
% g1 and g2 should be used in userf_plainStrain
   
   % Assume u is the solution: (u is displpacement)   
   u{1}=@(x,y)0.5.*exp(x).*sin(y) + 3.*y.^2;
   u{2}=@(x,y)x.^2 - 0.75.*exp(y).*sin(x);
   
   % Manufacture the Plain Strain PDE:
   syms x y nu E 
   
   C =(E/((1+nu)*(1-2*nu)))*[1-nu, nu, 0; nu,1-nu,0 ; 0 ,0 ,0.5*(1-2*nu)];   
   strain = [diff(u{1},x) ; diff(u{2},y) ; 0.5*(diff(u{2},x) + diff(u{1},y))];
   
   sigma = C * strain;
   
   g1_4_pde = -(diff(sigma(1),x) + diff(sigma(3),y)); 
   g2_4_pde = -(diff(sigma(3),x) + diff(sigma(2),y));
    
   g1=matlabFunction(g1_4_pde);
   g2=matlabFunction(g2_4_pde);
end