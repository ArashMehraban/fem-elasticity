function f = userdf_3d_elas(dlta_ue, grad_dlta_ue, xe) 
%USERF_3d_ELAS provides weak form of the linear 3D Elastisity problem to solve 
%
%  input:      ue: corresponding u for each element evalauted at quadrature points
%       : grad_ue: corresponding grad_u for each element evalauted at quadrature points
%       :      xe: quadrature points mapped to the reference elem
%
% output:  f0: any possible source from given problem
%       :  f1: any possible source from given problem



    x=xe(:,1);
    y=xe(:,2);
    z=xe(:,3);
    
   %grad_ue  structure:
   % 3D:          [D1*u1 | D1*u2 | D1*u3]   [du1/dx | du2/dx | du3/dx] 
   %    grad_ue = [D2*u1 | D2*u2 | D2*u3] = [du1/dx | du2/dy | du3/dy]
   %              [D3*u1 | D3*u2 | D3*u3]   [du1/dz | du2/dz | du3/dz] 
    
    
    % Young's modulus 
    E = 1;%2e11;
    % Poisson ratio 
    nu = 0.3;
        
    num_row = size(dlta_ue,1);
    
    %strain-stress matrix coefficient
    ss_coef = E/((1+nu)*(1-2*nu));
    
    %strain
    eps11 = grad_dlta_ue(1:num_row,1);
    eps22 = grad_dlta_ue(num_row+1:2*num_row,2);
    eps33 = grad_dlta_ue(2*num_row+1:end,3);
    eps12 = 0.5*(grad_dlta_ue(1:num_row,2)+grad_dlta_ue(num_row+1:2*num_row,1));
    eps13 = 0.5*(grad_dlta_ue(1:num_row,3)+grad_dlta_ue(2*num_row+1:end,1));
    eps23 = 0.5*(grad_dlta_ue(num_row+1:2*num_row,3)+grad_dlta_ue(2*num_row+1:end,2));
    
    sigma11 = ss_coef*((1-nu)*eps11 + nu*eps22 + nu*eps33);
    sigma22 = ss_coef*(nu*eps11 + (1-nu)*eps22 + nu*eps33);
    sigma33 = ss_coef*(nu*eps11 + nu*eps22 + (1-nu)*eps33);
    sigma12 = ss_coef*(0.5*(1-2*nu)*eps12);    
    sigma13 = ss_coef*(0.5*(1-2*nu)*eps13);  
    sigma23 = ss_coef*(0.5*(1-2*nu)*eps23);
   
    fu1 = 0*dlta_ue(:,1);
    fu2 = 0*dlta_ue(:,2);
    fu3 = 0*dlta_ue(:,3);
    
    f = [  fu1      fu2      fu3;
         sigma11  sigma12  sigma13;
         sigma12  sigma22  sigma23;
         sigma13  sigma23  sigma33];
    
end