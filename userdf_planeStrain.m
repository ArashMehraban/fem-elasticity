function f = userdf_planeStrain(dlta_ue, grad_dlta_ue, xe)
    
   % grad_dlta_ue structure:
   % 2D:             [D0*dlta_u1 | D0*dlta_u2] 
   %    grad_dlta_ue=[D1*dlta_u1 | D1*dlta_u2]
   
    % Young's modulus 
    E = 1;%2e11;
    % Poisson ratio 
    nu = 0.3;
        
    num_row = size(dlta_ue,1);
    
    %strain-stress matrix coefficient
    ss_coef = E/((1+nu)*(1-2*nu));
    
    %strain
    dlta_eps11 = grad_dlta_ue(1:num_row,1);
    dlta_eps12 = 0.5*(grad_dlta_ue(num_row+1:end,1)+grad_dlta_ue(1:num_row,2));
    dlta_eps22 = grad_dlta_ue(num_row+1:end,2);
    
    dlta_sigma11 = ss_coef*((1-nu)*dlta_eps11 + nu*dlta_eps22);
    dlta_sigma22 = ss_coef*(nu*dlta_eps11 +(1-nu)*dlta_eps22);
    dlta_sigma12 = ss_coef*(0.5*(1-2*nu)*dlta_eps12);    
      

    fu1 = 0*dlta_ue(:,1);
    fu2 = 0*dlta_ue(:,2);
    
    f = [fu1 fu2; 
         dlta_sigma11 dlta_sigma12;
         dlta_sigma12 dlta_sigma22 ];
end