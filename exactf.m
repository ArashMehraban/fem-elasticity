function exctSol = exactf(vtx_coords)
  x = vtx_coords(:,1);
  y = vtx_coords(:,2);
  
% % %   %u for Poisson problem
% % %   u=@(x,y)tanh(x).*exp(y)+sin(y); 
% % %   exctSol = u(x,y);
  
  % u1 and u2 for Plane_Strain problem
    u1=@(x,y)tanh(x).*exp(y)+sin(y);
    u2=@(x,y)tanh(x).*cos(y); 

    exctSol = [u1(x,y),u2(x,y)];
    
end