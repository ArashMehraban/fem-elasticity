function exctSol = exactf(vtx_coords)
  x = vtx_coords(:,1);
  y = vtx_coords(:,2);
  
  % u for Poisson problem
  %u=@(x,y)tanh(x).*exp(y)+sin(y); 
  
  % u1 and u2 for Plane_Strain problem
    u1=@(x,y)tanh(x).*exp(y)+sin(y);
    u2=@(x,y)tanh(x).*cos(y); 
    
    exctSol{1} = u1(x,y);
    exctSol{2} = u2(x,y);
  
  
%   u1=@(x,y)tanh(x).*exp(y)+exp(x);
%   u2=@(x,y)tanh(x).*exp(y);
%   u3=@(x,y)exp(x);
%   u4=@(x,y)tanh(x).*exp(y)+cos(y);

%   exctSol= [u(x,y),u1(x,y),u2(x,y),u3(x,y),u4(x,y)];
end