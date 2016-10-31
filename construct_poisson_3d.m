function g = construct_poisson_3d()
%CONSTRUCT_POISSON_3D manufactures the Poisson PDE for a given u:
%
%  -grad^2(u) = g . So it is written as :
%
%      \partial^2_u        \partial^2_u        \partial^2_u
%  -(---------------- +  ------------------ + -------------- )  = g  
%      \partial^2_x        \partial^2_y        \partial^2_z
%
% where u if a function of x, y and z
%
% g should be used in userf_poisson_3d
     
    %Assume u is sthe solution:
    u=@(x,y,z)tanh(x).*exp(y)+sin(y)+cos(z);
    
    %Manufacture the 3D Poisson PDE:
    syms x y z
    d2x = diff(u(x,y,z),x,2);
    d2y = diff(u(x,y,z),y,2);
    d2z = diff(u(x,y,z),z,2);
    g_4_pde = simplify(-(d2x+d2y+d2z));
    g = matlabFunction(g_4_pde);
end









