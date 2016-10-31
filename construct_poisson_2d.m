function g = construct_poisson_2d()
%CONSTRUCT_POISSON_2D manufactures the Poisson PDE for a given u:
%
%  -grad^2(u) = g . So it is written as :
%
%      \partial^2_u        \partial^2_u
%  -(---------------- +  ------------------)  = g  
%      \partial^2_x        \partial^2_y  
%
% where u if a function of x and y
%
% g should be used in userf_poisson_2d
    
    % Assume u is the solution:
    u=@(x,y)tanh(x).*exp(y)+sin(y);
    
    % Manufacture the Poisson PDE:
    syms x y
    d2x = diff(u(x,y),x,2);
    d2y = diff(u(x,y),y,2);
    g_4_pde = simplify(-(d2x+d2y));
    g = matlabFunction(g_4_pde);
end