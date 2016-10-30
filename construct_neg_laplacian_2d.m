function neg_laplacian_2d = construct_neg_laplacian_2d()
%CONSTRUCT_NEG_LAPLACIAN_2D returns the negative Laplacian of u to be used
%in get_userf function
% input: none (inside this function below)
% output: the negative laplacian of u below    
    u=@(x,y)tanh(x).*exp(y)+sin(y);
    
    syms x y
    d2x = diff(u(x,y),x,2);
    d2y = diff(u(x,y),y,2);
    poss2D = simplify(-(d2x+d2y));
    neg_laplacian_2d = matlabFunction(poss2D);
end