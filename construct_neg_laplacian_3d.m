function neg_laplacian_3d = construct_neg_laplacian_3d()
%CONSTRUCT_NEG_LAPLACIAN_3D returns the negative Laplacian of u to be used
%in get_userf function
% input: none (inside this function below)
% output: the negative laplacian of u below 
    u=@(x,y,z)tanh(x).*exp(y)+sin(y)+cos(z);
    
    syms x y z
    d2x = diff(u(x,y,z),x,2);
    d2y = diff(u(x,y,z),y,2);
    d2z = diff(u(x,y,z),z,2);
    poss3D = simplify(-(d2x+d2y+d2z));
    neg_laplacian_3d = matlabFunction(poss3D);
end