function neg_laplacian_3d = neg_laplacef_3d(u)
% input: givenF: user defined function
% output: the laplacian of givenF
    syms x y z
    d2x = diff(u(x,y,z),x,2);
    d2y = diff(u(x,y,z),y,2);
    d2z = diff(u(x,y,z),z,2);
    poss3D = simplify(-(d2x+d2y+d2z));
    neg_laplacian_3d = matlabFunction(poss3D);

end