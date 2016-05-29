function neg_laplacian_2d = neg_laplacef_2d(u)
% input: givenF: user defined function
% output: the laplacian of givenF
    syms x y
    d2x = diff(u(x,y),x,2);
    d2y = diff(u(x,y),y,2);
    poss2D = simplify(-(d2x+d2y));
    neg_laplacian_2d = matlabFunction(poss2D);

end