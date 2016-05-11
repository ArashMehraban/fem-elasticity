function neg_laplacian = neg_laplacef(u)
% input: givenF: user defined function
% output: the laplacian of givenF
    syms x y
    d2x = diff(u(x,y),x,2);
    d2y = diff(u(x,y),y,2);
    poss = simplify(-(d2x+d2y));
    neg_laplacian = matlabFunction(poss);

end