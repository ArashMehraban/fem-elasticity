%test cases for get_shape function: QUAD9
%test cases for Basis functions and their Derivatives, B and D's
clear
clc

%test for QUAD9 (9-noded 2D element) using Guassian Quadrature

f = @(x,y) [0.25.*(x.^2-x).*(y.^2-y), 0.5.*(1- x.^2).*(y.^2-y) ...
            0.25.*(x.^2+x).*(y.^2-y), 0.5.*(x.^2-x).*(1-y.^2) ...
            (1-x.^2).*(1-y.^2), 0.5.*(x.^2+x).*(1-y.^2) ...
            0.25.*(x.^2-x).*(y.^2+y),0.5.*(1-x.^2).*(y.^2+y)...
            0.25.*(x.^2+x).*(y.^2+y)];

syms x y 

dfx = diff(f,x);
dfy = diff(f,y);

fx = matlabFunction(dfx);
fy = matlabFunction(dfy);

[gs,~ ]= get_quadrature(3);
gs = abs(gs(1));
x = gs*[-1  0  1 -1  0  1 -1  0  1]';
y = gs*[-1 -1 -1  0  0  0  1  1  1]';

B = f(x,y);
D0= fx(x,y);
D1= fy(x,y);

[B1,D,~]=get_shape(9);
disp('QUAD9:')
disp('  Zero matrix below indicates same the values calculated numerically and analytically are equal.')
diffB = B -B1
diffD0 = D0 - D.D0
diffD1 = D1 - D.D1

%============================================================================%



