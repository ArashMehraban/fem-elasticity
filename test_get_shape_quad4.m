%test cases for get_shape function: QUAD4
%test cases for Basis functions and their Derivatives, B and D's
clear
clc

%test for QUAD4 (4-noded 2D element) using Guassian Quadrature

f = @(x,y) [0.25.*(1-x).*(1-y),0.25.*(1+x).*(1-y), 0.25.*(1-x).*(1+y), 0.25.*(1+x).*(1+y)];

syms x y 

dfx = diff(f,x);
dfy = diff(f,y);

fx = matlabFunction(dfx);
fy = matlabFunction(dfy);

[gs,~ ]= get_quadrature(2);
gs = abs(gs(1));
x = gs*[-1 1 -1 1]';
y = gs*[-1 -1 1 1]';

B = f(x,y);
D0= fx(y);
D1= fy(x);

[B1,D,~]=get_shape(4);
disp('QUAD4:')
disp('  Zero matrix below indicates same the values calculated numerically and analytically are equal.')
diffB = B -B1
diffD0 = D0 - D.D0
diffD1 = D1 - D.D1

clear

%============================================================================%







