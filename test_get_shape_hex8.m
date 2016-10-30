%test cases for get_shape function: Hex
%test cases for Basis functions and their Derivatives, B and D's
clear
clc

%test for QUAD9 (8-noded 3D element) using Guassian Quadrature

f = @(x,y,z) [0.1250.*(1-x).*(1-y).*(1-z),0.1250.*(1+x).*(1-y).*(1-z),...
            0.1250.*(1-x).*(1+y).*(1-z),0.1250.*(1+x).*(1+y).*(1-z),...
            0.1250.*(1-x).*(1-y).*(1+z),0.1250.*(1+x).*(1-y).*(1+z),...
            0.1250.*(1-x).*(1+y).*(1+z),0.1250.*(1+x).*(1+y).*(1+z)];

syms x y z

dfx = diff(f,x);
dfy = diff(f,y);
dfz = diff(f,z);

fx = matlabFunction(dfx);
fy = matlabFunction(dfy);
fz = matlabFunction(dfz);

[gs,~ ]= get_quadrature(2);
gs = abs(gs(1));
x = gs*[-1 1 -1 1 -1 1 -1 1]';
y = gs*[-1 -1 1 1 -1 -1 1 1]';
z = gs*[-1 -1 -1 -1 1 1 1 1]';

B = f(x,y,z);
D0= fx(y,z);
D1= fy(x,z);
D2= fz(x,y);

[B1,D,~]=get_shape(8);
disp('HEX8:')
disp('  Zero matrix below indicates same the values calculated numerically and analytically are equal.')
diffB = B -B1
diffD0 = D0 - D.D0
diffD1 = D1 - D.D1
diffD2 = D2 - D.D2

%============================================================================%



