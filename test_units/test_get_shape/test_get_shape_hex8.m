%test cases for get_shape function: HEX8
%test cases for Basis functions and their Derivatives, B and D_hat's
clear
clc

cd ../..

%test for HEX9 (8-noded 3D element) using Guassian Quadrature

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

num_quadr_pts_in_1d=2;
dimension=3;

[gs,~ ]= get_quadrature(num_quadr_pts_in_1d);
gs = abs(gs(1));
x = gs*[-1  1 -1  1 -1  1 -1 1]';
y = gs*[-1 -1  1  1 -1 -1  1 1]';
z = gs*[-1 -1 -1 -1  1  1  1 1]';

% analytical evaluation
B = f(x,y,z);
D0= fx(y,z);
D1= fy(x,z);
D2= fz(x,y);

% numerical evalution
[B1,D_hat,~]=get_shape(num_quadr_pts_in_1d, dimension);
num_gs_pts = size(D_hat,1)/dimension;


% difference between analytical and numerical evalutions of B and D_hat
diffB = B -B1;
diffD0 = D0 - D_hat(1:num_gs_pts,:);
diffD1 = D1 - D_hat(num_gs_pts+1:2*num_gs_pts,:);
diffD2 = D2 - D_hat(2*num_gs_pts+1:end,:);
if (norm(diffB,2) < 1.0e-14)
   tfB = 1;
else
    tfB = 0;
end
if (norm(diffD0,2) < 1.0e-14)
   tfD0 = 1;
else
    tfD0 = 0;
end
if (norm(diffD1,2) < 1.0e-14)
   tfD1 = 1;
else
    tfD1 = 0;
end
if (norm(diffD2,2) < 1.0e-14)
   tfD2 = 1;
else
    tfD2 = 0;
end

disp('HEX8:')
disp('1 = pass');
disp('0 = fail');
disp('     B     D0    D1    D2')
disp('     --------------------')
res = sprintf('     %d     %d     %d     %d\n',tfB, tfD0, tfD1,tfD2);
disp(res); 


cd test_units/test_get_shape


