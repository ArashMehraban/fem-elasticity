%test cases for get_shape function: QUAD9
%test cases for Basis functions and their Derivatives, B and D_hat's
clear
clc

cd ../..
%test for QUAD9 (9-noded 2D element) using Guassian Quadrature

f = @(x,y) [0.25.*(x.^2-x).*(y.^2-y), ...
            0.5.*(1- x.^2).*(y.^2-y), ...
            0.25.*(x.^2+x).*(y.^2-y), ...
            0.5.*(x.^2-x).*(1-y.^2), ...
            (1-x.^2).*(1-y.^2), ...
            0.5.*(x.^2+x).*(1-y.^2) ...
            0.25.*(x.^2-x).*(y.^2+y), ...
            0.5.*(1-x.^2).*(y.^2+y), ...
            0.25.*(x.^2+x).*(y.^2+y)];

syms x y 

dfx = diff(f,x);
dfy = diff(f,y);

fx = matlabFunction(dfx);
fy = matlabFunction(dfy);                     % (-1,1)     (0,1)    (1,1)
                                              %  o---------o---------o
num_quadr_pts_in_1d=3;                        %  |                   |
dimension=2;                                  %  o(-1,0)   o(0,0)    o(1,0)
                                              %  |                   |
[gs,~ ]= get_quadrature(num_quadr_pts_in_1d); %  o---------o---------o 
gs = abs(gs(1));                              % (-1,-1)    (0,-1)   (1,-1)
x = gs*[-1  0  1 -1  0  1 -1  0  1]';       
y = gs*[-1 -1 -1  0  0  0  1  1  1]';       

% analytical evaluation
B = f(x,y);
D0= fx(x,y);
D1= fy(x,y);

% numerical evalution
[B1,D_hat,~]=get_shape(num_quadr_pts_in_1d, dimension);
num_gs_pts = size(D_hat,1)/dimension;

% difference between analytical and numerical evalutions of B and D_hat
diffB = B -B1;
diffD0 = D0 - D_hat(1:num_gs_pts,:);
diffD1 = D1 - D_hat(num_gs_pts+1:end,:);
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

disp('QUAD9:')
disp('1 = pass');
disp('0 = fail');
disp('     B    D0   D1')
disp('     --------------')
res = sprintf('     %d     %d      %d\n',tfB, tfD0, tfD1);
disp(res); 

cd test_units/test_get_shape