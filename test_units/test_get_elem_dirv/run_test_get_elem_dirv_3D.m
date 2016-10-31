%3D
%test case for [~, D] = get_elem_dirv(element_1_vtx_coord, D_hat, dimension);
clear
clc

cd ../.. 

filename={'cylinder8_110e_us';'cylinder27_110e_us'};

folderName='mesh_Dirich_only';

addpath(fullfile(pwd,folderName));

cd test_units/test_get_elem_dirv


f={@(x,y,z)(0.*x + 0.*y + 0.*z + 5) ; @(x,y,z)(-3.*x + 2.*y - 4*z); @(x,y,z)(x.^2 + 0.5*y.^2 -2*z.^2)};
fx={@(x)(0.*x +0);                    @(x)(0.*x - 3);               @(x)(2.*x)};
fy={@(y)(0.*y +0);                    @(y)(0.*y + 2);               @(y)(y)};
fz={@(z)(0.*z +0);                    @(z)(0.*z - 4);               @(z)(-4.*z)};
df={fx;fy;fz};

num_quadr_pts_in_1d=[2,3];
dim=3;
ith_elem=[1, 10, 32 ,78, 110];

disp('Test for:');
disp('     [~, D] = get_elem_dirv(ith_element_vtx_coord, D_hat, dimension)');
disp(' ');
disp('        Pass: means numerical precision is AT LEAST within 1.0e-13 of analytical');
disp('       value: is the norm of the difference between numerical and analytical ');
disp(' ');
disp('  elem            func                   filename                 err_norm(D0)        err_norm(D1)      err_norm(D2) ');
disp('  -------------------------------------------------------------------------------------------------------------------');


for i=1:size(filename,1)
    for j=1:size(f,1)
        df={fx{j};fy{j};fz{j}};
        for k=1:size(ith_elem,2)
            [pfD0,pfD1,pfD2] = test_get_elem_dirv(filename{i},f{j},df,num_quadr_pts_in_1d(i),dim, ith_elem(k));
            res = sprintf('%5d  %30s  %15s   %15s  %15s  %15s\n',ith_elem(k), func2str(f{j}), filename{i}, pfD0, pfD1,pfD1);
            disp(res); 
        end
        disp('Next function:')
        disp('----------------------------------------------------------------------------------------------------------------')
    end
end

cd ../..
rmpath(fullfile(pwd,folderName));
cd test_units/test_get_elem_dirv