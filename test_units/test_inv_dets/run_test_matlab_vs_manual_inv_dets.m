%test case for inverses and determinants of an element Jacobian
clear
clc

format short

cd ../.. 

filename={'disk4_114e_us';'disk9_114e_us';'cylinder8_110e_us';'cylinder27_110e_us'};

folderName='mesh_Dirich_only';

addpath(fullfile(pwd,folderName));

cd test_units/test_inv_dets

ith_elem=[1, 10, 32 ,78, 110];
num_quadr_pts_in_1d = [2,3,2,3];
dim=[2,2,3,3];

disp('Test for the determinants and inverse Jacobian for an element in:');
disp('     get_elem_dirv(ith_element_vtx_coord, D_hat, dimension)');
disp(' ');
disp('       Pass: means vectorized inverse and determinants are equal to MATLAB "inv" and "det"');
disp('       value: means norm of difference between vertorized inverse and determinant and MATLAB "inv" and "det" ');
disp(' ');
disp('  elem       dets        inv               filename');
disp('  --------------------------------------------------------');

for i=1:size(filename,1)
    msh = get_mesh(filename{i},'exo','lex');
    elem_vtx_coords = msh.vtx_coords;
    conn = msh.conn;
    for j=1:size(ith_elem,2)
       ith_elem_vtx_coords = elem_vtx_coords(conn(ith_elem(j),:),:);
       [~,D_hat,~]=get_shape(num_quadr_pts_in_1d(i), dim(i));
       [tf_dets, tf_inv] = test_matlab_vs_manual_inv_det(ith_elem_vtx_coords, D_hat,dim(i));
       res = sprintf('%5d  %10s  %10s   %25s\n',ith_elem(j), tf_dets, tf_inv, filename{i});
       disp(res);
    end
    disp('  --------------------------------------------------------');    
end

cd ../..
rmpath(fullfile(pwd,folderName));
cd test_units/test_inv_dets
