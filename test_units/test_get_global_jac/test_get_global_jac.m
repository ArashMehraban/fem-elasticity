%Test for get_global_jac function with Plain Strain problem using unstructured mesh
clear
clc
format short

cd ../.. 


file_name_list{1} = {'disk4_114e_us'};
%file_name_list{1} = {'disk4_114e_us','disk4_274e_us'};
%file_name_list{2} = {'disk9_114e_us'};

folderName='mesh_Dirich_only';

addpath(fullfile(pwd,folderName));

files = cell(file_name_list);
for i=1:size(file_name_list,2)
    files{i} = get_files(folderName, file_name_list{i});
end

%Provide the number of quadrature points in 1D per element used per set of files
num_quadr_pts_in_1d=[2,3];
userf={@userf_plainStrain;@userf_plainStrain};
userdf={@userdf_plainStrain;@userdf_plainStrain};
max_iter_gmres = 80;
max_iter_nw = 10;
global_res_tol = 1.0e-8;
tol = 1.0e-9;
%To produce the Jacobian from GMRES
jac_flag =1;

for i=1:size(files,2)
    error_gmres = zeros(1,size(files{i},2));
    error_nw = zeros(1,size(files{i},2));
    h = zeros(1,size(files{i},2)); 
    for j=1:size(files{i},2)        
        filenames = files{i}{j};
        ext = 'exo';
        msh = get_mesh(filenames,ext,'lex');
        elem_type = msh.num_nodes_per_elem;
        h(j) = get_hsz(msh);
        
        sz_u_field = 2;
    
        %get all Dirichlet boundary node sets
        dir_bndry_nodes = get_all_dir_ns(msh);
    
        %NOTE: modify userf function according to given_u
        given_u{1}=@(x,y)tanh(x).*exp(y)+sin(y);
        given_u{2}=@(x,y)tanh(x).*cos(y); 


        %Construct manufactured solution:
        %========================================================================================%
        %get vertex coordinates from mesh                                     
        vtx_coords = msh.vtx_coords;                                          
        %get constructed dir_bndry_vals and exac Solutions on remaining nodes 
        [dir_bndry_val, exactSol] = get_exact_sol(vtx_coords,dir_bndry_nodes, given_u);
        %========================================================================================%
        
        % Solve the problem with GMRES:
        solver = {'gmres', max_iter_gmres,max_iter_nw,tol,global_res_tol, jac_flag};
        [fem_sol_gmres,global_res_norm_gmres,JACOB_gmres] =  get_fem_sol(msh, sz_u_field, dir_bndry_nodes, dir_bndry_val,num_quadr_pts_in_1d(i),userf{i},userdf{i},solver);      
        error_gmres(j) = norm(exactSol - fem_sol_gmres)/norm(fem_sol_gmres);
        gmres_cond = condest(JACOB_gmres);
        % Solve the problem with Newton:
        solver = {'newton',max_iter_nw, tol,global_res_tol};
        [fem_sol_nw,global_res_norm_nw,JACOB_nw] =  get_fem_sol(msh, sz_u_field, dir_bndry_nodes, dir_bndry_val,num_quadr_pts_in_1d(i),userf{i},userdf{i},solver);      
        error_nw(j) = norm(exactSol - fem_sol_nw)/norm(fem_sol_nw);
        nw_cond = condest(JACOB_nw);
        
        diffJac = JACOB_nw - JACOB_gmres;
                
        res0 = sprintf('filename: %s',filenames);
        res1= sprintf('   Consistent Tangent condition number at Solution (GMRES) %d:',gmres_cond);
        res2= sprintf('   Consistent Tangent condition number at Solution (fsolve) %d:',nw_cond);
        res3= sprintf('   Norm of the of the difference bewtween the Consistant Tangents from GMRES and fsolve %d:',norm(diffJac));
        res4= sprintf('   Norm of the difference between Manufatured solution vs FEM solution (GMRES): %d', error_gmres(j));
        res5= sprintf('   Norm of the difference between Manufatured solution vs FEM solution (fsolve): %d', error_nw(j));
        res6= sprintf('   Norm of global residual after each GMRES solve:');
        res7= sprintf('   Norm of global residual after fsolve:');
        disp(res0)
        disp(res1)
        disp(res2)
        disp(res3)
        disp(res4)
        disp(res5)
        disp(res6)
        disp(global_res_norm_gmres)
        disp(res7)
        disp(global_res_norm_nw)
        disp('-------------------------------------------------------------------------------------------------------------------')
    end
end

rmpath(fullfile(pwd,folderName));
cd test_units/test_get_global_jac


%% Results:
% filename: disk4_114e_us
%    Consistent Tangent condition number at Solution (GMRES) 8.615823e+01:
%    Consistent Tangent condition number at Solution (fsolve) 8.615823e+01:
%    Norm of the of the difference bewtween the Consistant Tangents from GMRES and fsolve 1.681775e-07:
%    Norm of the difference between Manufatured solution vs FEM solution (GMRES): 3.606818e-01
%    Norm of the difference between Manufatured solution vs FEM solution (fsolve): 3.606818e-01
%    Norm of global residual after each GMRES solve:
%    3.4088e-09
% 
%    Norm of global residual after fsolve:
%    8.5866e-09
%-------------------------------------------------------------------------------------------------------------------
%
% filename: disk4_274e_us
%    Consistent Tangent condition number at Solution (GMRES) 1.987300e+02:
%    Consistent Tangent condition number at Solution (fsolve) 1.987300e+02:
%    Norm of the of the difference bewtween the Consistant Tangents from GMRES and fsolve 1.571479e-07:
%    Norm of the difference between Manufatured solution vs FEM solution (GMRES): 3.139636e-01
%    Norm of the difference between Manufatured solution vs FEM solution (fsolve): 3.139636e-01
%    Norm of global residual after each GMRES solve:
%    2.8141e-09
% 
%    Norm of global residual after fsolve:
%    5.2183e-09
%-------------------------------------------------------------------------------------------------------------------
%
%filename: disk9_114e_us
%    Consistent Tangent condition number at Solution (GMRES) 7.330716e+02:
%    Consistent Tangent condition number at Solution (fsolve) 7.330716e+02:
%    Norm of the of the difference bewtween the Consistant Tangents from GMRES and fsolve 4.078288e-07:
%    Norm of the difference between Manufatured solution vs FEM solution (GMRES): 3.118932e-01
%    Norm of the difference between Manufatured solution vs FEM solution (fsolve): 3.118933e-01
%    Norm of global residual after each GMRES solve:
%    6.2060e-09
% 
%    Norm of global residual after fsolve:
%    2.0866e-08
%%
