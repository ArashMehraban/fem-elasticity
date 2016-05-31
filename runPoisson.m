%Solving Poisson problem using FEM with unstructured mesh
clear
clc
format short
  
%filenames{1} = 'disk4_37e_us'; 
%  filenames{1} = 'square9_4e_s';
%  filenames{2} = 'square9_16e_s';
%  filenames{3} = 'square9_36e_s';
% 
%     filenames{1} = 'disk9_4e';  
%   filenames{2} = 'disk9_118e';
%   filenames{3} = 'disk9_274e';
%  filenames{4} = 'disk9_641e';
 
 %filenames{1} = 'cylinder8_800e_s';
%  filenames{1} = 'cylinder8_110e_us';
%   filenames{2} = 'cylinder8_368e_us';
%   filenames{3} = 'cylinder8_1176e_us';  
%   filenames{4} = 'cylinder8_4180e_us';  

  %   filenames{1} = 'cube27_27e_s';
  
  filenames{1} = 'cylinder27_64e_us';
  filenames{2} = 'cylinder27_110e_us';
  %filenames{3} = 'cylinder27_1176e_us'; 

  ext = 'exo';

sz_u_field = 1;

err = zeros(1,size(filenames,2));

for i=1:size(filenames,2)
    
    msh = get_mesh(filenames{i},ext,'lex');
    %msh = create_mesh('QUAD9', 3,0,1, 3,0,1);
    %msh = create_mesh('QUAD9', 6,0,1, 6,0,1);
    
    %msh = create_mesh('HEX27', 10,0,1, 10,0,1, 5,0,1);
    
    %get all Dirichlet boundary node sets
    dir_bndry_nodes = get_all_dir_ns(msh);
    
    %NOTE: modify userf function according to given_u
    %given_u{1}=@(x,y)tanh(x).*exp(y)+sin(y);
    given_u{1}=@(x,y,z)tanh(x).*exp(y)+sin(y)+cos(z);
        
    %Construct manufactured solution:
    %========================================================================================%
    %get vertex coordinates from mesh                                     
    vtx_coords = msh.vtx_coords;                                          
    %get constructed dir_bndry_vals and exac Solutions on remaining nodes 
    [dir_bndry_val, exactSol] = get_exact_sol(vtx_coords,dir_bndry_nodes, given_u);
    %========================================================================================%
    
    fem_sol =  main(msh, sz_u_field, dir_bndry_nodes, dir_bndry_val);
    error = norm(exactSol - fem_sol)/norm(fem_sol);
    err(i) =error; 
    
end
    





