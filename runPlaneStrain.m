%+++++++++
% Old Plane Strain. Must be updated with new code!!
%++++++++

%Solving Plane Strain problem using FEM with unstructured mesh
clear
clc
format short


filenames{1} = 'square4';
% filenames{1} = 'disk9_4e';
% filenames{2} = 'disk9_118e';
% filenames{3} = 'disk9_274e';
% filenames{4} = 'disk9_641e';
ext = 'exo';

sz_u_field = 2;

err = zeros(1,size(filenames,2));

for i=1:size(filenames,2)
    
    msh = get_mesh(filenames{i},ext,'lex');
    co = msh.conn;
    a = co(1,:);
    con = fliplr(co);
    vtx = msh.vtx_coords;
    
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
    
    %dir_bndry_val = cellfun(@(x) x*2e11,dir_bndry_val,'un',0);
    
    fem_sol =  main(msh, sz_u_field, dir_bndry_nodes, dir_bndry_val);
    
    error = norm(exactSol - fem_sol)/norm(fem_sol);
    err(i) =error; 
    
end
    
