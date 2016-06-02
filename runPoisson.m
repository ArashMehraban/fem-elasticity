%Solving Poisson problem using FEM with unstructured mesh produced by Cubit
clear
clc
format short

file_name_list{1} = {'disk4_114e_us', 'disk4_169e_us', 'disk4_274e_us', 'disk4_641e_us','disk4_2598e_us' };
file_name_list{2} = {'disk9_114e_us', 'disk9_169e_us', 'disk9_274e_us', 'disk9_641e_us','disk9_2598e_us' };
file_name_list{3} = {'cylinder8_110e_us', 'cylinder8_368e_us', 'cylinder8_1176e_us', ...
                     'cylinder8_1635e_us','cylinder8_4180e_us' };
file_name_list{4} = {'cylinder27_64e_us', 'cylinder27_110e_us', 'cylinder27_368e_us', ...
                     'cylinder27_630e_us','cylinder27_1428e_us' };
folderName='mesh_Dirich_only';
files = cell(file_name_list);
for i=1:size(file_name_list,2)
    files{i} = get_files(folderName, file_name_list{i});
end

for i=1:size(files,2)
    error = zeros(1,size(files{i},2));
    h = zeros(1,size(files{i},2)); 
    for j=1:size(files{i},2)        
        filenames = files{i}{j};
        ext = 'exo';
        msh = get_mesh(filenames,ext,'lex');
        dims = msh.num_dims;
        elem_type = msh.num_nodes_per_elem;
        h(j) = get_hsz(msh);
         
        sz_u_field = 1;
        
        %get all Dirichlet boundary node sets
        dir_bndry_nodes = get_all_dir_ns(msh);

        %NOTE: modify userf function according to given_u
        if(dims==2)
            given_u{1}=@(x,y)tanh(x).*exp(y)+sin(y);
        elseif(dims==3)
            given_u{1}=@(x,y,z)tanh(x).*exp(y)+sin(y)+cos(z);
        end

        %Construct manufactured solution:
        %========================================================================================%
        %get vertex coordinates from mesh                                     
        vtx_coords = msh.vtx_coords;                                          
        %get constructed dir_bndry_vals and exac Solutions on remaining nodes 
        [dir_bndry_val, exactSol] = get_exact_sol(vtx_coords,dir_bndry_nodes, given_u);
        %========================================================================================%

        fem_sol =  get_sol(msh, sz_u_field, dir_bndry_nodes, dir_bndry_val);
        error(j) = norm(exactSol - fem_sol)/norm(fem_sol);
    end
    if(dims ==2)
        prb_title = 'Poisson: 2D';
        if(elem_type == 4)
            leg_enry_1 = 'FEM-QUAD4';
            lglg_factor_1 = 0.001;
            lglg_pwr_1 = h;
            lglg_factor_2 = 0.1;
            lglg_pwr_2 = h.^2;
            leg_enry_2 = 'O(h)';
            leg_enry_3 = 'O(h^2)';
        elseif(elem_type == 9)
            leg_enry_1 = 'FEM-QUAD9';
            lglg_factor_1 = 0.00001;
            lglg_pwr_1 = h.^2;
            lglg_factor_2 = 0.001;
            lglg_pwr_2 = h.^3;
            leg_enry_2 = 'O(h^2)';
            leg_enry_3 = 'O(h^3)';
        end
    elseif(dims ==3)
        prb_title = 'Poisson: 3D';
        if(elem_type == 8)
            leg_enry_1 = 'FEM-HEX8';
            lglg_factor_1 = 0.001;
            lglg_pwr_1 = h;
            lglg_factor_2 = 0.1;
            lglg_pwr_2 = h.^2;
            leg_enry_2 = 'O(h)';
            leg_enry_3 = 'O(h^2)';
        elseif(elem_type == 27)
            leg_enry_1 = 'FEM-HEX27';
            lglg_factor_1 = 0.0001;
            lglg_pwr_1 = h.^2;
            lglg_factor_2 = 0.001;
            lglg_pwr_2 = h.^3;        
            leg_enry_2 = 'O(h^2)';
            leg_enry_3 = 'O(h^3)';
        end
    end
    subplot(2,2,i);
    loglog(h,error,'r-o', h,lglg_factor_1*(lglg_pwr_1),'r:',h, lglg_factor_2*(lglg_pwr_2),'b--');
    legend(leg_enry_1,leg_enry_2,leg_enry_3,'Location','northwest')
    title(prb_title)
    xlabel('h')
    ylabel('error')    
end

    





