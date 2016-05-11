function [dir_bndry_val, exactSol] = get_exact_sol(vtx_coords,dir_bndry_nodes)
% Use GET_EXCAT_SOL function to construct Diricklet Boundary values
% at specifed Dirchlet boundary nodes passed to GET_EXCAT_SOL     

    %get the exact solution from constructed equations at each node/field
    exactSol = exactf(vtx_coords);
   
    % Allocate space for dir_bndry_vals
    dir_bndry_val = cellfun(@(x) x*0,dir_bndry_nodes,'un',0);
    
    % get the number of dir_bndry_nodes sets
    sz_dir_ns = size(dir_bndry_nodes,1);
    
    %populate the dir_bndry_vals and return exactSol arrays minus the dir_bndry_vals
    for i=1:sz_dir_ns
       tempExct = exactSol{i};
       temp_dir_bn_nods = dir_bndry_nodes{i};
       temp_dir_bd_idx = 1:size(temp_dir_bn_nods,1);
       dir_bndry_val{i} = tempExct(temp_dir_bn_nods(temp_dir_bd_idx));
       tempExct(temp_dir_bn_nods(temp_dir_bd_idx))=[];
       exactSol{i} = tempExct;
    end
end