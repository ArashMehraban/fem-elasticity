function [dir_bndry_val, exactSol] = get_exact_sol(vtx_coords,dir_bndry_nodes, given_u)
% Use GET_EXCAT_SOL function to construct Diricklet Boundary values
% at specifed Dirchlet boundary nodes passed to GET_EXCAT_SOL     

    %get the exact solution from constructed equations at each node/field
    exactSol = exactf(vtx_coords, given_u);
              
    dir_bndry_val = cellfun(@(x) 0*x,dir_bndry_nodes,'un',0);
    
    % get the number of dir_bndry_nodes sets
    sz_dir_ns = size(dir_bndry_nodes,1);
    
    %populate the dir_bndry_vals and return exactSol arrays minus the dir_bndry_vals
    for i=1:sz_dir_ns
       temp_dir_bn_nods = dir_bndry_nodes{i};
       temp_dir_bd_idx = 1:size(temp_dir_bn_nods,1);
       dir_bndry_val{i} = exactSol(temp_dir_bn_nods(temp_dir_bd_idx),:);
    end
    
    all_dir_bn_nods = cell2mat(dir_bndry_nodes);
    all_dir_bn_nods_idx = 1:size(all_dir_bn_nods,1);
    exactSol(all_dir_bn_nods(all_dir_bn_nods_idx),:)=[];
    
    exactSol = reshape(exactSol',[],1);
end
% helper function for coordinate evaluation
function exctSol = exactf(vtx_coords, given_u)
    
    x = vtx_coords(:,1);
    y = vtx_coords(:,2);
    
    num_dim = size(vtx_coords,2);
    if(num_dim>2)
        z = vtx_coords(:,3);        
    end
    
    sz_gvn_u = size(given_u,2);
    exctSol = zeros(size(vtx_coords,1),sz_gvn_u);    
    
    if(num_dim==2)
        for i=1:sz_gvn_u;
            exctSol(:,i) = given_u{i}(x,y);
        end
    else
        for i=1:sz_gvn_u;
            exctSol(:,i) = given_u{i}(x,y,z);
        end
    end

    
end