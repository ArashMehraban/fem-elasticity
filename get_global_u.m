function global_u =  get_global_u(u,dir_bndry_nodes,dir_bndry_val,global_idx_map)
%GET_GLOBAL_U returns a the global_u that includes Dirichlet Boundary values
%
% input:                u: solution/unknown u
%      : dir_bndry_nodes : Dirichlet boundary nodes
%      :   dir_bndry_val : Dirichlet boundary values
%      :  global_idx_map : global mapping of u
%
%output: global_u: global u that includes Dirichlet boundary values
     
    sz_global_u = size(global_idx_map);
    global_u =zeros(sz_global_u);
    for i=1:size(dir_bndry_nodes,1)
       tmp_nods =  dir_bndry_nodes{i};
       tmp_vals =  dir_bndry_val{i};
       global_u(tmp_nods,:) = tmp_vals;       
    end
     
    k=1;
    for i=1:size(global_idx_map,1)
        for j=1:size(global_idx_map,2)
            if(global_idx_map(i,j)>0)
                global_u(i,j) = u(k);
                k=k+1;
            end
        end
    end
    
end