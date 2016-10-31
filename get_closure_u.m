function closure_u =  get_closure_u(u,dir_bndry_nodes,dir_bndry_val,global_idx_map)
%GET_CLOSURE_U returns closure of u, a u vector that contains boundary values
%at boundary nodes
%
% input:                u: solution/unknown u
%      : dir_bndry_nodes : Dirichlet boundary nodes
%      :   dir_bndry_val : Dirichlet boundary values
%      :  global_idx_map : global mapping of u
%
% output: closure_u: global u that includes Dirichlet boundary values
     
    sz_global_u = size(global_idx_map);
    closure_u =zeros(sz_global_u);
    for i=1:size(dir_bndry_nodes,1)
       tmp_nods =  dir_bndry_nodes{i};
       tmp_vals =  dir_bndry_val{i};
       closure_u(tmp_nods,:) = tmp_vals;       
    end
     
    k=1;
    for i=1:size(global_idx_map,1)
        for j=1:size(global_idx_map,2)
            if(global_idx_map(i,j)>0)
                closure_u(i,j) = u(k);
                k=k+1;
            end
        end
    end
    
end