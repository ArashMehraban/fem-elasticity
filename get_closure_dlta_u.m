function closure_dlta_u =  get_closure_dlta_u(dlta_u,dir_bndry_nodes,global_idx_map)
%GET_CLOSURE_DELTA_U returns closure of delta_u, a delta_u vector that contains 0
%as variation of u on boundary nodes
%
%
% input:          dlta_u : variation of u
%      : dir_bndry_nodes : Dirichlet boundary nodes
%      :  global_idx_map : global mapping of u
%
% output: closure_dlta_u: global delta u (variation of u) that has 0 on Dirichlet
% boundary nodes
     
    sz_global_u = size(global_idx_map);
    closure_dlta_u =zeros(sz_global_u);
    for i=1:size(dir_bndry_nodes,1)
       tmp_nods =  dir_bndry_nodes{i};
       closure_dlta_u(tmp_nods,:) = 0;       
    end
     
    k=1;
    for i=1:size(global_idx_map,1)
        for j=1:size(global_idx_map,2)
            if(global_idx_map(i,j)>0)
                closure_dlta_u(i,j) = dlta_u(k);
                k=k+1;
            end
        end
    end
    
end