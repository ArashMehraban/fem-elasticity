function global_u =  get_global_u(u,dir_bndry_nodes,dir_bndry_val,global_idx_map)
%GET_GLOBAL_U returns a the global_u that includes Dirichlet Boundary values
     
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
    
%     s = global_idx_map(4);
%     i = 1:numel(global_idx_map);
%     a = find(global_idx_map(i)>0);
%     on_bndry = cell2mat(dir_bndry_nodes);
%     j = 1:size(on_bndry,1);
%     not_bndry = linspace(1,size(global_u,1),size(global_u,1));
%     not_bndry(on_bndry(j))=[];
%     global_u(not_bndry) = u;
end