function global_u =  setbdry(u,num_nodes, num_dims, dir_bndry_nodes,dir_bndry_val)
%SETBDRY_DIR returns a the global_u that includes Dichlet Boundary values

    global_u =zeros(num_nodes*num_dims,1);
    for i=1:size(dir_bndry_nodes,1)
       tmp_nods =  dir_bndry_nodes{i};
       tmp_vals =  dir_bndry_val{i};
       j=1:size(tmp_nods,1);
       global_u(tmp_nods(j)) = tmp_vals;       
    end
    
    on_bndry = cell2mat(dir_bndry_nodes);
    j = 1:size(on_bndry,1);
    not_bndry = linspace(1,size(global_u,1),size(global_u,1));
    not_bndry(on_bndry(j))=[];
    global_u(not_bndry) = u;
end