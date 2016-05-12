function [u, global_idx_map] = preproc(sz_vtx_coords,known_bd)
    
    field_sz = size(known_bd,1);
       
    u_section = ones(sz_vtx_coords,field_sz+1);
    u_section(:,1:3:end) =-1;
    for i=1:field_sz
        u_section(known_bd{i},1) = 0;
    end
    
    u_sz = -sum(sum(u_section(:,1)));
    u = zeros(u_sz*field_sz,1);
    
    for i=1:size(u_section,1)
        if(u_section(i,1) == 0)
            u_section(i,2:end) = -u_section(i,2:end);
        end        
    end
    
    global_idx_map=u_section(:,2:end);
    [r_u_sec, c_u_sec]=size(global_idx_map);
    
    global_idx_map=reshape(global_idx_map',[],1);
    
    non_bd_idx = find(global_idx_map > 0);
    
    sz_non_bd_idx = size(non_bd_idx,1);
    
    u_seq = 1:sz_non_bd_idx;
    
    global_idx_map(non_bd_idx) = u_seq;
    
    global_idx_map = reshape(global_idx_map,c_u_sec,r_u_sec)';
    
end