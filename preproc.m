function [u, global_idx_map] = preproc(sz_vtx_coords,known_bd)
    
    field_sz = size(known_bd,1);
       
    u_section = ones(sz_vtx_coords,field_sz*3);
    u_section(:,1:3:end) =-1;
    for i=1:field_sz
        u_section(known_bd{i},(i-1)*3+1) = 0;
    end
    
    u_sz = -sum(sum(u_section(:,1:3:end)));
    u = zeros(u_sz,1); 
    
    [r_u_section, c_u_section]=size(u_section);
    
    u_section=reshape(u_section',[],1);
           
    for i=5:3:size(u_section,1)
        if(u_section(i-1) == 0)
            u_section(i)= u_section(i-3);
            u_section(i+1)= -u_section(i);
        else
            u_section(i)= u_section(i-3)+1;
            u_section(i+1)= u_section(i);
        end
    end
    
    u_section = reshape(u_section,c_u_section,r_u_section)';
    
    global_idx_map = u_section(:,3:3:end);
      
end