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