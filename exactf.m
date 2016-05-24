function exctSol = exactf(vtx_coords, mesh_dim, given_u)
    
    x = vtx_coords(:,1);
    y = vtx_coords(:,2);
    if(mesh_dim>2)
        z = vtx_coords(:,3);        
    end
    
    exctSol = zeros(size(vtx_coords,1),mesh_dim);
    
    sz_gvn_u = size(given_u,2);
    if(mesh_dim==2)
        for i=1:sz_gvn_u;
            exctSol(:,i) = given_u{i}(x,y);
        end
    else
        for i=1:sz_gvn_u;
            exctSol(:,i) = given_u{i}(x,y,z);
        end
    end

    
end