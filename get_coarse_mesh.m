function coarse_mesh_conn = get_coarse_mesh(msh)
     
    elem_type = msh.num_nodes_per_elem;
    conn = msh.conn;
    i=1:size(conn,1);
    if(elem_type == 9)
        coarse = [1 3 7 9];
        coarse_mesh_conn = conn(i,coarse);
        
    end
    if(elem_type == 27) 
        coarse = [1 3 7 9 19 21 25 27];
        coarse_mesh_conn = conn(i,coarse);
    end

end