function max_h = get_hsz(msh)
%GET_HSZ returns the size of the largest element:
%     
%input: msh (mesh)
%
%output: max_h: 
%             2D: Square root of the largest element's area
%             3D: Cube root of the largest element's volume
    
    dims = msh.num_dims;
    vtx = msh.vtx_coords;
    conn = msh.conn;
    sz_conn = size(conn,1);
    h = zeros(1,size(conn,1));
    if(dims == 2)        
        for j=1:sz_conn
            min_x = min(vtx(conn(j,:),1));
            max_x = max(vtx(conn(j,:),1));
            min_y = min(vtx(conn(j,:),2));
            max_y = max(vtx(conn(j,:),2));
            h_x = max_x -min_x;
            h_y = max_y -min_y;
            h(j) = h_x *h_y;       
        end
        max_h = sqrt(max(h));
    else
        for j=1:sz_conn
            min_x = min(vtx(conn(j,:),1));
            max_x = max(vtx(conn(j,:),1));
            min_y = min(vtx(conn(j,:),2));
            max_y = max(vtx(conn(j,:),2));
            min_z = min(vtx(conn(j,:),3));
            max_z = max(vtx(conn(j,:),3));
            h_x = max_x -min_x;
            h_y = max_y -min_y;
            h_z = max_z -min_z;
            h(j) = h_x *h_y*h_z;       
        end
        max_h = nthroot(max(h),3);        
    end
end