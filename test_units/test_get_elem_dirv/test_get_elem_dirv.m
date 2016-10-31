function [tfD0,tfD1,tfD2] = test_get_elem_dirv(filename,f,df,num_quadr_pts_in_1d,dim, ith_elem)
  
    msh = get_mesh(filename,'exo','lex');
    elem_vtx_coords = msh.vtx_coords;
    conn = msh.conn;
    if(ith_elem > size(conn,1))
        msg = sprintf('Mesh has %d elements. %d is chosen',size(conn,1), ith_elem);
        disp(msg);
        error ('element does not exist!');
    end
    ith_elem_vtx_coord = elem_vtx_coords(conn(ith_elem,:),:);
    [B,D_hat,~]=get_shape(num_quadr_pts_in_1d, dim);
    [~, D] = get_elem_dirv(ith_elem_vtx_coord, D_hat, dim);
    num_gs_pts = size(D,1)/dim;
    xe = B*ith_elem_vtx_coord;
    
    if(dim == 2)
        evalf=f(ith_elem_vtx_coord(:,1),ith_elem_vtx_coord(:,2));
        df_dx_numerical = D(1:num_gs_pts,:)*evalf;
        df_dy_numerical = D(num_gs_pts+1:end,:)*evalf;
        fx = df{1};
        fy = df{2};
        fx_anal = fx(xe(:,1));
        fy_anal = fy(xe(:,2));
        diff_x= df_dx_numerical - fx_anal;
        diff_y= df_dy_numerical - fy_anal;
        if(norm(diff_x,2) < 1.0e-14)
            tfD0 = 'Pass';
        else
            tfD0 = num2str(norm(diff_x,2));
        end
        if(norm(diff_y,2) < 1.0e-14)
            tfD1 = 'Pass';
        else
            tfD1 = num2str(norm(diff_y,2));
        end
        
    end
    if(dim == 3)
        evalf=f(ith_elem_vtx_coord(:,1),ith_elem_vtx_coord(:,2),ith_elem_vtx_coord(:,3));
        df_dx_numerical = D(1:num_gs_pts,:)*evalf;
        df_dy_numerical = D(num_gs_pts+1:2*num_gs_pts,:)*evalf;
        df_dz_numerical = D(2*num_gs_pts+1:end,:)*evalf;
        fx = df{1};
        fy = df{2};
        fz = df{3};
        fx_anal = fx(xe(:,1));
        fy_anal = fy(xe(:,2));
        fz_anal = fz(xe(:,3));
        diff_x= df_dx_numerical - fx_anal;
        diff_y= df_dy_numerical - fy_anal;
        diff_z= df_dz_numerical - fz_anal;
        if(norm(diff_x,2) < 1.0e-13)
            tfD0 = 'Pass';
        else
            tfD0 = num2str(norm(diff_x,2));
        end
        if(norm(diff_y,2) < 1.0e-13)
            tfD1 = 'Pass';
        else
            tfD1 = norm(diff_y,2);
        end
        if(norm(diff_z,2) < 1.0e-13)
            tfD2 = 'Pass';
        else
            tfD2 = num2str(norm(diff_z,3));
        end
        
    end
end