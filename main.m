function u = main(msh, sz_u_field, dir_bndry_nodes, dir_bndry_val)
    
    %get number of nodes in mesh
    num_nodes = msh.num_nodes;  
    
    %get the known u vector and global to local mapping foe each u
    [u, global_idx_map] = preproc(num_nodes,dir_bndry_nodes,sz_u_field);
           
    global_res_norm=1;
    iter=1;
    max_iter = 3;
    gl_res_norm = zeros(max_iter,1);
    gl_res_norm_iter = 1;
    while((global_res_norm > 1.0e-9) && iter < max_iter)
        
        [global_res, jac] = eval_res(u, global_idx_map, msh, dir_bndry_val);          
                                       
        global_res_norm = norm(global_res);
        if(global_res_norm < 1.0e-9)
            break;
        end   
        gl_res_norm(gl_res_norm_iter) = global_res_norm;
        gl_res_norm_iter = gl_res_norm_iter+1; 
        if(gl_res_norm_iter >= 3)
            if(gl_res_norm(gl_res_norm_iter-1) > 100*gl_res_norm(gl_res_norm_iter-2) )
                error('The solution is diverging!!')
            end
        end
        
        fun = @(u)eval_res(u, global_idx_map, msh, dir_bndry_val);
               
        options = optimoptions(@fsolve,'Algorithm','trust-region-reflective','Jacobian','on');
        %options = optimoptions(@fsolve,'Algorithm','trust-region-reflective');
                
        u = fsolve(fun, global_res, options); 
        iter=iter+1;
    end  
      
 



