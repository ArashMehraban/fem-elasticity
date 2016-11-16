function [u,global_res_norm,JACOB__] = get_fem_sol(msh, sz_u_field, dir_bndry_nodes, dir_bndry_val, num_quadr_pts_in_1d,userf,userdf,solver)
%GET_FEM_SOL returns the solution to the PDE
%
%input :                msh: mesh object
%      :         sz_u_field: size of unknown field (eg. : Poisson 1, Plane Strain 2)
%      :    dir_bndry_nodes: Dirichlet Boundary nodes
%      :      dir_bndry_val: Dirichlet Boundary values
%      :num_quadr_pts_in_1d: number of quadrature points in 1 dimension
%      :              userf: any knowns physics resources as a function input for global residual 
%      :             userdf: any knowns physics resources as a function inout for consistent tangent 
%      :     max_iter_gmres: user specified maximum number of itertions before GMRES stops     
%      :        max_iter_nw: user specified maximum number newton steps before Newton stops
%      :     global_res_tol: user specified tolerance for norm of global_res before solver stops  
%
%output:              u: FEM Solution u
%           gl_res_norm: An array where each element contains the global residual norm after an iteration of GMRES
%               JACOB__: Jacobian of the system returned by the solver
%
% NOTE: return values that end with __ (e.g. JACOB__) are not intended to
% be used in the program per se. Be very careful before using them. 

    
    %get number of nodes in mesh
    num_nodes = msh.num_nodes;  
    
    %get the unknown u guess vector and global to local mapping for each u
    [u, global_idx_map] = get_global_map(num_nodes,dir_bndry_nodes,sz_u_field);
     dlta_u = u;
     dlta_u(2) =1;
    
    solverType = solver{1};
    
    iter=1;
    gl_res_norm = zeros(1);
    gl_res_norm_iter=0;
       
    if(strcmp(solverType,'gmres'))
        
        max_iter_gmres = solver{2};
        max_iter_nw = solver{3};
        tol = solver{4};
        global_res_tol = solver{5};
        
        while(true)
            
            % delme=get_global_Jv(dlta_u, global_idx_map, msh, num_quadr_pts_in_1d, userdf);
%             
             Je = get_global_jac(dlta_u, global_idx_map, msh,num_quadr_pts_in_1d,sz_u_field, userdf);
                    
            global_res = get_global_res(u, global_idx_map, msh, dir_bndry_val, num_quadr_pts_in_1d, userf); 
     
            global_res_norm = norm(global_res,inf);
        
            %store the norm of global_res after each iteration
            gl_res_norm(end+gl_res_norm_iter) = global_res_norm;
            gl_res_norm_iter = gl_res_norm_iter+1;
        
        
            if(global_res_norm < global_res_tol || iter > max_iter_nw )
                break;
            end           
        
            jac_fun = @(dlta_u)get_global_Jv(dlta_u, global_idx_map, msh, num_quadr_pts_in_1d, userdf);
            dlta_u= gmres(jac_fun, global_res,max_iter_gmres,tol);
            u=u-dlta_u;

            iter=iter+1;
        end
        if(size(solver,2) == 6)
            jac_flag = solver{6};
            if(jac_flag ==1)
                JACOB__ = zeros(size(u));
                for i=1:size(u,1)
                    Iden =eye(size(u,1));
                    du = Iden(:,i);
                    JACOB__(:,i)=get_global_Jv(du, global_idx_map, msh, num_quadr_pts_in_1d, userdf);
                end
            end
        else
            JACOB__ = 0;
        end
    end
    
    if(strcmp(solverType,'newton'))
        max_iter_nw = solver{2};
        tol = solver{3};
        global_res_tol = solver{4};
        
        while(true)
                            
            fun = @(u)get_global_res(u, global_idx_map, msh, dir_bndry_val, num_quadr_pts_in_1d, userf);        

            options = optimoptions(@fsolve,'Algorithm','trust-region-reflective', 'TolX',tol);%,'Jacobian','on');

            [u,~,~,~,JACOB__] = fsolve(fun, u,options); 

            global_res = get_global_res(u, global_idx_map, msh, dir_bndry_val, num_quadr_pts_in_1d, userf);

            global_res_norm = norm(global_res,inf);

            if(global_res_norm < global_res_tol || iter > max_iter_nw )
                    break;
            end
            
            iter=iter+1;
        end
        
    end
        
end