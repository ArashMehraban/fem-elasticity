function [u,global_res_norm,JACOB__] = get_fem_sol(vtk_dest_folder, vtk_filename, steps, origConn, msh, sz_u_field, dir_bndry_nodes, dir_bndry_val, num_quadr_pts_in_1d,userf,userdf,solver)
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
       
    solverType = solver{1};
    
    iter=1;
    gl_res_norm = zeros(1);
    gl_res_norm_iter=0;
       
    if(strcmp(solverType,'gmres'))
        
        if(solver{2}>size(u,1))
            max_iter_gmres = size(u,1);            
            gmresMsg = strcat('gmres_max_iter was reduced to row size of unknown matrix:', num2str(size(u,1)));
            warning(gmresMsg);
        else
            max_iter_gmres = solver{2};
        end
        max_iter_nw = solver{3};
        tol = solver{4};
        global_res_tol = solver{5};
        
        
       vtk_files = cell(steps,1);
       for s=1:steps+1
          vtk_files{s} = strcat(vtk_filename , num2str(s), '.vtk');    
       end
       zeroBd= dir_bndry_val;
       for b=1:size(dir_bndry_val,1)
           zeroBd{b} = 0* dir_bndry_val{b};
       end
       u_closure =  get_closure_u(u,dir_bndry_nodes,zeroBd,global_idx_map); 
       %graph at no boundary applied
       processVtk(vtk_files{1},origConn, msh,u_closure,sz_u_field);
       
       dt = linspace(0,1,steps+1);
        
 for m=1:steps 
     
     stepBndry = dir_bndry_val;     
     for b=1:size(dir_bndry_val,1)
           stepBndry{b} = dt(m+1)* dir_bndry_val{b};
     end
     
     
        tic
        global_jac = get_global_jac(global_idx_map, msh,num_quadr_pts_in_1d,sz_u_field, userdf);            
        elapsedTime=toc;
        msg = strcat('global Jacobian matrix assembly:', num2str(elapsedTime), ' seconds');
        disp(msg);
        
%         tic
%         [L,U] = lu(global_jac);
%         elapsedTime = toc;
%         msg = strcat('LU of global Jacobian matrix: ', num2str(elapsedTime), ' seconds');
%         disp(msg);
        
        
        while(true)
            
            global_res = get_global_res(u, global_idx_map, msh, stepBndry, num_quadr_pts_in_1d, userf); 
     
            global_res_norm = norm(global_res,inf);
        
            %store the norm of global_res after each iteration
            gl_res_norm(end+gl_res_norm_iter) = global_res_norm;
            gl_res_norm_iter = gl_res_norm_iter+1;        
        
            if(global_res_norm < global_res_tol || iter > max_iter_nw )
               break;
            end   
            
            jv_fun = @(dlta_u)get_global_Jv(dlta_u, global_idx_map, msh, num_quadr_pts_in_1d, userdf);
            %precond_jac_fun =@(u)(U\(L\u));  
            precond_jac_fun = @(u)(global_jac\u);
            tic
            dlta_u= gmres(jv_fun, global_res,max_iter_gmres,tol,[],precond_jac_fun);
            elapsedTime = toc;
            msg = strcat('gmres solve time: ', num2str(elapsedTime), ' seconds');
            disp(msg);
            u=u-dlta_u;
            
            iter=iter+1;
        end
        u_closure =  get_closure_u(u,dir_bndry_nodes,stepBndry,global_idx_map);       
        processVtk(vtk_files{m+1},origConn, msh,u_closure,sz_u_field);
 end
 addpath(fullfile(pwd,vtk_dest_folder));
 for n=1:steps+1
    movefile(vtk_files{n}, vtk_dest_folder);
 end
 rmpath(fullfile(pwd,vtk_dest_folder));
 
 
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
                            
            fun = @(u)get_global_res(u, global_idx_map, msh, stepBndry, num_quadr_pts_in_1d, userf);        

            options = optimoptions(@fsolve,'Algorithm','trust-region-reflective', 'TolX',tol);%,'Jacobian','on');

            [u,~,~,~,JACOB__] = fsolve(fun, u,options); 

            global_res = get_global_res(u, global_idx_map, msh, stepBndry, num_quadr_pts_in_1d, userf);

            global_res_norm = norm(global_res,inf);

            if(global_res_norm < global_res_tol || iter > max_iter_nw )
                    break;
            end
            
            iter=iter+1;
        end
        
    end
        
end