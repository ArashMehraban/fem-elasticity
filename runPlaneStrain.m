%Solving Poisson problem using FEM with unstructured mesh
clear
clc
format short

%mesh refinement
%vector of errors
j=1;
err=zeros(1,12);
h = zeros(1,12);
for i= 1:1 
    %disk9_4e
    filename = 'cylinder8';
    ext='exo';
   
    %get mesh information from exodus (.exo or .e) file
    msh = get_mesh(filename,ext,'lex');
    
    %get all Dirichlet boundary node_sets
    dir_bndry_nodes = get_all_dir_ns(msh);
    
    %get number of nodes in mesh
    num_nodes = msh.num_nodes;  
    
    %get the known u vector and global to local mapping foe each u
    [u, global_idx_map] = preproc(num_nodes,dir_bndry_nodes);
    
    %== Construct manufactured solution. Remove if dir_bndry_val given  ==%
    %get vertex coordinates from mesh                                     %
    vtx_coords = msh.vtx_coords;                                          %
    %get constructed dir_bndry_vals and exac Solutions on remaining nodes %
    [dir_bndry_val, exactSol] = get_exact_sol(vtx_coords,dir_bndry_nodes);%
    %=====================================================================%
        
    global_res_norm=1;
   
    iter=1;
    norm_iter=1;
    while((global_res_norm > 1.0e-9) && iter <3 )
        
        [global_res, jac] = eval_res(u, global_idx_map, msh, dir_bndry_val);
                                       
        global_res_norm = norm(global_res);
        if(global_res_norm < 1.0e-9)
            break;
        end         
        
        fun = @(u)eval_res(u, global_idx_map, msh, dir_bndry_val);
               
        options = optimoptions(@fsolve,'Algorithm','trust-region-reflective','Jacobian','on');
                
        u = fsolve(fun, global_res, options);        
        
        iter=iter+1;
    end  
      
    fem_sol = u;
       
    % element size (hsz)
    if (strcmp(elm_type,'Q1') || strcmp(elm_type,'Q2') )                         
        hsz = sqrt( ((x1-x0)/nx)^2 + ((y1-y0)/ny)^2);
    end
    if (strcmp(elm_type,'Q1H') || strcmp(elm_type,'Q2H') )                         
        hsz = sqrt(((x1-x0)/nx)^2 + ((y1-y0)/ny)^2 +((z1-z0)/nz)^2);
    end

    %calculate the error norm
    error = norm(exactSol - fem_sol)/norm(fem_sol);
    err(j) =error;
    h(j) = hsz;
    j=j+1;
    
end
figure
%Q1
loglog(h,err,'r-o', h,0.001*h,'r:',h, 0.01*(h.^2),'b--');
legend('FEM-Q1','O(h)','O(h^2)','Location','northwest')
title('Poisson: 25-4225 Q1 (2D) elements')

% % Q2
% loglog(h,err,'r-o', h,0.00001*(h.^3),'r:',h, 0.000001*(h.^2),'b--');
% legend('FEM-Q2','O(h^3)','O(h^2)','Location','northwest')
% title('Poisson: 100-2500 Q2 (2D) elements')

xlabel('h = sqrt(sum(element sides squared))')
ylabel('error')



