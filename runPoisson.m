%Solving Poisson problem using FEM
clear
clc
format short

%mesh refinement
%vector of errors
j=1;
err=zeros(1,12);
h = zeros(1,12);
for i= 4:4:4%6:5:11 %11:10:51 %66
    % User given parameters to generate mesh with quadliteral elements
    nx=i;
    ny=i;
    nz=i+2; 
%     x0 = 0;
%     y0 = 0;
%     z0 = -pi/2;
%     x1 = 8;
%     y1 = 8;
%     z1 = pi/2;
    x0 = -pi/2;
    y0 = -pi/2;
    z0 = -pi/2;
    x1 = pi/2;
    y1 = pi/2;
    z1 = pi/2;
    elm_type = 'Q1';
    
    % Array of of Number of dof(s) per Unknown Field:
      u_field = [2,3];   
        
    % generate 2D mesh based on element type
    if (strcmp(elm_type,'Q1') || strcmp(elm_type,'Q2'))
        [conn,vtx_coords,bndry_nodes,bndry_elems,geom] = create_mesh(elm_type, nx,x0,x1,ny, y0,y1); %l,r,t,b
    end
    % generate 3D mesh based on element type
    if (strcmp(elm_type,'Q1H') || strcmp(elm_type,'Q2H'))
        [conn,vtx_coords,bndry_nodes,bndry_elems] = create_mesh(elm_type, nx,x0,x1,ny, y0,y1,nz,z0,z1);
    end
        
%     dir_bndry_nodes = unique([b,l']);
%     nuemann_bndry_nodes = unique([t,r']);
    
    r = geom.right(2:end);
    t = geom.top(2:end);
    l = geom.left;
    b = geom.bottom;
    dir_bndry_u11 = unique([t,r]); %rt_top_nodes
    dir_bndry_u23 = dir_bndry_u11;
    dir_bndry_u12 = unique([b,l]); %bt_lft_nodes    
    dir_bndry_u21 = dir_bndry_u12;
    dir_bndry_u22 = dir_bndry_u12;
    
    dir_bndry_nodes{1} = dir_bndry_u11';
    dir_bndry_nodes{2} = dir_bndry_u12';
    dir_bndry_nodes{3} = dir_bndry_u21';
    dir_bndry_nodes{4} = dir_bndry_u22';
    dir_bndry_nodes{5} = dir_bndry_u23';
    
    [u, u_section] = preproc(size(vtx_coords,1),dir_bndry_nodes);
    
    dir_bndry_nodes = unique([b,l,t,r]);
         
    [dir_bndry_val, exactSol, all_nodes_exact] = get_exact_sol(vtx_coords,dir_bndry_nodes);
    %dir_bndry values on bt_lft_nodes
    dir_u12 = all_nodes_exact(dir_bndry_u12,2);
    dir_u21 = all_nodes_exact(dir_bndry_u12,3);
    dir_u22 = all_nodes_exact(dir_bndry_u12,4);
    
    %dir_bndry values on rt_top_nodes
    dir_u11 = all_nodes_exact(dir_bndry_u11,1);
    dir_u23 = all_nodes_exact(dir_bndry_u11,5);
    
    
    
    global_res_norm=1;
   
    iter=1;
    norm_iter=1;
    while((global_res_norm > 1.0e-9) && iter <3 )
        
        [global_res, jac] = eval_res(u, conn, vtx_coords, elm_type, u_section, bndry_nodes, dir_bndry_nodes,dir_bndry_val);
                                       
        global_res_norm = norm(global_res);
        if(global_res_norm < 1.0e-9)
            break;
        end         
        
        fun = @(u)eval_res(u, conn, vtx_coords, elm_type, u_section, bndry_nodes, dir_bndry_nodes, dir_bndry_val);
               
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


