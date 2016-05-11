%Solving the L2 projection problem using FEM
clear
clc
format short

%mesh refinement
%vector of errors
j=1;
err=zeros(1,12);
h = zeros(1,12);
for i=6:5:11 %6:5:66  %66
    % User given parameters to generate mesh with quadliteral elements
    nx=i;
    ny=i;
    nz=i+2; 
    x0 = -pi/2;
    y0 = -pi/2;
    z0 = -pi/2;
    x1 = pi/2;
    y1 = pi/2;
    z1 = pi/2;
    elm_type = 'Q1';
    
  
    % generate 2D mesh based on element type
    if (strcmp(elm_type,'Q1') || strcmp(elm_type,'Q2'))
        [conn,vtx_coords,~,~] = create_mesh(elm_type, nx,x0,x1,ny, y0,y1);
    end
    % generate 3D mesh based on element type
    if (strcmp(elm_type,'Q1H') || strcmp(elm_type,'Q2H'))
        [conn,vtx_coords,~,~] = create_mesh(elm_type, nx,x0,x1,ny, y0,y1,nz,z0,z1);
    end
    
    [conn,vtx_coords] = applybd2mesh(conn,vtx_coords,-1,-1);
    
    % number of unknowns/nodes on the mesh
    u_sz = size(vtx_coords,1);
    % solution vector
    u = zeros(1,u_sz)';
     
    norms = zeros(1,12);
    
    global_res_norm=1;
   
    iter=1;
    norm_iter=1;
    while((global_res_norm > 1.0e-9 || fdglobal_res_norm > 1.0e-9) && iter <3 )
        
        [global_res, jac] = eval_res(u, conn,vtx_coords, elm_type,1,1);
                
        global_res_norm = norm(global_res);
        if(global_res_norm < 1.0e-9)
            break;
        end

        norms(norm_iter)= global_res_norm;         
        
        fun = @(u)eval_res(u, conn, vtx_coords, elm_type,1,1);
               
        options = optimoptions(@fsolve,'Algorithm','trust-region-reflective','Jacobian','on');
        %options = optimset('Jacobian','on');
        sol = fsolve(fun, global_res, options);        
        u=sol;
        
        norm_iter=norm_iter+1;
        iter=iter+1;
    end  
      
    fem_sol = sol;
    
    exactSol = exactf(vtx_coords);
        
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
% Q1
loglog(h,err,'r-o', h,0.001*h,'r:',h, 0.01*(h.^2),'b--');
legend('FEM-Q1','O(h)','O(h^2)','Location','northwest')
title('L2-proj: 25-4225 Q1 (2D) elements')

% % Q2
% loglog(h,err,'r-o', h,0.001*(h.^3),'r:',h, 0.1*(h.^2),'b--');
% legend('FEM-Q2','O(h^3)','O(h^2)','Location','northwest')
% title('L2-proj: 25-4225 Q2 (2D) elements')

xlabel('h = sqrt(sum(rectangular element sides squared))')
ylabel('error')


