%L2projection: Q1-vs-Q2  problem using FEM
clear
clc
format short

%mesh refinement
%vector of errors
j=1;
errQ2=zeros(1,12);
hQ2 = zeros(1,12);
for i=6:5:26
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
    elm_type = 'Q2';
    
   
    % generate 2D mesh based on element type
    if (strcmp(elm_type,'Q1') || strcmp(elm_type,'Q2'))
        [conn,vtx_coords,bndry_nodes,bndry_elems] = create_mesh(elm_type, nx,x0,x1,ny, y0,y1);
    end
    % generate 3D mesh based on element type
    if (strcmp(elm_type,'Q1H') || strcmp(elm_type,'Q2H'))
        [conn,vtx_coords,bndry_nodes,bndry_elems] = create_mesh(elm_type, nx,x0,x1,ny, y0,y1,nz,z0,z1);
    end
    
    % number of unknowns/nodes on the mesh
    u_sz = size(vtx_coords,1);
    % solution vector
    u = zeros(1,u_sz)';
    
    % solution vector
     norms = zeros(2,10);
    
    global_res_norm=1;
    
   
    iter=1;
    norm_iter=1;
    while(global_res_norm > 1.0e-9  && iter <3 )
        
        [global_res, jac] = eval_res(u, conn,vtx_coords, elm_type);
                
        global_res_norm = norm(global_res);
        if(global_res_norm < 1.0e-9)
            break;
        end
        
         norms(1,norm_iter)= global_res_norm;
                 
        fun = @(u)eval_res(u, conn, vtx_coords, elm_type);
                
        options = optimset('Jacobian','on');
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
    
    errQ2(j) =error;
    
    hQ2(j) = hsz;
    j=j+1;
    
end


jj=1;
errQ1=zeros(1,12);
hQ1 = zeros(1,12);
for i=11:10:51
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
        [conn,vtx_coords,bndry_nodes,bndry_elems] = create_mesh(elm_type, nx,x0,x1,ny, y0,y1);
    end
    % generate 3D mesh based on element type
    if (strcmp(elm_type,'Q1H') || strcmp(elm_type,'Q2H'))
        [conn,vtx_coords,bndry_nodes,bndry_elems] = create_mesh(elm_type, nx,x0,x1,ny, y0,y1,nz,z0,z1);
    end
    
    % number of unknowns/nodes on the mesh
    u_sz = size(vtx_coords,1);
    % solution vector
    u = zeros(1,u_sz)';

    % solution vector
    norms = zeros(2,10);
    
    global_res_norm=1;
   
    iter=1;
    norm_iter=1;
    while(global_res_norm > 1.0e-9  && iter <3 )
        
        [global_res, jac] = eval_res(u, conn,vtx_coords, elm_type);
                
        global_res_norm = norm(global_res);
        if(global_res_norm < 1.0e-9)
            break;
        end
        
         norms(1,norm_iter)= global_res_norm;
                 
        fun = @(u)eval_res(u, conn, vtx_coords, elm_type);
                
        options = optimset('Jacobian','on');
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
    
    errQ1(jj) =error;
    
    hQ1(jj) = hsz;
    jj=jj+1;
    
end
figure
loglog(hQ2,errQ2,'r-o',hQ2,0.001*(hQ2.^3),'r--'  ,hQ1,errQ1,'b--*', hQ1,0.02*(hQ1.^2),'b:');
xlabel('h = sqrt(sum(element sides squared))')
ylabel('error')
legend('FEM-Q_2','O(h^3)', 'FEM-Q_1','O(h^2)','Location','northwest')
title('L_2-proj: 25-625 Q_2 VS 100-2500 Q_1 (2D) elements')
