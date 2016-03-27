%Solving the L2 projection problem using FEM
clear
clc
format short

%mesh refinement
%vector of errors
j=1;
err=zeros(1,12);
fd_err=zeros(1,12);
h = zeros(1,12);
for i=6:5:11
    % User given parameters to generate mesh with quadliteral elements
    nx=i;
    ny=i;
    nz=i; 
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
    
    fdu_sz = size(vtx_coords,1);
    % solution vector
    fdu = zeros(1,fdu_sz)';
    
    norms = zeros(2,10);
    
    global_res_norm=1;
    fdglobal_res_norm=1;
   
    iter=1;
    norm_iter=1;
    while((global_res_norm > 1.0e-9 || fdglobal_res_norm > 1.0e-9) && iter <3 )
        
        [global_res, jac] = eval_res(u, conn,vtx_coords, elm_type);
        [fdglobal_res, fdjac] = fdeval_res(fdu, conn,vtx_coords, elm_type);
        
        
        gl_res_diff =  fdglobal_res - global_res;       
        jac_diff = fdjac - jac;
                
        global_res_norm = norm(global_res);
        if(global_res_norm < 1.0e-9)
            fdglobal_res_norm = norm(fdglobal_res);
            break;
        end
        
        fdglobal_res_norm = norm(fdglobal_res);
        if(fdglobal_res_norm < 1.0e-9)
            break;
        end  
         norms(1,norm_iter)= global_res_norm;
         norms(2,norm_iter)= fdglobal_res_norm;

        
        fun = @(u)eval_res(u, conn, vtx_coords, elm_type);
        fdfun = @(fdu)fdeval_res(fdu, conn, vtx_coords, elm_type);
        
        options = optimoptions(@fsolve,'Algorithm','trust-region-reflective','Jacobian','on');
        %options = optimset('Jacobian','on');
        sol = fsolve(fun, global_res, options);
        
        fdoptions = optimset('Jacobian','on');
        fdsol =fsolve(fdfun, fdglobal_res, fdoptions);
        u=sol;
        fdu=fdsol;  
        
        udiff = fdsol-sol;
        norm_iter=norm_iter+1;
        iter=iter+1;
    end  
      
    fem_sol = sol;
    fd_fem_sol = fdsol;
    
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
    fd_error = norm(exactSol - fd_fem_sol)/norm(fd_fem_sol);
    err(j) =error;
    fd_err(j)=fd_error;
    h(j) = hsz;
    j=j+1;
    
end
figure
loglog(h,err,'r-o',h,fd_err,'b--*',h, 0.02*(h.^2),'b--');
xlabel('h = sqrt(sum(element sides squared))')
ylabel('error')
legend('FEM-Q_2-w/jac','FEM-Q_2-w/fd-jac','O(h^2)','Location','northwest')
title('L_2 Proj: 25-100 Q_2 (2D) elements')
