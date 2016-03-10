%Solving the L2 projection problem using FEM
clear
clc
format short

%mesh refinement
%vector of errors
j=1;
err=zeros(1,12);
h = zeros(1,12);
for i=5:5:10
    % User given parameters to generate mesh with quadliteral elements
    nx=i;
    ny=i+1;
    nz=i+2; 
    x0 = -pi/2;
    y0 = -pi/2;
    z0 = -pi/2;
    x1 = pi/2;
    y1 = pi/2;
    z1 = pi/2;
    elm_type = 'Q4';
    
   
    % generate 2D mesh based on element type
    if (strcmp(elm_type,'Q4') || strcmp(elm_type,'Q9'))
        [conn,vtx_coords,bndry_nodes,bndry_elems] = create_mesh(elm_type, nx,x0,x1,ny, y0,y1);
    end
    % generate 3D mesh based on element type
    if (strcmp(elm_type,'Hex8') || strcmp(elm_type,'Hex27'))
        [conn,vtx_coords,bndry_nodes,bndry_elems] = create_mesh(elm_type, nx,x0,x1,ny, y0,y1,nz,z0,z1);
    end
    
    % number of unknowns/nodes on the mesh
    u_sz = size(vtx_coords,1);
    % solution vector
    u = zeros(1,u_sz)';
    
    global_res_norm=1;
    while(global_res_norm > 1.0e-4)
        global_res = eval_res(u, conn,vtx_coords, elm_type);
        global_res_norm = norm(global_res);
        if(global_res_norm < 1.0e-4)
            break;
        end
        sol = fsolve(@(u)eval_res(u ,conn, vtx_coords, elm_type  ),global_res);
        u=sol;
    end  
     
    fem_sol = sol;
    
    exactSol = exactf(vtx_coords);
        
    % element size (hsz)
    if (strcmp(elm_type,'Q4') || strcmp(elm_type,'Q9') )                         
        hsz = sqrt(((y1-y0)/(i+2))^2+(x1-x0)/i^2);
    end
    if (strcmp(elm_type,'Hex8') || strcmp(elm_type,'Hex27') )                         
        hsz = sqrt(((y1-y0)/(i+2))^2+(x1-x0)/i^2+(z1-z0)/i^2);
    end

    %calculate the error norm
    error = norm(exactSol - fem_sol)/norm(fem_sol);
    err(j) =error;
    h(j) = hsz;
    j=j+1;
    
end
figure
loglog(h,err,'r-o',h,0.001*h,'b:',h, 0.1*(h.^2),'b--');
xlabel('h = sqrt(sum(element sides squared))')
ylabel('error')
legend('FEM-Hex8','O(h)','O(h^2)','Location','northwest')