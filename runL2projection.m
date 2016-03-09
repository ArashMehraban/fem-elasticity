%Solving the L2 projection problem using FEM
clear
clc
format short

%mesh refinement
%vector of errors
j=1;
err=zeros(1,12);
h = zeros(1,12);
for i=5:5:5
    % User given parameters to create 2D mesh with quadliteral elements
    nx=i;
    ny=i+1;
    nz=i+2; 
    x0 = 0;
    y0 = 0;
    z0 = 0;
    x1 = 1;
    y1 = 1;
    z1 = 1;
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
      
    % get global residual from FEM
    global_res = eval_res(conn,vtx_coords, u, elm_type);
      
      %while(global_res > 1.0e-4)
       % sol = fsolve(@eval_res,global_res);
      %end      
        
    % element size (hsz)
    if (strcmp(elm_type,'Q4') || strcmp(elm_type,'Q9') )                         
        hsz = sqrt(((y1-y0)/(i+2))^2+(x1-x0)/i^2);
    end
    if (strcmp(elm_type,'Q4') || strcmp(elm_type,'Q9') )                         
        hsz = sqrt(((y1-y0)/(i+2))^2+(x1-x0)/i^2+(z1-z0)/i^2);
    end

%     %calculate the error norm
%     error = norm(exact_f - sol)/norm(sol);
%     err(j) =error;
%     h(j) = hsize;
%     j=j+1;
    
end
% figure
% loglog(h,err,'r-o',h,0.001*h,'b:',h, 0.1*(h.^2),'b--');
% xlabel('h = sqrt(sum(element sides squared))')
% ylabel('error')
% legend('FEM','O(h)','O(h^2)','Location','northwest')