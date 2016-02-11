%Solving the L2 projection problem using FEM
clear
clc
format short
%         o---o---o---o (x1, y1)
%         |   |   |   |
%         o---o---o---o  
%         |   |   |   |
% (x0,x0) o---o---o---o

%mesh refinement
%vector of errors
j=1;
err=zeros(1,12);
h = zeros(1,12);
for i=5:5:60
    % User given parameters to create 2D mesh
    num_pts_y=i+2;
    num_pts_x=i;
    x0 = -pi/2;
    y0 = -pi/2;
    x1 = pi/2;
    y1 = pi/2;
    
    %user defined F
    givenF =@(x,y) tanh(x).*exp(y)+sin(y);
                                   
    %                                   o---o 
    % size for a rectangular element:   |  /|
    %                                   | / |  hsize: size of hypotnuse 
    % sqrt(sum(element sides squared))  o---o
    hsize = sqrt(((y1-y0)/(i+2))^2+(x1-x0)/i^2);
    
    % Create a rectangular 2D mesh with quad elements
    [conn,vtx_coords,dm] = create2D(num_pts_x,num_pts_y,x0,y0,x1,y1);
    
     %evaluate the given function at vtx_coords
     x = vtx_coords(:,1);
     y = vtx_coords(:,2);
     
     %calculat the user defined f on all vtx_coords
     exact_f= givenF(x,y);

    % Assemble stiffness matrix and right hand side F
      [K, F] = assembly(conn,vtx_coords,givenF);

    % Solve for u (the solution) in Ku=F 
      sol = K\F;   

    %calculate the error norm
    error = norm(exact_f - sol)/norm(sol);
    err(j) =error;
    h(j) = hsize;
    j=j+1;
    
end
figure
loglog(h,err,'r-o',h,0.001*h,'b:',h, 0.1*(h.^2),'b--');
xlabel('log(h) where h is sqrt(sum(element sides squared))')
ylabel('log(error)')
legend('FEM','y=x','y=2x','Location','northwest')

