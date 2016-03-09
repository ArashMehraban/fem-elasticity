function [conn,vtx_coords,bndry_nodes,bndry_elems] = create_mesh(elm_type, nx,x_init,x_final, ny,y_init,y_final,varargin)
%CREATE_MESH can genrates 2D or 3D rectangular or cubical LEXOGRAPHICALLY ORDERed mesh based
%on the elm_type input. 
%
%2D:
%[conn,vtx_coords,bndry_nodes,bndry_elems] = CREATE_MESH(elm_type, nx,x_init,x_final, ny,y_init,y_final)
%input: elm_type: Q4 (4-noded element) or Q9 (9-noded element)
%       nx: number of mesh points in x-direction
%       x_init:  lower-left x-coordiante of a rectangular mesh
%       x_final: upper-right x-coordiante of a rectangular mesh
%       ny: number of mesh points in y-direction  
%       y_init:  lower-left y-coordiante of a rectangular mesh
%       y_final: upper-right y-coordiante of a rectangular mesh
%
%output: conn: connectivity matrix for all elements in mesh
%      : vtx_coords: (x,y) coordinates for all nodes in mesh
%      : bndry_nodes: mesh nodes on Boundray
%      : bndry_elems: elements on Boundary
%      
%      2D mesh elements ordering: LEXOGRAPHICAL
%  
%      o---------o
%      |3       4|
%      |         |     <- 1 Q4 element with LEXOGRAPHICAL ORDERING
%      |1       2|
%      o---------o
%  e.g   
%  [conn, vtx_coords,bndry_nodes,bndry_elems] = create_mesh('Q4',4,0,1.5,3,0,0.8)
%  produces the mesh below:
%
%     
%         (0,0.8)   (0.5,0.8)   (1,0.8)   (1.5,0.8)  
%y_final=0.8 o----------o----------o----------o              __
%            |9         |10        |11        |12              |
%            |  elm_4   |  elm_5   |  elm_6   |                |  
%         (0,0.4)   (0.5,0.4)   (1,0.4)   (1.5,0.4)            |     
%            o----------o----------o----------o               ny=3
%            |5         |6         |7         |8               |
%            |   elm_1  |  elm_2   |  elm_3   |                | 
%  y_init=0  |1         |2         |3         |4               |              
%  x_init=0  o----------o----------o----------o x_final=1.5  __|
%          (0,0)     (0.5,0)     (1,0)     (1.5,0)
%
%           |_____________nx=4________________|
%
% conn =
% 
%      1     2     5     6
%      2     3     6     7
%      3     4     7     8
%      5     6     9    10
%      6     7    10    11
%      7     8    11    12
% 
% 
% vtx_coords =
% 
%          0         0
%     0.5000         0
%     1.0000         0
%     1.5000         0
%          0    0.4000
%     0.5000    0.4000
%     1.0000    0.4000
%     1.5000    0.4000
%          0    0.8000
%     0.5000    0.8000
%     1.0000    0.8000
%     1.5000    0.8000
% 
% bndry_nodes =
% 
%      1     2     3     4     5     8     9    10    11    12
% 
% 
% bndry_elems =
% 
%      1     2     3     4     5     6
%
%                   * * *              * * *                      * * *              * * *
%3D:
%[conn,vtx_coords,bndry_nodes,bndry_elems] = CREATE_MESH(elm_type, nx,x_init,x_final, ny,y_init,y_final,nz,z_init,z_final)
%input: elm_type: Hex8 (8-noded element) or Hex27 (27-noded element)
%       nx: number of mesh points in x-direction
%       x_init:  lower-left x-coordiante of a cubical mesh
%       x_final: upper-right x-coordiante of a cubical mesh
%       ny: number of mesh points in y-direction 
%       y_init:  lower-left y-coordiante of a cubical mesh
%       y_final: upper-right y-coordiante of a cubical mesh
%       nz: number of mesh points in z-direction
%       z_init:  lower-bound z-coordiante of a cubical mesh
%       z_final: upper-bound z-coordiante of a cubical mesh
%
%output: conn: connectivity matrix for all elements in mesh
%      : vtx_coords: (x,y,z) coordinates for all nodes in mesh
%      : bndry_nodes: mesh nodes on Boundray
%      : bndry_elems: elements on Boundary
%
%       3D mesh elements ordering: LEXOGRAPHICAL
%
%          7 o-------------o 8
%           /|            /|
%          / |           / |
%         /  |          /  |
%        /   o---------/---o 
%       /   / 3       /  4/
%    5 /   /       6 /   / 
%     o-------------o   /   <- 1 Hex8 element with LEXOGRAPHICAL ORDERING           
%     |  /          |  / 
%     | /           | /    
%     |/1          2|/
%     o-------------o

% Minimal error checking for user input:
if ~(strcmp(elm_type,'Q4') || strcmp(elm_type,'Q9') || strcmp(elm_type,'Hex8') || strcmp(elm_type,'Hex27'))
    error('elm_type not supported. Choose Q4, Q9, Hex8 or Hex27.')
end

if( (strcmp(elm_type,'Hex8') || strcmp(elm_type,'Hex27')) && nargin <10)
        error('Hex8 or Hex27 require z direction entries.')
end

if(strcmp(elm_type,'Q4') || strcmp(elm_type,'Q9') || strcmp(elm_type,'Hex8') || strcmp(elm_type,'Hex27'))
   if(x_final <= x_init || y_final <= y_init) 
       error('mesh final point value cannot be smaller than initial point value')
   end
   if(nx < 2 || ny < 2)
      error('For 2D mesh, number of points on x and y direction must be at least 2.') 
   end
end

if(strcmp(elm_type,'Hex8') || strcmp(elm_type,'Hex27'))
   if(varargin{3} <= varargin{2}) 
       error('mesh final point value cannot be smaller than initial point value')
   end
   if(varargin{1} < 2)
      error('For 3D mesh, number of points on z direction must be at least 2.') 
   end
end

%-----------------------Q4 element----------------------------------------%

  if(strcmp(elm_type,'Q4')) 
      A=linspace(1,nx*ny,nx*ny);
      msh=(reshape(A',[nx,ny])');        
     
      xx = linspace(x_init,x_final,nx);
      yy = linspace(y_init,y_final,ny)';       
      vtx_x = repmat(xx,1,ny)';
      vtx_y = repmat(yy,1,nx)';
      vtx_y = vtx_y(:);  
    
      % (x,y) coords for each Node in the mesh
      vtx_coords = [vtx_x, vtx_y];

      %Allocating space for connectivity matrix
      %Ignore this: (m-1)*(i-1)+j <-index of each element
      conn=zeros((ny-1)*(nx-1),4);
      iter=1;
      for i=1:ny-1
          for j=1:nx-1
              %populate connectivity matrix lexographically
              conn(iter,:) = [msh(i,j),msh(i,j+1), msh(i+1,j),msh(i+1,j+1)];
              iter=iter+1;
          end
      end
  end
  
%-----------------------Q9 element----------------------------------------%
   
  if(strcmp(elm_type,'Q9'))
           
      xx = linspace(x_init,x_final,nx);
      ave_x = mean([xx(1:end-1);xx(2:end)]);
      xtemp = zeros(1,size(xx,2)+size(ave_x,2));
      xtemp(1:2:end) = xx;
      xtemp(2:2:end) = ave_x;
        
      yy = linspace(y_init,y_final,ny);       
      ave_y = mean([yy(1:end-1);yy(2:end)]);
      ytemp = zeros(1,size(yy,2)+size(ave_y,2));
      ytemp(1:2:end) = yy;
      ytemp(2:2:end) = ave_y;
                
      vtx_x = repmat(xtemp,1,2*ny-1)';
      vtx_y = repmat(ytemp',1,2*nx-1)';  
      vtx_y = vtx_y(:); 

      % (x,y) coords for each Node in the mesh
      vtx_coords = [vtx_x, vtx_y];
      
      nx_q9 = size(xtemp,2);
      ny_q9 = size(ytemp,2);
      %mesh nodes
      A=linspace(1,nx_q9*ny_q9,nx_q9*ny_q9);
      msh=(reshape(A',[nx_q9,ny_q9])'); 

      %Allocating space for connectivity matrix 
      conn=zeros((ny-1)*(nx-1),9);
      iter=1;
      for i=1:2:ny_q9-1
          for j=1:2:nx_q9-1
              %populate connectivity matrix lexographically
              conn(iter,:) = [msh(i,j),msh(i,j+1), msh(i,j+2), msh(i+1,j), msh(i+1,j+1) ...
                              msh(i+1,j+2),msh(i+2,j), msh(i+2,j+1), msh(i+2,j+2)];
               iter=iter+1;
          end
      end
  end
  
  if(strcmp(elm_type,'Q4') || strcmp(elm_type,'Q9'))

  
%===============================%
%  Boundary Nodes (Q4) or (Q9)  %
%===============================%
                
      bottom_nodes = msh(1,1:end);
      top_nodes = msh(end,1:end);
      right_nodes = msh(:,end);
      left_nodes = msh(:,1);
      
      bndry_nodes = unique([bottom_nodes,right_nodes',left_nodes',top_nodes]);
    
%==================================%
% Boundary elements (Q4) or (Q9)   %
%==================================%
      
      B=linspace(1,(nx-1)*(ny-1),(nx-1)*(ny-1));
      bn_elem = reshape(B',[nx-1,ny-1])';
      
      bottom_elms = bn_elem(1,1:end);
      top_elms = bn_elem(end,1:end);
      right_elms = bn_elem(:,end);
      left_elms = bn_elem(:,1);
      
      bndry_elems = unique([bottom_elms,right_elms',left_elms',top_elms]);     
  end 
 
  
%3D  elements:

%-----------------------Hex8 element--------------------------------------%  
  if(strcmp(elm_type,'Hex8'))
      
      nz = varargin{1};
      z_init = varargin{2};
      z_final = varargin{3};
      
      xx = linspace(x_init,x_final,nx);
      yy = linspace(y_init,y_final,ny)';       
      vtx_x = repmat(xx,1,ny)';
      vtx_y = repmat(yy,1,nx)';
      vtx_y = vtx_y(:); 
           
      zz = linspace(z_init,z_final,nz)';
      vtx_3x = repmat(vtx_x,1,nz);
      vtx_3x = vtx_3x(:);
      vtx_3y = repmat(vtx_y,1,nz);
      vtx_3y = vtx_3y(:);
      vtx_3z = repmat(zz,1,size(vtx_x,1))';
      vtx_3z = vtx_3z(:);
      
      % (x,y,z) coords for each Node in the mesh
      vtx_coords = [vtx_3x, vtx_3y, vtx_3z];
     
      %mesh nodes
      A=linspace(1,ny*nx*nz,ny*nx*nz);
      msh=(reshape(A,[nx,ny,nz])); 
      
      %Allocating space for connectivity matrix
      conn=zeros((nx-1)*(ny-1)*(nz-1),8);
      iter=1;
      for k=1:nz-1
          for j = 1:ny-1
              for i = 1:nx-1
                  %populate connectivity matrix lexographically
                  conn(iter,:) = [msh(i,j,k), msh(i+1,j,k), msh(i,j+1,k), msh(i+1,j+1,k), ...
                                  msh(i,j,k+1), msh(i+1,j,k+1),msh(i,j+1,k+1) , msh(i+1,j+1,k+1)];
                  iter=iter+1;
              end
          end
      end
  end
  
%-----------------------Hex27 element-------------------------------------%  

  if(strcmp(elm_type,'Hex27'))
      
      nz = varargin{1};
      z_init = varargin{2};
      z_final = varargin{3};
      
      xx = linspace(x_init,x_final,nx);
      ave_x = mean([xx(1:end-1);xx(2:end)]);
      xtemp = zeros(1,size(xx,2)+size(ave_x,2));
      xtemp(1:2:end) = xx;
      xtemp(2:2:end) = ave_x;
      
      
      yy = linspace(y_init,y_final,ny);       
      ave_y = mean([yy(1:end-1);yy(2:end)]);
      ytemp = zeros(1,size(yy,2)+size(ave_y,2));
      ytemp(1:2:end) = yy;
      ytemp(2:2:end) = ave_y;
      
      vtx_x = repmat(xtemp,1,2*ny-1)';
      vtx_y = repmat(ytemp',1,2*nx-1)';  
      vtx_y = vtx_y(:); 
      
      zz = linspace(z_init,z_final,nz);
      ave_z = mean([zz(1:end-1);zz(2:end)]);
      ztemp = zeros(1,size(zz,2)+size(ave_z,2));
      ztemp(1:2:end) = zz;
      ztemp(2:2:end) = ave_z; 
      
      vtx_3x = repmat(vtx_x,1,size(ztemp,2));
      vtx_3x = vtx_3x(:);
      vtx_3y = repmat(vtx_y,1,size(ztemp,2));
      vtx_3y = vtx_3y(:);      
          
      vtx_3z = repmat(ztemp',1,size(vtx_x,1))';
      vtx_3z = vtx_3z(:);
      
      % (x,y,z) coords for each Node in the mesh
      vtx_coords = [vtx_3x, vtx_3y, vtx_3z];
      
            
      nx_q27 = size(xtemp,2);
      ny_q27 = size(ytemp,2);
      nz_q27 = size(ztemp,2);
      
      %mesh nodes
      A=linspace(1,nx_q27*ny_q27*nz_q27,nx_q27*ny_q27*nz_q27);
      msh=(reshape(A',[nx_q27,ny_q27,nz_q27])); 
      
      %Allocating space for connectivity matrix
      conn=zeros((nx-1)*(ny-1)*(nz-1),27); 
      iter=1;
      for k=1:nz-1
          for j = 1:ny-1
              for i = 1:nx-1
                  %populate connectivity matrix lexographically
                  conn(iter,:) = [msh(i,j,k), msh(i+1,j,k), msh(i+2,j,k), ...
                                  msh(i,j+1,k), msh(i+1,j+1,k), msh(i+2,j+1,k), ...
                                  msh(i,j+2,k), msh(i+1,j+2,k), msh(i+2,j+2,k), ...
                                  msh(i,j,k+1), msh(i+1,j,k+1), msh(i+2,j,k+1), ...
                                  msh(i,j+1,k+1), msh(i+1,j+1,k+1), msh(i+2,j+1,k+1), ...
                                  msh(i,j+2,k+1), msh(i+1,j+2,k+1), msh(i+2,j+2,k+1), ...                                  
                                  msh(i,j,k+2), msh(i+1,j,k+2), msh(i+2,j,k+2), ...
                                  msh(i,j+1,k+2), msh(i+1,j+1,k+2), msh(i+2,j+1,k+2), ...
                                  msh(i,j+2,k+2), msh(i+1,j+2,k+2), msh(i+2,j+2,k+2)];
                  iter=iter+1;
              end
          end
      end
  end
  
  if(strcmp(elm_type,'Hex8') || strcmp(elm_type,'Hex27'))
      
%===================================%
%  Boundary Nodes (Hex8) or (Hex27) %
%===================================%  
      
      btn = msh(:,:,1);
      bottom_nodes = btn(:)';
      tn = msh(:,:,end);
      top_nodes = tn(:)';
      ln = msh(:,1,1:end);
      left_nodes = ln(:)';
      rn = msh(:,end,1:end);
      right_nodes = rn(:)';
      fn = msh(end,:,1:end);
      front_nodes = fn(:)';
      bkn = msh(1,:,1:end);
      back_nodes = bkn(:)';
      
      bndry_nodes=unique([bottom_nodes,top_nodes,left_nodes,right_nodes,front_nodes,back_nodes]);
      
%=====================================%
% Boundary elements (Hex8) or (Hex27) %
%=====================================%
      
      B=linspace(1,(nx-1)*(ny-1)*(nz-1),(nx-1)*(ny-1)*(nz-1));
      bn_elem = reshape(B',[nx-1,ny-1, nz-1]);

      bte = bn_elem(:,:,1);
      bottom_elms = bte(:)';
      te = bn_elem(:,:,end);
      top_elms = te(:)';
      fe = bn_elem(end,:,:);
      front_elms = fe(:)';
      bke = bn_elem(1,:,:);
      back_elms = bke(:)';
      re = bn_elem(:,end,:);
      right_elms = re(:)';
      le = bn_elem(:,1,:);
      left_elms = le(:)';
      
      bndry_elems = unique([bottom_elms,right_elms,left_elms,front_elms,back_elms,top_elms]);
  end
   
end


