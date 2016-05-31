function msh = create_mesh(elm_type, nx,x_init,x_final, ny,y_init,y_final,varargin)

%CREATE_MESH can genrates 2D or 3D rectangular or cubical LEXOGRAPHICALLY ORDERed mesh based
%on the elm_type input. 
%
%2D:
%[conn,vtx_coords,bndry_nodes,bndry_elems] = CREATE_MESH(elm_type, nx,x_init,x_final, ny,y_init,y_final)
%input: elm_type: QUAD4 (4-noded element) or QUAD9 (9-noded element)
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
%      |         |     <- 1 QUAD4 element with LEXOGRAPHICAL ORDERING
%      |1       2|
%      o---------o
%  e.g   
%  msh = create_mesh('QUAD4',4,0,1.5,3,0,0.8)
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
%msh = CREATE_MESH(elm_type, nx,x_init,x_final, ny,y_init,y_final,nz,z_init,z_final)
%input: elm_type: HEX8 (8-noded element) or HEX27 (27-noded element)
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
%     o-------------o   /   <- 1 HEX8 element with LEXOGRAPHICAL ORDERING           
%     |  /          |  / 
%     | /           | /    
%     |/1          2|/
%     o-------------o

% Minimal error checking for user input:
if ~(strcmp(elm_type,'QUAD4') || strcmp(elm_type,'QUAD9') || strcmp(elm_type,'HEX8') || strcmp(elm_type,'HEX27'))
    error('elm_type not supported. Choose QUAD4, QUAD9, HEX8 or HEX27.')
end

if( (strcmp(elm_type,'HEX8') || strcmp(elm_type,'HEX27')) && nargin <10)
        error('HEX8 or HEX27 require z direction entries.')
end

if(strcmp(elm_type,'QUAD4') || strcmp(elm_type,'QUAD9') || strcmp(elm_type,'HEX8') || strcmp(elm_type,'HEX27'))
   if(x_final <= x_init || y_final <= y_init) 
       error('mesh final point value cannot be smaller than initial point value')
   end
   if(nx < 2 || ny < 2)
      error('For 2D mesh, number of points on x and y direction must be at least 2.') 
   end
end

if(strcmp(elm_type,'HEX8') || strcmp(elm_type,'HEX27'))
   if(varargin{3} <= varargin{2}) 
       error('mesh final point value cannot be smaller than initial point value')
   end
   if(varargin{1} < 2)
      error('For 3D mesh, number of points on z direction must be at least 2.') 
   end
end

%-----------------------QUAD4 element----------------------------------------%

  if(strcmp(elm_type,'QUAD4')) 
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
  
%-----------------------QUAD9 element----------------------------------------%
   
  if(strcmp(elm_type,'QUAD9'))
           
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
      
      nx_QUAD9 = size(xtemp,2);
      ny_QUAD9 = size(ytemp,2);
      %mesh nodes
      A=linspace(1,nx_QUAD9*ny_QUAD9,nx_QUAD9*ny_QUAD9);
      msh=(reshape(A',[nx_QUAD9,ny_QUAD9])'); 

      %Allocating space for connectivity matrix 
      conn=zeros((ny-1)*(nx-1),9);
      iter=1;
      for i=1:2:ny_QUAD9-1
          for j=1:2:nx_QUAD9-1
              %populate connectivity matrix lexographically
              conn(iter,:) = [msh(i,j),msh(i,j+1), msh(i,j+2), msh(i+1,j), msh(i+1,j+1) ...
                              msh(i+1,j+2),msh(i+2,j), msh(i+2,j+1), msh(i+2,j+2)];
               iter=iter+1;
          end
      end
  end
  
  if(strcmp(elm_type,'QUAD4') || strcmp(elm_type,'QUAD9'))

  
%=====================================%
%  Boundary Nodes (QUAD4) or (QUAD9)  %
%=====================================%
                
      bottom_nodes = msh(1,1:end)';
      top_nodes = msh(end,1:end)';
      right_nodes = msh(:,end);
      left_nodes = msh(:,1);
       
      if(strcmp(elm_type,'QUAD4'))
          num_node_per_elem = 4;
      else
          num_node_per_elem = 9;
      end

      num_elem = size(conn,1);        
      num_nods = size(vtx_coords,1);
      num_dim = 2;
      sz_ns_idx = 4;      
      
      field_nm = {'vtx_coords','conn','num_elem','num_nodes_per_elem','num_nodes','num_dims','num_node_sets', 'node_ns1','node_ns2','node_ns3','node_ns4'};
      field_vl = {vtx_coords,conn,num_elem,num_node_per_elem,num_nods,num_dim,sz_ns_idx,bottom_nodes,top_nodes,right_nodes,left_nodes };
      
      %return mesh structure
      msh=struct();
      for i=1:size(field_nm,2)
          msh.(field_nm{i}) = field_vl{i};
      end
    
    
%========================================%
% Boundary elements (QUAD4) or (QUAD9)   %
%========================================%
      
      B=linspace(1,(nx-1)*(ny-1),(nx-1)*(ny-1));
      bn_elem = reshape(B',[nx-1,ny-1])';
      
      bottom_elms = bn_elem(1,1:end);
      top_elms = bn_elem(end,1:end);
      right_elms = bn_elem(:,end);
      left_elms = bn_elem(:,1);
      
      bndry_elems = unique([bottom_elms,top_elms, right_elms',left_elms',]);     
  end 
 
  
%3D  elements:

%-----------------------HEX8 element--------------------------------------%  
  if(strcmp(elm_type,'HEX8'))
      
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
  
%-----------------------HEX27 element-------------------------------------%  

  if(strcmp(elm_type,'HEX27'))
      
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
      
            
      nx_QUAD97 = size(xtemp,2);
      ny_QUAD97 = size(ytemp,2);
      nz_QUAD97 = size(ztemp,2);
      
      %mesh nodes
      A=linspace(1,nx_QUAD97*ny_QUAD97*nz_QUAD97,nx_QUAD97*ny_QUAD97*nz_QUAD97);
      msh=(reshape(A',[nx_QUAD97,ny_QUAD97,nz_QUAD97])); 
      
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
  
  if(strcmp(elm_type,'HEX8') || strcmp(elm_type,'HEX27'))
      
%===================================%
%  Boundary Nodes (HEX8) or (HEX27) %
%===================================%  
      
      btn = msh(:,:,1);
      bottom_nodes = btn(:);
      tn = msh(:,:,end);
      top_nodes = tn(:);
      ln = msh(:,1,1:end);
      left_nodes = ln(:);
      rn = msh(:,end,1:end);
      right_nodes = rn(:);
      fn = msh(end,:,1:end);
      front_nodes = fn(:);
      bkn = msh(1,:,1:end);
      back_nodes = bkn(:);
      
     
      if(strcmp(elm_type,'HEX8'))
          num_node_per_elem = 8;
      else
          num_node_per_elem = 27;
      end

      num_elem = size(conn,1);        
      num_nods = size(vtx_coords,1);
      num_dim = 3;
      sz_ns_idx = 6;      
      
      field_nm = {'vtx_coords','conn','num_elem','num_nodes_per_elem','num_nodes','num_dims','num_node_sets', 'node_ns1','node_ns2','node_ns3','node_ns4','node_ns5','node_ns6'};
      field_vl = {vtx_coords,conn,num_elem,num_node_per_elem,num_nods,num_dim,sz_ns_idx,bottom_nodes,top_nodes,right_nodes,left_nodes,front_nodes,back_nodes };
      
      %return mesh structure
      msh=struct();
      for i=1:size(field_nm,2)
          msh.(field_nm{i}) = field_vl{i};
      end
      
%=====================================%
% Boundary elements (HEX8) or (HEX27) %
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
      
      bndry_elems = unique([bottom_elms,top_elms,right_elms,left_elms,front_elms,back_elms,]);
  end
   
end


