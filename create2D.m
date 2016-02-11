function [conn,vtx_coords,dispMsh] = create2D(m,n,x_init,y_init,x_final,y_final)
%      Developed by Arash Mehraban January, 2016
%
%input: m: number of nodes on x-direction of the mesh
%       n: number of nodes on y-direction of the mesh 
%       x_init:  lower-left x-coordiante of a rectangular mesh
%       y_init:  lower-left y-coordiante of a rectangular mesh
%       x_final: upper-right x-coordiante of a rectangular mesh
%       y_final: upper-right y-coordiante of a rectangular mesh
%
%output: LEXOGRAPHICAL Ordered Quad Mesh 
%      : conn: connectivity matrix for all elements in the mesh
%      : vtx_coords: x & y coordinates for each vertext of each element
%      : dispMsh: display the mesh created (not used in calculations)
%  
%      o---------o
%      |3       4|
%      |         |     <- 1 element with LEXOGRAPHICAL ORDERING
%      |1       2|
%      o---------o
%  e.g   
%  [conn, vtx_coords, dispMsh] = create2D(4,3,0,0,1.5,0.8)
%  produces the mesh below:
%
%         y=0.8   (0.5,0.8)    (1,0.8)   (1.5,0.8)  
%    (0,0.8)o----------o----------o----------o        
%           |9         |10        |11        |12
%           |          |          |          |
%           |      (0.5,0.4)   (1,0.4)   (1.5,0.4)    m=4, n=3     
%    (0,0.4)o----------o----------o----------o         
%           |5         |6         |7         |8
%           |          |          |          |
%           |1         |2         |3         |4                      
%           o----------o----------o----------o x=1.5
%      (0,0)    (0,0.5)      (0,1)      (0,1.5)
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
% 
% dispMsh =
% 
%      9    10    11    12
%      5     6     7     8
%      1     2     3     4

    A=linspace(1,n*m,n*m);
    msh =(reshape(A',[m,n])');
     
    %return vertex coordinate based on given x and y
    xx = linspace(x_init,x_final,m);
    yy = linspace(y_init,y_final,n)';
    vtx_x = repmat(xx,1,n);
    vtx_y = repmat(yy,1,m)';
    vtx_y = vtx_y(:);
              % Node x-coords y-coords
    vtx_coords = [vtx_x', vtx_y];

    %Allocating space for connectivity matrix
    %Ignore this: (m-1)*(i-1)+j 
    conn=zeros((n-1)*(m-1),4);
    k=1;
    for i=1:n-1
        for j=1:m-1
            %populate connectivity matrix
            conn(k,:) = [msh(i,j),msh(i,j+1),msh(i+1,j),msh(i+1,j+1)];
            k=k+1;
        end
    end
    
    % Displays the upside down grid points as normally drawn in the x-y coordinates
    dispMsh = flipud(msh);
    
end


