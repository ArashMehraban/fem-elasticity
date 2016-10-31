function  [detsJe, D, invJe__, permuted_jacs__] = get_elem_dirv(elem_vtx_coords, D_hat, dim)
%GET_ELEM_DIRV returns a matrix of an element Jacobian for all quadrature
%points and determinants. It also returns the inverse Jacobian and
%permuted jacobian for each element. The latter two variables must not be
%used in code. They are designed for testing puporses.
%
%input: elem_vtx: element vertex coordinates
%     : D_hat   : derivative of basis functions evaluated at quadrature points
%                   2D: D_hat = [D0 ; D1]      (xi and eta directions)
%                   3D: D_hat = [D0 ; D1 ; D2] (xi, eta and zeta directions)
%     : dim     : problem dimension (2D or 3D)
%
%output: detsJe: determinants for all quadrature points in element
%      :      D: element derivative for all quadratupre points 
%      :invJe__: invrse of Jacobian per element in vector format:
%                2D:
%                invJe__ = 1/dets .*[D0x -D1x -D0y D1y]
%                3D:
%                 Jacobian is matrix with 9 columns:
%                          1   2   3   4   5   6   7   8   9 
%                     Je =[D0x D1x D2x D0y D1y D2y D0z D1z D2z]
%
%    invJe__=
%       invJe__(:,1) = (Je(:,9).*Je(:,5)-Je(:,8).*Je(:,6)).*invdetJe; % 9*5-8*6
%       invJe__(:,2) = (Je(:,7).*Je(:,6)-Je(:,9).*Je(:,4)).*invdetJe; % 7*6-9*4
%       invJe__(:,3) = (Je(:,8).*Je(:,4)-Je(:,5).*Je(:,7)).*invdetJe; % 8*4-5*7
%       invJe__(:,4) = (Je(:,3).*Je(:,8)-Je(:,9).*Je(:,2)).*invdetJe; % 3*8-9*2 
%       invJe__(:,5) = (Je(:,9).*Je(:,1)-Je(:,3).*Je(:,7)).*invdetJe; % 1*9-3*7
%       invJe__(:,6) = (Je(:,2).*Je(:,7)-Je(:,8).*Je(:,1)).*invdetJe; % 2*7-8*1
%       invJe__(:,7) = (Je(:,6).*Je(:,2)-Je(:,3).*Je(:,5)).*invdetJe; % 6*2-3*5 
%       invJe__(:,8) = (Je(:,3).*Je(:,4)-Je(:,6).*Je(:,1)).*invdetJe; % 3*4-6*1
%       invJe__(:,9) = (Je(:,5).*Je(:,1)-Je(:,2).*Je(:,4)).*invdetJe; % 5*1-2*4
%  
%   NOTE: return values that end with __ (e.g. invJe__) are not intended to
%   be used in the program per se. Be very careful before using them.  
% -------------------------------------------------------------------------------
%   2D:  
  %
  % The general structure of Je per quadrature point:
  % Je = [partial_x/partial_xi  , partial_y/partial_xi]
  %      [partial_x/partial_eta , partial_y/partial_eta]
  %
  %-------------------------------------------------------------------------------
  % 3D:
  %
  % The general structure of Je per quadrature point:
  % Je = [partial_x/partial_xi ,  partial_y/partial_xi,   partial_z/partial_xi  ]
  %      [partial_x/partial_eta , partial_y/partial_eta,  partial_z/partial_eta ]
  %      [partial_x/partial_zeta, partial_y/partial_zeta, partial_z/partial_zeta]
  %
  %
  %---------------------------------------------------------------------------------
  
  %2D:                 
  %  permuted_jac = [D0x  D0y]
  %                 [D1x  D1y]
  %
  %3D:              [D0x D0y D0z]  
  %  permuted_jac = [D1x D1y D1z]
  %                 [D2x D2y D2z]
  permuted_jacs__ = D_hat*elem_vtx_coords;   
  num_gs_pts = size(D_hat,1)/dim;
  num_col_D_hat = size(D_hat,2);
  
  %2D:                 
  %  Je = [D0x D1x D0y D1y]
  %
  %3D:
  %  Je =[D0x D1x D2x D0y D1y D2y D0z D1z D2z] 
  Je=reshape(permuted_jacs__,num_gs_pts,[]);
  
  i=1:size(Je,1);
  detsJe=zeros(size(Je,1),1);
 
  %2D
  if(dim == 2) 
      %   dets = [ D0x  .*  D1y  -  D1x  .*  D0y]
      detsJe(i)= Je(:,1).*Je(:,4)-Je(:,2).*Je(:,3);
      
      tf = isempty(detsJe(detsJe < 0));
      if(tf~=1)
         error('Defective element! Negative determinant in element Jacobian');   
      end
      invdetJe = 1./detsJe;
      
      %inverse of Jacbian (invJe)
      % invJe__ = 1/dets .*[D0x -D1x -D0y D1y]
      invJe__(:,1) =  Je(:,1).*invdetJe; %  1/dets .* D0x   
      invJe__(:,2) = -Je(:,2).*invdetJe; % -1/dets .* D1x     
      invJe__(:,3) = -Je(:,3).*invdetJe; % -1/dets .* D0y       
      invJe__(:,4) =  Je(:,4).*invdetJe; %  1/dets .* D1y   
      
      
      
      %D =[Dx]
      %   [Dy]
      %Dx = D1y.*D0 + D0y.*D1   = partial_N\partial_x
      %Dy = D1x.*D0 + D0x.*D1   = partial_N\partial_y
      D1yD0 =zeros(num_gs_pts,num_col_D_hat);
      D0yD1 =zeros(num_gs_pts,num_col_D_hat);
      D1xD0 =zeros(num_gs_pts,num_col_D_hat);      
      D0xD1 =zeros(num_gs_pts,num_col_D_hat);
    
      for j=1:num_col_D_hat
          %  D1yD0  =  D1y     .*      D0    
          D1yD0(:,j)=invJe__(:,4).*D_hat(1:num_gs_pts,j);
          %  D0yD1  =  D0y     .*      D1
          D0yD1(:,j)=invJe__(:,3).*D_hat(num_gs_pts+1:end,j);
          %  D1xD0  =  D1x     .*      D0
          D1xD0(:,j)=invJe__(:,2).*D_hat(1:num_gs_pts,j);
          %  D0xD1  =  D0x     .*      D1 
          D0xD1(:,j)=invJe__(:,1).*D_hat(num_gs_pts+1:end,j);
       end
    
    D=[D1yD0 + D0yD1 ; D1xD0 + D0xD1];      
  end
  
  if(dim == 3)
  %3D:    
  %       1   2   3   4   5   6   7   8   9 
  %  Je =[D0x D1x D2x D0y D1y D2y D0z D1z D2z] 
  
      % element determinants
      detsJe(i)=(Je(:,1).*Je(:,5).*Je(:,9) + ...    %D0x.*D1y.*D2z  1-5-9  
                 Je(:,6).*Je(:,8).*Je(:,3) + ...    %D2y.*D1z.*D2x  6-8-3
                 Je(:,7).*Je(:,2).*Je(:,6) ) - ...  %D0z.*D1x.*D2y  7-2-6
                (Je(:,3).*Je(:,5).*Je(:,7) + ...    %D2x.*D1y.*D0z  3-5-7
                 Je(:,6).*Je(:,8).*Je(:,1) + ...    %D2y.*D1z.*D0x  6-8-1
                 Je(:,9).*Je(:,2).*Je(:,4));        %D2z.*D1x.*D0y  9-2-4
             
      tf = isempty(detsJe(detsJe < 0));
      if(tf~=1)
         error('Defective element! Negative determinant in element Jacobian');   
      end
      invdetJe = 1./detsJe;  
      
      % inverse of Jacbian (invJe__)
      invJe__(:,1) = (Je(:,9).*Je(:,5)-Je(:,8).*Je(:,6)).*invdetJe; % 9*5-8*6
      invJe__(:,2) = (Je(:,7).*Je(:,6)-Je(:,9).*Je(:,4)).*invdetJe; % 7*6-9*4
      invJe__(:,3) = (Je(:,8).*Je(:,4)-Je(:,5).*Je(:,7)).*invdetJe; % 8*4-5*7
      invJe__(:,4) = (Je(:,3).*Je(:,8)-Je(:,9).*Je(:,2)).*invdetJe; % 3*8-9*2 
      invJe__(:,5) = (Je(:,9).*Je(:,1)-Je(:,3).*Je(:,7)).*invdetJe; % 1*9-3*7
      invJe__(:,6) = (Je(:,2).*Je(:,7)-Je(:,8).*Je(:,1)).*invdetJe; % 2*7-8*1
      invJe__(:,7) = (Je(:,6).*Je(:,2)-Je(:,3).*Je(:,5)).*invdetJe; % 6*2-3*5 
      invJe__(:,8) = (Je(:,3).*Je(:,4)-Je(:,6).*Je(:,1)).*invdetJe; % 3*4-6*1
      invJe__(:,9) = (Je(:,5).*Je(:,1)-Je(:,2).*Je(:,4)).*invdetJe; % 5*1-2*4
      
      %             invJe              D_hat
      %               |                  |     
      %      [invJe1 invJe2 invJe3]    [D0]
      % D =  [invJe4 invJe5 invJe6]  * [D1] 
      %      [invJe7 invJe8 invJe9]    [D2] 
      %   [Dx]  
      %D =[Dy] 
      %   [Dz]
      %Dx=invJe1.*D0+invJe2.*D1+invJe3.*D2 = partial_N\partial_x
      %Dy=invJe4.*D0+invJe5.*D1+invJe6.*D2 = partial_N\partial_y
      %Dz=invJe7.*D0+invJe8.*D1+invJe9.*D2 = partial_N\partial_z
      invJe1 =zeros(num_gs_pts,num_col_D_hat);
      invJe2 =zeros(num_gs_pts,num_col_D_hat);
      invJe3 =zeros(num_gs_pts,num_col_D_hat);
      invJe4 =zeros(num_gs_pts,num_col_D_hat);
      invJe5 =zeros(num_gs_pts,num_col_D_hat);
      invJe6 =zeros(num_gs_pts,num_col_D_hat);
      invJe7 =zeros(num_gs_pts,num_col_D_hat);
      invJe8 =zeros(num_gs_pts,num_col_D_hat);
      invJe9 =zeros(num_gs_pts,num_col_D_hat);
      
      for j=1:num_col_D_hat
          invJe1(:,j)=invJe__(:,1).*D_hat(1:num_gs_pts,j);
          invJe2(:,j)=invJe__(:,2).*D_hat(num_gs_pts+1:2*num_gs_pts,j);
          invJe3(:,j)=invJe__(:,3).*D_hat(2*num_gs_pts+1:3*num_gs_pts,j);
          invJe4(:,j)=invJe__(:,4).*D_hat(1:num_gs_pts,j);
          invJe5(:,j)=invJe__(:,5).*D_hat(num_gs_pts+1:2*num_gs_pts,j);
          invJe6(:,j)=invJe__(:,6).*D_hat(2*num_gs_pts+1:3*num_gs_pts,j);
          invJe7(:,j)=invJe__(:,7).*D_hat(1:num_gs_pts,j);
          invJe8(:,j)=invJe__(:,8).*D_hat(num_gs_pts+1:2*num_gs_pts,j);
          invJe9(:,j)=invJe__(:,9).*D_hat(2*num_gs_pts+1:3*num_gs_pts,j);
       end
    
    D=[invJe1+invJe2+invJe3;invJe4+invJe5+invJe6;invJe7+invJe8+invJe9];
  end
  


  
  


