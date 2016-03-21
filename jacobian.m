function  [detsJe, invJe] = jacobian(elem_vtx_coords, D0, D1, varargin)
%   %2D
  if(nargin == 3 )  
      JD0 = D0*elem_vtx_coords;
      JD0=JD0(:);
      
      JD1 = D1*elem_vtx_coords;
      JD1=JD1(:);
      
      num_gs_pts = size(D0,1);
      
      invJe= cell(1,num_gs_pts);
      detsJe=zeros(1,num_gs_pts);
      
      for i=1:num_gs_pts
          invJe{i}=inv([JD0(i:size(D0,1):end),JD1(i:size(D1,1):end)]);
          detsJe(i)=det([JD0(i:size(D0,1):end),JD1(i:size(D1,1):end)]);
          if(detsJe < 0)
              error('Defective element! Negative determinant in Jacobian');
          end
      end
       
  end
  
  %3D
  if(nargin == 4) 
      D2 = varargin{1};
    
      JD0 = D0*elem_vtx_coords;
      JD0=JD0(:);
      
      JD1 = D1*elem_vtx_coords;
      JD1=JD1(:);
      
      JD2 = D2*elem_vtx_coords;
      JD2=JD2(:);
      
      num_gs_pts = size(D0,1);
      
      invJe= cell(1,num_gs_pts);
      detsJe=zeros(1,num_gs_pts);
      
      for i=1:num_gs_pts
          invJe{i}=inv([JD0(i:num_gs_pts:end),JD1(i:num_gs_pts:end),JD2(i:num_gs_pts:end)]);
          detsJe(i)=det([JD0(i:num_gs_pts:end),JD1(i:num_gs_pts:end),JD2(i:num_gs_pts:end)]);
          if(detsJe < 0)
              error('Defective element! Negative determinant in Jacobian');
          end
      end

  end

      
end  


