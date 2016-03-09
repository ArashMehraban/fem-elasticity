function  [all_dets, invJe] = jacobian(elem_vtx_coords, D0, D1, varargin)
  %2D
  if(nargin == 3 )  
      %interleaving D0 and D1 
      D0D1 = interleaveMat(D0,D1) ;

      %Calculate the Jacobian for all Guass points 
      Je= D0D1*elem_vtx_coords;
      
      %get number of rows for D0D1
      nrow = size(D0D1,1);
  end
  
  %3D
  if(nargin == 4) 
      D2 = varargin{1};
      %interleaving D0, D1 and D2
      D0D1D2 = interleaveMat(D0,D1,D2) ;
      
      %Calculate the Jacobian for all Guass points
      Je= D0D1D2*elem_vtx_coords;
      
      %get number of rows for D0D1
      nrow = size(D0D1D2,1);
  end

      %get number of determinants = number of Guass points. 
      num_gs_pts = size(D0,1);      
      
      blocksz = nrow/num_gs_pts;

      %Allocate space for all determinants
      all_dets = zeros(1,num_gs_pts);
      %Populate the matrix of all determinanats
      %Populate the matrix of all Jacobian inverse
      inJe = cell(1,num_gs_pts);
      
      j=1;
      for i=1:blocksz:nrow
          temp = Je(i:i+blocksz-1,:);
          jdet = det(temp);
          if(jdet < 0)
              error('Defective element! Negative determinant in Jacobian');
          end
          all_dets(j) = jdet;
          invTemp = inv(temp);
          inJe{j} = invTemp;
          j=j+1;
      end
      invJe = cell2mat(inJe);
end  


