function  [all_dets, Je_all_gauss_pts] = jacobian(elem_vtx_coords, D1, D0)

  %interleaving D1 and D0 
  nColumns = size(D1,2);
  D1D0 = [D1,D0]'; 
  D1D0 = reshape(D1D0(:),nColumns,[])'; 
  
  %Calculate the Jacobian for all Guass points by multiplying
  %each element_vtx_coords by derivatives wrt D1, D0 
  Je_all_gauss_pts= D1D0*elem_vtx_coords;
  
  %get the number of Guass points = number of determinants
  num_gs_pts = sqrt(size(D0,1));
  
  
  %Allocate space for all determinants
  all_dets = zeros(1,num_gs_pts);
  %Populate the matrix of all determinanats
  j=1;
  for i=1:num_gs_pts:size(D1D0,1)
      jdet = det(Je_all_gauss_pts(i:i+num_gs_pts-1,:));
      all_dets(j) = jdet;
      j=j+1;
  end

end

