function  [detsJe, invJe] = get_elem_jac(elem_vtx_coords, Ds)
%GET_ELEM_JAC returns a matrix of an element Jacobian for all quadrature points
%
%input: elem_vtx: element vertex coordinates
%     : Ds      : Ds: derivative of basis functions evaluated at quadrature points
%                 2D: Ds = D0  and D1 (xi and eta directions)
%                 3D: Ds = D0, D1 and D2 (xi, eta and zeta directions)
%
%output: invJe: A matrix of inverse Jacobian for all quadrature points in element
%             : size of invJe =  (num_derivatives*num_gs_pts) x num_derivatives        
%      : detJe: determinats of Jacobian for all element quadrature points

  D_sz = size(fieldnames(Ds),1);
  
  D0 = Ds.D0;
  D1 = Ds.D1;
  if(D_sz == 3)
     D2 = Ds.D2;
  end

  JD0 = D0*elem_vtx_coords;
  JD0=JD0(:);
      
  JD1 = D1*elem_vtx_coords;
  JD1=JD1(:);

  num_gs_pts = size(D0,1);

  invJe = zeros(D_sz*num_gs_pts,D_sz); 
  detsJe=zeros(num_gs_pts,1);
  
  %2D
   
  % structure of Je per quadrature point:
  % Je = [partial_x/partial_xi , partial_x/partial_eta]
  %      [partial_y/partial_xi , partial_y/partial_eta]
    
   if(D_sz == 2 )        
      idx = 2*(1:num_gs_pts)-1;
      for i=1:num_gs_pts
          tmpJe = [JD0(i:num_gs_pts:end),JD1(i:num_gs_pts:end)];
          invJe(idx(i):idx(i)+D_sz-1,:)=inv(tmpJe);
          detsJe(i)=det(tmpJe);
          if(detsJe(i) < 0)
              error('Defective element! Negative determinant in element Jacobian');
          end
      end       
   end
   
  %3D
  
  % structure of Je per quadrature point:
  % Je = [partial_x/partial_xi , partial_x/partial_eta, partial_x/partial_zeta]
  %      [partial_y/partial_xi , partial_y/partial_eta, partial_y/partial_zeta]
  %      [partial_z/partial_xi , partial_z/partial_eta, partial_z/partial_zeta]
  if(D_sz == 3)            
      JD2 = D2*elem_vtx_coords;
      JD2=JD2(:);
      
      idx = 3*(1:num_gs_pts)-2;
      for i=1:num_gs_pts
          tmpJe = [JD0(i:num_gs_pts:end),JD1(i:num_gs_pts:end),JD2(i:num_gs_pts:end)];
          invJe(idx(i):idx(i)+D_sz-1,:)=inv(tmpJe);
          detsJe(i)=det(tmpJe);
          if(detsJe(i) < 0)
              error('Defective element! Negative determinant in element Jacobian');
          end
      end

  end

      
end  


