function  [detsJe, invJe] = jacobian(elem_vtx_coords, Ds)

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

  invJe= cell(1,num_gs_pts);
  detsJe=zeros(1,num_gs_pts);
  
  %2D
  if(D_sz == 2 )        
      for i=1:num_gs_pts
          invJe{i}=inv([JD0(i:num_gs_pts:end),JD1(i:num_gs_pts:end)]);
          detsJe(i)=det([JD0(i:num_gs_pts:end),JD1(i:num_gs_pts:end)]);
          if(detsJe(i) < 0)
              error('Defective element! Negative determinant in Jacobian');
          end
      end
       
  end
  
  %3D
  if(D_sz == 3) 
           
      JD2 = D2*elem_vtx_coords;
      JD2=JD2(:);
      
      for i=1:num_gs_pts
          jac_tmp = [JD0(i:num_gs_pts:end),JD1(i:num_gs_pts:end),JD2(i:num_gs_pts:end)];
          invJe{i}=inv(jac_tmp);
          detsJe(i)=det(jac_tmp);
          if(detsJe(i) < 0)
              error('Defective element! Negative determinant in Jacobian');
          end
      end

  end

      
end  


