function [dets_mtlb, invJe_mtlb] = get_elem_jac(permuted_jacs,D_hat, dim)
    
     num_gs_pts = size(D_hat,1)/dim;
     [r,c]=size(permuted_jacs);  
     invJe_mtlb = zeros(r,c);  

     dets_mtlb=zeros(num_gs_pts,1);
     idx = dim*(1:num_gs_pts)-(dim-1);  

     for i=1:num_gs_pts
         tempJe = permuted_jacs(i:num_gs_pts:end,:);
         invJe_mtlb(idx(i):idx(i)+dim-1,:)=inv(tempJe);
         dets_mtlb(i)=det(tempJe);
             if(dets_mtlb(i) < 0)
                 error('Defective element! Negative determinant in element Jacobian');
             end
     end

end