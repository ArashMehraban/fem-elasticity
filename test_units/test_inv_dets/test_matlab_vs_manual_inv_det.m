function [tf_dets, tf_inv] = test_matlab_vs_manual_inv_det(ith_elem_vtx_coords, D_hat,dim)

     [dets, ~, invJe, permuted_jacs] = get_elem_dirv(ith_elem_vtx_coords, D_hat, dim);
     [dets_mtlb, invJe_mtlb] = get_elem_jac(permuted_jacs,D_hat, dim);
     
     if(dim == 2)
        tmp=invJe(:,1);
        invJe(:,1) = invJe(:,4);
        invJe(:,4) = tmp;
        tmp=invJe(:,2);
        invJe(:,2) = invJe(:,3);
        invJe(:,3) = tmp;
     end
      
      %reshape vectorized inverses to 2x2 (or 3x3) matrices
      invJe_cell=cell(size(invJe,1),1);
      for k=1:size(invJe,1)
         invJe_cell{k} = reshape(invJe(k,:),sqrt(size(invJe,2)),[])';
      end
      %convert the inverse cells into one matrix
      invJe_mat = cell2mat(invJe_cell);
      diffMat = invJe_mtlb - invJe_mat;
      diff_dets = dets - dets_mtlb;

      if(norm(diffMat,2) < 1.0e-14)
          tf_inv = 'Pass';
      else
          tf_inv = norm(diffMat,2);
      end
      if(norm(diff_dets,2) < 1.0e-14)
          tf_dets = 'Pass';
      else
          tf_dets = norm(diffMat,2);
      end

end