function [global_res, jac] = get_global_res(u, global_idx_map, msh, dir_bndry_val, problem_type)
% GET_GLOBAL_RES evaluates the global residual and the consistent tangent
%  input:              u: vector of unknowns 
%       : global_idx_map: global map of local u's
%       :            msh: mesh object (see get_mesh function)
%       :  dir_bndry_val: Dirchlet boundary values if any
%       :          userf: user supplied code for get_userf function
%
% output: global_res: global residual
%       :        jac: consistant tangent

     %get the number of elements
     num_elem = msh.num_elem; 
     
     %get connectivity matrix
     conn = msh.conn;
     
     %get vertex coordinates
     vtx_coords = msh.vtx_coords;
     
     %unknown size (size of dofs per unknown u)
     unknown_sz = size(u,1);
     
     %Allocate space for global residual for unkowns     
     global_res =zeros(unknown_sz,1);
     
     %get all dirichlet boundary node_sets
     dir_bndry_nodes = get_all_dir_ns(msh);
      
     global_u =  get_global_u(u,dir_bndry_nodes,dir_bndry_val,global_idx_map);     
     
     %get size of dofs per node
     sz_global_idx_map = size(global_idx_map,2);
     
     %get element type    
     num_quadr_pts = msh.num_nodes_per_elem;
     
     %get Weights, Basis (B) functions and their Derivatives (D0, D1 and D2)
     % D = partial_B/partial_x_i
     [B, Ds, W_hat] = get_shape(num_quadr_pts);
        
     %number of element derivative matrices
     num_elem_dirv = size(fieldnames(Ds),1);
     
     %Allocate space for globall Jacobian
     jac = zeros(unknown_sz,unknown_sz);
     
     %index for global Jacobian
     %row indices
     c = num_quadr_pts*(1:sz_global_idx_map)- (num_quadr_pts-1);
     r = (num_quadr_pts*(1:sz_global_idx_map));
     cidx = repmat(c,1,sz_global_idx_map)';
     ridx = repmat(r,1,sz_global_idx_map)';
     %col indices
     rr = (num_quadr_pts*(1:sz_global_idx_map)- (num_quadr_pts-1))';
     cc = (num_quadr_pts*(1:sz_global_idx_map))';
     jcidx = repmat(rr,1,sz_global_idx_map)';
     jcidx = jcidx(:);
     jridx = repmat(cc,1,sz_global_idx_map)';
     jridx = jridx(:);
     idx=[jcidx, jridx, cidx, ridx];     
               
     for i=1:num_elem
         
         %get number of dof per unknown u_i per element
         neldof_per_ui = size(conn(i,:),2);
         
         %get total number of dof per element
         neldof = neldof_per_ui*sz_global_idx_map;
          
         %get corresponding vertex coordinates for each element 
         element_vtx_coords = vtx_coords(conn(i,:),:);
         
         %get corresponding unknown/solution u for each element
         elem_u = global_u(conn(i,:),:);   
         
         %get mapping constituents from element jacobian 
         [dets, invJe] = get_elem_jac(element_vtx_coords, Ds);
         Di = get_elem_dirv(invJe, Ds); 
          
         % grad_ue structure:
         % 2D: [Di_1*u1, Di_2*u_1 |  Di_1*u2, Di_2*u_2]
         % 3D: [Di_1*u1, Di_2*u_1, Di_3*u1 | Di_1*u2, Di_2*u_2, Di_3*u2  | Di_1*u3, Di_2*u_3, Di_3*u3]
         grad_ue = zeros(size(elem_u,1),num_elem_dirv*size(elem_u,2));
         for  j=1:num_elem_dirv
             grad_ue(:,j:num_elem_dirv:end) = Di((j-1)*num_quadr_pts+1:j*num_quadr_pts,:)*elem_u;
         end
         
         % ue structure: 
         % 2D: [B*u1, B*u2]
         % 3D: [B*u1, B*u2, B*u3]
         ue = B*elem_u;             
 
         % mapped quadrature points to reference coordinate system
         mp_qd_pts= B*element_vtx_coords;
         
        [f0,f1,f00, f01, f10, f11] = get_userf(ue, grad_ue, mp_qd_pts, problem_type); 
                 
         %get Gauss Weights for the current element
         W = W_hat.*dets;
         
         
         wf1 = zeros(size(f1));
         for j=1:size(wf1,2);
             wf1(:,j) = W.*f1(:,j);
         end
         
         D_res = zeros(size(f1,1),sz_global_idx_map); 
         for j=1:num_elem_dirv 
             D_res = D_res + Di((j-1)*num_quadr_pts+1:j*num_quadr_pts,:)'*wf1(:,j:num_elem_dirv:end);
         end

         
         wf0 = zeros(size(f0));
         for j=1:size(f0,2)
             wf0(:,j) = W.*f0(:,j);
         end 
         
         % element residual evaluation
         res_e = B'*wf0 + D_res;
         
         
         %jac_e constituents:
         
%         f0u = B'*diag(W.*f00)*B;%
           
         f0u=zeros(neldof,neldof);
         for j=1:size(idx,1)
             f0u(idx(j,1):idx(j,2),idx(j,3):idx(j,4)) = B'*diag(W.*f00(:,j))*B;
         end
          
%          f01TD =zeros(num_quadr_pts,num_quadr_pts);
%          for j=1:num_elem_dirv 
%             f01TD = f01TD + diag(f01(:,j))*Di((j-1)*num_quadr_pts+1:j*num_quadr_pts,:);
%          end
%          f0gu = B'*diag(W)*f01TD;

         
         f01TD =zeros(num_quadr_pts,num_quadr_pts);
         f0gu=zeros(neldof,neldof);
         s=1;
         for j=1:size(idx,1)
             for k = 1:num_elem_dirv
                 f01TD = f01TD + diag(f01(:,s))*Di((k-1)*num_quadr_pts+1:k*num_quadr_pts,:);
                 s=s+1;
             end
            f0gu(idx(j,1):idx(j,2),idx(j,3):idx(j,4)) = B'*diag(W)*f01TD;
         end
         
%          
%          f1u = zeros(num_quadr_pts,num_quadr_pts);
%          for j=1:num_elem_dirv 
%             f1u = f1u + Di((j-1)*num_quadr_pts+1:j*num_quadr_pts,:)'*diag(W.*f10(:,j))*B;
%          end
         
         f1u = zeros(neldof,neldof);
         f1uTD = zeros(neldof_per_ui,neldof_per_ui);
         s=1;
         for j=1:size(idx,1) 
             for k = 1:num_elem_dirv
                 f1uTD = f1uTD + Di((k-1)*num_quadr_pts+1:k*num_quadr_pts,:)'*diag(W.*f10(:,s));
                 s=s+1;
             end
             f1u(idx(j,1):idx(j,2),idx(j,3):idx(j,4)) = f1uTD*B;             
         end

%          f1gu_old = zeros(num_quadr_pts,num_quadr_pts);
%          for j=1:num_elem_dirv
%             f1gu_old = f1gu_old + Di((j-1)*num_quadr_pts+1:j*num_quadr_pts,:)'* diag(W.*f11(:,j))*Di((j-1)*num_quadr_pts+1:j*num_quadr_pts,:); 
%          end

         f1gu = zeros(neldof,neldof);
         s=1;
         for j=1:size(idx,1)
             for k = 1:num_elem_dirv
                 f1gu(idx(j,1):idx(j,2),idx(j,3):idx(j,4)) = f1gu(idx(j,1):idx(j,2),idx(j,3):idx(j,4)) + Di((k-1)*num_quadr_pts+1:k*num_quadr_pts,:)'* diag(W.*f11(:,s))*Di((k-1)*num_quadr_pts+1:k*num_quadr_pts,:); 
                 s=s+1;
             end
             %f1gu(idx(j,1):idx(j,2),idx(j,3):idx(j,4)) = f1guTD*Di((k-1)*num_quadr_pts+1:k*num_quadr_pts,:);
         end
         
         
                
         % element consistant tnagent evaulation (jac_e)
         jac_e = f0u + f0gu + f1u + f1gu;           
         
         % global residual and jacobian assembly 
         temp=conn(i,:)';
         k=1:neldof_per_ui;
         kk=temp(k);
         in_glb = global_idx_map(kk,:);
         kk =in_glb(in_glb>=0);
         %global residual
         global_res(kk) = global_res(kk)+ res_e(in_glb>=0); 
         %global jacobian 
         jac(kk,kk)=jac(kk,kk)+jac_e(in_glb>=0, in_glb>=0);              
    
     end
end