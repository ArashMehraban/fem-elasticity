function [global_res , jac] = eval_res(u, global_idx_map, shape_bdw, msh, dir_bndry_val)
% EVAL_RES evaluates the global residual and the Jacobian
%  input: conn: mesh connectivity matrix 
%       : vtx_coords: mesh nodes veterx-coordinates 
%       : current solution (guess)
%       : elem_type = elment type
%
% output: global_res, jac

     %get the number of elements
     num_elem = msh.num_elem; 
     
     %get connectivity matrix
     conn = msh.conn;
     
     %get vertex coordinates
     vtx_coords = msh.vtx_coords;
     
     %unknown size (size of dofs per unknown u)
     %unknown_sz = sum(u_section((u_section(:,1)== -1),2));
     unknown_sz = size(u,1);
     
     %Allocate space for global residual for unkowns     
     global_res =zeros(unknown_sz,1);
     
     %get all dirichlet boundary node_sets
      dir_bndry_nodes = get_all_dir_ns(msh);
     
      num_nodes = msh.num_nodes;
      num_dims = msh.num_dims;
     % construct the global_u that contains the boundary values
      global_u = setbdry(u,num_nodes, num_dims, dir_bndry_nodes,dir_bndry_val);
     
     %Allocate space for globall Jacobian
     jac = zeros(unknown_sz,unknown_sz);
     
     sz_global_idx_map = size(global_idx_map,2);
           
     for i=1:num_elem
         
         %get number of dof for each element
         neldof = size(conn(i,:),2)*sz_global_idx_map;
          
         %get corresponding vertex coordinates for each element 
         element_vtx_coords = vtx_coords(conn(i,:),:);
         
         %get corresponding unknown/solution u for each element
         elem_u = global_u(conn(i,:));         
                 
         %get mapping constituents from jacobian
         if(strcmp(elm_type,'Q1') || strcmp(elm_type,'Q2') )
             [dets, invJe] = jacobian(element_vtx_coords, D0e, D1e);             
             Die = get_elem_dirv(invJe, D0e, D1e);
         end
         
         if(strcmp(elm_type,'Q1H') || strcmp(elm_type,'Q2H') )
            [dets, invJe] = jacobian(element_vtx_coords, D0e, D1e, D2e);
            Die = get_elem_dirv(invJe, D0e, D1e, D2e);
         end
         
         ue = Be*elem_u;
         
         grad_ue=cell(1,size(Die,2));
         for  j=1:size(Die,2)
             grad_ue{j}=Die{j}*elem_u;
         end
                  
         mp_qd_pts= Be*element_vtx_coords;
         
        [f0,f1,f00, f01, f10, f11] = userf(ue, grad_ue,mp_qd_pts); 
                 
         %get Gauss weights for the current element
         We = W_hat_e.*dets';
                  
         De_res = 0;
         for j=1:size(Die,2)
             De_res = De_res +(Die{j}'* (We.*f1{j}));
         end
         
         % element residual evaluation
         res_e = Be'*(We.*f0) + De_res;
                      
         %element jac_e constituents
         f0u = Be'*diag(We.*f00)*Be;
         
         f01TD =0;
         for j=1:size(Die,2) 
            f01TD = f01TD + diag(f01{j})*Die{j};
         end
         f0gu = Be'*diag(We)*f01TD;
         
         f1u = 0;
         for j=1:size(Die,2) 
            f1u = f1u + Die{j}'*diag(We.*f10{j})*Be;
         end
         
         f1gu = 0;
         for j=1:size(Die,2)
            f1gu = f1gu + Die{j}'* diag(We.*f11{j})*Die{j}; 
         end
         
         % element consistant tnagent evaulation (jac_e)
         jac_e = f0u + f0gu + f1u + f1gu;           
         
         % global residual and jacobian assembly 
         temp=conn(i,:)';
         k=1:neldof;
         kk=temp(k);
         in_glb = global_idx_map(kk);
         kk =in_glb(in_glb~=0);
         %global residual
         global_res(kk) = global_res(kk)+ res_e(in_glb~=0); 
         %global jacobian 
         jac(kk,kk)=jac(kk,kk)+jac_e(in_glb~=0, in_glb~=0);              
    
     end
end