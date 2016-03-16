function [global_res , jac] = eval_res(u ,conn, vtx_coords, elm_type) 
%
%  input: conn: mesh connectivity matrix 
%       : vtx_coords: mesh nodes veterx-coordinates 
%       : current solution (guess)
%       : elem_type = elment type
%
% output: global_res

    %get the number of elements (=row size of conn)
    nel=size(conn,1);
    
    %get number of dof for each element (=column size of conn)
    neldof = numel(conn(1,:));
    
    %get the number of nodes in mesh from vtx_coords
    num_nodes = size(vtx_coords,1);
    
    %Allocate space for global residual
    global_res =zeros(1,num_nodes)';
    jac = zeros(num_nodes,num_nodes);

    %for i=1:nel    
                 
         % Be: shape/basis functions at quadrature points for each element         
         % D0e: Derivative of shape functions with respect to xi for each element
         % D1e: Derivative of shape functions with respect to eta for each element
         % D2e: Derivative of shape functions with respect to zeta for each element
         % W_hat_e: Weights for each element
         
          if(strcmp(elm_type,'Q4') || strcmp(elm_type,'Q9'))
             [W_hat_e, Be, D0e, D1e] = get_shape(elm_type);  
          end
          if(strcmp(elm_type,'Hex8') || strcmp(elm_type,'Hex27'))
              [W_hat_e, Be, D0e, D1e, D2e] = get_shape(elm_type);
          end
         
     for i=1:nel  
          
         %get corresponding vertex coordinates for each element 
         element_vtx_coords = vtx_coords(conn(i,:),:);
         %get corresponding unknown/solution u for each element
         elem_u = u(conn(i,:));
                  
         %get mapping constituents from jacobian
         if (strcmp(elm_type,'Q4') || strcmp(elm_type,'Q9') )
             [dets, invJe] = jacobian(element_vtx_coords, D0e, D1e);             
             Die = get_elem_dirv(invJe, D0e, D1e);
         end
         
         if (strcmp(elm_type,'Hex8') || strcmp(elm_type,'Hex27') )
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
             De_res = De_res +(Die{j}'* We.*f1{j});
         end
         
               
         % element residual evaluation
         res_e = Be'*We.*f0 + De_res;
        
         
         % global residual assembly
         temp=conn(i,:)';
         for k=1:neldof
             kk=temp(k);
             if kk>0
                 global_res(kk)=global_res(kk)+res_e(k);
             end
         end
         
         %element jac_e constituents
         f0u = Be'*diag(We.*f00)*Be;
         
         f10TD =0;
         for j=size(Die,2) 
            f10TD = f10TD + diag(f10{j})*Die{j};
         end
         f0gu = Be'*diag(We)*f10TD;
         
         f1u = 0;
         for j=size(Die,2) 
            f1u = f1u + Die{j}'*diag(We.*f01{j})*Be;
         end
         
         f1gu = 0;
         for j=1:size(Die,2)
            f1gu = f1gu + Die{j}'* diag(We.*f11{j})*Die{j}; 
         end
         
         % element consistant tnagent evaulation (jac_e)
         jac_e = f0u + f0gu + f1u + f0gu;         
                         
         % global consistent tangent/jacobian assembly
         tempj=conn(i,:)';
         for k=1:neldof
             kk=tempj(k);
             if kk>0
                 for j=1:neldof
                     J=tempj(j);
                     if J>0
                         jac(kk,J)=jac(kk,J)+jac_e(k,j);
                     end
                 end
             end
         end          
     end
end