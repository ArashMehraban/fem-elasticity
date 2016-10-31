function global_res = get_global_res(u, global_idx_map, msh, dir_bndry_val,num_quadr_pts_in_1d,userf)
% GET_GLOBAL_RES evaluates the global residual and the consistent tangent
%  input:              u: vector of unknowns 
%       : global_idx_map: global map of local u's
%       :            msh: mesh object (see get_mesh function)
%       :  dir_bndry_val: Dirchlet boundary values if any
%       :          userf: user supplied code for get_userf function
%
% output: global_res: global residual
%       : 

     %get the number of elements
     num_elem = msh.num_elem; 
     
     %get connectivity matrix
     conn = msh.conn;
     
     %get vertex coordinates
     vtx_coords = msh.vtx_coords;
     
     %Allocate space for global residual for unkowns     
     global_res =zeros(size(u,1),1);
     
     %get all dirichlet boundary node_sets
     dir_bndry_nodes = get_all_dir_ns(msh);
     
     %get closure of u : a vector consisting the unknown and dirchlet boundary values
     u_closure =  get_closure_u(u,dir_bndry_nodes,dir_bndry_val,global_idx_map);     
     
     %get the dimension of the problem
     dim = msh.num_dims;
     
     %get Weights, Basis (B) functions and their Derivatives (D0, D1 and D2)
     %D_hat = partial_B/partial_x_i
     [B, D_hat, W_hat] = get_shape(num_quadr_pts_in_1d, dim);

     for i=1:num_elem
         
         %get number of dof per unknown u_i per element
         neldof = size(conn(i,:),2);
           
         %get corresponding vertex coordinates for each element 
         element_vtx_coords = vtx_coords(conn(i,:),:);
         
         %get corresponding unknown/solution u for each element
         elem_u = u_closure(conn(i,:),:);   
         
         %get mapping constituents from element derivative 
         % 2D: [D1] dN/dx
         %   D=[D2] dN/dy
         %
         % 3D: [D1] dN/dx
         %   D=[D2] dN/dy
         %     [D3] dN/dz
         [dets, D] = get_elem_dirv(element_vtx_coords, D_hat, dim);
         
         %get Gauss Weights for the current element
         W = W_hat.*dets;
          
         % grad_ue structure:
         % 2D:       [D1*u1 | D1*u_2]     [du1/dx | du2/dx] 
         %   grad_ue=[D2*u1 | D2*u_2]  =  [du1/dy | du2/dy]
         %
         % 3D:        [D1*u1 | D1*u2 | D1*u3]   [du1/dx | du2/dx | du3/dx] 
         %    grad_ue=[D2*u1 | D2*u2 | D2*u3] = [du1/dx | du2/dy | du3/dy]
         %            [D3*u1 | D3*u2 | D3*u3]   [du1/dz | du2/dz | du3/dz]      
         grad_ue=D*elem_u;
         
         % ue structure: 
         %   2D: ue=[B*u1 | B*u2]   
         %   3D: ue=[B*u1 | B*u2 | B*u3]
         ue = B*elem_u;  
         
         % vertex coordinates mapped to reference coordinate system
         xe= B*element_vtx_coords;
         
        
         [f0,f1] = userf(ue, grad_ue, xe); 

         % overwirte f0 by element-wise multiplication of its values by W
         for j=1:size(f0,2)
             f0(:,j)=W.*f0(:,j);
         end
         
         % reshape f1 for element-wise multiplication of its values by W
         tmp=reshape(f1,size(W,1),[]);         
         for j=1:size(tmp,2)
             tmp(:,j)=W.*tmp(:,j);
         end
                  
         % overwirte f1 by reshaped tmp
         f1=reshape(tmp,size(f1,1),[]);
         
         % element residual evaluation
         res_e = B'*f0 + D'*f1;
 
         
         % global residual and consistent tangent assembly 
         temp=conn(i,:)';
         k=1:neldof;
         kk=temp(k);
         in_glb = global_idx_map(kk,:);
         kk =in_glb(in_glb>=0);
         %global residual
         global_res(kk) = global_res(kk)+ res_e(in_glb>=0);    
     end
end