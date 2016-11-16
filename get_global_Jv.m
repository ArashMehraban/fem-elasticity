function Jv = get_global_Jv(dlta_u, global_idx_map, msh,num_quadr_pts_in_1d, userdf)
% GET_GLOBAL_MAT)FREEJ evaluates the action of Jaccobian (consistent
% tangent) on dlta_u. It produces a vector mfj=J*dlta_u
%  input:              dlta_u: vector of variation of unknowns 
%       :      global_idx_map: global map of local u's
%       :                 msh: mesh object 
%       : num_quadr_pts_in_1d: number of quadrature points in 1D
%       :              userdf: user supplied code for get_userdf function
%
% output: Jv: Jacobian-vector Product: the action of consistant tangent on dlta_u
%               

     %get the number of elements
     num_elem = msh.num_elem; 
     
     %get connectivity matrix
     conn = msh.conn;
     
     %get vertex coordinates
     vtx_coords = msh.vtx_coords;
     
     %Allocate space for global Matric Free Jacobian     
     Jv =zeros(size(dlta_u,1),1);
     
     %get all dirichlet boundary node_sets
     dir_bndry_nodes = get_all_dir_ns(msh);
      
     delta_u_closure = get_closure_dlta_u(dlta_u,dir_bndry_nodes,global_idx_map); 
     
     dim = msh.num_dims;                                            
     %get Weights, Basis (B) functions and   [D0]
     %  their Derivatives:           D_hat = [D1]
     %                                       [D2] <-- if 3D
     [B, D_hat, W_hat] = get_shape(num_quadr_pts_in_1d, dim);
  
               
     for i=1:num_elem
         
         %get number of dof per unknown u_i per element
         neldof = size(conn(i,:),2);
          
         %get corresponding vertex coordinates for each element 
         element_vtx_coords = vtx_coords(conn(i,:),:);
         
         %get corresponding unknown/solution u and du for each element
         elem_dlta_u = delta_u_closure(conn(i,:),:); 
         
         %get mapping constituents from element derivative 
         % 2D: [D1] du/dx
         %   D=[D2] du/dy
         %
         % 3D: [D1] du/dx
         %   D=[D2] du/dy
         %     [D3] du/dz
         [dets, D] = get_elem_dirv(element_vtx_coords, D_hat, dim);
         
         %get Gauss Weights for the current element
         W = W_hat.*dets;
         
         % grad_dlta_ue structure:
         % 2D:             [D1*dlta_u1 | D1*dlta_u2] 
         %    grad_dlta_ue=[D2*dlta_u1 | D2*dlta_u2]
         %
         % 3D:             [D1*dlta_u1 | D1*dlta_u2 | D1*dlta_u3]
         %    grad_dlta_ue=[D2*dlta_u1 | D2*dlta_u2 | D2*dlta_u3]
         %                 [D3*dlta_u1 | D3*dlta_u2 | D3*dlta_u3]
         grad_dlta_ue = D*elem_dlta_u;
         
         % dlta_ue structure: 
         %   2D: dlta_ue=[B*u1 | B*u2]   
         %   3D: dlta_ue=[B*u1 | B*u2 | B*u3]
         dlta_ue = B*elem_dlta_u;
                 
         % mapped quadrature points to reference coordinate system
         xe= B*element_vtx_coords;
         
         f = userdf(dlta_ue, grad_dlta_ue,xe);
         [r,c]=size(f);
         
         tmpf = reshape(f, size(W,1), []);
  
         for j=1:size(tmpf,2)
             tmpf(:,j) = W.*tmpf(:,j);
         end

         tmpf = reshape(tmpf,r,c);
         
         %Matrix Free Jacobian per element (This is a Vector)
         Jv_e = [B' D']*tmpf;
                 
         
         % global (the action of) jacobian assembly 
         temp=conn(i,:)';
         k=1:neldof;
         kk=temp(k);
         in_glb = global_idx_map(kk,:);
         kk =in_glb(in_glb>=0);
         Jv(kk) = Jv(kk)+ Jv_e(in_glb>=0);    
     end
end