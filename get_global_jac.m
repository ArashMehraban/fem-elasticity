function globa_jac = get_global_jac(dlta_u, global_idx_map, msh,num_quadr_pts_in_1d,sz_u_field, userdf)
% GET_GLOBAL_RES evaluates the actual Jacobian (consistent tangent)
%  input:              dlta_u: vector of variation of unknowns 
%       :      global_idx_map: global map of local u's
%       :                 msh: mesh object 
%       : num_quadr_pts_in_1d: number of quadrature points in 1D
%       :              userdf: user supplied code for get_userdf function
%
% output: mfj: Matrix Free Jacobian(mfJ) the action of consistant tangent on dlta_u
%               

     %get the number of elements
     num_elem = msh.num_elem; 
     
     %get connectivity matrix
     conn = msh.conn;
     
     %get vertex coordinates
     vtx_coords = msh.vtx_coords;
     
     %Allocate space for global Matric Free Jacobian     
     globa_jac =zeros(size(dlta_u,1),1);
     
%      %get all dirichlet boundary node_sets
      dir_bndry_nodes = get_all_dir_ns(msh);
%       
%      delta_u_closure = get_closure_dlta_u(dlta_u,dir_bndry_nodes,global_idx_map); 
     
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
         
%          %get corresponding unknown/solution u and du for each element
%          elem_dlta_u = delta_u_closure(conn(i,:),:); 
         
         %get mapping constituents from element derivative 
         % 2D: [D0] du/dx
         %   D=[D1] du/dy
         %
         % 3D: [D1] du/dx
         %   D=[D2] du/dy
         %     [D3] du/dz
         [dets, D] = get_elem_dirv(element_vtx_coords, D_hat, dim);
         
         %get Gauss Weights for the current element
         W = W_hat.*dets;
%          
%          % grad_dlta_ue structure:
%          % 2D:             [D0*dlta_u1 | D0*dlta_u2] 
%          %    grad_dlta_ue=[D1*dlta_u1 | D1*dlta_u2]
%          %
%          % 3D:             [D0*dlta_u1 | D0*dlta_u2 | D0*dlta_u3]
%          %    grad_dlta_ue=[D1*dlta_u1 | D1*dlta_u2 | D1*dlta_u3]
%          %                 [D2*dlta_u1 | D2*dlta_u2 | D2*dlta_u3]
%          grad_dlta_ue = D*elem_dlta_u;
%          
%          % dlta_ue structure: 
%          %   2D: dlta_ue=[B*u1 | B*u2]   
%          %   3D: dlta_ue=[B*u1 | B*u2 | B*u3]
%          dlta_ue = B*elem_dlta_u;
                 
         % mapped quadrature points to reference coordinate system
         xe= B*element_vtx_coords;
        
         jac_e = zeros(size(B,2)*sz_u_field);
         
         for k=1:size(B,2)
             f = userdf(kron(B(:,k),eye(sz_u_field)), kron(D(:,k),eye(sz_u_field)),xe);

             tmpf = zeros(size(f,1)/sz_u_field, sz_u_field^2);
             [r,c]=size(tmpf);
             
             s=1;
             for m=1:sz_u_field
                tmpf(:,1+(m-1)*sz_u_field:m*sz_u_field) = f(s:sz_u_field:end,:);                
                s=s+1;
             end
             
             tmpf = reshape(tmpf, size(W,1), []);

             for j=1:size(tmpf,2)
                 tmpf(:,j) = W.*tmpf(:,j);
             end

             tmpf = reshape(tmpf,r,c);
             
             %Jacobian Matrix per element 
             jac_e(:,1+(k-1)*sz_u_field:k*sz_u_field) = reshape(([B' D']*tmpf)',sz_u_field,[])';         
         end

         
         % global (the action of) jacobian assembly 
         temp=conn(i,:)';
         k=1:neldof;
         kk=temp(k);
         in_glb = global_idx_map(kk,:);
         kk =in_glb(in_glb>=0);
         globa_jac(kk) = globa_jac(kk)+ jac_e(in_glb>=0);    
     end
end