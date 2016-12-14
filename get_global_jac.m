function sparse_global_jac = get_global_jac(sz_u, global_idx_map, msh,num_quadr_pts_in_1d,sz_u_field, userdf)
% GET_GLOBAL_JAC evaluates the actual Jacobian (consistent tangent)
%  input:                sz_u: size of unknown vector 
%       :      global_idx_map: global map of local u's
%       :                 msh: mesh object 
%       : num_quadr_pts_in_1d: number of quadrature points in 1D
%       :          sz_u_field: dof per node
%       :              userdf: user supplied code for get_userdf function
%
% output:   sparse_global_jac: Jacobian Matrix (consistant tangent)
%               

     %get the number of elements
     num_elem = msh.num_elem; 
     
     %get connectivity matrix
     conn = msh.conn;
     
     %get vertex coordinates
     vtx_coords = msh.vtx_coords;
     
     sparse_global_jac= spalloc(sz_u,sz_u,sz_u*sz_u_field);
     
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
                 
         % mapped quadrature points to reference coordinate system
         xe= B*element_vtx_coords;
        
         jac_e = zeros(size(B,2)*sz_u_field);
         
         BI = kron(B,eye(sz_u_field));
         DI = kron(D,eye(sz_u_field));
         
         seq = 1:size(B,1)*sz_u_field;
         idx = reshape(reshape(seq',[],sz_u_field)',[],1);
         
         for k=1:size(BI,2)
             
             f = userdf(reshape(BI(:,k),sz_u_field,[])', reshape(DI(:,k),sz_u_field,[])',xe);
             [r,c]=size(f);

             tmpf = reshape(f, size(W,1), []);

             for j=1:size(tmpf,2)
                 tmpf(:,j) = W.*tmpf(:,j);
             end

             tmpf = reshape(tmpf,r,c);
             
             %Jacobian Matrix per element 
             jac_e(:,idx(k)) = reshape(([B' D']*tmpf),[],1);

         end

         
         % global (the action of) jacobian assembly 
         temp=conn(i,:)';
         k=1:neldof;
         kk=temp(k);
         in_glb = global_idx_map(kk,:);
         kk =in_glb(in_glb>=0);   
         %global jacobian  
         sparse_global_jac(kk,kk)=sparse_global_jac(kk,kk)+jac_e(in_glb>=0, in_glb>=0); 
     end
end