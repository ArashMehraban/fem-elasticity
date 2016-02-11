function [K, F] = assembly(conn,vtx_coords,givenF)
%
%  input: conn: connectivity matrix from mesh
%       : vtx_coords: veterx-coordinates for each node from the mesh
%       : givenF: user-defined F to be evalauted
% output: K: Global stiffness matrix
%       : F: rhs

    %get the number of elements (=row size of conn)
    nel=size(conn,1);
    
    %get number of dof for each element (=column size of conn)
    neldof = numel(conn(1,:));
    
    %get the number of nodes in mesh from vtx_coords
    num_nodes = size(vtx_coords,1);
    
    %Allocate space for stiffnes Matrix, K, and rhs vector, F
    K = zeros(num_nodes,num_nodes);
    F = zeros(num_nodes,1);

    for n=1:nel    
         %get 2-noded Gauss quadrature 
         [gx_pts, gs_w] = get_quadrature(2);
         
         %Evaluating shape functions and its derivatives wrt to the reference
         %coordinate system at quadrature points
         [B, D0, D1] = get_shapeF_dF_at_quadr_pts(gx_pts);
         
         %get Gauss weights for the current element
         W = kron(gs_w,gs_w);

         %get the corresponding vertex coordinates for each element from 
         %connectivity matrix (ScriptE matrix in Jed's paper)
         element_vtx_coords = vtx_coords(conn(n,:),:);
         
         %mapping using jacobian
         [dets, ~] = jacobian(element_vtx_coords, D1, D0);
         
         %element stiffness matrix
         k_e=B'*(B*diag(dets'));
         
         
         %rhs
         %f_e = sum_of_(Transpose_of_(B)*f*|J|)
         %
         % x's and y's must be mapped from global to local coordiante
         % system before evaluation by f (givenF)
         % mapped_x: x_i_(in_reference_coords)= sum_of_(B_i * x_i_(in_global_coords) )
         % mapped_y: y_i_(in_reference_coords)= sum_of_(B_i * y_i_(in_global_coords) )

         %mapping of x and y's:
         mapped_x = B*element_vtx_coords(:,1);
         mapped_y = B*element_vtx_coords(:,2);
         
         %calculating f_e
         f_e =  B'*givenF(mapped_x, mapped_y).*(dets'.*W);  
        
         % Global stiffness matrix, K, Assembley
         temp=conn(n,:)';
         for i=1:neldof,
             I=temp(i);
             if I>0
                 F(I)=F(I)+f_e(i);
                 for j=1:neldof,
                     J=temp(j);
                     if J>0
                         K(I,J)=K(I,J)+k_e(i,j);
                     end
                 end
             end
         end    
    end

end