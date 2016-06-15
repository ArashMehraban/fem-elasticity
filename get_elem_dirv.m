function Di = get_elem_dirv_mat(invJe, Ds)  
%   input: invJe: Jacobian Inverse
%      2D:
%        : D0: Drivative of shape/basis functions evaluated at quadrature
%              points in direction 1 = xi
%        : D1: Drivative of shape/basis functions evaluated at quadrature
%              points in direction 2 = eta
%      3D:
%        : D3: Drivative of shape/basis functions evaluated at quadrature
%              points in direction 3 = zeta
%
%  output: Di: A cell-array of element derivative in directions i (1,2 or 3)
%              Di_1_mat = partial_N/partial_x 
%              Di_2_mat = partial_N/partial_y
%              Di_3_mat = partial_N/partial_z
    % For calculation:2D (same trend for 3D)
    % rearrange invJe to perform block matrix mutiplication
    % eg.                     
    % invJe = [j1_11  j1_12]
    %         [j1_21  j1_22] 
    %         [j2_11  j2_12]
    %         [j2_21  j2_22]
    %         [j3_11  j3_12]
    %         [j3_21  j3_22]
    %         [j4_11  j4_12]
    %         [j4_21  j4_22]
    % block_invJe = [j1_11  j1_12  j1_21  j1_22]
    %               [j2_11  j2_12  j2_21  j2_22]
    %               [j3_11  j3_12  j3_21  j3_22]
    %               [j4_11  j4_12  j4_21  j4_22]

    D_sz = size(fieldnames(Ds),1);
  
    D0 = Ds.D0;
    D1 = Ds.D1;
    if(D_sz == 3)
       D2 = Ds.D2;
    end
    
    blocksz = size(invJe,1)/size(D0,1); 
    num_gs_pts = size(D0,1);
    Di = zeros(blocksz*num_gs_pts,num_gs_pts);
    
    %2D  
    
    % structure of invJe per quadrature point
    % invJe = [partial_xi/partial_x , partial_eta/partial_x]
    %         [partial_xi/partial_y , partial_eta/partial_y]
    % entries with \partial_xi must be multiplied by D0 (=partial_N/partial_xi)
    % entries with \partial_eta must be multiplied by D1 (=partial_N/partial_eta)
    
    % structure of Di per quadrature point
    %Di{1} = partial_N/partial_x 
    %Di{2} = partial_N/partial_y
    if(D_sz == 2)
        block_invJe = [invJe(1:blocksz:end,:) ,  invJe(2:blocksz:end,:)];
        Di(1:num_gs_pts,:) = diag(block_invJe(:,1))*D0 + diag(block_invJe(:,2))*D1;
        Di(num_gs_pts+1:end,:) = diag(block_invJe(:,3))*D0 + diag(block_invJe(:,4))*D1; 
    end
    
    %3D
    
    % structure of invJe per quadrature point
    % invJe = [partial_xi/partial_x , partial_eta/partial_x, partial_zeta/partial_x]
    %         [partial_xi/partial_y , partial_eta/partial_y, partial_zeta/partial_y]
    %         [partial_xi/partial_z , partial_eta/partial_z, partial_zeta/partial_z]
    % entries with \partial_x must be multiplied by D0
    % entries with \partial_y must be multiplied by D1
    % entries with \partial_z must be multiplied by D2
    if(D_sz == 3)
        block_invJe = [invJe(1:blocksz:end,:) , invJe(2:blocksz:end,:), invJe(3:blocksz:end,:)];
        Di(1:num_gs_pts,:) = diag(block_invJe(:,1))* D0 + diag(block_invJe(:,2))*D1 + diag(block_invJe(:,3))*D2;
        Di(num_gs_pts+1:2*num_gs_pts,:) = diag(block_invJe(:,4))* D0 + diag(block_invJe(:,5))*D1 + diag(block_invJe(:,6))*D2;
        Di(2*num_gs_pts+1:3*num_gs_pts,:) = diag(block_invJe(:,7))* D0 + diag(block_invJe(:,8))*D1 + diag(block_invJe(:,9))*D2;
    end

end
