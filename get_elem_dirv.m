function Di = get_elem_dirv(invJe, Ds)  
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


    % rearrange invJe to perform block matrix mutiplication
    % eg.                     
    % invJe = [j1_11  j1_12 | j2_11  j2_12 | j3_11  j3_12 | j4_11 j4_12]
    %         [j1_21  j1_22 | j2_21  j2_22 | j3_21  j3_22 | j4_21 j4_22] 
    %
    % block_invJe = [j1_11  j1_12  j1_21  j1_22]
    %               [j2_11  j2_12  j2_21  j2_22]
    %               [j3_11  j3_12  j3_21  j3_22]
    %               [j4_11  j4_12  j4_21  j4_22]
    
    % rearrange invJe to perform block matrix mutiplication
    % eg.
    % invJe = [j1_11  j1_12  j1_13 | j2_11  j2_12  j2_13 | j3_11  j3_12  j3_13 | j4_11 j4_12 j4_13]
    %         [j1_21  j1_22  j1_23 | j2_21  j2_22  j2_23 | j3_21  j3_22  j3_21 | j4_21 j4_22 j4_23] 
    %         [j1_31  j1_32  j1_33 | j2_31  j2_32  j2_33 | j3_31  j3_32  j3_31 | j4_31 j4_32 j4_33]
    %
    % block_invJe = [j1_11  j1_12  j1_13  j1_21  j1_22  j1_23  j1_31  j1_32  j1_33]
    %               [j2_11  j2_12  j2_13  j2_21  j2_22  j2_23  j2_31  j2_32  j2_33]
    %               [j3_11  j3_12  j3_13  j3_21  j3_22  j3_23  j3_31  j3_32  j3_33]
    %               [j4_11  j4_12  j4_13  j4_21  j4_22  j4_23  j4_31  j4_32  j4_33]
    
    D_sz = size(fieldnames(Ds),1);
  
    D0 = Ds.D0;
    D1 = Ds.D1;
    if(D_sz == 3)
       D2 = Ds.D2;
    end
    
    invJe = cell2mat(invJe);
    blocksz = size(invJe,2)/size(D0,1); 
    
    %2D  
    if(D_sz == 2)
        block_invJe = [invJe(:, 1:blocksz:end) ; invJe(:, 2:blocksz:end)]';
        Di{1} = diag(block_invJe(:,1))*D0 + diag(block_invJe(:,2))*D1;
        Di{2} = diag(block_invJe(:,3))*D0 + diag(block_invJe(:,4))*D1; 
    end
    
    
    
    %3D
    if(D_sz == 3)
        block_invJe = [invJe(:, 1:blocksz:end) ; invJe(:, 2:blocksz:end); invJe(:, 3:blocksz:end)]';
        Di{1} = diag(block_invJe(:,1))* D0 + diag(block_invJe(:,2))*D1 + diag(block_invJe(:,3))*D2;
        Di{2} = diag(block_invJe(:,4))* D0 + diag(block_invJe(:,5))*D1 + diag(block_invJe(:,6))*D2;
        Di{3} = diag(block_invJe(:,7))* D0 + diag(block_invJe(:,8))*D1 + diag(block_invJe(:,9))*D2;
    end
    

end
