function Di = get_elem_dirv(invJe, D0, D1, varargin)  
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
    % invJe = [j1_11  j1_12  j2_11  j2_12  j3_11  j3_12  j4_11 j4_12]
    %         [j1_21  j1_22  j2_21  j2_22  j3_21  j3_22  j4_21 j4_22] 
    %
    % block_invJe = [j1_11  j1_12  j1_21  j1_22]
    %               [j2_11  j2_12  j2_21  j2_22]
    %               [j3_11  j3_12  j3_21  j3_22]
    %               [j4_11  j4_12  j4_21  j4_22]
    
    invJe = cell2mat(invJe);

        blocksz = size(invJe,2)/size(D0,1);
        block_invJe = [invJe(:, 1:blocksz:end) ; invJe(:, 2:blocksz:end)]';
   
    
    %2D     
         Di{1} = (diag(block_invJe(:,1)) + diag(block_invJe(:,2)))*D0;
         Di{2} = (diag(block_invJe(:,3)) + diag(block_invJe(:,4)))*D1;    
    
    %3D
     if(nargin > 3)
        D2 = varargin{1};
        Di{3} = diag(block_invJe(:,5))* D2 + diag(block_invJe(:,6))*D2;                
    end


end
