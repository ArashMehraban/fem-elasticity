function [B, D_hat, W_hat] = get_shape(num_quadr_pts_in_1d, dim)
% GET_SHAPE returns shape/basis functions and their deivatievs evaluated at
% quadrature points in xi, eta and (zeta if 3D problem) directions:
% 
%        (+ shows the positive direction)
%                                     +
%          +                        zeta  +   
%         eta                         |  /eta
%          |                          | /
%   2D:    |______ xi + or       3D:  |/__ __ __ xi +
%                                 
%
% output: 
%       : W_hat: Weights for the elements
%       : B (N in most FEM contexts): 
%               basis/shape functions evaluted at quadrature points 
%       : D_hat (partial_N/partial_xi, partial_N/partial_eta, partial_N/partial_zeta) 
%               derivative of basis/shape functions evaluated at 
%               quadrature points in xi, eta and (zeta if 3D) directions.
%         D_hat consists of D0, D1 and (D2 if 3D problem) where
%           D0 (d_N/d_xi):   derivative of basis/shape functions wrt xi 
%                            evaluted at quadrature points 
%           D1 (d_N/d_eta):  derivative of basis/shape functions wrt eta 
%                            evaluted at quadrature points
%           D2 (d_N/d_zeta): derivative of basis/shape functions wrt zeta 
%                            evaluted at quadrature points
%
%           Structure of D_hat returned: 
%                 |D0|            |D0|
%           D_hat=|D1|   ,  D_hat=|D1|
%                                 |D2|      
%           To compute D0, D1 and D2:    
%
%               bHat: Shape/Basis functions in 1D
%               dHat: Derivative of shape/basis functions in 1D
%               Kronecker Product = (K)
%                  D0 = bHat (k) dHat     D0 = bHat (K) bHat (K) dHat   
%                  D1 = dHat (K) bHat     D1 = bHat (K) dHat (K) bHat
%                                         D2 = dHat (K) bHat (K) bHat
%
       % x: Guass points    w: Gauss weights
       [x, w] = get_quadrature(num_quadr_pts_in_1d);
       
       if(num_quadr_pts_in_1d == 2)
           bHat = [(1-x)/2, (1+x)/2];
           dHat = [-1/2+0*x, 1/2+0*x];
       end
       if(num_quadr_pts_in_1d == 3)
           bHat = [(x.^2 - x)/2, (1-x.^2), (x.^2+x)/2];
           dHat = [x-1/2, -2*x, x+1/2];
       end
       
       %2D
       if(dim == 2)
           % Basis/Shape functions (B)
           B = kron(bHat,bHat); 
           % Derivative of Basis/Shapefunctions (D)
           D_hat = [kron(bHat,dHat);kron(dHat,bHat)];           
           % weights (W_hat)
           W_hat = kron(w,w);
       end

       %3D
       if(dim == 3)
           % Basis/Shape functions (B)
           B = kron(kron(bHat,bHat),bHat); 
           % Derivative of Basis/Shapefunctions (D_hat)
           D_hat = [kron(kron(bHat,bHat),dHat); kron(kron(bHat,dHat),bHat);kron(kron(dHat,bHat),bHat)];
           % weights (W_hat)
           W_hat = kron(kron(w,w), w);
       end    
end