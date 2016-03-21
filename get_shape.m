function [W_hat, B, D0, D1, D2] = get_shape(elm_type)
% GET_SHAPE returns shape/basis functions for following Element Types:
%     2D elements: Q1 , Q2,
%     3D elements: Q1H , Q2H 
%
% input: elm_type: Element type
%
%         eta                     eta 
%          |                       |
%   2D:    |______ xi   or    3D:  |_____ xi
%                                 /
%                                /  
%                             zeta
% 2D:
% output: B: Basis/Shape functions evaluted at quadrature points.
%       : D0: derivative of basis/shape functions wrt xi evaluted at quadrature points
%       : D1: derivative of basis/shape functions wrt eta evaluted at quadrature points
%       : W_hat: Weights for the elements
% 3D:
% output: B: Basis/Shape functions evaluted at quadrature points.
%       : D0 = Derivative of basis/shape functions wrt xi evaluted at quadrature points
%       : D1 = Derivative of basis/shape functions wrt eta evaluted at quadrature points
%       : D2 = Derivative of basis/shape functions wrt zeta evaluted at quadrature points
%       : W_hat: Weights for the elements
%
%    bHat: Shape/Basis functions in 1D
%    dHat: Derivative of shape/basis functions in 1D
%        Kronecker Product = (K)
% Note: Tensor product order changes form 2D to 3D cases.
%         2D                       3D
%   ------------------     --------------------------       
%   D0 = bHat (k) dHat     D0 = bHat (K) bHat (K) dHat   
%   D1 = dHat (K) bHat     D1 = bHat (K) dHat (K) bHat
%                          D2 = dHat (K) bHat (K) bHat
%
%
       if (strcmp(elm_type,'Q1') || strcmp(elm_type,'Q1H'))
            n_gs_pts = 2;
            % x: Guass points    w: Gauss weights
            [x, w] = get_quadrature(n_gs_pts);  
            bHat = [(1-x)/2, (1+x)/2];
            dHat = [-1/2+0*x, 1/2+0*x];
       end
       if(strcmp(elm_type, 'Q2') || strcmp(elm_type,'Q2H'))
            n_gs_pts = 3;
            % x: Guass points    w: Gauss weights
            [x, w] = get_quadrature(n_gs_pts); 
            bHat = [(x.^2 - x)/2, (1-x.^2), (x.^2+x)/2];
            dHat = [x-1/2, -2*x, x+1/2];
       end
      

%   2D 
    if (nargout == 4) 
        % Basis/Shape functions
        B = kron(bHat,bHat);  
        % xi
        D0 = kron(bHat,dHat);
        % eta
        D1 = kron(dHat,bHat);
        % weights
        W_hat = kron(w,w);

    end
    
%   3D
    if (nargout == 5) 
        % Basis/Shape functions 
        B = kron(kron(bHat,bHat),bHat);  
        % xi 
        D0 = kron(kron(bHat,bHat),dHat);
        % eta 
        D1 = kron(kron(bHat,dHat),bHat);
        % zeta
        D2 = kron(kron(dHat,bHat),bHat);
        % weights
        W_hat = kron(kron(w,w), w);
    end
end