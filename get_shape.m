function [B, Ds, W_hat] = get_shape(elm_type)
% GET_SHAPE returns shape/basis functions for following Element Types:
%     2D elements: 4: QUAD4  or  9: QUAD9
%     3D elements: 8: HEX8   or 27: HEX27 
%
% input: elm_type: Element type 
%        (+ shows the positive direction)
%                                  +
%          +                     zeta  +   
%         eta                      |  /eta
%          |                       | /
%   2D:    |______ xi + or    3D:  |/__ __ __ xi +
%                                 
% 2D:
% output: shape object that contains:
%       : B: Basis/Shape functions evaluted at quadrature points.
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
%         2D                       3D
%   ------------------     --------------------------       
%   D0 = bHat (k) dHat     D0 = bHat (K) bHat (K) dHat   
%   D1 = dHat (K) bHat     D1 = bHat (K) dHat (K) bHat
%                          D2 = dHat (K) bHat (K) bHat
%
%
       if(elm_type == 4 || elm_type == 8)
            n_gs_pts = 2;
            % x: Guass points    w: Gauss weights
            [x, w] = get_quadrature(n_gs_pts);  
            bHat = [(1-x)/2, (1+x)/2];
            dHat = [-1/2+0*x, 1/2+0*x];
       end
       if(elm_type == 9 || elm_type == 27)
            n_gs_pts = 3;
            % x: Guass points    w: Gauss weights
            [x, w] = get_quadrature(n_gs_pts); 
            bHat = [(x.^2 - x)/2, (1-x.^2), (x.^2+x)/2];
            dHat = [x-1/2, -2*x, x+1/2];
       end
      

%   2D 
    if(elm_type == 4 || elm_type == 9)         
        
        field_names = {'D0', 'D1'};
        
        % Basis/Shape functions (B)
        B = kron(bHat,bHat);  
        % xi (D0)
        field_vals{1} = kron(bHat,dHat);
        % eta (D1)
        field_vals{2} = kron(dHat,bHat);
        % weights (W_hat)
        W_hat = kron(w,w);
    end
    
%   3D
    if(elm_type == 8 || elm_type == 27) 
        
        field_names = {'D0', 'D1' ,'D2'};
        
        % Basis/Shape functions (B)
        B = kron(kron(bHat,bHat),bHat);  
        % xi (D0)
        field_vals{1} = kron(kron(bHat,bHat),dHat);
        % eta (D1)
        field_vals{2} = kron(kron(bHat,dHat),bHat);
        % zeta (D2)
        field_vals{3} = kron(kron(dHat,bHat),bHat);
        % weights (W_hat)
        W_hat = kron(kron(w,w), w);
    end
    
    Ds=struct();
    for i=1:size(field_names,2)
        Ds.(field_names{i}) = field_vals{i};
    end
    
    
end