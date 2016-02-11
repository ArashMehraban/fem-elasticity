function [B, D0, D1] = get_shapeF_dF_at_quadr_pts(quadrature_pts)
%input: quadrature points
%output: B: shape functions evaluted at quadrature points.
%      : D0 = derivative of shape functions wrt xi evaluted at quadrature points
%      : D1 = derivative of shape functions wrt eta evaluted at quadrature points
    x = quadrature_pts;  
    bHat = [(1-x)/2, (1+x)/2];
    %Shape functions evaluated at quadrature points
    B = kron(bHat,bHat);
    
    dHat = [-1/2+0*x, 1/2+0*x];
    %Derivative of shape functions with respect to xi
    D0 = kron(dHat,bHat);
    %Derivative of shape functions with respect to eta
    D1 = kron(bHat,dHat);

end