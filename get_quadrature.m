function [x,w] = get_quadrature(n)
%input: n: Number of quadrature points (Gauss)
%output:x: Gauss quadrature points
%       w: Gauss weights
% Golub-Welsch algorithm: (Brute force version by Trefethen-Bau)
% to calculate Gauss points and weights using Legendre weight function 
%
    beta = 0.5./sqrt(1-(2*(1:n-1)).^(-2));
    [Q,D]=eig(diag(beta,1)+diag(beta,-1));
    [x,i]=sort(diag(D)); 
    w=2*Q(1,i).^2';
end