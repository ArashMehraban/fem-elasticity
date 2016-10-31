function [f0,f1] = userf_L2_2d(ue, grad_ue, xe) 
% userf_L2_2d provides weak form of the L2-Projection problem to solve 
%  integral[v(u-g)] = 0 (v: test/weight function) 
%  
%    Weak Form: 
%      integral(v*f0(ue,grad_ue) + grad_v : f1(ue,grad_ue)) = 0 
% compare to:
%      integral[v(u-g)] = 0 
% Therfore,
%      f0 = u-g
%      f1 = 0   (as there is no grad_u in L2 projection problem)
%
%
%  input:      ue: corresponding u for each element evalauted at quadrature points
%       : grad_ue: corresponding grad_u for each element evalauted at quadrature points
%       :      xe: quadrature points mapped to the reference elem
%
% output:  f0: any possible source from given problem : u-g
%       :  f1: any possible source from given problem : 0



    x=xe(:,1);
    y=xe(:,2);

   % RHS 
   g = tanh(x).*exp(y)+sin(y);

   f0 = ue - g; 
   f1 = 0*grad_ue;
end