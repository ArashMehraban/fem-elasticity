function f = userdf_L2_3d(dlta_ue, grad_dlta_ue, xe) 
%USERF_L2_3d provides weak form of the L2-Projection problem to solve 
%
%  input:      ue: corresponding u for each element evalauted at quadrature points
%       : grad_ue: corresponding grad_u for each element evalauted at quadrature points
%       :      xe: quadrature points mapped to the reference elem
%
% output:  f0: any possible source from given problem
%       :  f1: any possible source from given problem
%       : f00: partial of f0 wrt u (algebriac operations only)
%       : f01: partial of f0 wrt grad_u (algebriac operations only)
%       : f10: partial of f1 wrt u (algebriac operations only)
%       : f11: partial of f1 wrt grad_u (algebriac operations only)

   % Weak Form for L2-Projection problem (v: test/weight function):  
   %   integral(v*f0(ue,grad_ue) + grad_v : f1(ue,grad_ue)) = 0   
 
   f = [dlta_ue ; 0*grad_dlta_ue];

end