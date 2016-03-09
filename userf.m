function [f0, f1] = userf(ue, grad_ue, mp_qd_pts)
%USERF provides weak form of the problem to solve 

%  input: ue: corresponding u for each element evalauted at quadrature points
%       : grad_ue: corresponding grad_u for each element evalauted at quadrature points (cell not matrix)
%       : mp_qd_pts: quadrature points mapped in the reference elem
%
% output: f0: any possible source from given problem
%       : f1: any possible source from given problem
     x = mp_qd_pts(:,1);
     y = mp_qd_pts(:,2);
     
     %user defined rhs
     g = tanh(x).*exp(y)+sin(y);   
     
     % Weak Form: (v: test/weight function)
     % integral(v*f0(ue,grad_ue) + grad_v : f1(ue,grad_ue)) = 0 
     % for L2-Projection problem:
     f0 = ue - g; 
     
     f1=cell(1,size(grad_ue,2));
     for  j=1:size(grad_ue,2)
          f1{j} = 0*grad_ue{j};
     end
end


