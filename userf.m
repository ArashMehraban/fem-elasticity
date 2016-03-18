function [f0, f1,f00, f01, f10, f11] = userf(ue, grad_ue, xe)
%USERF provides weak form of the problem to solve 

%  input: ue: corresponding u for each element evalauted at quadrature points
%       : grad_ue: corresponding grad_u for each element evalauted at quadrature points (cell not matrix)
%       : xe: quadrature points mapped in the reference elem
%
% output: f0: any possible source from given problem
%       : f1: any possible source from given problem
     x=xe(:,1);
     y=xe(:,2);
     %user defined rhs
     g = tanh(x).*exp(y)+sin(y); 
     
     % Weak Form: (v: test/weight function)
     % integral(v*f0(ue,grad_ue) + grad_v : f1(ue,grad_ue)) = 0 
     % for L2-Projection problem:
       
       %==== L2 ======%
       f0 = ue - g; 
      
       f1=cell(1,size(grad_ue,2));
       for  j=1:size(grad_ue,2)
           f1{j} = 0*grad_ue{j};
       end
       %==============%
       
       %f_0,i = partial(f_i)/partial(ue)
       f00 = ones(size(f0,1),1);
       f01{1} = zeros(size(f1{1},1),1); 
       f01{2} = zeros(size(f1{1},1),1);
       
       
       %f_1,i = partial(f_i)/partial(gradue)
       f10{1} = zeros(size(f0,1),1);
       f10{2} = zeros(size(f0,1),1);
       f11{1} = zeros(size(f1{1},1),1);
       f11{2} = zeros(size(f1{1},1),1);       
       
%      %==== -\nabla^2 (u) + u = g  with 0 Dirichlet B.C. ====%
%      f0 = ue - g; 
%      f1 = grad_ue;
%      %=====================================================%
%        
%      %f0 = @(ue, grad_ue) ue-g;
%      %f1 = @(ue, grad_ue) grad_ue;

       %f_i,0 = partial(f_i)/partial(ue)
%      f00 = ones(size(f0,1),1);
%      f10{1} = zeros(size(f1{1},1)); 
%      f10{2} = zeros(size(f1{1},1)); 
%      
%      f_i,1 = partial(f_i)/partial(gradue)
%      f01 = zeros(size(f0,1),1);        
%      f11{1} = ones(size(f1{1},1));
%      f11{2} = ones(size(f1{1},1));
       
             
end


