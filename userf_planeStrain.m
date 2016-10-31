function [f0,f1] = userf_planeStrain(ue, grad_ue, xe) 
%USERF_PLANESTRAIN provides weak form of the Plain Strain problem to solve 
%
%  input:      ue: corresponding u for each element evalauted at quadrature points
%       : grad_ue: corresponding grad_u for each element evalauted at quadrature points
%       :      xe: quadrature points mapped to the reference elem
%
% output:  f0: any possible source from given problem
%       :  f1: any possible source from given problem



    x=xe(:,1);
    y=xe(:,2);
    
    
    % grad_ue structure:
    % 2D:       [D1*u1 | D1*u_2]     [du1/dx | du2/dx] 
    %   grad_ue=[D2*u1 | D2*u_2]  =  [du1/dy | du2/dy]
    
    
    % Young's modulus 
    E = 1;%2e11;
    % Poisson ratio 
    nu = 0.3;
        
    num_row = size(ue,1);
    
    %strain-stress matrix coefficient
    ss_coef = E/((1+nu)*(1-2*nu));
    
    %strain
    eps11 = grad_ue(1:num_row,1);
    eps12 = 0.5*(grad_ue(num_row+1:end,1)+grad_ue(1:num_row,2));
    eps22 = grad_ue(num_row+1:end,2);
    
    sigma11 = ss_coef*((1-nu)*eps11 + nu*eps22);
    sigma22 = ss_coef*(nu*eps11 +(1-nu)*eps22);  
    sigma12 = ss_coef*(0.5*(1-2*nu)*eps12);
   

    % RHS: manufactured from construct_pln_strn_rhs.m
    %g1= -(E.*(nu-1.0./2.0).*(sin(y).*-2.0+exp(y).*tanh(x).*2.0+sin(y).*(tanh(x).^2-1.0).*2.0))./((nu.*2.0-1.0).*(nu+1.0))+(E.*nu.*sin(y).*(tanh(x).^2-1.0))./((nu.*2.0-1.0).*(nu+1.0))-(E.*exp(y).*tanh(x).*(tanh(x).^2-1.0).*(nu-1.0).*2.0)./((nu.*2.0-1.0).*(nu+1.0));
    %g2= (E.*(exp(y).*(tanh(x).^2-1.0).*2.0-cos(y).*tanh(x).*(tanh(x).^2-1.0).*4.0).*(nu-1.0./2.0))./((nu.*2.0-1.0).*(nu+1.0))-(E.*nu.*exp(y).*(tanh(x).^2-1.0))./((nu.*2.0-1.0).*(nu+1.0))+(E.*cos(y).*tanh(x).*(nu-1.0))./((nu.*2.0-1.0).*(nu+1.0));
     
    %const g
    %g1=(E.*(nu-1.0).*-2.0)./((nu.*2.0-1.0).*(nu+1.0))-(E.*(nu-1.0./2.0))./((nu.*2.0-1.0).*(nu+1.0))+0*x;
    %g2=(E.*(nu-1.0).*1.0)./((nu.*2.0-1.0).*(nu+1.0))-(E.*(nu-1.0./2.0).*2.0)./((nu.*2.0-1.0).*(nu+1.0))+0*y;
    
    %linear g
    %g1=(E.*x.*(nu-1.0).*-2.0)./((nu.*2.0-1.0).*(nu+1.0))-(E.*y.*(nu-1.0./2.0).*6.0)./((nu.*2.0-1.0).*(nu+1.0));
    %g2=(E.*x.*(nu-1.0./2.0).*-6.0)./((nu.*2.0-1.0).*(nu+1.0))+(E.*y.*(nu-1.0).*3.0)./((nu.*2.0-1.0).*(nu+1.0));
    
    %quadratic g
    %g1=(E.*x.^2.*(nu-1.0).*-3.0)./((nu.*2.0-1.0).*(nu+1.0))-(E.*y.^2.*(nu-1.0./2.0).*1.2e1)./((nu.*2.0-1.0).*(nu+1.0));
    %g2=(E.*x.^2.*(nu-1.0./2.0).*-1.8e1)./((nu.*2.0-1.0).*(nu+1.0))+(E.*y.^2.*(nu-1.0).*9.0)./((nu.*2.0-1.0).*(nu+1.0));
    
    %g1=(E.*(nu-1.0./2.0).*(cos(x).*3.75e-1+x.*sin(y).*2.5e-1-3.0))./((nu.*2.0-1.0).*(nu+1.0))-(E.*nu.*cos(x).*7.5e-1)./((nu.*2.0-1.0).*(nu+1.0));
    %g2=-(E.*(nu-1.0./2.0).*(cos(y).*2.5e-1+y.*sin(x).*3.75e-1+1.0))./((nu.*2.0-1.0).*(nu+1.0))+(E.*nu.*cos(y).*5.0e-1)./((nu.*2.0-1.0).*(nu+1.0));
    
    g1=(E.*(nu-1.0./2.0).*(exp(y).*cos(x).*3.75e-1+exp(x).*sin(y).*2.5e-1-3.0))./((nu.*2.0-1.0).*(nu+1.0))-(E.*nu.*exp(y).*cos(x).*7.5e-1)./((nu.*2.0-1.0).*(nu+1.0))-(E.*exp(x).*sin(y).*(nu-1.0).*5.0e-1)./((nu.*2.0-1.0).*(nu+1.0));
    g2=-(E.*(nu-1.0./2.0).*(exp(x).*cos(y).*2.5e-1+exp(y).*sin(x).*3.75e-1+1.0))./((nu.*2.0-1.0).*(nu+1.0))+(E.*nu.*exp(x).*cos(y).*5.0e-1)./((nu.*2.0-1.0).*(nu+1.0))+(E.*exp(y).*sin(x).*(nu-1.0).*7.5e-1)./((nu.*2.0-1.0).*(nu+1.0));
    
    f0 = [-g1, -g2]; 
    f1(:,1) = [sigma11;sigma12];
    f1(:,2) = [sigma12;sigma22]; 
end