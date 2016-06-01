function [f0, f1,f00, f01, f10, f11] = get_userf(ue, grad_ue, xe)
%USERF provides weak form of the problem to solve 

%  input:      ue: corresponding u for each element evalauted at quadrature points
%       : grad_ue: corresponding grad_u for each element evalauted at quadrature points (cell not matrix)
%       :      xe: quadrature points mapped to the reference elem
%
% output: f0: any possible source from given problem
%       : f1: any possible source from given problem

    x=xe(:,1);
    y=xe(:,2);
    sz_dim = size(xe,2);
    if(sz_dim == 3)
        z = xe(:,3);
    end
    
     
     
%         % For L2-Projection problem:
%         % Weak Form: (v: test/weight function)
%         % integral(v*f0(ue,grad_ue) + grad_v : f1(ue,grad_ue)) = 0   
%       
%        %==== L2 ======%
%        %user defined rhs
%        g = tanh(x).*exp(y)+sin(y);
%        f0 = ue - g; 
%       
%        f1=cell(1,size(grad_ue,2));
%        for  j=1:size(grad_ue,2)
%            f1{j} = 0*grad_ue{j};
%        end
%        %==============%
%        
%        %f_0,i = partial(f_i)/partial(ue)
%        f00 = ones(size(f0,1),1);
%        f01{1} = zeros(size(f1{1},1),1); 
%        f01{2} = zeros(size(f1{1},1),1);
%        
%        
%        %f_1,i = partial(f_i)/partial(gradue)
%        f10{1} = zeros(size(f0,1),1);
%        f10{2} = zeros(size(f0,1),1);
%        f11{1} = zeros(size(f1{1},1),1);
%        f11{2} = zeros(size(f1{1},1),1);   
     
   if(sz_dim == 2)
     %2D  
     %Poisson eqn:
     g=sin(y)+exp(y).*tanh(x)-exp(y).*tanh(x).^3.*2.0;
     %==== -\nabla^2 (u) = g  with 0 Dirichlet B.C. ====%
     f0{1} =  0*ue -g; 
     f1{1} = grad_ue{1};
     f1{2} = grad_ue{2};
     %=====================================================%
       
     %f0 = @(ue, grad_ue) ue-g;
     %f1 = @(ue, grad_ue) grad_ue;

     %f_i,0 = partial(f_i)/partial(ue)
     f00{1} = zeros(size(f0{1},1),1);
     f10{1} = zeros(size(f1{1},1),1); 
     f10{2} = zeros(size(f1{1},1),1); 
     
     %f_i,1 = partial(f_i)/partial(gradue)
     f01{1} = zeros(size(f0{1},1),1);   
     f01{2} = zeros(size(f1{1},1),1);
     f11{1} = ones(size(f1{1},1),1);
     f11{2} = ones(size(f1{1},1),1);
   end
   
   if(sz_dim == 3)
     %2D  
     %Poisson eqn:
     g=cos(z)+sin(y)+exp(y).*tanh(x)-exp(y).*tanh(x).^3.*2.0;
     %==== -\nabla^2 (u) = g  with 0 Dirichlet B.C. ====%
     f0{1} =  0*ue -g; 
     f1{1} = grad_ue{1};
     f1{2} = grad_ue{2};
     f1{3} = grad_ue{3};
     %=====================================================%
       
     %f0 = @(ue, grad_ue) ue-g;
     %f1 = @(ue, grad_ue) grad_ue;

     %f_i,0 = partial(f_i)/partial(ue)
     f00{1} = zeros(size(f0{1},1),1);
     f10{1} = zeros(size(f1{1},1),1); 
     f10{2} = zeros(size(f1{2},1),1); 
     f10{3} = zeros(size(f1{3},1),1);
     
     %f_i,1 = partial(f_i)/partial(gradue)
     f01{1} = zeros(size(f1{1},1),1);   
     f01{2} = zeros(size(f1{2},1),1);
     f01{3} = zeros(size(f1{3},1),1);
     f11{1} = ones(size(f1{1},1),1);
     f11{2} = ones(size(f1{2},1),1);
     f11{3} = ones(size(f1{3},1),1);
   end
   
   
   
   
       
% %      %==== -\nabla^2 (u) + u = g  with 0 Dirichlet B.C. ====%
% %      f0 = ue - g; 
% %      f1 = grad_ue;
% %      %=====================================================%
% %        
% %      %f0 = @(ue, grad_ue) ue-g;
% %      %f1 = @(ue, grad_ue) grad_ue;
% % 
% %      %f_i,0 = partial(f_i)/partial(ue)
% %      f00 = ones(size(f0,1),1);
% %      f10{1} = zeros(size(f1{1},1)); 
% %      f10{2} = zeros(size(f1{1},1)); 
% %      
% %      %f_i,1 = partial(f_i)/partial(gradue)
% %      f01 = zeros(size(f0,1),1);        
% %      f11{1} = ones(size(f1{1},1));
% %      f11{2} = ones(size(f1{1},1));
   
% % % 
% % %    
% % %    %=== Plane Strain Problem =====%
% % %    
% % %    %Young's modulus 
% % %    E = 2;
% % %    % Poisson ratio 
% % %    nu = 0.3;
% % %    %strain/stress matrix
% % %    C =(E/((1+nu)*(1-2*nu)))*[1-nu,nu,0; nu,1-nu,0; 0,0,(1-2*nu)/2];
% % %     
% % %    %grad_ue{1} = [partial_u1/partial_x , partial_u2/partial_x ]
% % %    %grad_ue{2} = [partial_u1/partial_y , partial_u2/partial_y ]
% % %    
% % %    %strain = [partial_u1/partial_x, partial_u2/partial_y, 1/2*(partial_u2/partial_x + partial_u1/partial_y)]   
% % %    strain = [grad_ue{1}(:,1), grad_ue{2}(:,2), 0.5*(grad_ue{1}(:,2)+grad_ue{2}(:,1))];     
% % %     
% % %    %stress (sigma) = (strain/stress matrix)*strain 
% % %    sigma = strain*C';
% % %    
% % %    g1= -(E.*(nu-1.0./2.0).*(sin(y).*(-1.0./2.0)+exp(y).*tanh(x).*(1.0./2.0)+sin(y).*(tanh(x).^2-1.0).*(1.0./2.0)))./((nu.*2.0-1.0).*(nu+1.0))+(E.*nu.*sin(y).*(tanh(x).^2-1.0))./((nu.*2.0-1.0).*(nu+1.0))-(E.*exp(y).*tanh(x).*(tanh(x).^2-1.0).*(nu-1.0).*2.0)./((nu.*2.0-1.0).*(nu+1.0));
% % %    g2= (E.*(exp(y).*(tanh(x).^2-1.0).*(1.0./2.0)-cos(y).*tanh(x).*(tanh(x).^2-1.0)).*(nu-1.0./2.0))./((nu.*2.0-1.0).*(nu+1.0))-(E.*nu.*exp(y).*(tanh(x).^2-1.0))./((nu.*2.0-1.0).*(nu+1.0))+(E.*cos(y).*tanh(x).*(nu-1.0))./((nu.*2.0-1.0).*(nu+1.0));
% % %    
% % % 
% % %    f0{1} = g1;
% % %    f0{2} = g2;
% % %     
% % %    f1{1} = [sigma(:,1),sigma(:,3)] ;
% % %    f1{2} = [sigma(:,3),sigma(:,2)];
% % %    
% % %    %f00{1} = [partial_u1/partial_x, partial_u1/partial_y]
% % %    %f00{2} = [partial_u2/partial_x, partial_u2/partial_y]
% % %    f00{1} = zeros(size(f0{1},1),2);
% % %    f00{2} = zeros(size(f0{2},1),2);   
% % %    
% % %    f01{1} = zeros(size(f0{1},1),2);
% % %    f01{2} = zeros(size(f0{2},1),2);
% % %    
% % %    f10{1} = zeros(size(f1{2},1),2);
% % %    f10{2} = zeros(size(f1{2},1),2);
% % %    
% % %    one = ones(size(f1{1},1),3);
% % %    f1_gu = one*C';
% % %    f11{1}=[f1_gu(:,1),f1_gu(:,3)];
% % %    f11{2}=[f1_gu(:,3),f1_gu(:,2)];       
             
end


