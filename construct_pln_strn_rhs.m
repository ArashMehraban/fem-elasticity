function [g1,g2, g1_dx, g1_dy,g2_dx, g2_dy] = construct_pln_strn_rhs()
% g1 and g2 should be used in userf as constructed RHS
      
   u{1}=@(x,y)tanh(x).*exp(y)+sin(y);
   u{2}=@(x,y)tanh(x).*cos(y); 
   syms x y nu E
   strain =[diff(u{1},x); diff(u{1},y); diff(u{2},x)+diff(u{1},y)];
   C =(E/((1+nu)*(1-2*nu)))*[1-nu, nu, 0; nu,1-nu,0 ; 0 ,0 ,0.5*(1-2*nu)];
   sigma = C*strain;
   g1_rhs = -(diff(sigma(1),x)+diff(sigma(3),y));
   g2_rhs = -(diff(sigma(3),x)+diff(sigma(2),y));
   g1_userf_dx  = diff(g1_rhs,x);
   g1_userf_dy  = diff(g1_rhs,y);
   g2_userf_dx  = diff(g2_rhs,x);
   g2_userf_dy  = diff(g2_rhs,y);   
   
   g1=matlabFunction(g1_rhs);
   g2=matlabFunction(g2_rhs);
   g1_dx = matlabFunction(g1_userf_dx);
   g1_dy = matlabFunction(g1_userf_dy);
   g2_dx = matlabFunction(g2_userf_dx);
   g2_dy = matlabFunction(g2_userf_dy);   
end