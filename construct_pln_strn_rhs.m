function [g1,g2] = construct_pln_strn_rhs()
% g1 and g2 should be used in userf as constructed RHS
      
   u{1}=@(x,y)tanh(x).*exp(y)+sin(y);
   u{2}=@(x,y)tanh(x).*cos(y); 
   syms x y nu E
   strain =[diff(u{1},x); diff(u{1},y); diff(u{2},x)+diff(u{1},y)];
   C =(E/((1+nu)*(1-2*nu)))*[1-nu, nu, 0;
                             nu,1-nu,0 ;
                             0 ,0 ,(1-2*nu)/2];
   sigma = C*strain;
   g1=matlabFunction(-(diff(sigma(1),x)+diff(sigma(3),y)));
   g2=matlabFunction(-(diff(sigma(3),x)+diff(sigma(2),y)));
end