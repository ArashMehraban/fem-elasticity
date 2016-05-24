function [g1,g2] = construct_pln_strn_rhs()
% g1 and g2 should be used in userf as constructed RHS
   
   % u is displpacement
   u{1}=@(x,y)tanh(x).*exp(y)+sin(y);
   u{2}=@(x,y)tanh(x).*cos(y); 
 
   syms x y nu E 
   
   C =(E/((1+nu)*(1-2*nu)))*[1-nu, nu, 0; nu,1-nu,0 ; 0 ,0 ,0.5*(1-2*nu)];   
   strain = [diff(u{1},x) ; diff(u{2},y) ; 0.5*(diff(u{2},x) + diff(u{1},y))];
   
   sigma = C * strain;
   
   g1_rhs = -(diff(sigma(1),x) + diff(sigma(3),y)); 
   g2_rhs = -(diff(sigma(3),x) + diff(sigma(2),y));
    
   g1=matlabFunction(g1_rhs);
   g2=matlabFunction(g2_rhs);
end