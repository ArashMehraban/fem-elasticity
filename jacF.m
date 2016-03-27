function J=jacF(func,x)
% JACF computes the Jacobian of a function using Finite differences
    
    n=length(x);
    % Allocate space for teh Jacobian
    J =zeros(n,n);
    fx=feval(func,x);
    eps=1.e-8; 
    xperturb=x;
    for i=1:n
        xperturb(i)=xperturb(i)+eps;
        J(:,i)=(feval(func,xperturb)-fx)/eps;
        xperturb(i)=x(i);
    end
    
end