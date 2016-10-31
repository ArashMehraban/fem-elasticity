function J=jacF(func,x)
% computes the Jacobian of a function
    n=length(x);
    J =zeros(n,n);
    fx=feval(func,x);
    eps=1.e-8;  % could be made better
    xperturb=x;
    for i=1:n
        xperturb(i)=xperturb(i)+eps;
        J(:,i)=(feval(func,xperturb)-fx)/eps;
        xperturb(i)=x(i);
    end
    
end 