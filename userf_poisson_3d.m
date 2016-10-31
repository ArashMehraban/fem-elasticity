function [f0,f1] = userf_poisson_3d(ue, grad_ue, xe) 

    x=xe(:,1);
    y=xe(:,2);
    z=xe(:,3);

    %RHS: manufactured from construct_neg_laplacian_rhs_3d.m 
    g=cos(z)+sin(y)+exp(y).*tanh(x)-exp(y).*tanh(x).^3.*2.0;

    f0 = 0*ue -g; 
    f1 = grad_ue;
end