%test Di's

clear
clc

u=@(x,y) y.*sin(x)+ x.*cos(y);
syms x y
u_dx = matlabFunction(diff(u,x));
u_dy = matlabFunction(diff(u,y));

vtx =[-1 -1; 1 ,-1; -1 1; 1 1];

elm_type = 4;
[B, Ds, W_hat] = get_shape(elm_type);
[dets, invJe] = jacobian(vtx, Ds);         
Di = get_elem_dirv(invJe, Ds);

x = vtx(:,1);
y = vtx(:,2);
fval = u(x,y);

Dx = Di{1}*fval;
Dy = Di{2}*fval;

n_gs_pts = 2;
[x, ~] = get_quadrature(n_gs_pts);
gs_pts = [x(1), -x(2); -x(1), -x(2);x(1), x(2);-x(1), x(2) ];
x = gs_pts(:,1);
y = gs_pts(:,2);

ux = u_dx(x,y);
uy = u_dy(x,y);
diff_x = Dx - ux;
diff_y = Dy - uy;







