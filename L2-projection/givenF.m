syms x y

f = tanh(x)*exp(y)+sin(y);
d2x = diff(f,x,2);
d2y = diff(f,y,2);
l = simplify(d2x+d2y);