function exctSol = exactf(vtx_coords)
  x = vtx_coords(:,1);
  y = vtx_coords(:,2);
  g=@(x,y)tanh(x).*exp(y)+sin(y);
  
  exctSol= g(x,y);

end