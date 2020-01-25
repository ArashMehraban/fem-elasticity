function vtk_file = processVtk(vtk_file, conn, mesh,u)
  fid = fopen(vtk_file, 'w');
  fprintf(fid, '# vtk DataFile Version 2.0\n');
  fprintf(fid, vtk_file);
  fprintf(fid, '\nASCII\n');
  fprintf(fid, 'DATASET UNSTRUCTURED_GRID\n\n');
  fprintf(fid, 'POINTS %d float\n',mesh.num_nodes);
  if(mesh.num_dims == 2)
      u = [u,zeros(size(u,1),1)];
  end
  myu = u';
  myu = myu(:);
  myu = reshape(myu,[],3);
  
  if(mesh.num_dims == 2)
      vtx = [(mesh.vtx_coords),zeros(size(mesh.vtx_coords,1),1)]';
  else
      vtx = (mesh.vtx_coords)';
  end
  vtx = vtx(:);
  vtx = reshape(vtx,[],3);  
  vtx = vtx +myu;
  f = '%f ';
  fs = repmat(f,1,3);
  fs = strcat(fs,'\n');
  fprintf(fid, fs ,vtx);
  fprintf(fid, '\nCELLS %d %d\n',mesh.num_elem,mesh.num_elem*(mesh.num_nodes_per_elem+1));
  vecNumElems = uint64(mesh.num_nodes_per_elem * ones(mesh.num_elem,1));
  theConn = [vecNumElems, conn-1];
  conny = theConn';
  conny = conny(:);
  conny = reshape(conny,[],mesh.num_nodes_per_elem +1);
  d = '%d ';
  ds = repmat(d,1,mesh.num_nodes_per_elem +1);
  r = strcat(ds,'\n');
  fprintf(fid, r,conny);
  fprintf(fid, '\nCELL_TYPES %d\n',mesh.num_elem);
  if (mesh.num_dims == 2 && mesh.num_nodes_per_elem == 4)
      cellType = 9;
  elseif (mesh.num_dims == 2 && mesh.num_nodes_per_elem == 9) 
      cellType = 23;
  elseif (mesh.num_dims == 3 && mesh.num_nodes_per_elem == 8)
      cellType = 12;
  elseif (mesh.num_dims == 3 && mesh.num_nodes_per_elem == 27)
      cellType = 25;
  else
      disp('I do NOT know this mesh to prepare for a vtk file');
  end
  
  ctype = uint8(cellType * ones(mesh.num_elem,1));
  fprintf(fid,'%d\n', ctype);
  
  fprintf(fid, '\nPOINT_DATA %d\n',mesh.num_nodes);
  fprintf(fid,'VECTORS displacement float\n');
  fprintf(fid,fs, myu);
  
  fclose(fid);

end


