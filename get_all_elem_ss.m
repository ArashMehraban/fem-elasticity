function elem_ss = get_all_elem_ss(msh)
%GET_ALL_ELEM_SS returns all elements corresponding to sidesets in one cell array.
%input: mesh
%output: All elements for sidesets
   
    fld_names = fieldnames(msh);
    elem_names = fld_names((strncmp(fld_names,'elem_ss',7)));
    sz_elem_names = size(elem_names,1);
    sz_delem_ss = zeros(sz_elem_names,1);
    for i=1:sz_elem_names
       sz_delem_ss(i) = size(msh.(elem_names{i}),1);
    end
    elem_ss_vec = zeros(sum(sz_delem_ss),1);
    elem_ss = mat2cell(elem_ss_vec, sz_delem_ss);
      
    for i=1:sz_elem_names
        elem_ss{i} = msh.(elem_names{i});
    end
   
end