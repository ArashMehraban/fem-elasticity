function dir_ns = get_all_dir_ns(msh)
%GET_ALL_DIR_NS returns all of nodesets in one cell array.
%input: mesh
%output: All nodestes
   
    fld_names = fieldnames(msh);
    ns_names = fld_names((strncmp(fld_names,'node_ns',7)));
    sz_ns_names = size(ns_names,1);
    sz_dir_ns = zeros(sz_ns_names,1);
    for i=1:sz_ns_names
       sz_dir_ns(i) = size(msh.(ns_names{i}),1);
    end
    dir_ns_vec = zeros(sum(sz_dir_ns),1);
    dir_ns = mat2cell(dir_ns_vec, sz_dir_ns);
      
    for i=1:sz_ns_names
        dir_ns{i} = msh.(ns_names{i});
    end
   
end