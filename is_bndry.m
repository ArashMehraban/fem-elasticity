function [is_bn, not_bndry] = is_bndry(sz_vtx_coords,bndry_nodes)
%SETBDRY_DIR returns a vector indicating Diricklet Boundray Values on each
%Node
%  0 indicates the node is on boundary
% -1 indicates the node is NOT on boundary

    is_bn =-ones(sz_vtx_coords,1);
    i = 1:size(bndry_nodes,2);
    is_bn(bndry_nodes(i)) = 0;
    
    not_bn_sz=sz_vtx_coords- size(bndry_nodes,2);
    i=1:not_bn_sz;
    not_bndry(i) = find(is_bn~= 0);

end