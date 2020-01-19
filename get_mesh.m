function [origConn, msh] = get_mesh(filename, ext, varargin)
%GET_MESH function receives an exodus format file (.exo or .e) as input and
%returns as output a mesh object that can be used in MATLAB with finite
%element code. The output mesh object contains:
%
% vtx_coords        : vertex coordinates for all nodes in the mesh
% conn              : connectivity matrix for vtx_coords
% num_elem          : number of elements
% num_nodes_per_elem: number of nodes per element
% num_nodes         : total number of nodes
% num_dims          : number of dimensions for mesh (2 = 2D or 3 = 3D)
% num_node_sets     : number of node sets
% num_side_sets     : number of side sets
% node_ns(i)        : (i) indicates multiple sets. node_ns1 means node set 1
% side_ss(i)        : (i) indicates multiple sets. side_ss1 means side set 1
% elem_ss(i)        : (i) indicates multiple sets. elem_ss1 means elements set 1
%                   : elem_ss corresponds to side_ss. Example: elem_ss1
%                     contains the elements that have side_ss1 
% side_ss(i)_nodes  : (i) indicates multiple sets. side_ss_nodes returns the
%                     nodes that are on each entry in side_ss 
%
% Note: This code handles 1 block element type only!! For more info about blocks, 
% refer to cubit manuals (keyword: num_el_blk).
%
% WARNING: This code is desinged (and tested) to handle exodus files 
% produced by cubit 12.1 only. It is expected to work with higher versions of
% cubit. This code can handle the following element types output by cubit:
% 1)  4-noded QAUD element (QUAD4 in cubit: 2D element)
% 2)  9-noded QAUD element (QUAD9 in cubit: 2D element)
% 3)  8-noded HEX element  ( HEX8 in cubit: 3D element)
% 4) 27-noded HEX element  (HEX27 in cubit: 3D element)
%
%
% Mesh ordering may be reordered lexicographically by the user if desired using
% the keyword 'lex' as an optional input parameter described below:
% 
% inputs: (necessary): filename 
%       : (necessary): ext (file extension) 
%       : (optional):  lex (optional for lexicographical ordering) 
%
% output: mesh object (struct)
%   
% example on how to use the function: 
%    filename = 'disk';
%    ext = 'exo';
%    msh = get_mesh(filename,ext);  
%    or
%    msh = get_mesh(filename,ext,'lex');  (for lexicographical ordering)
%
% Developed by Arash Mehraban, University of Colorado, Boulder (April 2016)


    %open exodus file
    if(~(strcmp(ext,'exo') || strcmp(ext,'e')))
        error('Wrong file format input. get_mesh can read .exo and .e files only!!')
    end
    
    %get file handle
    ncid = netcdf.open(strcat(filename,'.',ext));

    % netCDF files have three header sections:
    % 1) dimensions
    % 2) variables
    % 3) attributes (global) (Not needed/extracted for get_mesh function)

    %get the number of dimensions and variables 
    [numdims, numvars, ~, ~] = netcdf.inq(ncid);

    %Allocate space for netCDF file dimensions 
    dimen_names = cell(numdims,1);
    dimen_lens = zeros(numdims,1);
    dimen_ID = zeros(numdims,1);

    %get dimensions information
    for i=0:numdims-1    
        [dim_names, dimlen] = netcdf.inqDim(ncid,i);
        dimen_ID(i+1) = netcdf.inqDimID(ncid,dim_names);
        dimen_names{i+1} = dim_names;
        dimen_lens(i+1) = dimlen;     
    end

    %Allocate space for netCDF file variables 
    var_names = cell(numvars,1);
    var_ID = zeros(numvars,1);
    
    %Allocate space for data extracted from netCDF variables
    mesh_data = cell(numvars,1);

    %get variables information
    for i=0:numvars-1   
        var_names{i+1} = netcdf.inqVar(ncid,i);
        var_ID(i+1) = netcdf.inqVarID(ncid,var_names{i+1});
        mesh_data{i+1} = netcdf.getVar(ncid, var_ID(i+1));
    end
    
    %build a map of all var_names
    M = containers.Map(var_names,zeros(size(var_names,1),1));
    
    %get number of nodes in element
    num_node_per_elem = dimen_lens((strncmp(dimen_names,'num_nod_per_el',11)));
    %get connectivity matrix
    conn = mesh_data{strncmp(var_names,'connect',7)}';
    origConn = conn;
    
    %get number of elelments
    num_elem = dimen_lens(strcmp(dimen_names,'num_elem')); 
    %get number of nodes in mesh
    num_nods = dimen_lens(strcmp(dimen_names,'num_nodes'));
    %get the problem dimension (2D, 3D)
    num_dim = dimen_lens(strcmp(dimen_names,'num_dim'));
    %get vertx coordinates
    TF = isKey(M,'coord');
    if(TF==1)
        coords = mesh_data{strcmp(var_names,'coord')};
    elseif(num_dim == 1)
        coords = mesh_data{strcmp(var_names,'coordx')};        
    elseif(num_dim == 2)
        coords = [mesh_data{strcmp(var_names,'coordx')},mesh_data{strcmp(var_names,'coordy')}];
    elseif(num_dim == 3)
        coords = [mesh_data{strcmp(var_names,'coordx')},mesh_data{strcmp(var_names,'coordy')}, mesh_data{strcmp(var_names,'coordz')}];        
    end
    
    
    %get the array indecies of the node sets in netCDF variables
    ns_idx= find(strncmp(var_names,'node_ns',7));
    %get the array indecies of the side sets in netCDF variables
    ss_idx= find(strncmp(var_names,'side_ss',7));
    %get the array indecies of the element sets in netCDF variables
    elem_idx= find(strncmp(var_names,'elem_ss',7));
    
    %get size of node, side and element sets
    sz_ns_idx = size(ns_idx,1);
    sz_ss_idx = size(ss_idx,1);
    sz_elem_idx =size(elem_idx,1);
    
    %Allocate space for each node set
    node_set = cell(1,sz_ns_idx);
    %populate each node set
    for i =1:sz_ns_idx
        node_set{i} = mesh_data{ns_idx(i)};
    end
     
    %Allocate space for side sets sizes
    side_set_szs = zeros(1,sz_ss_idx);
    %Allocate space for each side set
    side_set = cell(1,sz_ss_idx);
    %populate side sets and array of side set sizes
    for i =1:sz_ss_idx
        side_set{i} = mesh_data{ss_idx(i)};
        side_set_szs(i) = size(side_set{i},1);
    end
    
    %Allocate space for each element set (correspondin to side sets)
    elem_side_set = cell(1,sz_elem_idx);
    %populate each element set (correspondin to side sets)
    for i =1:sz_elem_idx
        elem_side_set{i} = mesh_data{elem_idx(i)};
    end
    
    %get the number of nodes in all side sets combined
    ss_nod_sz = sum(side_set_szs);
    %Allocate space for number of nodes in all side sets combined
    side_set_nodes=cell(1,ss_nod_sz);    
    
    % check if lexicographical (LEX) ordering mesh is requested
    if(nargin >2)
        lex = varargin{1};
        if(~strcmp(lex,'lex'))
            error('Ordering was not recognized. Choose lex for lexicographical ordering!')
        end
        
        %node mapping for LEX ordering           
        if(num_node_per_elem == 4)
            lex_conn = [1 2 4 3];            
            % Side Set Node Ordering (indecies)
            QUAD = [1 2;
                    2 3;
                    3 4;
                    4 1];
        elseif((num_node_per_elem == 9))
            lex_conn = [1 5 2 8 9 6 4 7 3]; %used originally
            % Side Set Node Ordering (indecies)
            QUAD = [1 5 2;
                    2 6 3;
                    3 7 4;
                    4 8 1];
        elseif((num_node_per_elem == 8))
            lex_conn = [1 2 4 3 5 6 8 7];
            % Side Set Node Ordering (indecies)
            HEX = [1 2 5 6;
                   2 3 6 7;
                   3 4 7 8;
                   1 5 4 8;
                   1 4 2 3;
                   5 6 8 7];
        elseif((num_node_per_elem == 27))
            lex_conn = [1 9 2 12 22 10 4 11 3 13 26 14 24 21 25 16 27 15 5 17 6 20 23 18 8 19 7];
            % Side Set Node Ordering (indecies)
            HEX = [1 9  2 13 26 14 5 17 6;
                   2 10 3 14 25 15 6 18 7;
                   3 11 4 15 27 16 7 19 8;
                   1 13 5 12 24 20 4 16 8;
                   1 12 4  9 22 11 2 10 3;
                   5 17 6 20 23 18 8 19 7];
        end
        
        % permute conectivity matrix to LEX ordering 
        i=1:num_elem;
        conn(i,:) = conn(i,lex_conn);
        
        %get nodes on side sets (in cell format per element)
        k=1;
        for i =1:sz_ss_idx
            side_set_sides = side_set{i};
            cur_elem_side_set = elem_side_set{i};
            for j = 1:size(side_set_sides,1)
                elem_side_set_nodes_idx = side_set_sides(j);
                if(num_node_per_elem == 9)
                    local_elem_side_nodes_idx = QUAD(elem_side_set_nodes_idx,:);
                elseif(num_node_per_elem == 4)
                    local_elem_side_nodes_idx = QUAD(elem_side_set_nodes_idx,:);
                elseif(num_node_per_elem == 27)
                    local_elem_side_nodes_idx = HEX(elem_side_set_nodes_idx,:);
                 elseif(num_node_per_elem == 8)
                     local_elem_side_nodes_idx = HEX(elem_side_set_nodes_idx,:);
                end
                side_set_nodes{k} = conn(cur_elem_side_set(j),local_elem_side_nodes_idx);
                k=k+1;
            end   
        end
          
    else
        % Side Set Node Ordering (indecies)
        QUAD =[1, 2, 5; 
               2, 3, 6;
               3, 4, 7; 
               4, 1, 8];
        HEX = [1, 2, 6, 5, 9,  14, 17, 13, 26;
               2, 3, 7, 6, 10, 15, 18, 14, 25;
               3, 4, 8, 7, 11, 16, 19, 15, 27;
               1, 5, 8, 4, 13, 20, 16, 12, 24;
               1, 4, 3, 2, 12, 11, 10, 9,  22;
               5, 6, 7, 8, 17, 18, 19, 20, 23];       


        %%get connectivity matrix mapped to geometry in cubit screen:
        %node_num_map = mesh_data{strcmp(var_names,'node_num_map')}';
        % connect = mesh_data{connect_idx};
        % connect_reshape = connect(:);
        % i=1:size(connect_reshape,1);
        % conn_reshape = node_num_map(connect(i));
        % [r_connect, c_connect] = size(connect);
        % conn_transpose = reshape(conn_reshape,r_connect,c_connect);
        % conn = conn_transpose';

        %get nodes on side sets (in cell format per element)
        k=1;
        for i =1:sz_ss_idx
            side_set_sides = side_set{i};
            cur_elem_side_set = elem_side_set{i};
            for j = 1:size(side_set_sides,1)
                elem_side_set_nodes_idx = side_set_sides(j);
                if(num_node_per_elem == 9)
                    local_elem_side_nodes_idx = QUAD(elem_side_set_nodes_idx,:);
                elseif(num_node_per_elem == 4)
                    local_elem_side_nodes_idx = QUAD(elem_side_set_nodes_idx,1:end-1);
                elseif(num_node_per_elem == 27)
                    local_elem_side_nodes_idx = HEX(elem_side_set_nodes_idx,:);
                 elseif(num_node_per_elem == 8)
                     local_elem_side_nodes_idx = HEX(elem_side_set_nodes_idx,1:4);
                end
                side_set_nodes{k} = conn(cur_elem_side_set(j),local_elem_side_nodes_idx);
                k=k+1;
            end   
        end
    
    end
    % convert sides sets to matrix format
    side_set_nodes = cell2mat(side_set_nodes)';   
   
   % get node set variable names from netCDF variables
   ns_names=cell(1,sz_ns_idx);
   for i=1:sz_ns_idx
       ns_names{i}= var_names{ns_idx(i)};
   end
   
   % get element variable names from netCDF variables
   elem_names=cell(1,sz_elem_idx);
   for i=1:sz_elem_idx
       elem_names{i}= var_names{elem_idx(i)};
   end
   
   % get side set variable names from netCDF variables
   ss_names=cell(1,sz_ss_idx);
   for i=1:sz_ss_idx
       ss_names{i}= var_names{ss_idx(i)};
   end
   
   ss_node_names=cell(1,sz_ss_idx);
   for i=1:sz_ss_idx
       ss_node_names{i}= strcat(ss_names{i},'_nodes');
   end
   
   % seting common field names for mesh object (partial names)
   field_names={'vtx_coords','conn','num_elem','num_nodes_per_elem','num_nodes','num_dims','num_node_sets','num_side_sets'};
   
   %size of common field names 
   partial_fld_sz = size(field_names,2);
   
   %create a cell containing node sets, element sets and side sets names
   rest_of_field_names = {ns_names,elem_names,ss_names,ss_node_names};   
   
   %size of entities in rest_of_field_name 
   sz_entity_in_rest_field = size(rest_of_field_names,2);
   
   % setting rest of fields names after common field names in field_names
   i=partial_fld_sz+1;
   for k=1:sz_entity_in_rest_field
       for j=1:size(rest_of_field_names{k},2)
           field_names{i}= cell2mat(rest_of_field_names{k}(j));
           i=i+1;
       end
   end
   
   % seting common field values for mesh object (partial values)
   field_vals={coords,conn,num_elem,num_node_per_elem,num_nods,num_dim,sz_ns_idx,sz_ss_idx};
   
   %setting node set field values in field_vals after common values
   for i=1:sz_ns_idx
       field_vals{i+partial_fld_sz} = node_set{i};
   end
   
   %setting element set field values in field_vals after common values and
   %node set values
   for i=1:sz_elem_idx
       field_vals{i+partial_fld_sz+sz_ns_idx} = elem_side_set{i};
   end
   
   %setting side set field values in field_vals after common values,
   %node set values and element set values
   for i=1:sz_ss_idx
       field_vals{i+partial_fld_sz+sz_ns_idx+sz_elem_idx} = side_set{i};
   end
   
   % factor: number of entities in each row of side_set_nodes
   if(num_node_per_elem == 4)
       factor = 2;
   elseif(num_node_per_elem == 9)
       factor = 3;
   elseif(num_node_per_elem == 8)
       factor = 4;
   elseif(num_node_per_elem == 27)
       factor = 9;
   end  
   
   %setting side set node field values in field_vals after common values,
   %node set values, element values and side set values. Reshaped by conversion factor above
   k=1;   
   for j=1:size(side_set_szs,2)
       if(j==1)
            field_vals{partial_fld_sz+sz_elem_idx+sz_ns_idx+sz_ss_idx+j} = reshape(side_set_nodes(k:factor*side_set_szs(j)),factor,[])';
       else
           field_vals{partial_fld_sz+sz_elem_idx+sz_ns_idx+sz_ss_idx+j} = reshape(side_set_nodes(k:factor*side_set_szs(j-1)+factor*side_set_szs(j)),factor,[])';
       end
       k=factor*side_set_szs(j)+1;
   end
   
    %return mesh structure
    msh=struct();
    for i=1:size(field_names,2)
        msh.(field_names{i}) = field_vals{i};
    end

    %close the file
    netcdf.close(ncid);
end