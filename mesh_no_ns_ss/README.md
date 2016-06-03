The mesh included in this folder have been gnerated by cubit.

The files are exodus format. 

These mesh files contain no boundary nodesets or sideset. Good for running 
L2 Projection problem.

The file names include info about the mesh:

nameX_YYYYe_us.exo or nameX_YYYYe_s.exo
filename[elem_type]_[num_elem]e_[structured/unstructured]

elem_type: 8 & 27 are 3D and 4 & 9 are 2D

us: unstructured mesh

s: structured mesh

Example:

cube8_110e_us.exo

filename: cube

element type: 8 

num_elems: 110

us: unstructured mesh

