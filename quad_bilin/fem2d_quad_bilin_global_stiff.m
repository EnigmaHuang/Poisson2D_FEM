function K = fem2d_quad_bilin_global_stiff(coords, ien)
% Generate the global stiffness matrix for $-\nabla^2 u$ in 2D case
% using bilinear basis functions and quadrilateral elements
% [IN]  coords : n * 2 matrix, n is the number of grid points,
%                each row has the x and y coordinate of a grid point
% [IN]  ien    : m * 4 matrix, m is the number of quadrilateral elements,
%                each row has 4 vertex point ids of the element, in 
%                counter clockwise order
% [OUT] K      : n * n Global stiffness matrix
	
	n = size(coords, 1);
	m = size(ien, 1);
	
	K = sparse(n, n);
	
	for i_elem = 1 : m
		elem_vertex_ids = ien(i_elem, :);
		vertex_coords = coords(elem_vertex_ids, :);
		
		% Transpose vertex_coords to use in fem2d_quad_bilin_local_stiff
		k = fem2d_quad_bilin_local_stiff(vertex_coords');
		
		% Assemble to global stiffness matrix
		K(elem_vertex_ids, elem_vertex_ids) = K(elem_vertex_ids, elem_vertex_ids) + k;
	end
end