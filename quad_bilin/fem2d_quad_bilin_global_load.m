function b = fem2d_quad_bilin_global_load(coords, ien)
% Generate the right hand side for $-\nabla^2 u = f$ in 2D case
% using bilinear basis functions and quadrilateral elements
% [IN]  coords : n * 2 matrix, n is the number of grid points,
%                each row has the x and y coordinate of a grid point
% [IN]  ien    : m * 4 matrix, m is the number of quadrilateral elements,
%                each row has 4 vertex point ids of the element, in 
%                counter clockwise order 
% [OUT] n      : n * 1 RHS vector
	
	n = size(coords, 1);
	m = size(ien, 1);
	
	b = zeros(n, 1);
	
	for i_elem = 1 : m
		elem_vertex_ids = ien(i_elem, :);
		vertex_coords = coords(elem_vertex_ids, :)';
		
		ub = fem2d_quad_bilin_unit_load(vertex_coords);
		b(elem_vertex_ids) = b(elem_vertex_ids) + ub;
	end
end