function b = fem2d_quad_bilin_rhs(coords, ien, gpie)
% Generate the right hand side for $-\nabla^2 u = f$ in 2D case
% using bilinear basis functions and quadrilateral elements
% [IN]  coords : n * 2 matrix, n is the number of grid points,
%                each row has the x and y coordinate of a grid point
% [IN]  ien    : m * 4 matrix, m is the number of quadrilateral elements,
%                each row has 4 vertex point ids of the element, in 
%                counter clockwise order
% [IN]  gpie   : Grid Poinr In Element, each row is the information of a 
%                grid point: gpie(i, 1) is the number of elements this grid
%                point belongs to, gpie(i, 2:end) are the element ids 
% [OUT] n      : n * 1 RHS vector
	
	n = size(coords, 1);
	m = size(ien, 2);
	
	b = zeros(n, 1);
	
	for gpid = 1 : n
		n_elem = gpie(gpid, 1);
		for i_elem = 1 : n_elem
			elem_id = gpie(gpid, 1 + i_elem);
			elem_vertex_ids = ien(elem_id, :);
			for phi_id = 1 : 4
				if (elem_vertex_ids(phi_id) == gpid)
				end
			end
			
			vertex_coords = coords(elem_vertex_ids, :)';
			b(gpid) = b(gpid) + fem2d_quad_bilin_int_f(vertex_coords, phi_id);
		end
	end
end