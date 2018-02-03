function b = fem2d_tri_lin_rhs(coords, ien, gpie)
% Generate the right hand side for $-\nabla^2 u = f$ in 2D case
% using linear basis functions and quadrilateral elements
% [IN]  coords : n * 2 matrix, n is the number of grid points,
%                each row has the x and y coordinate of a grid point
% [IN]  ien    : m * 3 matrix, m is the number of quadrilateral elements,
%                each row has 3 vertex point ids of the element, in 
%                counter clockwise order
% [IN]  gpie   : Grid Poinr In Element, each row is the information of a 
%                grid point: gpie(i, 1) is the number of elements this grid
%                point belongs to, gpie(i, 2:end) are the element ids 
% [OUT] n      : n * 1 RHS vector
	
	n = size(coords, 1);
	m = size(ien, 1);
	
	b = zeros(n, 1);
	
	% $\int_{\Omega} f v ds$
	for gpid = 1 : n
		n_elem = gpie(gpid, 1);
		for i_elem = 1 : n_elem
			elem_id = gpie(gpid, 1 + i_elem);
			elem_vertex_ids = ien(elem_id, :);
			for phi_id = 1 : 3
				if (elem_vertex_ids(phi_id) == gpid)
				end
			end
			
			vertex_coords = coords(elem_vertex_ids, :)';
			b(gpid) = b(gpid) + fem2d_tri_lin_int_f(vertex_coords, phi_id);
		end
	end
	
	
	b2 = zeros(n, 1);
	% $\int_{\Omega^e} g v ds$
	for i_elem = 1 : m
		elem_vertex_ids = ien(i_elem, :);
		vertex_coords = coords(elem_vertex_ids, :)';
		
		if (( (vertex_coords(1, 1) == 0) && (vertex_coords(1, 2) == 0) ) ...  % on x == 0
		 || ( (vertex_coords(1, 1) == 1) && (vertex_coords(1, 2) == 1) ) ...  % on x == 1
		 || ( (vertex_coords(2, 1) == 0) && (vertex_coords(2, 2) == 0) ) ...  % on y == 0
		 || ( (vertex_coords(2, 1) == 1) && (vertex_coords(2, 2) == 1) ))     % on y == 1
			b2(elem_vertex_ids) = b2(elem_vertex_ids) + fem2d_tri_lin_int_g(vertex_coords, 1, 2);
		end
		
		if (( (vertex_coords(1, 2) == 0) && (vertex_coords(1, 3) == 0) ) ...  % on x == 0
		 || ( (vertex_coords(1, 2) == 1) && (vertex_coords(1, 3) == 1) ) ...  % on x == 1
		 || ( (vertex_coords(2, 2) == 0) && (vertex_coords(2, 3) == 0) ) ...  % on y == 0
		 || ( (vertex_coords(2, 2) == 1) && (vertex_coords(2, 3) == 1) ))     % on y == 1 
			b2(elem_vertex_ids) = b2(elem_vertex_ids) + fem2d_tri_lin_int_g(vertex_coords, 2, 3);
		end
		
		if (( (vertex_coords(1, 3) == 0) && (vertex_coords(1, 1) == 0) ) ...  % on x == 0
		 || ( (vertex_coords(1, 3) == 1) && (vertex_coords(1, 1) == 1) ) ...  % on x == 1
		 || ( (vertex_coords(2, 3) == 0) && (vertex_coords(2, 1) == 0) ) ...  % on y == 0
		 || ( (vertex_coords(2, 3) == 1) && (vertex_coords(2, 1) == 1) ))     % on y == 1
			b2(elem_vertex_ids) = b2(elem_vertex_ids) + fem2d_tri_lin_int_g(vertex_coords, 3, 1);
		end
	end
	
	b = b + b2;
end