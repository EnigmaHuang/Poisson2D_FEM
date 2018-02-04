function k = fem2d_quad_bilin_unit_stiff(x)
% Generate the unit stiffness matrix for a quadrilateral element 
% using bilinear basis functions
% [IN]  x : 2 * 4 matrix, the geometric coordinates of the element,
%           first row is x coordinates, second row is y coordinates,
%           points should be in counter clockwise order
% [OUT] k : 4 * 4 unit stiffness matrix, k(i, j) = \int_{\Omega^e}
%           phi_{i}^{'}(x, y) phi_{j}^{'}(d, y) dx dy
%           phi_{1, 2, 3, 4} are in counter clockwise order
	
	k = zeros(4, 4);
	
	% Gauss quadrature points and weight function
	q = [-1 / sqrt(3), 1 / sqrt(3)];
	w = [1, 1];
	n_quadrature = size(q, 2);
	
	% In-element 2D integral for $\phi_i^{'} * \phi_j^{'}$
	% Numerical integral using Gauss quadrature
	for ix = 1 : n_quadrature
		for iy = 1 : n_quadrature
			[sh, dtm] = fem2d_quad_bilin_shape(q(ix), q(iy), x, 0);
			d_N = sh(1 : 2, :);
			k = k + d_N' * d_N * dtm * w(ix) * w(iy);
		end
	end
	
	% Boundary line integral for alpha(x, y)
	k2 = zeros(4, 4);
	if (( (x(1, 1) == 0) && (x(1, 2) == 0) ) ...  % on x == 0
	 || ( (x(1, 1) == 1) && (x(1, 2) == 1) ) ...  % on x == 1
	 || ( (x(2, 1) == 0) && (x(2, 2) == 0) ) ...  % on y == 0
	 || ( (x(2, 1) == 1) && (x(2, 2) == 1) ))     % on y == 1
		k2 = k2 + fem2d_quad_bilin_int_alpha(x, 1, 2);
	end
	
	if (( (x(1, 2) == 0) && (x(1, 3) == 0) ) ...  % on x == 0
	 || ( (x(1, 2) == 1) && (x(1, 3) == 1) ) ...  % on x == 1
	 || ( (x(2, 2) == 0) && (x(2, 3) == 0) ) ...  % on y == 0
	 || ( (x(2, 2) == 1) && (x(2, 3) == 1) ))     % on y == 1
		k2 = k2 + fem2d_quad_bilin_int_alpha(x, 2, 3);
	end
	
	if (( (x(1, 3) == 0) && (x(1, 4) == 0) ) ...  % on x == 0
	 || ( (x(1, 3) == 1) && (x(1, 4) == 1) ) ...  % on x == 1
	 || ( (x(2, 3) == 0) && (x(2, 4) == 0) ) ...  % on y == 0
	 || ( (x(2, 3) == 1) && (x(2, 4) == 1) ))     % on y == 1
		k2 = k2 + fem2d_quad_bilin_int_alpha(x, 3, 4);
	end
	
	if (( (x(1, 4) == 0) && (x(1, 1) == 0) ) ...  % on x == 0
	 || ( (x(1, 4) == 1) && (x(1, 1) == 1) ) ...  % on x == 1
	 || ( (x(2, 4) == 0) && (x(2, 1) == 0) ) ...  % on y == 0
	 || ( (x(2, 4) == 1) && (x(2, 1) == 1) ))     % on y == 1
		k2 = k2 + fem2d_quad_bilin_int_alpha(x, 4, 1);
	end
	
	k = k + k2;
end