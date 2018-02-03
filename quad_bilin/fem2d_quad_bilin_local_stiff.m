function k = fem2d_quad_bilin_local_stiff(x)
% Generate the local stiffness matrix for a quadrilateral element 
% using bilinear basis functions
% [IN]  x : 2 * 4 matrix, the geometric coordinates of the element,
%           first row is x coordinates, second row is y coordinates,
%           points should be in counter clockwise order
% [OUT] k : 4 * 4 local stiffness matrix, k(i, j) = \int_{\Omega^e}
%           phi_{i}^{'}(x, y) phi_{j}^{'}(d, y) dx dy
%           phi_{1, 2, 3, 4} are in counter clockwise order
	
	k = zeros(4, 4);
	
	% Gauss quadrature points and weight function
	q = [-1 / sqrt(3), 1 / sqrt(3)];
	w = [1, 1];
	n_quadrature = size(q, 2);
	
	% Numerical integral using Gauss quadrature
	for ix = 1 : n_quadrature
		for iy = 1 : n_quadrature
			[sh, dtm] = fem2d_quad_bilin_shape(q(ix), q(iy), x, 0);
			d_N = sh(1 : 2, :);
			k = k + d_N' * d_N * dtm * w(ix) * w(iy);
		end
	end
end