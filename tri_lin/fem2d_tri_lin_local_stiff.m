function k = fem2d_tri_lin_global_stiff(x)
% Generate the local stiffness matrix for a triangle element 
% using linear basis functions
% [IN]  x : 2 * 3 matrix, the geometric coordinates of the element,
%           first row is x coordinates, second row is y coordinates,
%           points should be in counter clockwise order
% [OUT] k : 3 * 3 local stiffness matrix, k(i, j) = \int_{\Omega^e}
%           phi_{i}^{'}(x, y) phi_{j}^{'}(d, y) dx dy
%           phi_{1, 2, 3} are in counter clockwise order
	
	k = zeros(3, 3);
	
	% Gauss quadrature points and weight function
	qx = [1.0/3.0, 0.6, 0.2, 0.2];
	qy = [1.0/3.0, 0.2, 0.2, 0.6];
	w  = [-27.0/96.0, 25.0/96.0, 25.0/96.0, 25.0/96.0];
	n_quadrature = size(w, 2);
	
	% Numerical integral using Gauss quadrature
	for iq = 1 : n_quadrature
		[sh, dtm] = fem2d_tri_lin_shape(qx(iq), qy(iq), x);
		d_N = sh(1 : 2, :);
		k = k + d_N' * d_N * dtm * w(iq);
	end
end