function bb = fem2d_tri_lin_int_g(x, i, j)
% Perform line integral $\int_{\partial \Omega} g * \phi_j d s$ on 
% the boundary edge x(i)-->x(j) of a given triangular element
% [IN]  x    : 2 * 3 matrix, the geometric coordinates of the element's nodes,
%              first row is x coordinates, second row is y coordinates,
%              points should be in counter clockwise order
% [IN]  i, j : The index of the nodes on the boundary edge
% [OUT] k    : 3 * 3 matrix, contains the result of line integral
	
	qx = [-0.3399810435848563, 0.3399810435848563, -0.8611363115940526, 0.8611363115940526];
	qw = [ 0.6521451548625461, 0.6521451548625461,  0.3478548451374538, 0.3478548451374538];
	n_quadrature = size(qx, 2);
	
	dx = x(1, j) - x(1, i);
	dy = x(2, j) - x(2, i);
	l  = sqrt(dx * dx + dy * dy);
	semi_bma = l / 2;
	semi_bpa = l / 2;
	
	bb = zeros(3, 1);
	bi = 0; bj = 0;
	for iq = 1 : n_quadrature
		% Convert from [-1, 1] to [0, l] to evaluate linear basis function
		t = semi_bma * qx(iq) + semi_bpa;
		t_over_l = t / l;
		
		% Convert to coordinate on line node(i)-->node(j) to evaluate g(x,y)
		g_x = x(1, i) + t_over_l * dx;
		g_y = x(2, i) + t_over_l * dy;
		gxy = poisson2d_robin_bc_g(g_x, g_y);
		
		bi = bi + ((1 - t_over_l) * gxy) * qw(iq);
		bj = bj + (t_over_l       * gxy) * qw(iq);
	end
	
	bb(i) = semi_bma * bi;
	bb(j) = semi_bpa * bj;
end