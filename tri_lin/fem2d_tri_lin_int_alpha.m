function k2 = fem2d_tri_lin_int_alpha(x, i, j)
% Perform line integral $\int_{\partial \Omega} \alpha * \phi_i * \phi_j d s$ on 
% the boundary edge x(i)-->x(j) of a given triangle element
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
	
	k2 = zeros(3, 3);
	fii = 0; fij = 0; fjj = 0;
	% Use Gauss quadrature to perform 1D integral
	for iq = 1 : n_quadrature
		% Convert from [-1, 1] to [0, l] to evaluate linear basis function
		t = semi_bma * qx(iq) + semi_bpa;
		t_over_l = t / l;
		
		% Convert to coordinate on line node(i)-->node(j) to evaluate alpha
		alpha_x = x(1, i) + t_over_l * dx;
		alpha_y = x(2, i) + t_over_l * dy;
		alpha   = poisson2d_robin_bc_alpha(alpha_x, alpha_y);
		
		fii = fii + ((1 - t_over_l) * (1 - t_over_l) * alpha) * qw(iq);
		fij = fij + ((1 - t_over_l) * t_over_l       * alpha) * qw(iq);
		fjj = fjj + (t_over_l       * t_over_l       * alpha) * qw(iq);
	end
	
	k2(i, i) = semi_bma * fii;
	k2(i, j) = semi_bma * fij;
	k2(j, i) = semi_bma * fij;
	k2(j, j) = semi_bma * fjj;
end