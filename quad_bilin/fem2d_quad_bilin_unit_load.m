function ub = fem2d_quad_bilin_unit_load(x)
% Use variable substitution to transform the f(x, y) and phi_i(x, y) on a
% element to [-1, 1] * [-1, 1], then use Gauss quadrature to integrate.
% [IN]  x  : 2 * 4 matrix, the geometric coordinates of the element's vertexes
%            1st row is x coordinates, 2nd row is y coordinates, 
%            vertexes shoule be in counter clockwise order
% [OUT] ub : 4 * 1 vector, integral result
	
	w = [ 0.6521451548625461, 0.6521451548625461,  0.3478548451374538, 0.3478548451374538];
	q = [-0.3399810435848563, 0.3399810435848563, -0.8611363115940526, 0.8611363115940526];
	
	ub = zeros(4, 1);
	
	% $\int_{\Omega^e} f * \phi ds$
	for ix = 1 : 4
		for iy = 1 : 4
			[f_x, f_y] = fem2d_quad_xi_eta_to_x_y(q(ix), q(iy), x);
			f_xy = poisson2d_rhs_f(f_x, f_y);
			
			dtm = fem2d_quad_bilin_det(q(ix), q(iy), x);
			
			N(1) = 0.25 * (1.0 - q(ix)) * (1.0 - q(iy));
			N(2) = 0.25 * (1.0 + q(ix)) * (1.0 - q(iy));
			N(3) = 0.25 * (1.0 + q(ix)) * (1.0 + q(iy));
			N(4) = 0.25 * (1.0 - q(ix)) * (1.0 + q(iy));
			
			ub(1) = ub(1) + dtm * f_xy * N(1) * w(ix) * w(iy);
			ub(2) = ub(2) + dtm * f_xy * N(2) * w(ix) * w(iy);
			ub(3) = ub(3) + dtm * f_xy * N(3) * w(ix) * w(iy);
			ub(4) = ub(4) + dtm * f_xy * N(4) * w(ix) * w(iy);
		end
	end
	
	% $\int_{\partial \Omega^e} g * \phi dt $
	if (( (x(1, 1) == 0) && (x(1, 2) == 0) ) ...  % on x == 0
	 || ( (x(1, 1) == 1) && (x(1, 2) == 1) ) ...  % on x == 1
	 || ( (x(2, 1) == 0) && (x(2, 2) == 0) ) ...  % on y == 0
	 || ( (x(2, 1) == 1) && (x(2, 2) == 1) ))     % on y == 1
		ub = ub + fem2d_quad_bilin_int_g(x, 1, 2);
	end
	
	if (( (x(1, 2) == 0) && (x(1, 3) == 0) ) ...  % on x == 0
	 || ( (x(1, 2) == 1) && (x(1, 3) == 1) ) ...  % on x == 1
	 || ( (x(2, 2) == 0) && (x(2, 3) == 0) ) ...  % on y == 0
	 || ( (x(2, 2) == 1) && (x(2, 3) == 1) ))     % on y == 1 
		ub = ub + fem2d_quad_bilin_int_g(x, 2, 3);
	end
	
	if (( (x(1, 3) == 0) && (x(1, 4) == 0) ) ...  % on x == 0
	 || ( (x(1, 3) == 1) && (x(1, 4) == 1) ) ...  % on x == 1
	 || ( (x(2, 3) == 0) && (x(2, 4) == 0) ) ...  % on y == 0
	 || ( (x(2, 3) == 1) && (x(2, 4) == 1) ))     % on y == 1
		ub = ub + fem2d_quad_bilin_int_g(x, 3, 4);
	end
	
	if (( (x(1, 4) == 0) && (x(1, 1) == 0) ) ...  % on x == 0
	 || ( (x(1, 4) == 1) && (x(1, 1) == 1) ) ...  % on x == 1
	 || ( (x(2, 4) == 0) && (x(2, 1) == 0) ) ...  % on y == 0
	 || ( (x(2, 4) == 1) && (x(2, 1) == 1) ))     % on y == 1
		ub = ub + fem2d_quad_bilin_int_g(x, 4, 1);
	end
end