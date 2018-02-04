function res = fem2d_quad_bilin_int_f(geo_coord, phi_id)
% Use variable substitution to transform the f(x, y) and phi_i(x, y) on a
% element to [-1, 1] * [-1, 1], then use Gauss quadrature to integrate.
% [IN]  geo_coord : 2 * 4 matrix, the geometric coordinates of the element's vertexes
%                   1st row is x coordinates, 2nd row is y coordinates, 
%                   vertexes shoule be in counter clockwise order
% [IN]  phi_id    : Which bilinear basis function are to be integrated, 1 to 4
% [OUT] res       : Integral result
	
	w = [ 0.6521451548625461, 0.6521451548625461,  0.3478548451374538, 0.3478548451374538];
	q = [-0.3399810435848563, 0.3399810435848563, -0.8611363115940526, 0.8611363115940526];
	
	res = 0;
	for ix = 1 : 4
		for iy = 1 : 4
			[f_x, f_y] = fem2d_quad_xi_eta_to_x_y(q(ix), q(iy), geo_coord);
			f_xy = poisson2d_rhs_f(f_x, f_y);
			
			dtm = fem2d_quad_bilin_det(q(ix), q(iy), geo_coord);
			
			N(1) = 0.25 * (1.0 - q(ix)) * (1.0 - q(iy));
			N(2) = 0.25 * (1.0 + q(ix)) * (1.0 - q(iy));
			N(3) = 0.25 * (1.0 + q(ix)) * (1.0 + q(iy));
			N(4) = 0.25 * (1.0 - q(ix)) * (1.0 + q(iy));
			
			res = res + dtm * f_xy * N(phi_id) * w(ix) * w(iy);
		end
	end
end