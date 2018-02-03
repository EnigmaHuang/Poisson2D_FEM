function res = fem2d_int_lin_f(geo_coord, phi_id)
% Use variable substitution to transform the f(x, y) and phi_i(x, y) on a
% element to triangle (0,0)-(1,0)-(0,1), then use Gauss quadrature to integrate
% [IN]  geo_coord : 2 * 3 matrix, the geometric coordinates of the element's vertexes
%                   1st row is x coordinates, 2nd row is y coordinates, 
%                   vertexes shoule be in counter clockwise order
% [IN]  phi_id    : Which bilinear basis function are to be integrated, 1 to 3
% [OUT] res       : Integral result
	
	qx = [1.0/3.0, 0.6, 0.2, 0.2];
	qy = [1.0/3.0, 0.2, 0.2, 0.6];
	w  = [-27.0/96.0, 25.0/96.0, 25.0/96.0, 25.0/96.0];
	n_quadrature = size(w, 2);
	
	res = 0;
	for iq = 1 : n_quadrature
		[f_x, f_y] = fem2d_tri_xi_eta_to_x_y(qx(iq), qy(iq), geo_coord);
		f_xy = poisson2d_rhs_f(f_x, f_y);
		
		dtm = fem2d_tri_lin_det(qy(iq), qy(iq), geo_coord);
		
		N(1) = 1.0 - qx(iq) - qy(iq);
		N(2) = qx(iq);
		N(3) = qy(iq);
		
		res = res + dtm * f_xy * N(phi_id) * w(iq);
	end
end