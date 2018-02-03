function dtm = fem2d_tri_lin_det(qx, qy, x)
% Generate the determinant of the Jacobian on the given quadrature point
% in a triangle element (a simplified version of fem2d_tri_lin_shape)
% [IN]  qx, qy   : The given quadrature point (to evaluate the function value)
% [IN]  x        : 2 * 3 matrix, the geometric coordinates of the element,
%                  first row is x coordinates, second row is y coordinates,
%                  points should be in counter clockwise order
% [OUT] dtm      : The determinant of the Jacobian
	
	% Compute the partial derivatives between x, y and xi, eta
	d_x_d_xi  = x(1, 2) - x(1, 1);
	d_x_d_eta = x(1, 3) - x(1, 1);
	d_y_d_xi  = x(2, 2) - x(2, 1);
	d_y_d_eta = x(2, 3) - x(2, 1); 
	
	% Compute the determinant of Jacobian
	dtm = d_x_d_xi * d_y_d_eta - d_y_d_xi * d_x_d_eta;
end