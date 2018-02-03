function [sh, dtm] = fem2d_quad_bilin_shape(qx, qy, x, dtm_only)
% Generate 4 bilinear basis functions' information on the given quadrature point
% in a quadrilateral element
% [IN]  qx, qy   : The given quadrature point (to evaluate the function value)
% [IN]  x        : 2 * 4 matrix, the geometric coordinates of the element,
%                  first row is x coordinates, second row is y coordinates,
%                  points should be in counter clockwise order
% [IN]  dtm_only : Only compute Jacobian
% [OUT] sh       : 3 * 4 matrix, the bilinear basis functions' information on the quadrature
%                  point, first row is dN/dx, second row is dN/dy, third row is N
% [OUT] dtm      : The determinant of the Jacobian
	
	% Evaluate the partial derivatives of basis functions w.r.t. xi and eta at quadrature point
	d_N1_d_xi  = -0.25 * (1.0 - qy);
	d_N1_d_eta = -0.25 * (1.0 - qx);
	d_N2_d_xi  =  0.25 * (1.0 - qy);
	d_N2_d_eta = -0.25 * (1.0 + qx);
	d_N3_d_xi  =  0.25 * (1.0 + qy);
	d_N3_d_eta =  0.25 * (1.0 + qx);
	d_N4_d_xi  = -0.25 * (1.0 + qy);
	d_N4_d_eta =  0.25 * (1.0 - qx);
	
	% Compute the partial derivatives between x, y and xi, eta
	d_x_d_xi  = d_N1_d_xi  * x(1, 1) + d_N2_d_xi  * x(1, 2) + d_N3_d_xi  * x(1, 3) + d_N4_d_xi  * x(1, 4);
	d_x_d_eta = d_N1_d_eta * x(1, 1) + d_N2_d_eta * x(1, 2) + d_N3_d_eta * x(1, 3) + d_N4_d_eta * x(1, 4);
	d_y_d_xi  = d_N1_d_xi  * x(2, 1) + d_N2_d_xi  * x(2, 2) + d_N3_d_xi  * x(2, 3) + d_N4_d_xi  * x(2, 4);
	d_y_d_eta = d_N1_d_eta * x(2, 1) + d_N2_d_eta * x(2, 2) + d_N3_d_eta * x(2, 3) + d_N4_d_eta * x(2, 4);
	
	% Compute the determinant of Jacobian
	dtm = d_x_d_xi * d_y_d_eta - d_y_d_xi * d_x_d_eta;
	
	if (dtm_only == 1)
		sh = 0;
		return;
	end
	
	% Evaluate the partial derivatives of basis functions w.r.t. x and y at quadrature point
	sh(1, 1) = (d_N1_d_xi *  d_y_d_eta + d_N1_d_eta * -d_y_d_xi) / dtm;
	sh(2, 1) = (d_N1_d_xi * -d_x_d_eta + d_N1_d_eta *  d_x_d_xi) / dtm;
	sh(1, 2) = (d_N2_d_xi *  d_y_d_eta + d_N2_d_eta * -d_y_d_xi) / dtm;
	sh(2, 2) = (d_N2_d_xi * -d_x_d_eta + d_N2_d_eta *  d_x_d_xi) / dtm; 
	sh(1, 3) = (d_N3_d_xi *  d_y_d_eta + d_N3_d_eta * -d_y_d_xi) / dtm;
	sh(2, 3) = (d_N3_d_xi * -d_x_d_eta + d_N3_d_eta *  d_x_d_xi) / dtm;
	sh(1, 4) = (d_N4_d_xi *  d_y_d_eta + d_N4_d_eta * -d_y_d_xi) / dtm;
	sh(2, 4) = (d_N4_d_xi * -d_x_d_eta + d_N4_d_eta *  d_x_d_xi) / dtm;
	
	% Evaluate the basis function at quadrature point
	sh(3, 1) = 0.25 * (1.0 - qx) * (1.0 - qy);
	sh(3, 2) = 0.25 * (1.0 + qx) * (1.0 - qy);
	sh(3, 3) = 0.25 * (1.0 + qx) * (1.0 + qy);
	sh(3, 4) = 0.25 * (1.0 - qx) * (1.0 + qy);
end