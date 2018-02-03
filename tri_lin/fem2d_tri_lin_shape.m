function [sh, dtm] = fem2d_tri_lin_shape(qx, qy, x)
% Generate 3 linear basis functions' information on the given quadrature point
% in a triangle element
% [IN]  qx, qy   : The given quadrature point (to evaluate the function value)
% [IN]  x        : 2 * 3 matrix, the geometric coordinates of the element,
%                  first row is x coordinates, second row is y coordinates,
%                  points should be in counter clockwise order
% [IN]  dtm_only : Only compute Jacobian
% [OUT] sh       : 3 * 3 matrix, the bilinear basis functions' information on the quadrature
%                  point, first row is dN/dx, second row is dN/dy, third row is N
% [OUT] dtm      : The determinant of the Jacobian
	
	% Evaluate the partial derivatives of basis functions w.r.t. xi and eta at quadrature point
	d_N1_d_xi  = -1.0;
	d_N1_d_eta = -1.0;
	d_N2_d_xi  =  1.0;
	d_N2_d_eta =  0.0;
	d_N3_d_xi  =  0.0;
	d_N3_d_eta =  1.0;
	
	% Compute the partial derivatives between x, y and xi, eta
	d_x_d_xi  = d_N1_d_xi  * x(1, 1) + d_N2_d_xi  * x(1, 2) + d_N3_d_xi  * x(1, 3);
	d_x_d_eta = d_N1_d_eta * x(1, 1) + d_N2_d_eta * x(1, 2) + d_N3_d_eta * x(1, 3);
	d_y_d_xi  = d_N1_d_xi  * x(2, 1) + d_N2_d_xi  * x(2, 2) + d_N3_d_xi  * x(2, 3);
	d_y_d_eta = d_N1_d_eta * x(2, 1) + d_N2_d_eta * x(2, 2) + d_N3_d_eta * x(2, 3);
	
	% Compute the determinant of Jacobian
	dtm = d_x_d_xi * d_y_d_eta - d_y_d_xi * d_x_d_eta;
	
	% Evaluate the partial derivatives of basis functions w.r.t. x and y at quadrature point
	sh(1, 1) = (d_N1_d_xi *  d_y_d_eta + d_N1_d_eta * -d_y_d_xi) / dtm;
	sh(2, 1) = (d_N1_d_xi * -d_x_d_eta + d_N1_d_eta *  d_x_d_xi) / dtm;
	sh(1, 2) = (d_N2_d_xi *  d_y_d_eta + d_N2_d_eta * -d_y_d_xi) / dtm;
	sh(2, 2) = (d_N2_d_xi * -d_x_d_eta + d_N2_d_eta *  d_x_d_xi) / dtm; 
	sh(1, 3) = (d_N3_d_xi *  d_y_d_eta + d_N3_d_eta * -d_y_d_xi) / dtm;
	sh(2, 3) = (d_N3_d_xi * -d_x_d_eta + d_N3_d_eta *  d_x_d_xi) / dtm;
	
	% Evaluate the basis function at quadrature point
	sh(3, 1) = 1.0 - qx - qy;
	sh(3, 2) = qx;
	sh(3, 3) = qy;
end