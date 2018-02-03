function [x, y] = xi_eta_to_xy(xi, eta, geo_coord)
% Transform local coordinate (xi, eta) to the geometric coordinate in a given element
% x(xi, eta) = \sum x_{1,2,3,4}e * N{1,2,3,4}(xi, eta), same for y(xi, eta)
% geo_coord are the geometric coordinates of the element's vertex, first row is 
% x coordinates, second row is y coordinates, should be in the counter clockwise order
	x =     geo_coord(1, 1) * 0.25 * (1.0 - xi) * (1.0 - eta);
	x = x + geo_coord(1, 2) * 0.25 * (1.0 + xi) * (1.0 - eta);
	x = x + geo_coord(1, 3) * 0.25 * (1.0 + xi) * (1.0 + eta);
	x = x + geo_coord(1, 4) * 0.25 * (1.0 - xi) * (1.0 + eta);
	y =     geo_coord(2, 1) * 0.25 * (1.0 - xi) * (1.0 - eta);
	y = y + geo_coord(2, 2) * 0.25 * (1.0 + xi) * (1.0 - eta);
	y = y + geo_coord(2, 3) * 0.25 * (1.0 + xi) * (1.0 + eta);
	y = y + geo_coord(2, 4) * 0.25 * (1.0 - xi) * (1.0 + eta);
end