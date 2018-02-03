function [x, y] = fem2d_tri_xi_to_x_y(xi, eta, geo_coord)
% Transform local coordinate (xi, eta) to the geometric coordinate in a given element
% x(xi, eta) = \sum x_{1,2,3}e * N{1,2,3}(xi, eta), same for y(xi, eta)
% geo_coord are the geometric coordinates of the element's vertex, first row is 
% x coordinates, second row is y coordinates, should be in the counter clockwise order

	x = geo_coord(1, 1) * (1.0 - xi - eta) + geo_coord(1, 2) * xi + geo_coord(1, 3) * eta;
	y = geo_coord(2, 1) * (1.0 - xi - eta) + geo_coord(2, 2) * xi + geo_coord(2, 3) * eta;
end