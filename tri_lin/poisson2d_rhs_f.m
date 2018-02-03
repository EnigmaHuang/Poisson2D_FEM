function res = poisson2d_rhs_f(x, y)
% Right hand side function f(x, y) in 2D Poisson equation
	res = 5 * x * y;
	%res = 6 * y - 3 * x;
end