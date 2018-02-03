function res = rhs_f(x, y)
% Right hand side in $\nabla^2 u = f$
	res = (x + y);
	%res = 6 * y - 3 * x;
end