function res = bv_g(x, y)
% Dirichlet boundary condition

	if (abs(y - 0) < eps) 
		res = x;          % g(x, 0) = x
		%res = x;          % g(x, 0) = x
	end
	
	if (abs(y - 1) < eps) 
		res = 0;          % g(x, 1) = 0
		%res = 1 + x;      % g(x, 1) = 1 + x
	end
	
	if (abs(x - 0) < eps) 
		res = 0;          % g(0, y) = 0
		%res = y;          % g(0, y) = y
	end
	
	if (abs(x - 1) < eps) 
		res = 1 - y;      % g(1, y) = 1 - y
		%res = 1 + y;      % g(1, y) = 1 + y
	end
end