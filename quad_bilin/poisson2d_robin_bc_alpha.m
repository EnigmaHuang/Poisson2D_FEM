function res = poisson_robin_bc_alpha(x, y)
% Function $\alpha(x, y)$ in Robin boundary condition 
% $\frac{\partial u}{\partial \vec{n}} + \alpha * u = g$
% If you want to use Dirichlet boundary condition, set res = 0
	res = 1;
	%res = 0;
end