function u = test_fem2d(h)
% Test my 2D FEM solver, h is the grid point interval
	
	if (h > 0)
		x = 0 : h : 1; x(end) = 1;
		y = x;
	else
		x = [0 : 0.05 : 0.09, 0.1 : 0.1 : 0.8, 0.9 : 0.05 : 1.0];
		y = [0 : 0.05 : 0.19, 0.2 : 0.1 : 0.7, 0.8 : 0.05 : 1.0];
	end
	
	nx = size(x, 2);
	ny = size(y, 2);
	
	% Generate the square mesh
	[coords, ien, gpie, bgp] = fem2d_tri_mesh(x, y);
	
	% Generate the global stiffness matrix
	K = fem2d_tri_lin_global_stiff(coords, ien);
	
	% Generate the right hand side
	b = fem2d_tri_lin_rhs(coords, ien, gpie);
	
	% Apply the Dirichlet boundary condition
	% If we want to use Dirichlet boundary condition, set use_dirichlet_bc = 1,
	% and set poisson2d_robin_bc_alpha() = 0, poisson2d_robin_bc_g() = 0
	use_dirichlet_bc = 0;
	if (use_dirichlet_bc == 1)
		n_bgp = max(size(bgp));
		for i = 1 : n_bgp
			bgp_id = bgp(i);
			K(bgp_id, :) = 0;
			K(bgp_id, bgp_id) = 1;
			bv_x = coords(bgp_id, 1);
			bv_y = coords(bgp_id, 2);
			bv_val = poisson2d_dirichlet_bc_g(bv_x, bv_y);
			b(bgp_id) = bv_val;
		end
	end
	
	% Solve the linear system
	u = K \ b;
	
	% Plot the result
	xm = reshape(coords(:, 1), nx, ny);
	ym = reshape(coords(:, 2), nx, ny);
	zm = reshape(u, nx, ny);
	
	subplot(1, 2, 1);
	contourf(xm, ym, zm, 'ShowText','on');
	subplot(1, 2, 2);
	trisurf(ien, xm, ym, zm);
end