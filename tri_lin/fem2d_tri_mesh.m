function [coords, ien, gpie, bgp] = fem2d_tri_mesh(x, y)
% Generate triangle mesh on [0, 1] * [0, 1] for 2D FEM
% [IN]  x, y   : Grid point on x and y direction
% [OUT] coords : N * 2 matrix, N is the number of grid points,
%                each row has the x and y coordinate of a grid point
% [OUT] ien    : M * 3 matrix, M is the number of quadrilateral elements,
%                each row has 3 vertex point ids of the element, in 
%                counter clockwise order
% [OUT] gpie   : N * 7 matrix, each row is the information of a 
%                grid point: gpie(i, 1) is the number of elements this grid
%                point belongs to, gpie(i, 2:1+gpie(i,1)) are the element ids 
% [OUT] bgp    : 2 * (nx + ny - 2) vector, the ids of the grid points on the boundary
	
	nx = max(size(x));
	ny = max(size(y));
	
	N = nx * ny;
	M = 2 * (nx - 1) * (ny - 1);
	
	coords = zeros(N, 2);
	ien    = zeros(M, 3);
	gpie   = zeros(N, 7);
	bgp    = zeros(2 * (nx + ny - 2), 1);
	
	ibgp   = 0;
	for iy = 1 : ny
		for ix = 1 : nx
			% Record the coordinates of this grid point
			icoord = (iy - 1) * nx + ix;
			coords(icoord, :) = [x(ix) y(iy)];
			
			% Handle the element that on the top right of the current grid point
			if ((iy < ny) && (ix < nx))
				% The id of the element
				ielem  = (iy - 1) * (nx - 1) + ix;
				ielem1 = 2 * ielem - 1;
				ielem2 = 2 * ielem;
				
				% The ids of 4 vertexes in the element, counter clockwise order
				icoord1 = icoord;
				icoord2 = icoord + 1;
				icoord3 = icoord + 1 + nx;
				icoord4 = icoord + nx;
				
				% Record the vertex ids of the element
				ien(ielem1, 1) = icoord1;
				ien(ielem1, 2) = icoord2;
				ien(ielem1, 3) = icoord3;
				
				ien(ielem2, 1) = icoord1;
				ien(ielem2, 2) = icoord3;
				ien(ielem2, 3) = icoord4;
				
				
				% Mark that this element contains these vertexes
				gpie(icoord1, 1) = gpie(icoord1, 1) + 1;
				gpie(icoord1, 1 + gpie(icoord1, 1)) = ielem1;
				
				gpie(icoord2, 1) = gpie(icoord2, 1) + 1;
				gpie(icoord2, 1 + gpie(icoord2, 1)) = ielem1;
				
				gpie(icoord3, 1) = gpie(icoord3, 1) + 1;
				gpie(icoord3, 1 + gpie(icoord3, 1)) = ielem1;
				
				gpie(icoord1, 1) = gpie(icoord1, 1) + 1;
				gpie(icoord1, 1 + gpie(icoord1, 1)) = ielem2;
				
				gpie(icoord3, 1) = gpie(icoord3, 1) + 1;
				gpie(icoord3, 1 + gpie(icoord3, 1)) = ielem2;
				
				gpie(icoord4, 1) = gpie(icoord4, 1) + 1;
				gpie(icoord4, 1 + gpie(icoord4, 1)) = ielem2;
			end
			
			% Record boundary grid point
			if ((iy == 1) || (iy == ny) || (ix == 1) || (ix == nx))
				ibgp = ibgp + 1;
				bgp(ibgp) = icoord;
			end
		end
	end
end