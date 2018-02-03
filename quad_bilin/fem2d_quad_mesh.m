function [coords, ien, gpie, bgp] = fem2d_rect_mesh(x, y)
% [IN]  x, y   : Grid point on x and y direction
% [OUT] coords : N * 2 matrix, N is the number of grid points,
%                each row has the x and y coordinate of a grid point
% [OUT] ien    : M * 4 matrix, M is the number of quadrilateral elements,
%                each row has 4 vertex point ids of the element, in 
%                counter clockwise order
% [OUT] gpie   : N * 5 matrix, each row is the information of a 
%                grid point: gpie(i, 1) is the number of elements this grid
%                point belongs to, gpie(i, 2:1+gpie(i,1)) are the element ids 
% [OUT] bgp    : 2 * (nx + ny - 2) vector, the ids of the grid points on the boundary
	
	nx = max(size(x));
	ny = max(size(y));
	
	N = nx * ny;
	M = (nx - 1) * (ny - 1);
	
	coords = zeros(N, 2);
	ien    = zeros(M, 4);
	gpie   = zeros(N, 5);
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
				ielem = (iy - 1) * (nx - 1) + ix;
				
				% The ids of 4 vertexes in the element, counter clockwise order
				icoord1 = icoord;
				icoord2 = icoord + 1;
				icoord3 = icoord + 1 + nx;
				icoord4 = icoord + nx;
				
				% Record the vertex ids of the element
				ien(ielem, 1) = icoord1;
				ien(ielem, 2) = icoord2;
				ien(ielem, 3) = icoord3;
				ien(ielem, 4) = icoord4;
				
				% Mark that this element contains these vertexes
				gpie(icoord1, 1) = gpie(icoord1, 1) + 1;
				gpie(icoord1, 1 + gpie(icoord1, 1)) = ielem;
				
				gpie(icoord2, 1) = gpie(icoord2, 1) + 1;
				gpie(icoord2, 1 + gpie(icoord2, 1)) = ielem;
				
				gpie(icoord3, 1) = gpie(icoord3, 1) + 1;
				gpie(icoord3, 1 + gpie(icoord3, 1)) = ielem;
				
				gpie(icoord4, 1) = gpie(icoord4, 1) + 1;
				gpie(icoord4, 1 + gpie(icoord4, 1)) = ielem;
			end
			
			% Record boundary grid point
			if ((iy == 1) || (iy == ny) || (ix == 1) || (ix == nx))
				ibgp = ibgp + 1;
				bgp(ibgp) = icoord;
			end
		end
	end
end