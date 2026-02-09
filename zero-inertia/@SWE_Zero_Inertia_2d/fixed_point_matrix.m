% Mon 15 Sep 14:37:36 CEST 2025
% Karl Kastner, Berlin
%
% This program is free software: you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation, either version 3 of the License, or
% (at your option) any later version.
% 
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
% 
% You should have received a copy of the GNU General Public License
% along with this program.  If not, see <https://www.gnu.org/licenses/>.
% note: the function returns dh/dt, is is named dz_dt to be comptible
% h(t+dt) = -int dq/dx dt 
% h(t+dt) = -int duh/dx dt 
% h(t+dt)^(m+1) = -int d u^m h^(m+1)/dx dt 
function A = fixed_point_matrix(obj,t,h,arg)
	obj.interface_values(t,h,false,false,true);

	%buf3 = [    Di(1:end-1);
	%	    (-Di(1:end-1) - Di(2:end));
	%	        Di(2:end);
	%           ]/(dx*dx);

	A =       [ flat(obj.aux.Di(1:end-1,:)), ...
		    flat(obj.aux.Di(2:end,:)), ...
		    flat(obj.aux.Dj(:,1:end-1)), ...
		    flat(obj.aux.Dj(:,2:end)), ...
		    flat( (   (-obj.aux.Di(1:end-1,:) - obj.aux.Di(2:end,:)) ...
			    + (-obj.aux.Dj(:,1:end-1) - obj.aux.Dj(:,2:end)) ...
			  ) ...
			);
	           ];

	if (obj.opt.returnmat)
		% TODO store precomputed indices
		A = diags2mat2d(t,A,obj.nx,obj.boundary_condition{3}, ...
					obj.boundary_condition{4}, ...
					obj.boundary_condition{1}, ...
					obj.boundary_condition{2}, ...
					obj.dx);
		%J  = diags2mat2d(obj.aux.Jd,obj.n,true);
	end % if returnmat
	if (0)
	buf1 = repmat((1:nx)',3,1);
	buf2 = [(0:nx-1)';
                (1:nx)';
                (2:nx+1)'];


	% apply boundary conditions
	buf2(1)   = nx;
	buf2(end) = 1;
	b31 = buf3(1);
	b3n = buf3(end);

	bcl        = obj.boundary_condition{1};
	bcr        = obj.boundary_condition{2};
	buf3(1)    = b31*((2*bcl(3)*dx)./(bcl(1)*dx - 2*bcl(2) + 2*bcl(3)*dx));
	buf3(nx+1) = buf3(nx+1) + b31.*((-2*bcl(2) - bcl(1)*dx))./(bcl(1)*dx - 2*bcl(2) + 2*bcl(3)*dx);
	buf3(3*nx) = b3n*((2*bcr(3)*dx)./(bcr(1)*dx + 2*bcr(2) + 2*bcr(3)*dx));
	buf3(2*nx) = buf3(2*nx) + b3n.*(( 2*bcr(2) - bcr(1)*dx))./(bcr(1)*dx + 2*bcr(2) + 2*bcr(3)*dx);

	A = sparse(buf1,buf2,buf3,nx,nx);
	end
end

