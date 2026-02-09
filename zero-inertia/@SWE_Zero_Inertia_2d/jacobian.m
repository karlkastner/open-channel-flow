% Fri  4 Oct 11:41:53 CEST 2024
% Tue  8 Oct 21:18:13 CEST 2024
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
%
%
% J = d/dh(dh/dt)
function [J] = jacobian(obj,t,h,varargin)
	if (obj.opt.analytic_jacobian)
	if (obj.nx(1)<=3 || obj.nx(2)<=3)
		error('not correctly computed for overlapping boundaries');
	end
		% note: this approximates the jacobian by a 5-point stencil
		% and ignores the influence of the direction term S/|S| in the jacobian
		h = obj.input_factor*h;
		obj.interface_values(t,h,false,true,false);

		dx = obj.dx;
		J = -[   -flat(obj.aux.dqxi_dhl(1:end-1,:))'/dx(1); ...
                          flat(obj.aux.dqxi_dhc(2:end,:))'/dx(1);
                         -flat(obj.aux.dqyj_dhl(:,1:end-1))'/dx(2); ...
                          flat(obj.aux.dqyj_dhc(:,2:end))'/dx(2);
		          flat(  ( -obj.aux.dqxi_dhc(1:end-1,:) + obj.aux.dqxi_dhl(2:end,:)) / dx(1) ...
			       + ( -obj.aux.dqyj_dhc(:,1:end-1) + obj.aux.dqyj_dhl(:,2:end)) / dx(2) ...
			      )';
			];
		if (~isempty(obj.infiltration_rate_linear))
			% TODO only for 5-point kernel
			J(5,:) = J(5,:) + rvec(obj.infiltration_rate_linear(:));
		end
		% TODO the corresponding columns have to be zeroed as well
		J(:,obj.aux.ri_limited) = 0;
		% estimate the max eigenvalue
		% the condition number will then be: c ~ 1 + dt*eJmax
		% by gershgorins circle theorem
		% emax = max(abs(aii) + sum(|aij|,j!=i))
		%      = 2*max(abs(aii)), as rowsum is zero
		%      = 2*abs*min(aii)), as aii<0
		obj.aux.max_eig_J = 2*abs(min(J(5,:)));

		nx = obj.nx;
		J = reshape(J',prod(obj.nx),5);
		if (obj.opt.returnmat)
			% TODO store precomputed indices
			J = diags2mat2d(t,J,obj.nx,obj.boundary_condition{3}, ...
						obj.boundary_condition{4}, ...
						obj.boundary_condition{1}, ...
						obj.boundary_condition{2}, ...
						obj.dx);
		end % if returnmat
		J = obj.output_factor*J;
	else
		hmin = 0;
		%h = reshape(h,obj.nx);
		%J = jacobian_regular_grid_2d(@obj.dh_dt,0,h,obj.nx,obj.opt.epsJ,obj.opt.twosided,hmin,varargin{:});
		J = jacobian_numeric_unstructured(@(h) obj.dh_dt(t,h(:)),h,true);
	end
end

