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


% J = d/dh(dh/dt)
function [J, bl] = jacobian_y(obj,t,h,p,varargin)
	if (obj.opt.analytic_jacobian)
		h = obj.input_factor*h;
		obj.interface_values_y(t,h,false,true);

		dx = obj.dx;	
		z=zeros(prod(obj.nx),1);
		J  = -[   z, ... %-0.*flat(obj.aux.dqi_dhl(1:end-1,:))'/dx(1); ...
                          z, ... % 0.*flat(obj.aux.dqi_dhc(2:end,:))'/dx(1);
                          -flat(obj.aux.dqyj_dhl(:,1:end-1))/dx(2),... ...
                           flat(obj.aux.dqyj_dhc(:,2:end))/dx(2),...
		           flat(  ... %( -obj.aux.dqi_dhc(1:end-1,:) + obj.aux.dqi_dhl(2:end,:)) / dx(1) ...
			       ( -obj.aux.dqyj_dhc(:,1:end-1) + obj.aux.dqyj_dhl(:,2:end)) / dx(2) ...
			       )
			];

		if (~isempty(obj.infiltration_rate_linear))
			% TODO only for 5-point kernel
			J(:,5) = J(:,5) + rvec(obj.infiltration_rate_linear);
		end

		nx = obj.nx;
		%J = reshape(J',prod(obj.nx),5);
		if (obj.opt.returnmat)
			% TODO store precomputed indices
			J = diags2mat2d(t,J,obj.nx,obj.boundary_condition{3}, ...
						obj.boundary_condition{4}, ...
						obj.boundary_condition{1}, ...
						obj.boundary_condition{2}, ...
						obj.dx);
			%J  = diags2mat2d(obj.aux.Jd,obj.n,true);
		end % if returnmat
		J = obj.output_factor*J;
	else
		h = reshape(h,obj.nx);
		J = jacobian_regular_grid_2d_j(@obj.dh_dt_j,0,h,obj.n,obj.eps,obj.twosided);
	end
end % jacobian_y

