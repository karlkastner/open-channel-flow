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
%
% h(t+dt) = -int dq/dx dt 
% h(t+dt) = -int duh/dx dt 
% h(t+dt)^(m+1) = -int d u^m h^(m+1)/dx dt 
function A = fixed_point_matrix(obj,t,h,arg)
	obj.interface_values(t,h,false,false,true);
	A = [    obj.aux.Di(1:end-1), ...
	         (-obj.aux.Di(1:end-1) - obj.aux.Di(2:end)), ...
		 obj.aux.Di(2:end);
	    ];

	%if (~isempty(obj.infiltration_rate_linear))
	%	buf(:,2) = J(2,:) + rvec(obj.infiltration_rate_linear);
	%end

	if (obj.opt.returnmat)
		nx = obj.nx;
		% TODO store precomputed indices
		A = diags2mat1d(t,A,obj.nx,obj.boundary_condition{1},...
					   obj.boundary_condition{2},obj.dx);
	end % if returnmat
end

