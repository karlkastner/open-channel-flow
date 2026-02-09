% 2025-09-09 17:43:40.530158718 +0200
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
function ri = infiltration_rate(obj,h,dh_dt)
	if (~isscalar(obj.infiltration_rate_linear))
		obj.infiltration_rate_linear = reshape(obj.infiltration_rate_linear,size(h));
	end
	if (~isscalar(obj.infiltration_rate_constant))
		obj.infiltration_rate_constant = reshape(obj.infiltration_rate_constant,size(h));
	end
	ri = obj.infiltration_rate_linear.*h + obj.infiltration_rate_constant;
	if (~obj.opt.step_split_source)
		h0    = obj.aux.zold;
		siz   = size(dh_dt);
		rimax = max(0,reshape(h0,siz)/obj.aux.dt + dh_dt);
		obj.aux.ri_limited = (ri>rimax);
		ri = min(ri,rimax);
	else
		obj.aux.ri_limited = [];
	end
end

