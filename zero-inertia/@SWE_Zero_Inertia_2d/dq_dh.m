% 2025-10-28 18:03:46.680763713 +0100
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
function [dqxi_dhl,dqxi_dhc,dqyj_dhl,dqyj_dhc] = dqxi_dh(obj,t,h0)
	h0 = reshape(h0,obj.nx);
	[qi0,qj0] = obj.interface_values(t,h0);
	epsh = sqrt(eps);
%	epsh = 1e-4;
%	epsh = 0.1;
	dqxi_dhl = zeros(obj.nx(1)+1,obj.nx(2));
	dqxi_dhc = zeros(obj.nx(1)+1,obj.nx(2));
	dqyj_dhl = zeros(obj.nx(1),obj.nx(2)+1);
	dqyj_dhc = zeros(obj.nx(1),obj.nx(2)+1);
	for idx=1:3
		id = (idx:3:obj.nx(2))'*ones(1,obj.nx(2)/3);
		for jdx=1:3
			jd = ones(obj.nx(1)/3,1)*(jdx:3:obj.nx(2));
			h = h0;
			ih = sub2ind(size(h),id,jd);
			% perturb
%			h(ih) = h(ih) + epsh
			h(id,jd) = h(id,jd) + epsh;
			[qi,qj] = obj.interface_values(t,h);
%qi
			dqi_dh = (qi - qi0)/epsh;
			dqj_dh = (qj - qj0)/epsh;
			dqxi_dhc(id,jd)   = dqi_dh(id,jd);
			dqxi_dhl(id+1,jd) = dqi_dh(id+1,jd);
			dqyj_dhc(id,jd) = dqj_dh(id,jd);
			dqyj_dhl(id,jd+1) = dqj_dh(id,jd+1);
%pause
		end
	end
end

