% 2025-12-08 15:47:26.669228873 +0100
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
function flag = output_event(obj,t,z,dt)
	flag = false;
%	if (obj.opt.output.on_output_event)
%	rms_zold = rms(obj.aux.zold);
%	rms_z = rms(z);
%	% wet to dry
%	if (rms_zold >= obj.opt.heps & rms_z < obj.opt.heps)
%		flag = true;
%	end
%	% dry to wet (first wet step)
%	% TODO push the old event	if (rms_zold < obj.opt.heps & rms_z >= obj.opt.heps)
%	if (rms_z < obj.opt.heps && (t+dt >= obj.aux.tevent(obj.aux.edx+1)))
%		flag = true;
%	end
%	end
end


