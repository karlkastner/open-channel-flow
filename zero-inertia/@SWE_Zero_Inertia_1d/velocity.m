% 2025-09-08 10:59:34.396667299 +0200
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
function [ul,a,d] = velocity(obj,h)
	C = obj.C;
	h = rvec(h);
	h = real(h);
	%h = max(h,0);
	zb = rvec(obj.zb);
	dx = obj.L./(obj.nx-1);
	% h = max(h,0);
	% dh/dt = -d/x(u*h) 
	%       = -d/dx(C*h^(3/2)*sqrt(d/dx(h+zb)))
	%                   = -3/2*C*h^(1/2)*sqrt(S)*dh/dx -1/2*C*h^(3/2)/sqrt(S)*d2/dx2(h+zb)
	%                   = -3/2*u*dh/dx -1/2*h*u/S*d2/dx2(h+zb)
	% surface slope at interface
	Sl = (h+zb - left(h+zb,1))/dx;
	%Sl(Sl == 0) = sqrt(eps);
	% velocity at left interface
	ul = -C*sign(Sl).*sqrt(0.5*(left(h,1)+h).*abs(Sl));
end

