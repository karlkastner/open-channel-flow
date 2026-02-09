% Mon 25 Aug 14:50:03 CEST 2025
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
% function c = flow_combined_resistance_diffusion(h,S,c1,c2)
function [Ds,Dp] = flow_combined_resistance_diffusion(h,S,c1,c2)
	if (~issym(h))
		g = Constant.gravity;
	else
		syms g;
	end
	den = sqrt(c1.*c1 + 4*g*abs(S).*c2.*h.*h.*h);
	Ds = (g*h.*h.*h)./den;
	Dp = (den-c1)./(2*c2.*abs(S));
end

