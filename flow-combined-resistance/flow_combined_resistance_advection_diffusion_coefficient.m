% 2025-05-13 12:00:39.611392059 +0200
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
function [a,d] = flow_combined_resistance_advection(h,S,c1,c2)
	if (~issym(h))
	g = 9.81;
	else
		syms g;
	end
	%u  = 0.5*sign(S).*(C*(C*lcd - sqrt(C*C*lcd*lcd + 4*abs(S).*h*g^2)))/g;
	den = sqrt(c1.*c2 + 4*g*c2.*abs(S).*h.*h.*h);

	a = - (3*g*S.*h.*h)./den;

	% diffusion of the disharge wave, diffusion coefficient at the interface
	% had sign S^2 in it
	d =  (g.*h.^3)./den;
end

