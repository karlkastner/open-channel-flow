% Thu  3 Oct 17:14:25 CEST 2024
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
function [dh_dt, J] = linear_diffusion_2d_dhdt(h,z,a,d,L,n,mode,returnmat)

	dx  = L./n;
	dx2 = dx.*dx;

	isvector_h = isvector(h);
	h = reshape(h,n);

	if (size(h,3)>1)
		z = h(:,:,2);
	else
		if (~isempty(z))
			z = reshape(z,n);
		else
			z = 0;
		end
	end

	% mesh peclet number
	Pel = -(0.5*dx(1))*a(1)./d(1);
	Peu = -(0.5*dx(2))*a(2)./d(2);

	% weights of difference terms at left and down interface
	switch (mode)
	case {'central'}
		bl = zeros(n);
		bu = zeros(n);
	case{'upwind'}
		bl = 2*(al<0)-1;
		bu = 2*(au<0)-1;
	case{'optimal'}
		% magical coefficient
		bl = coth(Pel) - 1./Pel;
		bu = coth(Ped) - 1./Ped;
	case {'optimal-linearized'}
		% linearized magical coefficient
		bl = max(-1,min(1,Pel/3));
		bu = max(-1,min(1,Peu/3));
	end
	br = right(bl,1);
	bd = down(bu,1);

	% TODO diffusion
	dh_dt = ( -(   0.5*      a(1).*( (1-br).*h          + (1+br).*right(h,1)) ...
	             - 0.5*      a(1).*( (1-bl).*left(h,1)  + (1+bl).*h) ...
		   )/dx(1) ...
	          -(   0.5*      a(2).*( (1-bd).*h          + (1+bd).*down(h,1)) ...
	             - 0.5*      a(2).*( (1-bu).*up(h,1)    + (1+bu).*   h   ) ...
		   )/dx(2) ...
		 );
	if (isvector_h)
		dh_dt = reshape(dh_dt,n(1)*n(2),1);
	end		
	if (nargout()>0)
		J = swe_simplified_2d_jacobian(h,z,K,L,n,mode,returnmat);
	end
end % swe_linearized_2d_dhdt

