% Mon  7 Oct 21:12:44 CEST 2024
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
function [dh_dt,bl,bu,vl,vu,Pel,dl,Sl] = swe_zero_inertia_2d_dhdt(h,z,C,L,n,mode,returnmat)
	%h  = rvec(h);
	isvector_h = isvector(h);
	e = sqrt(eps);
	h = reshape(h,n);
	z  = reshape(z,n);
	dx = L./n;


	% g h dzs/dx + cd V u = 0
	% h = max(h,0);
	% dh/dt = -d/x(u*h) 
	%       = -d/dx(C*h^(3/2)*sqrt(d/dx(h+z)))
	%                   = -3/2*C*h^(1/2)*sqrt(S)*dh/dx -1/2*C*h^(3/2)/sqrt(S)*d2/dx2(h+z)
	%                   = -3/2*u*dh/dx -1/2*h*u/S*d2/dx2(h+z)
	% surface slope in x-direction at left interface
	Slx   = (h + z - left(h+z,1))/dx(1);
	% surface slope in y-direction at up-interface
	Suy = (h+z - up(h+z,1))/dx(2);

	% surface level above and below left interface midpoint
	zsl_u = 0.5*(left(  up(h,1),1) +   up(h,1));
	zsl_d = 0.5*(left(down(h,1),1) + down(h,1));

	% surface slope in y-direction at left interface
	Sly   = 0.5*(zsl_d - zsl_u)/dx(2);




	% surface level left and right of above interface midpoint
	zsu_l = 0.5*( left(up(h,1),1) +  left(h,1));
	zsu_r = 0.5*(right(up(h,1),1) + right(h,1));
	Sux = 0.5*(zsu_r - zsu_l)/dx(1);

	% slope magnitude at left interface
	Sl = sqrt(e + Slx.^2 + Sly.^2);
	% slope in x-direction at up-interface

	% slope magnitude at up interface
	Su = sqrt(e + Sux.^2 + Suy.^2);

	% velocity magnitude at left interface
	Vl = C*sqrt(0.5*(left(h,1)+h).*Sl);
	% velocity-magnitude at up-interface
	Vu = C*sqrt(0.5*(up(h,1)+h).*Su);

	% velocity in x-direction at left interface 
	vl = -Vl.*Slx./Sl;

	% velocity in y-direcion at up-interface
	vu = -Vu.*Suy./Su;

	% celerity in x-directon, advection coefficient at left interface
	al = 1.5*vl;
	% celerity in y-directon, advection coefficient at left interface
	au = 1.5*vu;

	% diffusion in x-direction, coefficient at left interface
	dl = -0.25*(h+left(h,1)).*Vl./Sl;
	du = -0.25*(h+  up(h,1)).*Vu./Su;

	% mesh peclet number
	s = -1;
	Pel = -s*(0.5*dx(1))*al./dl;
	Peu = -s*(0.5*dx(2))*au./du;

	% weights of difference terms at left interface
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
		bu = coth(Peu) - 1./Peu;
	case {'optimal-linearized'}
		% linearized magical coefficient
		bl = max(-1,min(1,Pel/3));
		bu = max(-1,min(1,Peu/3));
	end
	fdx = isnan(bl);
	bl(fdx) = 2*(al(fdx)<0)-1;
	br = right(bl,1);
	fdx = isnan(bu);
	bu(fdx) = 2*(au(fdx)<0)-1;
	bd = down(bu,1);
 
	dh_dt = -( ( ...
		       0.5*right(vl,1).*( (1-br).*h         + (1+br).*right(h,1)) ...
	             - 0.5*      vl   .*( (1-bl).*left(h,1) + (1+bl).*h) ...
		   )/dx(1) ...
		 + ( ...
		       0.5*down(vu,1).*( (1-bd).*h          + (1+bd).*down(h,1)) ...
	             - 0.5*     vu   .*( (1-bu).*up(h,1)    + (1+bu).*h) ...
		   )/dx(2) ...
		 );
	if (any(h(:)<0))
		'honk'
		dh_dt = 1/sqrt(eps)*ones(size(dh_dt));
	end

	if (isvector_h)
		dh_dt = reshape(dh_dt,n(1)*n(2),1);
	end
end % swe_zero_inertia_2d_dhdt

