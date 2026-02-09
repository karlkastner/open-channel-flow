% Fri  4 Oct 14:13:58 CEST 2024
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
function [qxi] = interface_values_x(obj,t,h,convert,compute_derivatives)
	if (nargin()>3 && ~isempty(convert) && convert)
		h = obj.input_factor*h;
	end
	if (nargin()<5)
		compute_derivatives = false;
	end
	
	h = reshape(h,obj.nx);

	g   = obj.g;
	cd1 = obj.p.cd1;
	h   = real(h);
	zb  = obj.zb;
	dx  = obj.dx;

	siz_zb = size(zb);
	siz_h  = size(h);
	if (siz_zb(1) ~= siz_h(1)+2 || siz_zb(2) ~= siz_h(2)+2)
		error('zb must contain exterior values');
	end
	
	% h = max(h,0);
	% dh/dt = -d/x(u*h) 
	%       = -d/dx(C*h^(3/2)*sqrt(d/dx(h+zb)))
	%                   = -3/2*C*h^(1/2)*sqrt(S)*dh/dx -1/2*C*h^(3/2)/sqrt(S)*d2/dx2(h+zb)
	%                   = -3/2*u*dh/dx -1/2*h*u/S*d2/dx2(h+zb)

	% padd values extrapolated across the boundary, vector length increases to 2+1
	% TODO avoid copying by applying it only to interior grid cells and treating boundaries separately
	he  = padd_boundary_2d(t, h,obj.boundary_condition{3} ...
				,obj.boundary_condition{4} ...
				,obj.boundary_condition{1} ...
				,obj.boundary_condition{2} ...
				,dx);

	if (isscalar(obj.p.cd2))
		cd2i = obj.p.cd2;
	else
		siz_cd2  = size(obj.p.cd2);
		if (( siz_cd2(1) ~= (siz_h(1)+2) || siz_cd2(2) ~= (siz_h(2)+2) ) )
			error('cd2 must contain exterior values');
		end
		cd2i = 0.5*(obj.p.cd2(1:end-1,:)+obj.p.cd2(2:end,:));
	end

	% surface elevation
	zs = he + zb;

	% depth at interfaces
	hi = 0.5*(he(1:end-1,2:end-1)+he(2:end,2:end-1));

	% surface elevation at interfaces
	zsi = 0.5*(zs(1:end-1,:)+zs(2:end,:));

	% surface slope components perpendicular to interfaces
	Sxi = (zs(2:end,2:end-1) - zs(1:end-1,2:end-1))/dx(1);

	% surface slope parallel to interfaces
	Syi = (zsi(:,3:end) - zsi(:,1:end-2))/(2*dx(2));

	% surface slope magnitude at interfaces
	Smi = hypot(Sxi, Syi);

	% dh/dt + d(uh)/dx = 0
	% dh/dt + ai*dh/dx + di*d2h/dx2 = d2zb/dx2-term
	deni = sqrt(cd1.*cd1 + (4*g)*cd2i.*Smi.*hi.*hi.*hi);
	% velocity magnitude at interfaces
	umi  = -(cd1 - deni)./(2*cd2i*hi);
	umi(hi<obj.opt.heps) = 0;
	% velocity across interfaces
	uxi = -umi.*Sxi./Smi;
	uxi(Smi < obj.opt.Seps) = 0;
		
	% dq_dx = - (4*cd2*g*sign(S)*h^3*dS_dx + 12*cd2*g*abs(s)*h^2*dh_dx)/(4*cd2*sign(S)*sqrt(cd1^2 + 4*cd2*g*abs(S)*h^3)) 
	% celerity of the discharge wave = advection coefficient at interfaces
	ai = - (3*g*Sxi.*hi.*hi)./deni;
	% diffusion of the disharge wave, diffusion coefficient at the interface
	di =  (g*hi.*hi.*hi)./deni;

	% mesh peclet number
	Pei = -(0.5*dx(1))*ai./di;

	% weights of difference terms at left interface
	switch (obj.opt.spatial_discretization)
	case {'central'}
		bi = zeros(obj.nx(1)+1,obj.nx(2));
	case{'upwind'}
		bi = 2*(ai<0)-1;
	case{'optimal'}
		% magical coefficient
		bi = coth(Pei) - 1./Pei;
	case {'optimal-linearized'}
		% linearized magical coefficient
		bi = max(-1,min(1,Pei/3));
	otherwise
		error('here');
	end
	fdx     = ~isfinite(Pei) | ~isfinite(bi); % dl
	bi(fdx) = sign(uxi(fdx));

	% discharge across interfaces
	hiup = 0.5*( (1-bi).*he(1:end-1,2:end-1)  + (1+bi).*he(2:end,2:end-1));
	qxi    = uxi.*hiup;

	if (compute_derivatives)
		duxi_div_2hi       = uxi./(2*hi);
		duxi_div_2hi(hi<obj.opt.heps) = 0;
		duxi_dhl =  Syi.^2./Smi.^3.*(-1/dx(1)).*-umi - duxi_div_2hi - g*((1.5*Sxi.*hi) - (hi.*hi).*(Sxi./Smi).^2/dx(1))./deni;
		duxi_dhc =  Syi.^2./Smi.^3.*(+1/dx(2)).*-umi - duxi_div_2hi - g*((1.5*Sxi.*hi) + (hi.*hi).*(Sxi./Smi).^2/dx(1))./deni;
		%duxi_dhl(Smi<obj.opt.Seps) = 0;
		%duxi_dhc(Smi<obj.opt.Seps) = 0;
		obj.aux.dqxi_dhl = (duxi_dhl.*hiup + 0.5*uxi.*(1-bi));
		obj.aux.dqxi_dhc = (duxi_dhc.*hiup + 0.5*uxi.*(1+bi));
		fdx= (Smi<obj.opt.Seps);
		obj.aux.dqxi_dhl(fdx) = +g*(hi(fdx).*hi(fdx).*hi(fdx))/(cd1*dx(1));
		obj.aux.dqxi_dhc(fdx) = -g*(hi(fdx).*hi(fdx).*hi(fdx))/(cd1*dx(1));
	end
end % interface_values_x

