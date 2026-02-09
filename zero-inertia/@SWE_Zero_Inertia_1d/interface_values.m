% Fri  4 Oct 14:13:58 CEST 2024
% Karl Kastner, Berlin
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
% Cary plover
% a*((1-b)/2dx*(hr - hl) + bp/dx*(hc - hl) + bn/dx*(hr - h))
% w : w=0, switches off h in the advection term for testing
% d/dx(h*d/dx(  h + z)) =   h*d2h/dx2 +   d/dx(h+z)*dh/dx + h*d2z/dx2
% d/dx(h*d/dx(w*h + z)) = w*h*d2h/dx2 + d/dx(w*h+z)*dh/dx + h*d2z/dx2
% function [qxi] = interface_values(obj,t,h,convert,comute_dir,compute_fp)
function qxi = interface_values(obj,t,h,convert,compute_derivatives,compute_fp_matrix)
	if (nargin()>3 && convert)
		h = obj.input_factor*h;
	end
	if (nargin()<5)
		compute_derivatives = false;
	end
	if (nargin()<6)
		compute_fp_matrix = false;
	end

	g   = obj.g;
	h   = rvec(h);
	h   = real(h);
	h = max(h,0);
	zb  = rvec(obj.zb);
	dx  = obj.dx;

	if (length(zb) ~= (length(h)+2))
		error('zb must contain exterior values');
	end
	
	% h = max(h,0);
	% dh/dt = -d/x(u*h) 
	%       = -d/dx(C*h^(3/2)*sqrt(d/dx(h+zb)))
	%                   = -3/2*C*h^(1/2)*sqrt(S)*dh/dx -1/2*C*h^(3/2)/sqrt(S)*d2/dx2(h+zb)
	%                   = -3/2*u*dh/dx -1/2*h*u/S*d2/dx2(h+zb)

	% padd values extrapolated across the boundary, vector length increases to 2+1
	% TODO padd boundary is slow, optimize
	he  = padd_boundary_1d(t,h,obj.boundary_condition{1},...
					obj.boundary_condition{2},dx);
	if (isscalar(obj.p.cd1))
		cd1i = obj.p.cd1;
	else
		if (length(obj.p.cd1) ~= (length(h)+2))
			error('cd1 must contain exterior values');
		end
		cd1i = 0.5*rvec(obj.p.cd1(1:end-1)+obj.p.cd1(2:end));
	end
	if (isscalar(obj.p.cd2))
		cd2i = obj.p.cd2;
	else
		if (length(obj.p.cd2) ~= (length(h)+2))
			error('cd2 must contain exterior values');
		end
		cd2i = 0.5*rvec(obj.p.cd2(1:end-1)+obj.p.cd2(2:end));
	end

	% surface elevation
	zs = he + zb;

	% depth at interfaces
	hi = 0.5*(he(1:end-1)+he(2:end));

	% surface slope across interfaces
	Sxi = (zs(2:end) - zs(1:end-1))/dx(1);

	% dh/dt + d(uh)/dx = 0
	% dh/dt + ai*dh/dx + di*(d2h/dx2 - d2zb/dx2) = 0
	% note that x.*x is faster than x^2, even for cubes
	deni = sqrt(cd1i.*cd1i + (4*g)*cd2i.*abs(Sxi).*hi.*hi.*hi);
	% velocity of interfaces
	if (cd2i==0)
		uxi = -g*(hi.*hi.*Sxi)./cd1i;
	else
		uxi  = sign(Sxi).*(cd1i - deni)./(2*cd2i.*hi);
	end
	uxi(hi < obj.opt.heps) = 0;

	% celerity of the discharge wave = advection coefficient at interfaces
	ai = - (3*g*Sxi.*hi.*hi)./deni;
	% diffusion of the disharge wave, diffusion coefficient at the interface
	di =  (g*hi.*hi.*hi)./deni;

	% mesh peclet number
	Pei = -(0.5*dx(1))*ai./di;

	% weights of difference terms at left interface
	% c.f. Donea 1983 eq 8, Carey 1983                                             
	switch (obj.opt.spatial_discretization)
	case {'central'}
		bi = zeros(1,obj.nx+1);
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
	% note that uxi can get undefined when hi is zero 
	fdx     = ~isfinite(Pei) | ~isfinite(bi);
	bi(fdx) = sign(uxi(fdx));

	% discharge across interfaces
	hiup = 0.5*( (1-bi).*he(1:end-1) + (1+bi).*he(2:end) );
	qxi  = uxi.*hiup;

	if (compute_derivatives)
		% qxi = uxi*hiup = Sxi/Smi*(cd1 - deni)/(2*cd2i.*hi)*hiup
		% qxi = uxi*hiup = Six/Si*(cd1 - sqrt(cd1^2 + 4*g*cd2i*|Sxi|*hi^3)/(2*cd2i.*hi))*hiup
	
		duxi_div_2hi = uxi./(2*hi);
		duxi_div_2hi(hi<obj.opt.heps) = 0;
		s = 1;
		duxi_dhl = - duxi_div_2hi - g*((1.5*Sxi.*hi) - s*(hi.*hi)/dx(1))./deni;
		duxi_dhc = - duxi_div_2hi - g*((1.5*Sxi.*hi) + s*(hi.*hi)/dx(1))./deni;
		obj.aux.dqxi_dhl = (duxi_dhl.*hiup + 0.5*uxi.*(1-bi));
		obj.aux.dqxi_dhc = (duxi_dhc.*hiup + 0.5*uxi.*(1+bi));
	end

	if (compute_fp_matrix)
%	qx = -|q/Sx| Sx
%	qx = -|q/Sx/dx| dzs
%          = -|ux hiup/Sx| Sx
%          = -|sign(Sxi).*(cd1 - deni)./(2*cd2i.*hi) hiup/Sx| Sx
%          = -|(cd1 - deni)./(2*cd2i.*hi) hiup/Sx| Sx
%          = -|Sxi/Smi.*(cd1 - deni)./(2*cd2i.*hi) hiup| Sx
%	   = Di Sx
%	(g hup Sx + cd1 qx + cd2 qm qx == 0)
%	 qx = g hup Sx/(cd1 + cd2 qm)

		%bi=-bi;
		Di = g*(hi.*hi.*hiup)./(obj.p.cd1i + obj.p.cd2i.*abs(qxi));	
		%Di = g*((0.5*((1-bi).*he(1:end-1) + (1+bi).*he(2:end))).^3)./(obj.p.cd1 + obj.p.cd2*abs(qxi));
		obj.aux.Di = cvec(Di)/(dx(1)*dx(1));
	end

	if (nargin()>3 && ~isempty(convert) && convert)
		qxi = obj.output_factor*qxi/obj.input_factor;
		hi = hi/obj.input_factor;	
		he = he/obj.input_factor;	
		uxi = obj.output_factor*uxi;
		ai = obj.output_factor*ai;
		di = obj.output_factor*di;
	end

	if (obj.opt.output.store_interface_values)
	obj.aux.qxi = qxi;
	obj.aux.he = he;
	obj.aux.hi = hi;
	obj.aux.hiup = hiup;
	obj.aux.deni = deni;
	obj.aux.uxi = uxi;
	obj.aux.bi = bi;
	obj.aux.Sxi = Sxi;	
	obj.aux.ai = ai;
	obj.aux.di = di;
	obj.aux.Pei = Pei;
	end
end % interface_values

