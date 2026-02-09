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
%
% note: the function returns dh/dt, is is named dz_dt to be comptible
% TODO merge with interface values x, y
function [qxi,qyj] = interface_values(obj,t,h,convert,compute_derivatives,compute_fp_matrix)
	if (nargin()>3 && ~isempty(convert) && convert)
		h = obj.input_factor*h;
	end
	if (nargin()<5)
		compute_derivatives = false;
	end
	if (nargin()<6)
		compute_fp_matrix = false;
	end

	h = reshape(h,obj.nx);

	g   = obj.g;
	cd1 = obj.p.cd1;
	h   = real(h);
	h = max(h,0);
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
	he  = padd_boundary_2d(t, h, obj.boundary_condition{3} ...
				, obj.boundary_condition{4} ...
				, obj.boundary_condition{1} ...
				, obj.boundary_condition{2} ...
				, dx);

	if (isscalar(obj.p.cd1))
		cd1i = obj.p.cd1;
		cd1j = obj.p.cd1;
	else
		siz_cd1  = size(obj.p.cd1);
		if (( siz_cd1(1) ~= (siz_h(1)+2) || siz_cd1(2) ~= (siz_h(2)+2) ) )
			error('cd1 must contain exterior values');
		end
		cd1i = 0.5*(obj.p.cd1(1:end-1,2:end-1)+obj.p.cd1(2:end,2:end-1));
		cd1j = 0.5*(obj.p.cd1(2:end-1,1:end-1)+obj.p.cd1(2:end-1,2:end));
	end
	if (isscalar(obj.p.cd2))
		cd2i = obj.p.cd2;
		cd2j = obj.p.cd2;
	else
		siz_cd2  = size(obj.p.cd2);
		if (( siz_cd2(1) ~= (siz_h(1)+2) || siz_cd2(2) ~= (siz_h(2)+2) ) )
			error('cd2 must contain exterior values');
		end
		cd2i = 0.5*(obj.p.cd2(1:end-1,2:end-1)+obj.p.cd2(2:end,2:end-1));
		cd2j = 0.5*(obj.p.cd2(2:end-1,1:end-1)+obj.p.cd2(2:end-1,2:end));
	end

	% surface elevation
	zs = he + zb;

	% qi1 is flow between cell beyond boundary and first cell
	% depth at interfaces
	hi = 0.5*(he(1:end-1,2:end-1)+he(2:end,2:end-1));
	hj = 0.5*(he(2:end-1,1:end-1)+he(2:end-1,2:end));

	% surface elevation at interfaces
	zsi = 0.5*(zs(1:end-1,:)+zs(2:end,:));
	zsj = 0.5*(zs(:,1:end-1)+zs(:,2:end));

	% surface slope components perpendicular to interfaces
	Sxi = (zs(2:end  ,2:end-1) - zs(1:end-1,2:end-1))/dx(1);
	Syj = (zs(2:end-1,2:end  ) - zs(2:end-1,1:end-1))/dx(2);

	% surface slope parallel to interfaces
	Syi = (zsi(:,3:end) - zsi(:,1:end-2))/(2*dx(2));
	Sxj = (zsj(3:end,:) - zsj(1:end-2,:))/(2*dx(1));

	% surface slope magnitude at interfaces
	Smi = hypot(Sxi, Syi);
	Smj = hypot(Sxj, Syj);

	% dh/dt + d(uh)/dx = 0
	% dh/dt + ai*dh/dx + di*d2h/dx2 = d2zb/dx2-term
	deni = sqrt(cd1i.*cd1i + (4*g)*cd2i.*Smi.*hi.*hi.*hi);
	denj = sqrt(cd1j.*cd1j + (4*g)*cd2j.*Smj.*hj.*hj.*hj);
	% velocity magnitude at interfaces
	umi  = -(cd1i - deni)./(2*cd2i.*hi);
	if (cd2i == 0)
	%fdx = (c2 == 0);
	%u = -g*(h.*h.*S)./c1;
		umi = g*(hi.*hi.*Smi)./cd1i;
	end
	umi(hi<obj.opt.heps) = 0;
	umj  = -(cd1j - denj)./(2*cd2j.*hj);

	umj(hj<obj.opt.heps) = 0;
	if (cd2j == 0)                                                         
        %fdx = (c2 == 0);                                                        
        %u(fdx) = -g*(h(fdx).*h(fdx).*S(fdx))./c1(fdx);                          
        umj = g*(hj.*hj.*Smj)./cd1j;                          
        end 
	% velocity across interfaces
	uxi = -umi.*Sxi./Smi;
	% TODO flag = sqrt(eps)*cd1.*cd1 > (4*g)*cd2i.*Smi.*hi.*hi.*hi
	uxi(Smi < obj.opt.Seps) = 0;
	uyj = -umj.*Syj./Smj;
	uyj(Smj < obj.opt.Seps) = 0;

	% TODO capture the case cd2 == 0
	
	% ai = (uxi - uxi./hi  + 3*g*Si.*hi./sqrt(cd1.^2 + 4*cd2*g*abs(Si).*hi.^3));
	% dq_dx = - (4*cd2*g*sign(S)*h^3*dS_dx + 12*cd2*g*abs(s)*h^2*dh_dx)/(4*cd2*sign(S)*sqrt(cd1^2 + 4*cd2*g*abs(S)*h^3)) 
	% celerity of the discharge wave = advection coefficient at interfaces
	ai = - (3*g*Sxi.*hi.*hi)./deni;
	aj = - (3*g*Syj.*hj.*hj)./denj;

	% diffusion of the disharge wave, diffusion coefficient at the interface
	di =  (g*hi.*hi.*hi)./deni;
	dj =  (g*hj.*hj.*hj)./denj;
			
	% mesh peclet number
	% TODO this treats the dimensions separately, I am not sure if this is the
	% appropriate extension of the 1D scheme to 2D,
	% maybe it is better to compute Pe in the direction of the slope,
	% combined with central differences in the direction perpendicular to the slope,
	% and than rotate the factors back to the grid
	Pei = -(0.5*dx(1))*ai./di;
	Pej = -(0.5*dx(2))*aj./dj;

	% weights of difference terms at left interface
	switch (obj.opt.spatial_discretization)
	case {'central'}
		bi = zeros(obj.nx(1)+1,obj.nx(2));
		bj = zeros(obj.nx(1),obj.nx(2)+1);
	case{'upwind'}
		bi = 2*(ai<0)-1;
		bj = 2*(aj<0)-1;
	case{'optimal'}
		% magical coefficient
		bi = coth(Pei) - 1./Pei;
		bj = coth(Pej) - 1./Pej;
	case {'optimal-linearized'}
		% linearized magical coefficient
		bi = max(-1,min(1,Pei/3));
		bj = max(-1,min(1,Pej/3));
	otherwise
		error('here');
	end
	fdx     = ~isfinite(Pei) | ~isfinite(bi); % dl
	bi(fdx) = sign(uxi(fdx));
	fdx     = ~isfinite(Pej) | ~isfinite(bj); % dl
	bj(fdx) = sign(uyj(fdx));

	% discharge across interfaces
	hiup = 0.5*( (1-bi).*he(1:end-1,2:end-1)  + (1+bi).*he(2:end,2:end-1));
	qxi   = uxi.*hiup;
	hjup = 0.5*( (1-bj).*he(2:end-1,1:end-1)  + (1+bj).*he(2:end-1,2:end));
	qyj   = uyj.*hjup;

	if (compute_derivatives)
		% Jij = d/dhj dqx_i/du + d/dhj qy_i/dy
		% d/dh dq/dx = d/dh dq/dx( udir*umag*hup )
		% note that the subscripts in Jij refer to cell indices
		% but subscripts in hi and hj refer to values at horizontal and vertical interfaces, respectively

% error here should go dzs_dxi !!!
% Si = sqrt((hc - hl + zbc - zbl)^2/dx^2 + Sij^2)
% d/dhl Sim = 1/2*2*(-1/dx)*Sxi/Si = -Sxi/(Sim*dx) not 1/dx !!!!
% d/dhc Sim = 2/dx*1/2*Sxi/Si = +Sxi/(Si*dx) not 1/dx !!!!
% 
% deni = sqrt(cd1.*cd1 + (4*g)*cd2i.*abs(Si).*hi.*hi.*hi);
% uxi  = (cd1 - deni)./(2*cd2i.*hi)*Sxi/Si;
% qxi = uxi*hiup = Sxi/Si*(cd1 - deni)/(2*cd2i.*hi)*hiup
% qxi = uxi*hiup = Sxi/Si*(cd1 - sqrt(cd1^2 + 4*g*cd2i*|Si|*hi^3)/(2*cd2i.*hi))*hiup
% qxi = uxi*hiup 
%    = Sxi/Si*(cd1 - sqrt(cd1^2 + 4*g*cd2i*sqrt(Sxi^3 + Syi^2)*hi^3)/(2*cd2i.*hi))*hiup
%    TODO add derivative for d/dhl Sxi/Si = Syi^2/(Sxi^2 + Syi^2)^(3/2)*dSi/dhl
if (obj.opt.analytic_flux_derivatives)

		duxi_div_2hi       = uxi./(2*hi);
		duxi_div_2hi(hi<obj.opt.heps) = 0;
		duxi_dhl =  Syi.^2./Smi.^3.*(-1/dx(1)).*-umi - duxi_div_2hi - g.*((1.5*Sxi.*hi) - (hi.*hi).*(Sxi./Smi).^2/dx(1))./deni;
		duxi_dhc =  Syi.^2./Smi.^3.*(+1/dx(2)).*-umi - duxi_div_2hi - g.*((1.5*Sxi.*hi) + (hi.*hi).*(Sxi./Smi).^2/dx(1))./deni;
		obj.aux.dqxi_dhl = (duxi_dhl.*hiup + 0.5*uxi.*(1-bi));
		obj.aux.dqxi_dhc = (duxi_dhc.*hiup + 0.5*uxi.*(1+bi));
		fdx= (Smi<obj.opt.Seps);
		if (isscalar(cd1i))
		obj.aux.dqxi_dhl(fdx) = +g*(hi(fdx).*hi(fdx).*hi(fdx))./(cd1i*dx(1));
		obj.aux.dqxi_dhc(fdx) = -g*(hi(fdx).*hi(fdx).*hi(fdx))./(cd1i*dx(1));
		else
		obj.aux.dqxi_dhl(fdx) = +g*(hi(fdx).*hi(fdx).*hi(fdx))./(cd1i(fdx)*dx(1));
		obj.aux.dqxi_dhc(fdx) = -g*(hi(fdx).*hi(fdx).*hi(fdx))./(cd1i(fdx)*dx(1));
		end
		
		duyj_div_2hj       = uyj./(2*hj);
		duyj_div_2hj(hj<obj.opt.heps) = 0;
		duyj_dhl =  Sxj.^2./Smj.^3.*(-1/dx(2)).*-umj - duyj_div_2hj - g.*((1.5*Syj.*hj) - (hj.*hj).*(Syj./Smj).^2/dx(2))./denj;
		duyj_dhc =  Sxj.^2./Smj.^3.*(+1/dx(2)).*-umj - duyj_div_2hj - g.*((1.5*Syj.*hj) + (hj.*hj).*(Syj./Smj).^2/dx(2))./denj;
		obj.aux.dqyj_dhl = (duyj_dhl.*hjup + 0.5*uyj.*(1-bj));
		obj.aux.dqyj_dhc = (duyj_dhc.*hjup + 0.5*uyj.*(1+bj));
		fdx = (Smj<obj.opt.Seps);
		if (isscalar(cd1j))
			obj.aux.dqyj_dhl(fdx) = +g*(hj(fdx).*hj(fdx).*hj(fdx))./(cd1j*dx(2));
			obj.aux.dqyj_dhc(fdx) = -g*(hj(fdx).*hj(fdx).*hj(fdx))./(cd1j*dx(2));
		else
		obj.aux.dqyj_dhl(fdx) = +g*(hj(fdx).*hj(fdx).*hj(fdx))./(cd1j(fdx)*dx(2));
		obj.aux.dqyj_dhc(fdx) = -g*(hj(fdx).*hj(fdx).*hj(fdx))./(cd1j(fdx)*dx(2));
		end
else
		% note: this is for testing purposes and does not compute values at the boundary
		[obj.aux.dqxi_dhl,obj.aux.dqxi_dhc,obj.aux.dqyj_dhl,obj.aux.dqyj_dhc] = obj.dq_dh(t,h); 
end

	
	end

	if (compute_fp_matrix)
	%	qx = -|q/S| Sx
	%	   = Di Sx

		%bi=-bi;
		Di = g*(hi.*hi.*hiup)./(obj.p.cd1 + obj.p.cd2.*(umi.*hiup));	
		Dj = g*(hj.*hj.*hjup)./(obj.p.cd1 + obj.p.cd2.*(umj.*hjup));	
		% 1/dx for dq/dx and 1/dx for Sx = dh/dx
		obj.aux.Di = Di/(dx(1)*dx(2));
		obj.aux.Dj = Dj/(dx(2)*dx(2));
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
		obj.aux.qxi  = qxi;
		obj.aux.qyj  = qyj;
		obj.aux.he  = he;
		obj.aux.hi  = hi;
		obj.aux.hiup = hiup;
		obj.aux.hjup = hjup;
		obj.aux.hj  = hj;
		obj.aux.deni = deni;
		obj.aux.denj = denj;
		obj.aux.uxi  = uxi;
		obj.aux.umi  = umi;
		obj.aux.umj  = umj;
		obj.aux.uyj  = uyj;
		obj.aux.bi  = bi;
		obj.aux.bj  = bj;
		obj.aux.Sxi  = Sxi;
		obj.aux.Syi  = Syi;
		obj.aux.Sxj  = Sxj;
		obj.aux.Syj  = Syj;
		obj.aux.Smi  = Smi;
		obj.aux.Smj  = Smj;
		obj.aux.ai  = ai;
		obj.aux.aj  = aj;
		obj.aux.di  = di;
		obj.aux.dj  = dj;
		obj.aux.Pei = Pei;
		obj.aux.Pej = Pej;
	end
	
end % interface_values

