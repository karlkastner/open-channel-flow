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
function [qyj] = interface_values_j(obj,t,h,convert,compute_derivatives)
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
				,obj.boundary_condition{4} ...
				,obj.boundary_condition{1} ...
				,obj.boundary_condition{2} ...
				,dx);

	if (isscalar(obj.p.cd2))
		cd2j = obj.p.cd2;
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
	hj = 0.5*(he(2:end-1,1:end-1)+he(2:end-1,2:end));

	% surface elevation at interfaces
	zsj = 0.5*(zs(:,1:end-1)+zs(:,2:end));

	% surface slope components perpendicular to interfaces
	Syj = (zs(2:end-1,2:end)  - zs(2:end-1,1:end-1))/dx(2);

	% surface slope parallel to interfaces
	Sxj = (zsj(3:end,:) - zsj(1:end-2,:))/(2*dx(1));

	% surface slope magnitude at interfaces
	Smj = hypot(Sxj, Syj);

	% dh/dt + d(uh)/dx = 0
	% dh/dt + ai*dh/dx + di*d2h/dx2 = d2zb/dx2-term
	denj = sqrt(cd1.*cd1 + (4*g)*cd2j.*Smj.*hj.*hj.*hj);
	% velocity magnitude at interfaces
	umj  = -(cd1 - denj)./(2*cd2j.*hj);
	umj(hj<obj.opt.heps) = 0;
	% velocity across interfaces
	uyj = -umj.*Syj./Smj;
	uyj(Smj < obj.opt.Seps) = 0;
		
	% dq_dx = - (4*cd2*g*sign(S)*h^3*dS_dx + 12*cd2*g*abs(s)*h^2*dh_dx)/(4*cd2*sign(S)*sqrt(cd1^2 + 4*cd2*g*abs(S)*h^3)) 
	% celerity of the discharge wave = advection coefficient at interfaces
	aj = - (3*g*Syj.*hj.*hj)./denj;
	% diffusion of the disharge wave, diffusion coefficient at the interface
	dj =  (g*hj.*hj.*hj)./denj;

	% mesh peclet number
	Pej = -(0.5*dx(2))*aj./dj;

	% weights of difference terms at left interface
	switch (obj.opt.spatial_discretization)
	case {'central'}
		bj = zeros(obj.nx(1),obj.nx(2)+1);
	case{'upwind'}
		bj = 2*(aj<0)-1;
	case{'optimal'}
		% magical coefficient
		bj = coth(Pej) - 1./Pej;
	case {'optimal-linearized'}
		% linearized magical coefficient
		bj = max(-1,min(1,Pej/3));
	otherwise
		error('here');
	end
	fdx     = ~isfinite(Pej) | ~isfinite(bj); % dl
	bj(fdx) = sign(uyj(fdx));

	% discharge across interfaces
	hjup = 0.5*( (1-bj).*he(2:end-1,1:end-1)  + (1+bj).*he(2:end-1,2:end));
	qyj    = uyj.*hjup;

	if (compute_derivatives)
		duyj_div_2hj       = uyj./(2*hj);
		duyj_div_2hj(hj<obj.opt.heps) = 0;
		duyj_dhl =  Sxj.^2./Smj.^3.*(-1/dx(2)).*-umj - duyj_div_2hj - g*((1.5*Syj.*hj) - (hj.*hj).*(Syj./Smj).^2/dx(2))./denj;
		duyj_dhc =  Sxj.^2./Smj.^3.*(+1/dx(2)).*-umj - duyj_div_2hj - g*((1.5*Syj.*hj) + (hj.*hj).*(Syj./Smj).^2/dx(2))./denj;
		obj.aux.dqyj_dhl = (duyj_dhl.*hjup + 0.5*uyj.*(1-bj));
		obj.aux.dqyj_dhc = (duyj_dhc.*hjup + 0.5*uyj.*(1+bj));
		fdx = (Smj<obj.opt.Seps);
		obj.aux.dqyj_dhl(fdx) = +g*(hj(fdx).*hj(fdx).*hj(fdx))/(cd1*dx(1));
		obj.aux.dqyj_dhc(fdx) = -g*(hj(fdx).*hj(fdx).*hj(fdx))/(cd1*dx(1));
	end
end % interface_values_y

