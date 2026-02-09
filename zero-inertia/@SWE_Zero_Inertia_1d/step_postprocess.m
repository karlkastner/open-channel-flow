% 2025-09-09 12:48:27.831974714 +0200
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
function step_postprocess(obj,hnew)
	dt = obj.aux.dt;
	% call parent class function
	step_postprocess@RAD_Model(obj,hnew);
	%obj.aux.told,hnew,dt);
	if (obj.opt.output.store_fluxes)
		odx = obj.aux.odx;
		c   = obj.aux.butcher_table.node;
		w   = obj.aux.butcher_table.weight;
		% TODO this is only correct for theta-schemes and midpoint	
		% for other schemes the stage values have to be stored
		store_iv = obj.opt.output.store_interface_values;
		obj.opt.output.store_interface_values = true;
		rp    = obj.precipitation_rate(obj.aux.told,dt);
		%sc = 0;
		obj.out.duration(odx)         = obj.out.duration(odx) + dt*(rp>0);
		for idx=1:length(c)
			t    = obj.aux.told+c(idx)*dt;
			h    = (1-c(idx))*obj.aux.zold(:) + c(idx)*hnew(:);
			%dqdx = obj.dh_dt(t,h);
			dqdx = obj.dq_dx(t,h);
			obj.out.h(odx,:) = obj.out.h(odx,:) + (dt*w(idx))*rvec(h);
			obj.out.flow_x(odx,:) = obj.out.flow_x(odx,:) + (dt*w(idx))*rvec(obj.aux.qxi(:));
			obj.out.precipitation(odx,1) = obj.out.precipitation(odx,1) + (dt*w(idx))*rp;
			% dh_dt without infiltration
			dh_dt_woi                   = dqdx + rp;
			ri                          = obj.infiltration_rate(h,dh_dt_woi(:));
			obj.out.infiltration(odx,:) = obj.out.infiltration(odx,:) + (dt*w(idx)).*ri';
			obj.out.celerity_x(odx,:)   = obj.out.celerity_x(odx,:)   + (dt*w(idx)).*rvec(obj.aux.ai(:));
			obj.out.diffusion_x(odx,:)  = obj.out.diffusion_x(odx,:)  + (dt*w(idx)).*rvec(obj.aux.di(:));
			obj.out.twet(odx,:)         = obj.out.twet(odx,:)         + (dt*w(idx)).*(rvec(h>obj.opt.heps));
		%	sc = dt*w(idx)*(cvec(dqdx) + rp - cvec(ri));
		%	dqdx_ = -diff(obj.aux.qxi(:)',[],2)/obj.dx;
		if (1 == idx)
			% time step estimate
			dx = obj.dx;
			if (obj.ndim == 1)
			dt_zie = 1./(abs(mid(obj.aux.ai))/dx(1)  + 2*mid(obj.aux.di)/dx(1)^2);
			else
			dt_zie = 1./(  abs(mid(obj.aux.ai))/dx(1)   + 2*mid(obj.aux.di)/dx(1)^2 ...
			             + abs(mid(obj.aux.aj,2))/dx(2) + 2*mid(obj.aux.dj,2)/dx(2)^2);
			end
			fdx    = (h>obj.opt.heps);
			dt_zie = min(dt_zie(fdx),[],'all');
			if (~isempty(dt_zie))
			n_step_zie = dt./dt_zie;
			obj.out.n_step_zie(odx) = obj.out.n_step_zie(odx) + n_step_zie;
			end
			% note: the eigenvalues of the swe are
			% ux +/- sqrt(gh), and uv +/- sqrt(gh)
			% dimensional splitting will require 
			hi = obj.aux.hiup; %inner2outer(reshape(h,obj.nx));
			dt_swe = dx(1)./(abs(obj.aux.qxi)./hi + sqrt(obj.g*hi));
			fdx = (hi>obj.opt.heps);
			dt_swe = min(dt_swe(fdx),[],'all');
			if (obj.ndim == 2)
				%hj = inner2outer(reshape(h,obj.nx),2);
				hj = obj.aux.hjup; %inner2outer(reshape(h,obj.nx));
				dt_swe_ = dx(2)./(abs(obj.aux.qyj)./hj + sqrt(obj.g*hj));
				fdx = (hj>obj.opt.heps);
				dt_swe_ = min(dt_swe_(fdx),[],'all');
				if (~isempty(dt_swe_))
					dt_swe = min(dt_swe,dt_swe_);
				end
			end
			if (~isempty(dt_swe))
			n_step_swe = dt./dt_swe;
			obj.out.n_step_swe(odx) = obj.out.n_step_swe(odx) + n_step_swe;
			end
		end
		end % for length c
		%rms((hnew-obj.aux.zold) - cvec(sc))/rms(hnew-obj.aux.zold)
		obj.opt.output.store_interface_values = store_iv;
	end % if store_fluxes
end % step_postporcess

