% Tue 20 May 20:09:28 CEST 2025
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
function compute_derivatives(obj,h)
	flag = 0;
	%h    = rvec(h);
	h    = real(h);
	zb   = obj.zb;
	g    = obj.g;

	dx = obj.dx;
	if (true) %obj.analytic_derivatives)
		% TODO trap cases where S = 0
		[qi,qj,hi,hj,he,ui,uj] = obj.interface_values(h);

		umi = obj.aux.umi;
		umj = obj.aux.umj;
		bi = obj.aux.bi;
		bj = obj.aux.bj;
		Si = obj.aux.Si;
		Sj = obj.aux.Sj;
		dzs_dx_i = obj.aux.dzs_dx_i;
		dzs_dy_i = obj.aux.dzs_dy_i;
		dzs_dx_j = obj.aux.dzs_dx_j;
		dzs_dy_j = obj.aux.dzs_dy_j;
		deni = obj.aux.deni;
		denj = obj.aux.denj;
		ai = 0.25*(g*hi.^2.*dzs_dx_i./deni + ui).*dzs_dy_i./(dx(2)*Si.^2) ...
			.*( (1-bi).*he(1:end-1,2:end-1)  + (1+bi).*he(2:end,2:end-1));
		fdx = (Si == 0);
		ai_ = 0.25*(g*hi.^2./deni)./(dx(2)) ...
			.*( (1-bi).*he(1:end-1,2:end-1)  + (1+bi).*he(2:end,2:end-1));
		ai(fdx) = ai_(fdx);
		
		dqi_dhlu_ = -ai;
		
				%deni = sqrt(c1^2 + 4*c2*g*Si.*hi.^3);
				%ui  = (c1 - deni)./(2*c2*hi);
				% velocity across interfaces
				%ui = ui.*dzs_dx_i./Si;
		
		dqi_dhlc_ = 0.5*( - ui./(dzs_dx_i*dx(1)) ...
			    - ui./(2*hi) ...
			    + ui.*dzs_dx_i./(dx(1)*Si.^2) ...
			    - g*hi.^2.*dzs_dx_i./deni.*(1.5./hi - dzs_dx_i./(dx(1)*Si.^2))  ...
			    ).*( (1-bi).*he(1:end-1,2:end-1)  + (1+bi).*he(2:end,2:end-1)) ...
			   + 0.5*ui.*( (1-bi) );
		dqi_dhlc_(ui == 0) = 0;
		fdx_ = (dzs_dx_i == 0);
		[sum(fdx(:)),sum(fdx_(:))]
		dqi_dhlc__ = 0.5*( - umi./(Si*dx(1))).*( (1-bi).*he(1:end-1,2:end-1)  + (1+bi).*he(2:end,2:end-1));
		dqi_dhlc_(fdx_) = dqi_dhlc__(fdx_);
		dqi_dhlc__ =  0.5*( - g*hi.^2.*1./deni.*( - 1./(dx(1))) ... 
				.*( (1-bi).*he(1:end-1,2:end-1)  + (1+bi).*he(2:end,2:end-1))) ;
		dqi_dhlc_(fdx) = dqi_dhlc__(fdx);
		
		dqi_dhld_ =   ai;
		dqi_dhcu_ =  -ai;
		% TODO ui/dzs_dx_i can be non-zero 
		dqi_dhcc_ =   0.5*(   ui./(dx(1)*dzs_dx_i) ...
				- ui./(2*hi) ...
				- ui.*dzs_dx_i./(dx(1)*Si.^2) ...
				- g*hi.^2.*dzs_dx_i./deni.*(1.5./hi + dzs_dx_i./(dx(1)*Si.^2)) ...
				).*( (1-bi).*he(1:end-1,2:end-1)  + (1+bi).*he(2:end,2:end-1)) ...
				+ 0.5*ui.*( (1+bi) );
		dqi_dhcc_(ui == 0) = 0;
		dqi_dhcc__ = 0.5*( + umi./(Si*dx(1))).*( (1-bi).*he(1:end-1,2:end-1)  + (1+bi).*he(2:end,2:end-1));
		dqi_dhcc_(fdx_) = dqi_dhcc__(fdx_);
		dqi_dhcc__ =  0.5*( - g*hi.^2.*1./deni.*( + 1./(dx(1))) ) ... 
				.*( (1-bi).*he(1:end-1,2:end-1)  + (1+bi).*he(2:end,2:end-1));
		dqi_dhcc_(fdx) = dqi_dhcc__(fdx);
		
		dqi_dhcd_ = ai;
		fdx = (Sj == 0);
		fdx_ = (dzs_dy_j == 0);
		[sum(fdx(:)),sum(fdx_(:))]
		aj =  0.25*(g*hj.^2.*dzs_dy_j./denj + uj).*dzs_dx_j./(dx(1)*Sj.^2) ...
			.*( (1-bj).*he(2:end-1,1:end-1)  + (1+bj).*he(2:end-1,2:end));
		aj_ =  0.25*(g*hj.^2./denj)./(dx(1)) ...
			.*( (1-bj).*he(2:end-1,1:end-1)  + (1+bj).*he(2:end-1,2:end));
		aj(fdx) = aj_(fdx);
		
		dqj_dhlu_ = aj;
		dqj_dhcu_ = 0.5*(	-  uj./(dx(2)*dzs_dy_j) ...
				-  uj./(2*hj) ...
				+  uj.*dzs_dy_j./(dx(2)*Sj.^2) ...
				- g*hj.^2.*dzs_dy_j./denj.*(1.5./hj - dzs_dy_j./(dx(2)*Sj.^2))...
			        ).*( (1-bj).*he(2:end-1,1:end-1)  + (1+bj).*he(2:end-1,2:end)) ...
			+ 0.5*uj.*( (1-bj) );
		dqj_dhcu_(uj == 0) = 0;
		dqj_dhcu__ = 0.5*( -  umj./(dx(2)*Sj) ).*( (1-bj).*he(2:end-1,1:end-1)  + (1+bj).*he(2:end-1,2:end));
		dqj_dhcu_(fdx_) = dqj_dhcu__(fdx_);
		
		dqj_dhcu__ = 0.5*(	- g*hj.^2.*1./denj.*(- 1./(dx(2))) ...
			  .*( (1-bj).*he(2:end-1,1:end-1)  + (1+bj).*he(2:end-1,2:end)));
		dqj_dhcu_(fdx) = dqj_dhcu__(fdx);
		%max(abs(dqj_dhcu__(:)))
		%max(abs(hj(:)))
		%max(abs(denj(:)))
		%min(abs(denj(:)))
		%pause
		
		dqj_dhru_ = -aj;
		dqj_dhlc_ = aj;
		dqj_dhcc_ = 0.5*( uj./(dx(2)*dzs_dy_j) ... 
			  - uj./(2*hj) ...
			  - uj.*dzs_dy_j./(dx(2)*Sj.^2) ...
			  - g*hj.^2.*dzs_dy_j./denj.*(1.5./hj + dzs_dy_j./(dx(2)*Sj.^2))  ... 
			  ).*( (1-bj).*he(2:end-1,1:end-1)  + (1+bj).*he(2:end-1,2:end)) ...
			  + 0.5*uj.*( (1+bj) );
		dqj_dhcc_(uj == 0) = 0;
		dqj_dhcu__ = 0.5*(  +  umj./(dx(2)*Sj) ).*( (1-bj).*he(2:end-1,1:end-1)  + (1+bj).*he(2:end-1,2:end));
		dqj_dhcu_(fdx_) = dqj_dhcu__(fdx_);
		dqj_dhcc__ = 0.5*(	- g*hj.^2.*1./denj.*(+ 1./(dx(2)))) ...
			  .*( (1-bj).*he(2:end-1,1:end-1)  + (1+bj).*he(2:end-1,2:end));
		dqj_dhcc_(fdx) = dqj_dhcc__(fdx);
		dqj_dhrc_ = -aj;

		obj.aux.dqi_dhlu = dqi_dhlu;
		obj.aux.dqi_dhlc = dqi_dhlc;
		obj.aux.dqi_dhld = dqi_dhld;
		obj.aux.dqi_dhcu = dqi_dhcu;
		obj.aux.dqi_dhcc = dqi_dhcc;
		obj.aux.dqi_dhcd = dqi_dhcd;

		obj.aux.dqj_dhlu = dqj_dhlu;
		obj.aux.dqj_dhlc = dqj_dhlc;
		obj.aux.dqj_dhld = dqj_dhld;
		obj.aux.dqj_dhcu = dqj_dhcu;
		obj.aux.dqj_dhcc = dqj_dhcc;
		obj.aux.dqj_dhcd = dqj_dhcd;


	else
		% qi-derivatives
		for i1=1:obj.n(1)
		 for i2=1:obj.n(2)
			% lu
			id_ = left(up(id,1),1);
			hl  = h;
			hr(id_(i1,i2)) = hr(id_(i1,i2)) + e;
			[qil,qj] = obj.interface_values(h_);
			hl(id_(i1,i2)) = hl(id_(i1,i2)) - e;
			[qir,qj] = obj.interface_values(h_);
			dqi_dhlu = (qir(id)-qil(id))/(2*e);
			% lc
			id_ = left(id,1);
			hl  = h;
			hr(id_(i1,i2)) = hr(id_(i1,i2)) + e;
			[qil,qj] = obj.interface_values(h_);
			hl(id_(i1,i2)) = hl(id_(i1,i2)) - e;
			[qir,qj] = obj.interface_values(h_);
			dqi_dhlc = (qir(id)-qil(id))/(2*e);
			% ld
			id_ = left(down(id,1),1);
			hl  = h;
			hr(id_(i1,i2)) = hr(id_(i1,i2)) + e;
			[qil,qj] = obj.interface_values(h_);
			hl(id_(i1,i2)) = hl(id_(i1,i2)) - e;
			[qir,qj] = obj.interface_values(h_);
			dqi_dhld = (qir(id)-qil(id))/(2*e);
			% cu
			id_ = up(id,1);
			hl  = h;
			hr(id_(i1,i2)) = hr(id_(i1,i2)) + e;
			[qil,qj] = obj.interface_values(h_);
			hl(id_(i1,i2)) = hl(id_(i1,i2)) - e;
			[qir,qj] = obj.interface_values(h_);
			dqi_dhcu = (qir(id)-qil(id))/(2*e);
			% cc
			id_ = (id);
			hl  = h;
			hr(id_(i1,i2)) = hr(id_(i1,i2)) + e;
			[qil,qj] = obj.interface_values(h_);
			hl(id_(i1,i2)) = hl(id_(i1,i2)) - e;
			[qir,qj] = obj.interface_values(h_);
			dqi_dhcc = (qir(id)-qil(id))/(2*e);
			% cd
			id_ = down(id,1);
			hl  = h;
			hr(id_(i1,i2)) = hr(id_(i1,i2)) + e;
			[qil,qj] = obj.interface_values(h_);
			hl(id_(i1,i2)) = hl(id_(i1,i2)) - e;
			[qir,qj] = obj.interface_values(h_);
			dqi_dhcd = (qir(id)-qil(id))/(2*e);
		 end
		end
		
	end
end

