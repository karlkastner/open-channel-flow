% Fri 18 May 14:02:51 CEST 2018
% ikeda, 1976
%% transverse velocity profile in a meander bend
%
% u    : depth averaged velocity
% baru : depth averaged velocity at the centreline
% 
function v    = bend_transverse_velocity(r,h,u,baru,C)
	barF1 = h.*(log(2.56e2)-5.0).*(1.0./3.0);
	barF2 = h.*(-5.0./2.7e1)+(1.0./2.0)./h+1.0./h.^2.*(1.0./6.0) ...
		-h.*(1.0./h+1.0).^3.*(log(1.0./h+1.0)-1.0./3.0).*(1.0./3.0) ...
		+1.0./h.^2.*(h+1.0).^3.*(log(1.0./h+1.0).*-6.0+log(1.0./h+1.0).^2.*9.0+2.0).*(1.0./2.7e1) ...
		+4.0./2.7e1;

	kappa = Constant.KAPPA;
	F     = u/baru;
	v     = F.^2.*h./r*1./kappa*(barF1 - 1/kappa*sqrt(g)/C*barF2);
end

