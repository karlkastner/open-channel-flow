% 2024-10-07 20:17:54.410620707 +0200
% run for 10 time steps 
% dh/dt = K*d/dx(h*dh/dx)
%
% implicit euler:
% h(t+dt) - h(t) = dt*K*d/dx*(h_next*(dh_next/dx - dzb/dx))
% (I - dt*K*D*(h_n*D)*h_n - dt*K*D*(h*D*zb/dx) = h
%
% advection part:	Dz*(Dh) + d^2z/dx^2*h
%
% res  = (h - (I - dt*K*D*h_n*D)*h_n)
% fobj = 1/2*res'*res
% gauss newton
%
% fobj = fobj(0) + Delta h*dfobj/dh + 1/2(Delta h)^2*d^2fobj/dh^2 + hot
% dfobj/dDelta h = 0
% - df/dh = d^2f/dh^2*Delta h
% Delta h = -(d^2f/dh^2)^-1 df/dh

% dfobj/dh     = (dres/dh)'*res = 0
% d^2fobj/dh^2 = (dres/dh)'*(dres/dh) + _(d^2res/dh^2)*res_neglegtec

% Delta h = -((dres/dh)'*dresh/dh)^-1*(dres/dh)'*res
% Delta h = -(dres/dh)^-1*res

% central differences are not stable for sloping ground, even though they are
% suggested by gilad:
%  d(h(dh/dx - dzb/dx)
%  = (dh/dx)^2 + h d2h/dx^2 - dh/dx*dzb/dx - h/dzb/dx^2
%  = h d2h/dx^2 + (dh/dx-dzb/dx)*dh/dx - h*d2zb/dx^2
% for h small, and d2zb/dx^2 = 0, this is the advection equation
% -> -dzb/dx*dh/dx

%   (hr-hc)*(hr + hc) - (hc - hl)*(hl + hc)
%   hr^2 - hc^2 - hc^2 + hl^2
% vs
%   hc*(hl - 2 hc + hr) + (hr - hl)^2
%   
% (dzb/dx) > 0 (hr-hc)*(zr - zl)


syms h L I dt K hl hr hu hd dx dy zu zd zr zl z br bl


	  al  = -K*(h+z-hl-zl)/dx;
	  ar  = -K*(hr+zr-h-z)/dx;

	  dh_dt = -(   0.5*ar.*( (1-br).*h   + (1+br).*hr ) ...
	             - 0.5*al.*( (1-bl).*hl  + (1+bl).*h  ) ...
		   )/dx;

h_ = [hl,hr,h]
clear df
for idx=1:length(h_)
	h_(idx)
	df(idx,1) = (diff(dh_dt,h_(idx)))
	%dfz(idx) = diff(resz,h_(idx))
end


