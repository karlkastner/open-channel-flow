% Sat 23 May 21:31:44 +08 2020
% T = 2
function [xy,uv,cc] = meander_centreline(Rcm,alpha,T,nt);
	if (length(T)<2)
		T = [0,T];
	end
	t   = linspace(T(1),T(2),nt)';
	cm  = 1./Rcm;
	cc  = (cm*sin(pi*t));
	um  = -pi*(alpha+pi/2)/cm;
	t_  = (cm*um*cos(pi*t))/pi;
	uv  = [um*cos(t_), um*sin(t_)];
	xy=[0,0;
            cumsum(mid(uv).*(t(2)-t(1)))];
end

