% Mon  4 Sep 18:03:16 CEST 2017
%
%% analytical solution to the gradually varied flow equation (backwater equation)
%% c.f. Chow, Bresse
function f = gvf_x_chow(u,yc,yn,C,S0)
	if (issym(u))
		syms g a;
		yc_div_yn_3 = a^(1/3);
	else
		g = 9.81;
		yc_div_yn_3 = C^2*S0/g;
	end

	% this is numerically unstable when u -> 1

	%f3 = 1/6*(log((u.^2+u+1)./(u-1).^2) + 1/sqrt(3)*atan((2*u+1)/sqrt(3)) - 1/sqrt(3)*atan(1/sqrt(3)));
	%f3 = 1/6*(log((u.^2+u+1)./(u-1).^2) - 1/sqrt(3)*acot((2*u+1)/sqrt(3)));
	f3 = 1/6*log((u.^2+u+1)./(u-1).^2) - 1/sqrt(3)*acot((2*u+1)/sqrt(3));
	%(2*u+1)/sqrt(3)) - 1/sqrt(3)*atan(1/sqrt(3)));
	f = -yn/S0*(u - (1-yc_div_yn_3)*f3);
end

