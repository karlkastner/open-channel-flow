% Fri  4 Dec 13:59:40 +08 2020
function dzs_dr = transverse_slope_rozovskii(v,R,C)
% 1.40 : parabolic, m
% 1.59 : log, kappa, C
% 1.50 : elliptic, P
% 1.67 : power law (exp), n
	g = 9.81;
	kappa = 0.41;

% 1.78
a0 = 1 + g./(kappa.^2.*C.^2)
% 1.38 :  a0 = 1.03
dzs_dr = a0.*v.^2./(g*R)
 
end

