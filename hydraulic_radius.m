% 2020-01-06 17:18:52.703667971 +0800 hydraulic_radius.m
% vanoni-brooks
function R = hydraulic_radius(d,b)
	R     = (b.*d)./(b+2*d);
end
