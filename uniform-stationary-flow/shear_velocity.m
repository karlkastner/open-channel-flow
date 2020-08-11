% Tue 10 Dec 09:26:46 +08 2019
%function us=shear_velocity(U,C)
function us=shear_velocity(U,C)
	if (~issym(U))
	g  = Constant.gravity; 
	else
		syms g
	end
	us = sqrt(g)./C.*U;
end

