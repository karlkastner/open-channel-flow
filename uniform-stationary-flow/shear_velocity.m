% Tue 10 Dec 09:26:46 +08 2019
function us=shear_velocity(U,C)
	g  = Constant.gravity; 
	us = sqrt(g)./C.*U;
end

