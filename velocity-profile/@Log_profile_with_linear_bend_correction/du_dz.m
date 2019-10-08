% 2018-09-21 15:25:30.396877971 +0800
%% velocity shear along vertical
function du_dz = du_dz(obj,S,h)
	c = obj.param;
	z = S*h;
	du_dz = (c(3) + c(1)./z + c(3)*log(z) + (c(2)*c(3))/c(1));
end
