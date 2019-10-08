% 2018-09-21 15:32:41.721250642 +0800
%% velocity shear
function du_dz = du_dz(obj,S,h)
	c = obj.param;
	z = S*h;
	du_dz = c(1)./z - c(3)*((6*z)/h^2 - 2/h);
end

